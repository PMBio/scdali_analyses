"""Utility functions operating on anndata objects."""
import numpy as np
import pandas as pd
import scipy.sparse
from scipy.stats import fisher_exact


def binarize(adata):
    """Binarize adata.X."""
    adata.X = (adata.X > 0).astype(np.int_)


def filter_chromosomes(adata, chrs=['chrX', 'chrY']):
    """Filter out chromosomes in place.

    Args:
        adata (anndata.AnnData): AnnData object.
        chrs (list): List of chromosomes to filter.
    """
    filter_chrs = ~ adata.var.chr.str.contains('|'.join(chrs))
    adata._inplace_subset_var(filter_chrs)


def apply_to_layer(adata, layer, func):
    """Apply function to adata, treating adata.layers[layer] as adata.X.

    Args:
        adata (anndata.AnnData): AnnData object.
        layer (str): adata.layers key
        func (function): Function that modifies adata inplace.

    Returns:
        adata (anndata.AnnData): Modified AnnData object.
    """
    # save adata.X and set layer element as new X
    adata = adata.copy()
    if 'X_tmp' in adata.layers.keys():
        raise ValueError('X_tmp must not be a key in adata.layers.')
    adata.layers['X_tmp'] = adata.X
    adata.X = adata.layers[layer]
    # apply function
    func(adata)
    # swap X and layers[layer]
    adata.layers[layer] = adata.X
    adata.X = adata.layers['X_tmp']

    del adata.layers['X_tmp']
    return adata


def drop_obsm_cols(adata, key='X_pca', comps=[]):
    """Moves columns to obs (in place).

    Remove columns from obsm and place them in obs

    Args:
        adata (anndata.AnnData): Annotated data matrix.
        key (str): Key in adata.obsm.
        comps (list): Columns to drop.
    """
    if not isinstance(comps, list):
        comps = [comps]
    ids = [i for i in range(adata.obsm[key].shape[1]) if i not in comps]
    comps_to_remove = adata.obsm[key][:, comps]
    adata.obsm[key] = adata.obsm[key][:, ids]
    adata.obs[['X_pca_' + str(i) for i in comps]] = pd.DataFrame(comps_to_remove)


def calculate_tfidf(adata):
    """Compute tf idf transformation in place."""
    from sklearn.feature_extraction.text import TfidfTransformer
    adata.X = TfidfTransformer().fit_transform(adata.X)


def compute_obsm_count_corr(adata, key='X_pca'):
    """Compute Pearson correlation with total counts (in place).

    Compute correlation with total counts and number of detected features in
    each cell. Requires to run calculate_qc_metrics first.

    Args:
        adata (anndata.AnnData): Annotated data matrix.
        key (str): Key in adata.obsm.

    Raises:
        ValueError: Requires running calculate_qc_metrics first.
        ValueError: Could not find key.
    """
    if ('n_genes_by_counts' not in adata.obs_keys() or
        'total_counts' not in adata.obs_keys()):
        raise ValueError('Could not find general qc metrics.')
    if key not in adata.obsm.keys():
        raise ValueError('Could not find pca components.')

    uns_key = key.split('_')[1]

    adata.uns[uns_key]['corr_n_genes_by_counts'] = _corr_coef(
        adata.obs['n_genes_by_counts'].values, adata.obsm[key])
    adata.uns[uns_key]['corr_total_counts'] = _corr_coef(
            adata.obs['total_counts'].values, adata.obsm[key])


def _corr_coef(X, v):
    """Computes correlation between columns of X and v."""
    if v.ndim == 1:
        v = v.reshape(-1, 1)

    X_mX = X - X.mean(0)
    v_mv = v - v.mean(0)

    ssX = (X_mX**2).sum(0)
    ssv = (v_mv**2).sum(0)

    return np.dot(v_mv.T, X_mX).squeeze() / np.sqrt(ssv * ssX)


def compute_enrichment(
        adata,
        terms,
        top_frac=1.,
        score_key='pvals_adj',
        score_threshold=None,
        n_terms=10,
        keep_unknown=False):
    """Perform enrichment analysis based on Fisher exact test.

    Args:
        adata (anndata.AnnData): AnnData object with rank_genes_groups.
        terms (pandas.Series): Series with ';'-separated annotation terms,
            where the index is a subset of adata.var_names.
        top_frac (float): Fraction of top markers used for positive set.
            Ignored when score_threshold is set.
        score_key (str): Key in adata.uns['rank_genes_groups'] used for
            ranking features.
        score_threshold (float): Cutoff for top results.
        n_terms (int): Number of top terms to report.
        keep_unknown (bool): Keep peaks with no known annotation.
    """
    markers = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
    markers_scores = pd.DataFrame(adata.uns['rank_genes_groups'][score_key])

    if top_frac < 0 or top_frac > 1:
        raise ValueError('top_frac has to be between 0 and 1.')

    if score_threshold is None:
        markers = markers.head(int(markers.shape[0] * top_frac))

    if not 'rank_genes_groups' in adata.uns.keys():
        raise ValueError(
            'Please identify markers first using rank_genes_groups.')

    if adata.uns['rank_genes_groups']['params']['use_raw']:
        all_peaks = adata.raw.var_names
    else:
        all_peaks = adata.var_names

    all_terms = set(_terms_to_list(terms))
    n_terms = min(len(all_terms), n_terms)

    results = pd.DataFrame(1., index=all_terms, columns=markers.columns)
    for group in markers.columns:
        group_markers = markers[group].tolist()

        if score_threshold is not None:
            n_markers = (markers_scores[group] < score_threshold).sum()
            group_markers = group_markers[:n_markers]

        background = all_peaks.difference(group_markers).tolist()

        n_markers = len(group_markers)
        n_background = len(background)
        group_markers = terms.index.intersection(group_markers).tolist()
        background = terms.index.intersection(background).tolist()

        if not keep_unknown:
            n_markers = len(group_markers)
            n_background = len(background)

        # run test for each term
        for term in all_terms:
            # count term occurences in markers and background
            term_count_marker = terms[group_markers].str.contains(term).sum()
            term_count_bg = terms[background].str.contains(term).sum()
            other_count_marker = n_markers - term_count_marker
            other_count_bg = n_background - term_count_bg

            table = np.array([
                [term_count_marker, term_count_bg],
                [other_count_marker, other_count_bg]])
            if table[0, 0] > 0:
                _, pval = fisher_exact(table, alternative='greater')
                results.loc[term, group] = pval

        adata.uns['enrichment'] = dict()
        adata.uns['enrichment']['terms'] = results.apply(lambda x: x.nsmallest(
            n_terms).index).to_records(index=False, column_dtypes='U75')
        adata.uns['enrichment']['pvals'] = results.apply(lambda x: x.nsmallest(
            n_terms).values).to_records(index=False, column_dtypes='float64')


def _terms_to_list(terms):
    """Transforms pandas.Series of ';'-separated terms into
    list of substrings.

    Strips whitespaces and removes empty strings.

    Args:
        terms (pandas.Series): Series of ';'-separated strings.

    Returns:
        list: List of substrings.
    """
    return [x.strip() for x in filter(None, terms.str.split(';').sum())]
