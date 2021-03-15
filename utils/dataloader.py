"""
================================================================================
Functions for loading data and associated sample/feature annotation.
================================================================================
"""


import os
import re

import anndata
import pandas as pd
import scipy.sparse

import utils.settings as settings


CACHE_DIR = os.path.join(settings.DATA_DIR, '.cache')


def clear_cache():
    import shutil
    shutil.rmtree(CACHE_DIR)


#===============================================================================
# Load peak annotation
#===============================================================================
def load_gene_annotation(timepoint):
    """Load gene annotation

    Creates dataframe with peak index and gene / term annotation

    Args:
        timepoint (str): Timepoint for subsetting.

    Raises:
        FileNotFoundError: If closest-gene file was not computed.

    Returns:
        pd.DataFrame: Annotation dataframe.
    """
    fname = os.path.join(settings.RESOURCES_DIR, 'genes', 'peak_to_gene.txt')
    if not os.path.exists(fname):
        msg = 'Could not find %s.' % fname
        msg += '\nPlease run:'
        msg += '\nbedtools closest -d -a %s ' % settings.TOTAL_COUNTS_PEAKS
        msg += '-b %s ' % os.path.join(settings.RESOURCES_DIR, 'genes', 'expr.bodycoord.all.txt')
        msg += '| cut -f1,2,3,4,5,6,7,9 > %s' % fname
        raise FileNotFoundError(msg)
    peak_to_gene = pd.read_csv(
        fname, sep='\t', header=None).drop([3,4,5,7], axis=1)
    peak_to_gene.columns = ['chr', 'start', 'end', 'gene']

    stage_annotation = list()
    for stage in settings.TIMEPOINT_TO_STAGE[timepoint]:
        fname = os.path.join(
            settings.RESOURCES_DIR,
            'genes/expr.bodycoord.' + stage + '.txt')
        if os.path.exists(fname):
            anno = pd.read_csv(fname, sep='\t', header=None)
            anno.columns = ['chr', 'start', 'end', 'gene', 'term']
            anno = anno.set_index('gene')
            anno['term'] = anno['term'].str.split(',').str.join(';')
            stage_annotation.append(anno)

    # match peaks
    peak_to_gene['peak']= generate_var_data(peak_to_gene[['chr', 'start', 'end']]).index
    df = peak_to_gene[['peak', 'gene']]

    # replace gene name with term annotation
    def _get_terms(gene):
        terms = ''
        for anno in stage_annotation:
            if gene in anno.index:
                if terms != '':
                    terms += ';'
                terms += anno.loc[gene, 'term']
        if terms != '':
            terms = ';'.join(set(terms.split(';')))
        return terms

    df['term'] = df['gene'].apply(_get_terms)

    # merge multiple entries for each peak
    df = df.groupby('peak', as_index=False).agg(
        lambda x: ';'.join(set(x.astype(str))))
    df['term'] = df['term'].str.strip(';')

    # remove peaks without annotation
    df = df[df['term'] != ''].set_index('peak')
    return df


def load_enhancer_annotation(timepoint=None):
    """Load enhancer annotation.

    Creates dataframe with peak index and enhancer / term annotation.

    Args:
        timepoint (str, optional): Timepoint for subsetting
            annotation. Defaults to None.

    Raises:
        FileNotFoundError: If peak-enhancer overlap was not computed.

    Returns:
        pd.DataFrame: Annotation dataframe.
    """
    fname = os.path.join(settings.RESOURCES_DIR, 'enhancers', 'peak_overlaps.txt')
    if not os.path.exists(fname):
        msg = 'Could not find %s.' % fname
        msg += '\nPlease run:'
        msg += 'bedtools intersect '
        msg += '-a %s ' % settings.TOTAL_COUNTS_PEAKS
        msg += '-b  ' % os.path.join(settings.RESOURCES_DIR, 'enhancers', 'CAD4_plus_vienna_dm6.bed')
        msg += '-wb > %s' % fname
        raise FileNotFoundError(msg)
    peak_overlaps = pd.read_csv(
        fname, sep='\t', header=None).drop([7, 8], axis=1)

    enhancer_anno = pd.read_csv(
        os.path.join(
            settings.RESOURCES_DIR,
            'enhancers',
            'CAD4_plus_vienna_flat_annotations.txt'),
        sep='\t', encoding='ISO-8859-1')

    # replacement table to deal with encoding issues ...
    REPLACEMENT_TABLE = {
        'Mef2 IEd5': 'Mef2_IEd5',
        'eIF-4a': 'eIF-4a_p5′-436eIF4awt-luc/lacZ',
        'eya_meso': 'eya_meso_eya(1)',
        'eya_ecto': 'eya_ecto_eya(2)',
        '1421/': 'ush_−1421/−25'
    }
    for pat, repl in REPLACEMENT_TABLE.items():
        enhancer_anno['elementName'] = enhancer_anno['elementName'].apply(
            lambda x: repl if pat in x else x)
    enhancer_anno = enhancer_anno.set_index('elementName')

    # restrict to desired stage
    if timepoint is None or timepoint == 'all':
        stages = ['allTerms']
    else:
        stages = settings.TIMEPOINT_TO_STAGE[timepoint]
    enhancer_anno = enhancer_anno[stages]
    enhancer_anno[enhancer_anno.isna()] = ''

    # merge entries from all matched stages
    enhancer_anno = enhancer_anno.apply(
        lambda x: ';'.join(set(x.astype(str))), 1)

    peak_overlaps['peak'] = (
        peak_overlaps[0].map(str)
        + '_' + peak_overlaps[1].map(str)
        + '_' + peak_overlaps[2].map(str))

    # lookup peak annotation using enhancer label
    df = pd.concat(
        [peak_overlaps['peak'], enhancer_anno.loc[peak_overlaps[6]].reset_index()],
        axis=1)
    df.columns = ['peak', 'enhancer', 'term']

    # merge multiple entries for each peak
    df = df.groupby('peak', as_index=False).agg(
        lambda x: ';'.join(set(x.astype(str))))
    df['term'] = df['term'].str.strip(';')

    # remove peaks without annotation
    df = df[df['term'] != ''].set_index('peak')
    return df


def load_variant_info(
        exp_id,
        indels=False,
        filtering='stringent'):
    """Load variant info file.

    Args:
        exp_id (str): Experiment id.
        indels (bool, optional): Include indels. Defaults to False.
        filter (str, optional): Filtering, stringent or lenient.

    Returns:
        pandas.DataFrame: Variant info.
    """
    if filtering == 'lenient':
        df = pd.read_csv(settings.VARIANT_INFO_LENIENT(exp_id), sep='\t')
    elif filtering == 'stringent':
        df = pd.read_csv(settings.VARIANT_INFO_STRINGENT(exp_id), sep='\t')
    else:
        raise ValueError('Unkown filtering level %s' % filtering)

    if not indels:
        mask = (df['allele1'].str.len() == 1) & (df['allele2'].str.len() == 1)
        df = df.loc[mask]
    return df


#===============================================================================
# Load count data
#===============================================================================
def load_total_count_data(
        exp_ids='all',
        timepoints='all',
        cache=True):
    """Load counts for given experiment ids and timepoints.

    Args:
        exp_ids (list): List of experiment ids.
        timepoints (list): List of timepoints.
        cache (bool, optional): Cache files. Defaults to True.

    Returns:
        anndata.AnnData: AnnData count files.
    """
    exp_ids = parse_experiment_ids(exp_ids)
    timepoints = parse_timepoints(timepoints)

    adatas = []
    for exp_id in exp_ids:
        if exp_id == 'cusanovich': continue
        adatas.append(read_counts_f1(
            counts_file=settings.TOTAL_COUNTS(exp_id),
            region_file=settings.TOTAL_COUNTS_PEAKS,
            cell_index=settings.CELL_INDEX(exp_id),
            cache=cache))

    if 'cusanovich' in exp_ids:
        adatas.append(read_counts_cusanovich(cache=cache))

    adata = adatas[0].concatenate(
        adatas[1:],
        batch_key='exp_id',
        batch_categories=exp_ids)

    # make categorical
    adata.obs['exp_id'] = adata.obs['exp_id'].astype('category')
    adata.obs['timepoint'] = adata.obs['timepoint'].astype('category')
    adata.var['chr'] = adata.var['chr'].astype('category')

    return adata[adata.obs['timepoint'].isin(timepoints)].copy()


def load_allelic_count_data(
        exp_ids='f1',
        timepoints='all',
        wasp_corrected=True,
        regions='peaks',
        cache=True):
    """Load allelic counts for given experiment ids and timepoints.

    Args:
        exp_ids (list): List of experiment ids.
        timepoints (list): List of timepoints.
        wasp_corrected (bool): Use wasp corrected counts.
        regions (str): Region descriptor.
        cache (bool, optional): Cache files. Defaults to True.

    Returns:
        anndata.AnnData: AnnData file with allele1 counts and layer for
            allelic total.
    """
    exp_ids = parse_experiment_ids(exp_ids)
    timepoints = parse_timepoints(timepoints)

    if 'cusanovich' in exp_ids:
        raise ValueError('No allelic counts for cusanovich data.')

    adatas = []
    for exp_id in exp_ids:
        counts_dict = settings.COUNTS_DICT(exp_id, regions, wasp_corrected)
        adatas.append(read_allelic_counts_f1(
            allele1_file=settings.ALLELE1_COUNTS(exp_id),
            allele2_file=settings.ALLELE2_COUNTS(exp_id),
            region_file=settings.ALLELIC_COUNTS_PEAKS,
            cell_index=settings.CELL_INDEX(exp_id),
            cache=cache))

    adata = adatas[0].concatenate(
        adatas[1:],
        batch_key='exp_id',
        batch_categories=exp_ids)

    # make categorical
    adata.obs['exp_id'] = adata.obs['exp_id'].astype('category')
    adata.obs['timepoint'] = adata.obs['timepoint'].astype('category')
    adata.var['chr'] = adata.var['chr'].astype('category')

    return adata[adata.obs['timepoint'].isin(timepoints)].copy()


#===============================================================================
# Files to anndata
#===============================================================================
def read_counts_cusanovich(cache=True):
    """Load Cusanovich et al. data."""
    cache_file = CACHE_DIR + '/cusanovich.h5ad'
    if cache:
        print('Checking if file exists: ' + str(cache_file) + ' ... ', end='')
        if os.path.isfile(cache_file):
            print('found.\n -> loading cached anndata')
            return anndata.read(cache_file)
        else:
            print('not found.')
    print(' -> downloading ... ')

    adatas = []
    for timepoint in settings.TIMEPOINTS:
        url = settings.CUSANOVICH_SUMMIT_MATRIX(timepoint)
        X = pd.read_csv(url, sep='\t')

        obs_data = pd.DataFrame(X.columns[4:])
        obs_data['timepoint'] = timepoint
        obs_data = generate_obs_data(obs_data)

        region_file = settings.CUSANOVICH_PEAKS
        var_data = pd.read_csv(
            region_file,
            sep='\t', header=None)
        var_data = generate_var_data(var_data)

        if var_data.shape[0] != X.shape[0]:
            raise ValueError(
                'dm6 peak file does not have '
                'the same dimensions as ' + url)

        X = X[obs_data.index]
        X = scipy.sparse.csr_matrix(X.values.transpose())

        adatas.append(anndata.AnnData(X=X, obs=obs_data, var=var_data))

    adata = adatas[0].concatenate(
        adatas[1:],
        batch_key='timepoint',
        batch_categories=settings.TIMEPOINTS)

    adata.obs['timepoint'] = adata.obs['timepoint'].astype('category')
    adata.var['chr'] = adata.var['chr'].astype('category')

    if cache and not os.path.isfile(cache_file):
        print(' -> caching')
        adata.write(cache_file)
    return adata


def read_counts_f1(
        counts_file,
        region_file,
        cell_index,
        cache=True):
    """Create anndata from count matrix & annotation.

    Args:
        counts_file (str): Tab-seperated file with count data in COO format.
        region_file (str): File with regions (chr, start, end)
        cell_index (str): Cell index file with timepoints.
        cache (bool, optional): Cache anndata. Defaults to True.

    Returns:
        anndata.AnnData: AnnData file with count data.
    """
    file_prefix = re.sub('.txt|.gz', '', os.path.basename(counts_file))
    cache_file = CACHE_DIR + '/' + file_prefix + '.h5ad'
    if cache:
        print('Checking if file exists: ' + cache_file + ' ... ', end='')
        if os.path.isfile(cache_file):
            print('found.\n -> loading cached anndata ...')
            return anndata.read(cache_file)
        else:
            print('not found.')
    print(' -> creating anndata ...')

    X = pd.read_csv(counts_file, sep='\t', header=None)

    obs_data = pd.read_csv(
        cell_index,
        sep='\t', header=None)
    obs_data = generate_obs_data(obs_data)

    var_data = pd.read_csv(
        region_file,
        sep='\t', header=None)
    var_data = generate_var_data(var_data)

    n_cells = obs_data.shape[0]
    n_regions = var_data.shape[0]

    X = scipy.sparse.coo_matrix(
        (X[2], (X[1], X[0])),
        shape=(n_cells, n_regions)).tocsr()

    adata = anndata.AnnData(X=X, obs=obs_data, var=var_data)

    if cache and not os.path.isfile(cache_file):
        print(' -> caching')
        adata.write(cache_file)
    return adata


def read_allelic_counts_f1(
        allele1_file,
        allele2_file,
        region_file,
        cell_index,
        cache=True):
    """Create anndata for allelic counts files.

    Args:
        allele1_file (str): Tab-seperated file with count data in COO format.
        allele2_file (str): Tab-seperated file with count data in COO format.
        region_file (str): File with regions (chr, start, end)
        cell_index (str): Cell index file with timepoints.
        cache (bool, optional): Cache anndata. Defaults to True.

    Returns:
        anndata.AnnData: AnnData file with allele1 counts and layer for
            allelic total.
    """
    allele1 = read_counts_f1(allele1_file, region_file, cell_index, cache)
    allele2 = read_counts_f1(allele2_file, region_file, cell_index, cache)

    adata = allele1
    adata.layers['allelic_total'] = allele1.X + allele2.X
    return adata


#===============================================================================
# Utility functions
#===============================================================================
def parse_timepoints(timepoints):
    """Parse timepoints into desired format (list)."""
    if not isinstance(timepoints, list):
        timepoints = [timepoints]
    if len(timepoints) == 1 and timepoints[0] == 'all':
        timepoints = settings.TIMEPOINTS
    else:
        timepoints = list(map(format_timepoint, timepoints))
    return timepoints


def format_timepoint(timepoint):
    """Set to default timepoint format.

    Args:
        timepoint (str): Timepoint string, e.g. '2-4hr'.

    Returns:
        str: Adjusted timepoint string.
    """
    return 'to'.join(re.findall('[0-9]+', timepoint))


def parse_experiment_ids(exp_ids):
    """Parse experiment ids into desired format (list)."""
    if not isinstance(exp_ids, list):
        exp_ids = [exp_ids]

    if len(exp_ids) == 1 and exp_ids[0] == 'all':
        exp_ids = settings.F1_EXP_IDS + ['cusanovich']
    if len(exp_ids) == 1 and exp_ids[0] == 'f1':
        exp_ids = settings.F1_EXP_IDS
    return exp_ids


def generate_var_data(region_data):
    """Generate variable info.

    Args:
        region_data (pandas.DataFrame): Dataframe with 3 columns corresponding
        to chromosome, start and end.

    Returns:
        pandas.DataFrame: Formatted variable info with region length.
    """
    region_data.columns = ['chr', 'start', 'end']
    region_data['chr'] = 'chr' + \
        region_data['chr'].str.replace(r'Chr|chr', '')
    region_data['chr'] = region_data['chr'].astype('category')
    region_data.index = region_data.apply(
        lambda x: '_'.join(x.map(str)), 1)

    region_data[['start', 'end']] = region_data[[
        'start', 'end']].apply(pd.to_numeric)
    region_data['length'] = region_data['end'] - region_data['start'] + 1
    return region_data


def generate_obs_data(cell_index):
    """Generate observation info.

    Args:
        cell_index (pandas.DataFrame): Dataframe with 3 columns corresponding
            to cell id and timepoint.

    Returns:
        pandas.DataFrame: Formatted observation info.
    """
    cell_index.columns = ['cell', 'timepoint']
    cell_index['timepoint'] = cell_index['timepoint'].apply( format_timepoint)
    cell_index['timepoint'] = cell_index['timepoint'].astype('category')
    cell_index['timepoint'] = cell_index['timepoint'].cat.set_categories(settings.TIMEPOINTS)

    cell_index.set_index('cell', inplace=True)
    return cell_index

