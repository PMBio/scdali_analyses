import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def plot_rate_umap(
        adata,
        peak,
        plot_covered_only=False,
        show_colorbar=True,
        ax=None):
    """Visualize GP posterior means on top of UMAP.

    Requires that adata contains the following fields:

    Args:
        adata: Anndata object with the following fields:
            - adata.obsm['X_umap']
            - adata.layers['gp_post_mean']
            - adata.layers['allelic_total']
        peak: Element from adata.var_names.
        plot_covered_only: Only plot cells with nonzero
            total counts. If False, covered cells are
            indicated by black dots.
        show_colorbar: Show colorbar.
        ax: matplotlib.axes.Axes.

    """
    rasterized = True

    if ax is None:
        ax = plt.gca()

    adata = adata[:, peak]
    ids_covered = adata.layers['allelic_total'] > 0
    if plot_covered_only:
        # restrict to covered cells
        adata = adata[ids_covered, :]

    p = ax.scatter(
        adata.obsm['X_umap'][:, 0],
        adata.obsm['X_umap'][:, 1],
        c=adata.layers['gp_post_mean'],
        vmin=0, vmax=1, rasterized=rasterized,
        s=plt.rcParams['lines.markersize'] ** 2 * 2,
        cmap=plt.cm.coolwarm)

    if not plot_covered_only:
        ax.scatter(
            adata[ids_covered, :].obsm['X_umap'][:, 0],
            adata[ids_covered, :].obsm['X_umap'][:, 1],
            rasterized=rasterized,
            c='black')
    ax.set(xlabel='UMAP1', ylabel='UMAP2', xticks=[], yticks=[])
    if show_colorbar:
        plt.colorbar(p, ax=ax, label='Estimated rate', orientation='horizontal')
    sns.despine(ax=ax)
    return ax


def plot_rate_time(
        adata,
        peak,
        plot_covered_only=True,
        fit_reg=False,
        c='black',
        ax=None):
    """Visualize GP posterior means across time.

    Requires that adata contains the following fields:

    Args:
        adata: Anndata object with the following fields:
            - adata.obs['time_vae']
            - adata.layers['gp_post_mean']
            - adata.layers['allelic_total']
        peak: Element from adata.var_names.
        plot_covered_only: Only plot cells with nonzero
            total counts. If False, covered cells are
            indicated by black dots.
        add_regplot: Add regplot.
        ax: matplotlib.axes.Axes.

    """
    rasterized = True

    if ax is None:
        ax = plt.gca()

    adata = adata[:, peak]
    ids_covered = adata.layers['allelic_total'] > 0
    if plot_covered_only:
        # restrict to covered cells
        adata = adata[ids_covered, :]

    sns.regplot(
        x=adata.obs['time_vae'],
        y=adata.layers['gp_post_mean'],
        fit_reg=fit_reg,
        order=3, scatter=True,
        n_boot=500, color='orange',
        scatter_kws={'color': c, 'rasterized': rasterized},
        ax=ax)
    ax.set(
        xlabel='VAE Pseudo-time',
        ylabel='Estimated rate',
        ylim=(0, 1))
    sns.despine(ax=ax)
    return ax


def alpha_threshold_closure(levels, alphas):
    if (len(levels) - 1) != len(alphas):
        raise ValueError('Number of intervals needs to match number of alphas.')
    def alpha_by_nobs(nobs):
        for i in range(len(alphas)):
            if (nobs >= levels[i]) and (nobs < levels[i+1]):
                return alphas[i]
        raise ValueError('Invalid nobs.')
    return alpha_by_nobs


def plot_rate_by_group(
        adata,
        peak,
        group_key,
        group_colors,
        plot_covered_only=True,
        alpha_by_nobs=alpha_threshold_closure(
            [0, 20, 50, np.inf], [0.05, 0.3, 1]),
        ax=None):
    """Visualize GP posterior means by group.

    Requires that adata contains the following fields:

    Args:
        adata: Anndata object with the following fields:
            - adata.obs[group_key]
            - adata.layers['gp_post_mean']
            - adata.layers['allelic_total']
            - addata.var['gp_mean']
        peak: Element from adata.var_names.
        group_key: Key in adata.obs.
        group_colors: Dict mapping elements from adata.obs[group_key]
            to colors.
        plot_covered_only: Only plot cells with nonzero
            total counts. If False, covered cells are
            indicated by black dots.
        alpha_by_nobs: Function mapping integers to values in (0, 1)
        ax: matplotlib.axes.Axes.

    """
    rasterized = True

    if ax is None:
        ax = plt.gca()

    adata = adata[:, peak]
    ids_covered = adata.layers['allelic_total'] > 0
    if plot_covered_only:
        # restrict to covered cells
        adata = adata[ids_covered, :]

    df = pd.DataFrame(adata.obs[group_key])
    df['gp_posterior'] = adata.layers['gp_post_mean']

    ax.axhline(
        adata.var['gp_mean'].item(),
        c='red',
        linestyle='--',
        alpha=.5)
    ax.axhline(
        .5,
        c='grey',
        alpha=.5)
    sns.violinplot(
        x=group_key,
        y='gp_posterior',
        data=df,
        palette=group_colors,
        ax=ax)

    nobs = df.groupby(group_key)['gp_posterior'].size().values
    for g, violin in enumerate(ax.collections[::2]):
        violin.set_alpha(alpha_by_nobs(nobs[g]))
    sns.stripplot(
        x=group_key,
        y='gp_posterior',
        data=df,
        color='black',
        size=plt.rcParams['lines.markersize'],
        rasterized=rasterized,
        ax=ax)
    ax.set(
        xticklabels=[s.get_text().replace(' ', '\n') for s in ax.get_xticklabels()],
        ylabel='Estimated rate',
        ylim=(0, 1))
    ax.tick_params(axis='x', which='both', length=0)
    sns.despine(ax=ax, bottom=True)
    return ax


def plot_rate_quantiles(
        adata,
        peak,
        group_key,
        group_colors,
        plot_covered_only=True,
        annotate=False,
        q=.05):
    """Quantile plot of GP posterior means.

    Requires that adata contains the following fields:

    Args:
        adata: Anndata object with the following fields:
            - adata.obs[group_key]
            - adata.layers['gp_post_mean']
            - adata.layers['allelic_total']
            - addata.var['gp_mean']
        peak: Element from adata.var_names.
        group_key: Key in adata.obs.
        group_colors: Dict mapping elements from adata.obs[group_key]
            to colors.
        plot_covered_only: Only plot cells with nonzero
            total counts. If False, covered cells are
            indicated by black dots.
        annotate: Annotate plot.
        q: Quantile to indiciate.
    """
    adata = adata[:, peak]
    ids_covered = adata.layers['allelic_total'] > 0
    if plot_covered_only:
        # restrict to covered cells
        adata = adata[ids_covered, :]

    df = pd.DataFrame(adata.obs[group_key])
    df['gp_posterior'] = adata.layers['gp_post_mean']
    df = df.sort_values('gp_posterior').reset_index()
    df['order'] = range(df.shape[0])

    g = sns.JointGrid(height=4, ylim=(0, 1), ratio=2)
    sns.lineplot(
        x='order',
        y='gp_posterior',
        color='dimgray',
        linewidth = 1.5 * plt.rcParams['lines.linewidth'],
        data=df,
        ax=g.ax_joint)

    qp = np.quantile(df['gp_posterior'], q)
    qq = np.quantile(df['gp_posterior'], 1-q)
    g.ax_joint.axvline(q * df.shape[0], linestyle='--', color='gray')
    g.ax_joint.axvline((1-q) * df.shape[0], linestyle='--', color='gray')
    g.ax_joint.axhline(qp, linestyle='--', color='gray')
    g.ax_joint.axhline(qq, linestyle='--', color='gray')

    sns.kdeplot(
        x='order',
        hue=group_key,
        clip=(0, df.shape[0]),
        palette=group_colors,
        common_norm=False,
        fill=True,
        data=df,
        ax=g.ax_marg_x)
    # modify legend
    l = g.ax_marg_x.get_children()[-2]
    for t in l.get_texts():
        t.set_text(t.get_text().replace(' ', '\n'))
    l.set_bbox_to_anchor((1.5, 1, 0, 0))
    l.set_title(group_key)

    if annotate:
        qdiffx = .95 * g.ax_joint.get_xlim()[1]
        g.ax_joint.annotate(
            s='', xy=(qdiffx, qp), xytext=(qdiffx, qq),
            arrowprops=dict(arrowstyle='<->', color='red'))
        g.ax_joint.annotate(
            s='%d%% quantile \ndifference \n(effect size)' % (100 * q),
            xy=(qdiffx+25, qp + .8*(qq - qp)/2))

    g.ax_marg_y.remove()
    g.set_axis_labels(xlabel='Cell order', ylabel='Estimated rate')

    g.fig.set_figwidth(7)
    g.fig.set_figheight(4)
    return g


def format_peak(peak):
    peak_split = peak.split('_')
    return peak_split[0] + ':' + '-'.join(peak_split[1:])


def lighten_color(color, amount=0.5):
    """From https://gist.github.com/ihincks.

    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])
