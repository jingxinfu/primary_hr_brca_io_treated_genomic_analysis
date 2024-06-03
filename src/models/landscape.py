from typing import List
import pandas as pd
import click
from comut import comut
import matplotlib
import matplotlib.pyplot as plt
from ..dataset import bulkDataset
from ..settings import init_logger
from ..settings import COLOR_PAlETTE, PlotStyle, FigureDir

logger = init_logger(__name__)
plt.style.use(PlotStyle)


def plot_comut(
    mutation_data: pd.DataFrame,
    metadata: pd.DataFrame,
    visual_category_columns: List[str],
    visual_continuous_columns:List[str],
    mutation_signature_columns: List[str],
    sample_order_by_columns: List[str] = [
        'RCB',
        'Treatment_Arm',
        'Patient',
        'Timepoint',
    ],
    n_genes=25,
    figsize=(10, 10),
    wspace=0.3,
    x_padding=0.04,
    y_padding=0.04,
    tri_padding=0.03,
):
    # select overlapped samples
    overlap_samples = metadata.index.intersection(
        mutation_data['sample']
    ).drop_duplicates()

    mutation_data = mutation_data.loc[mutation_data['sample'].isin(overlap_samples), :]
    metadata = metadata.loc[overlap_samples, :]

    # Order samples
    samples = metadata.sort_values(sample_order_by_columns).index

    # reset index
    metadata.index.name = 'sample'
    metadata.reset_index(inplace=True)
    indicator_df = metadata[['sample', 'Patient']].rename(columns={'Patient': 'group'})
    indicator_df['group'] = indicator_df['group'].astype('category').cat.codes

    # Rank and select genes
    # genes = (
    #     mutSigCV.sort_values(by='q', ascending=False)
    #     .index.intersection(mutation_data['category'])
    #     .tolist()
    # )
    mut_freqs = (
        mutation_data[['sample','category']]
        .drop_duplicates()['category']
        .value_counts()
        .rename('Mutated samples')
        .reset_index()
        )
    mut_freqs.columns = ['category','Mutated samples']
    genes = mut_freqs['category'].tolist()
    genes = genes[:n_genes] if len(genes) > n_genes else genes
    genes.reverse()
    mutation_data = mutation_data.loc[mutation_data['category'].isin(genes),:]
    mut_freqs = mut_freqs.loc[mut_freqs['category'].isin(genes),:]
    # --------------------Comut Plotting
    comut_plot = comut.CoMut()
    comut_plot.samples = samples

    #add indicators first, since they will be at the bottom
    indicator_kwargs = {
        'color': 'black',
        'marker': 'o',
        'linewidth': 1,
        'markersize': 5,
    }
    comut_plot.add_sample_indicators(
        indicator_df.loc[indicator_df['sample'].isin(samples), :],
        name='Same patient',
        plot_kwargs=indicator_kwargs,
    )

    # if 'facets_wgd_bool' in metadata.columns:
    #     name = 'facets_wgd_bool'
    #     cat_df = metadata.melt(
    #         id_vars=['sample'], value_vars=[name], var_name='category'
    #     )
    #     comut_plot.add_categorical_data(
    #         cat_df, name='Whole Genome Duplication', mapping=COLOR_PAlETTE[name]
    #     )

    # then add somatic mutation information
    name = 'Mutation type'
    comut_plot.add_categorical_data(
        mutation_data,
        name=name,
        category_order=genes,
        mapping=COLOR_PAlETTE[name],
        tick_style='italic',
    )
    # Add mutation frequency
    side_mapping = {'Mutated samples': 'darkgrey'}
    bar_kwargs = {'height': 0.8}
    comut_plot.add_side_bar_data(
        mut_freqs,
        paired_name='Mutation type',
        name='Mutated samples',
        xlabel='# of Mutated samples',
        position='left',
        bar_kwargs=bar_kwargs,
        mapping=side_mapping,
    )
    
    # add other category rows on top of the somatic mutation plot
    for col_name in visual_category_columns:
        cat_df = metadata.melt(
            id_vars=['sample'], value_vars=[col_name], var_name='category'
        )
        cat_df['category'] = cat_df['category'].map(
            lambda x:x.replace('bulkRNA_','').replace('WES_',"").replace('_',' ').replace('facets wgd bool','Whole Genome Duplication').replace('er status','ER Status')
            )
        comut_plot.add_categorical_data(
            cat_df, name=col_name.replace('bulkRNA_','').replace('WES_',"").replace('_',' ').replace('facets wgd bool','Whole Genome Duplication').replace('er status','ER Status'), mapping=COLOR_PAlETTE[col_name]
        )
        
    # add other continuous rows on top of the somatic mutation plot
    for col_name in visual_continuous_columns:
        cont_df = metadata.melt(
            id_vars=['sample'], value_vars=[col_name], var_name='category'
        )
        cont_df['category'] =cont_df['category'].map(lambda x:x.replace('bulkRNA_','').replace('WES_',"").replace('_',' ').replace('absolute purity','Tumor Purity'))
        comut_plot.add_continuous_data(
            cont_df, name=col_name.replace('bulkRNA_','').replace('WES_',"").replace('_',' ').replace('absolute purity','Tumor Purity'), mapping=COLOR_PAlETTE[col_name]['color'],
            cat_mapping = {'Absent':{'facecolor':COLOR_PAlETTE[col_name]['N/A']}}
        )

    ## add mutation signature
    if len(mutation_signature_columns) != 0:
        mutation_signature_df = metadata[['sample'] + mutation_signature_columns]
        mutation_signature_df.columns = mutation_signature_df.columns.str.replace('WES_','')

        comut_plot.add_bar_data(
            mutation_signature_df,
            name='Mutational signatures',
            mapping=COLOR_PAlETTE['Mutation_Signature'],
            stacked=True,
            ylabel='Mutational signatures',
        )
    bar_kwargs = {'width': 0.8, 'edgecolor': 'black'}

    ## add TMB
    tmb_df = metadata[
        ['sample'] + [k for k in COLOR_PAlETTE['Mutation clonality'] if k != 'N/A']
    ]
    comut_plot.add_bar_data(
        tmb_df,
        name='Mutation clonality',
        mapping=COLOR_PAlETTE['Mutation clonality'],
        stacked=True,
        bar_kwargs=bar_kwargs,
        ylabel='Muts/Mb',
    )

    # plot comut and add unified legend
    comut_plot.plot_comut(
        x_padding=x_padding,
        y_padding=y_padding,
        tri_padding=tri_padding,
        figsize=figsize,
        hspace=0.03,
        wspace=wspace,
    )
    percentages = (mut_freqs['Mutated samples']/len(samples)*100).round(1).astype(str) + '%'
    # # set location of yticks
    for i,s in enumerate(percentages[::-1]):
        comut_plot.axes['Mutated samples'].text(x=0,y=i+.5,s=s,va='center',ha='right')

    # comut_plot.axes['MutsigQ'].axvline(
    #     1, color='black', linestyle='dotted', linewidth=2
    # )
    comut_plot.axes['Mutation clonality'].axhline(
        10, color='black', linestyle='dotted', linewidth=2
    )
    comut_plot.axes['Mutation clonality'].text(x=2,
        y=10.5, s='TMB = 10'
    )
        
    comut_plot.axes['Mutation clonality'].spines.left.set_visible(True)
    # For the white color label.
    borders = ['N/A']
    border_white = ['False']
    comut_plot.axes['Same patient'].set_xticklabels(metadata.set_index('sample').loc[samples,'Patient'].tolist())

    # comut_plot.axes['Mutation type'].set_xticklabels([])
    comut_plot.add_unified_legend(
        bbox_to_anchor=(1.3, 1.2), frameon=False, border_white=border_white,ignored_values=['Same patient']
    )  # ,borders=borders)

    if len(mutation_signature_columns) != 0:
        # delete the ytick labels for the mutational signature plot. It's implied the scale is 0 - 1. Also add a ylabel
        comut_plot.axes['Mutational signatures'].tick_params(
            axis='y', which='both', length=0, labelleft=False
        )
        comut_plot.axes['Mutational signatures'].set_ylabel(
            'Mutational\nsignatures', rotation='horizontal', ha='right', va='center'
        )
        
    # color bars must be added manually based on figure coordinates - [left, bottom, width, height]
    pfs_ax =comut_plot.figure.add_axes([.85, 0, 0.08, 0.014])
    # purity ranges 
    norm = matplotlib.colors.Normalize(vmin=metadata[visual_continuous_columns[0]].min(), 
                                       vmax=metadata[visual_continuous_columns[0]].max())

    # create the colorbar with colormap used when the continuous data was added (purp_7)
    pfs_colorbar = comut_plot.figure.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap='Purples'),
                                                    cax=pfs_ax, orientation='horizontal')

    # remove tick marks and move tick labels slightly inwards. Also remove black border
    pfs_colorbar.ax.tick_params(size=0)
    pfs_colorbar.outline.set_visible(False)
    # set title of colorbar to line up with other legend elements
    pfs_colorbar.set_label('Tumor Purity',labelpad = -40, x = .5, fontsize = 12)
    
    return comut_plot


@click.command()
@click.argument('data_processed_folder', type=click.Path(exists=True))
def main(data_processed_folder):
    dataset = bulkDataset(data_processed_folder)
    metadata=dataset.metadata
    metadata = metadata.loc[(metadata.WES_combined_TiN<=0.3)&
                            (metadata.WES_combined_fracContam<=0.04)&
                            (metadata.WES_Profile==True)#&(metadata.Timepoint=='Baseline')
                            ,:]
    timepoint ='Baseline'
    metadata.rename(columns={'er_status':'HR status','BluePrint':'Tumor type'},inplace=True)
    metadata['HR status'] = metadata['HR status'].replace(
        {
        "Weakly Positive (1-10% cell staining)":"HR-low positive",
        "Positive (>10% cell staining)":"HR positive"
        }
    )
    comut_plot = plot_comut(
        mutation_data=dataset.load('somatic_mutation'),
        metadata=metadata.loc[metadata.Timepoint==timepoint,:],
        visual_category_columns=[
            'Treatment_Arm',
            'RCB',
            'Tumor type',
            'HR status',
            'WES_facets_wgd_bool'
        ],
        visual_continuous_columns=[
            'WES_absolute_purity'
            ],
        mutation_signature_columns=[
            'WES_'+x for x in COLOR_PAlETTE['Mutation_Signature'].keys() if x != 'N/A'
        ],
        figsize=(25, 10),
        wspace=0.1,
        n_genes=15,
    )
    n_samples = (metadata.Timepoint==timepoint).sum()
    comut_plot.figure.suptitle(f"{timepoint} (N={ n_samples})")
    comut_plot.figure.savefig(
        FigureDir / f'Somatic_Mutation_CoMutPlot_{timepoint}.svg', bbox_inches='tight', dpi=300
    )