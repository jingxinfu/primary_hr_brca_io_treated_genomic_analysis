from typing import List, Union
import click
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import distinctipy
from .utils import bulkAnalysis, getDEGs
from ..dataset import bulkDataset
from ..settings import init_logger, PlotStyle, FigureDir

plt.style.use(PlotStyle)
logger = init_logger(__name__)


class rnaAnalysis(bulkAnalysis):
    def __init__(
        self,
        count: pd.DataFrame,
        tpm: pd.DataFrame,
        metadata: pd.DataFrame,
        *args,
        **kwargs,
    ):
        # transpose the gene expression matrix since they are often being written as gene by sample matrix.
        self.tpm = tpm.T
        super().__init__(
            feature_mtrx=count.T,
            metadata=metadata,
            min_feature_variance=None,
            significance_column='padj',
            statistics_column='log2FC',
            statistics_threshold=np.log2(1.5),
            *args,
            **kwargs,
        )

    def test_function(
        self,
        X: pd.DataFrame,
        y: pd.Series,
        covars: Union[pd.Series, pd.DataFrame, None],
       *args,
       **kwargs,
    ) -> pd.DataFrame:
        """
        Find the features in the X that associates with the

        Parameters
        ----------
        X : pd.DataFrame
            independent variables (sample,feature)
        y : pd.Series
            dependent variables (sample,)
        covars : Union[pd.Series,pd.DataFrame,None]
            confounding variables (sample,confounding variable)

        Returns
        -------
        pd.DataFrame
            association analysis results for each feature
        """
        X_tpm = self.tpm.loc[X.index, :]
        result = getDEGs(count=X.T, tpm=X_tpm.T, y=y)
        return result

    def analyze_batch_effect(self, batch: List[str]):
        from sklearn.preprocessing import StandardScaler
        from sklearn.decomposition import PCA

        # It is necessary to normalize data before performing PCA.
        # The PCA calculates a new projection of your data set.
        # And the new axis are based on the standard deviation of your variables.
        # So a variable with a high standard deviatin will have a higer weight for the
        # calculation of axis than a variable with a low standard deviation.
        # if your normalize your data, all variables have the same standard deviation,
        # thus all variables have the same weight and your PCA calculates relevant axis.
        scaler = StandardScaler()
        X = scaler.fit_transform(self.tpm)
        pca = PCA(n_components=2)
        pca.fit(X)
        X_pca = pca.transform(X)
        pca_names = [
            f'PC{i+1}({pca.explained_variance_ratio_[i]*100:.1f}%)' for i in range(2)
        ]
        X_pca = pd.DataFrame(X_pca, columns=pca_names, index=self.tpm.index)
        X_pca = X_pca.merge(
            self.metadata, left_index=True, right_index=True, how='left'
        )

        fig, axs = plt.subplots(len(batch), 2, figsize=(8, 4 * len(batch)))
        for i, b in enumerate(batch):
            ax = axs[i, 0]
            data = X_pca[pca_names + [b]]
            if X_pca[b].dtype == 'O':
                N = data[b].nunique()
                colors = distinctipy.get_colors(N)
                sns.scatterplot(
                    data=data,
                    x=pca_names[0],
                    y=pca_names[1],
                    hue=b,
                    ax=ax,
                    palette=colors,
                )
            else:
                sns.scatterplot(data=data, x=pca_names[0], y=pca_names[1], hue=b, ax=ax)
            ax.legend(loc=(-0.8, 0))

            ## quantative visualization
            ax = axs[i, 1]
            data = data.melt(id_vars=b, value_vars=pca_names)
            if X_pca[b].dtype == 'O':
                sns.boxplot(
                    data=data, x='variable', y='value', hue=b, ax=ax, palette=colors
                )
                ax.legend(loc=(1.05, 0.5), title=b, ncol=int(np.ceil(N / 10)))
                ax.set(xlabel='', ylabel='')
            else:
                from scipy.stats import pearsonr

                sns.scatterplot(data=data, x=b, y='value', hue='variable', ax=ax)
                ax.legend(loc=(1.05, 0), title=b)
                ax.set(xlabel='', ylabel='')
                cor_p = (
                    data.dropna()
                    .groupby('variable')
                    .apply(lambda v: pearsonr(v['value'], v[b]))
                    .apply(pd.Series)
                )
                cor_p.columns = ['cor', 'p']
                for i, (name, row) in enumerate(cor_p.iterrows()):
                    text = '{}:Rho={:.2f},Pvalue={:.3f}'.format(
                        name, row['cor'], row['p']
                    )
                    color = 'red' if row['p'] < 0.05 else 'black'

                    ax.text(1, 1 - i * 0.1, text, color=color, transform=ax.transAxes)
        plt.subplots_adjust(hspace=0.3)
        return fig


@click.command()
@click.argument('data_processed_folder', type=click.Path(exists=True))
def main(data_processed_folder):
    dataset = bulkDataset(data_processed_folder)
    count = dataset.load('gene_counts', index_col=0)
    tpm = dataset.load('TPM', index_col=0)
    metadata = dataset.metadata
    analyzer = rnaAnalysis(
        count=count,
        tpm=tpm,
        metadata=metadata.loc[metadata.BulkRNA_Profile, :],
        confounders=None,
    )
    output_prefix = f"{data_processed_folder}/bulkRNA"

    logger.info('-------------------- Evaluate the batch effect --------------------')
    fig = analyzer.analyze_batch_effect(
        batch=[
            'bulkRNA_Sequence_Batch',
            'Patient',
            'Treatment_Arm',
        ]
    )
    fig.savefig(f'{FigureDir}/bulkRNA_batch.pdf')
    plt.close('all')

    logger.info(
        '--------------------Calculate clinical benefits of gene expression--------------------'
    )
    analyzer.analyze_arm_response_difference(
        output_prefix=output_prefix, pathway_analysis=True
    )

    # Additional analysis on samples without AC treatments.
    logger.info(
        '-------------------Additional analysis on samples without AC treatments.--------------------'
        )
    analyzer = rnaAnalysis(
        count=count,
        tpm=tpm,
        metadata=metadata.loc[(metadata.BulkRNA_Profile)&
                              (metadata.AC_Treatment=='No'), :],
        confounders=None,
        timepoints_of_interested= ['Baseline', 'W3D1', 'W7D1'],
    )
    logger.info(
        '--------------------Calculate clinical benefits of gene expression--------------------'
    )
    analyzer.analyze_arm_response_difference(
        output_prefix=output_prefix+'_NoAC', pathway_analysis=True
    )