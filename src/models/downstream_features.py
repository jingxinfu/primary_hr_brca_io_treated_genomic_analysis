from typing import Union
import click
import numpy as np
import pandas as pd
from .utils import bulkAnalysis, fdrAdjust, wilcoxonTest, contigencyTableTest
from ..dataset import bulkDataset
from ..settings import init_logger

logger = init_logger(__name__)


class featureAnalysis(bulkAnalysis):
    def __init__(
        self,
        feature_mtrx: pd.DataFrame,
        metadata: pd.DataFrame,
        *args,
        **kwargs,
    ):
        super().__init__(
            feature_mtrx=feature_mtrx,
            metadata=metadata,
            min_feature_variance=None,
            significance_column='adj',
            statistics_column='Zscore',
            statistics_threshold=0,
            *args,
            **kwargs,
        )

    def test_function(
        self,
        X: pd.DataFrame,
        y: pd.Series,
        paired:Union[pd.Series,None]=None,
        covars: Union[pd.Series, pd.DataFrame, None] = None,
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
        result = X.apply(
            lambda c: wilcoxonTest(c, y)
            if pd.api.types.is_numeric_dtype(c)
            else contigencyTableTest(x=c, y=y),
            axis=0,
        ).T
        result = result.apply(pd.Series)
        result.columns = ['Zscore', 'Pvalue', 'Method']
        result['padj'] = np.nan
        result.loc[~result.Pvalue.isna(), 'padj'] = fdrAdjust(
            p_value=result['Pvalue'][~result.Pvalue.isna()]
        )

        return result.reset_index()


@click.command()
@click.argument('data_processed_folder', type=click.Path(exists=True))
def main(data_processed_folder):
    dataset = bulkDataset(data_processed_folder)
    metadata = dataset.metadata
    selected_features = metadata.columns[metadata.columns.str.contains('Cytokine') | 
                                         metadata.columns.str.contains('CIBERSORT-ABS') | 
                                         metadata.columns.str.contains('ssGSEA')
                                         ].tolist()
    feature_mtrx = metadata.loc[
        metadata['WES_Profile'] & metadata['BulkRNA_Profile'],
        metadata.columns.str.startswith('WES')
        | metadata.columns.str.startswith('bulkRNA'),
    ]
    feature_mtrx = feature_mtrx.loc[
        :, feature_mtrx.apply(lambda v: v.nunique() > 1, axis=0)
    ]
    feature_mtrx = feature_mtrx.apply(
        lambda v: v.map({True: '1', False: '0'}) if v.dtypes.name == 'bool' else v,
        axis=0,
    )
    to_dummies_cols = [
        c
        for c in feature_mtrx.columns
        if feature_mtrx[c].dtype == "O" and feature_mtrx[c].nunique() > 2
    ]
    feature_mtrx = pd.get_dummies(feature_mtrx, columns=to_dummies_cols, dtype=str)
    feature_mtrx = feature_mtrx.loc[:, ~feature_mtrx.columns.str.contains('N/A')]
    feature_mtrx.columns.name = 'Feature'
    feature_mtrx = feature_mtrx[selected_features]

    analyzer = featureAnalysis(
        feature_mtrx=feature_mtrx, metadata=metadata, confounders=None
    )
    output_prefix = f"{data_processed_folder}/DownstreamFeatures"
    logger.info(
        '--------------------Calculate clinical benefits of features --------------------'
    )
    analyzer.analyze_arm_response_difference(output_prefix=output_prefix)
    logger.info(
        '--------------------Detect enriched highly expressed features after treatment exposure--------------------'
    )
    analyzer.analyze_treatment_impact(output_prefix=output_prefix)
