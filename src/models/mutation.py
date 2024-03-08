from typing import List
import click
import pandas as pd
from .utils import bulkAnalysis
from ..dataset import bulkDataset
from ..settings import init_logger

logger = init_logger(__name__)

def calculate_mutation_frequency(mutation_data:pd.DataFrame,metadata:pd.DataFrame,groupby: List[str])->pd.DataFrame:
    """
    Calulcate the mutation frequency based on the groupby values

    Parameters
    ----------
    mutation_data : pd.DataFrame
        mutation load with three columns: sample, category, value
    metadata : pd.DataFrame
        metadata
    groupby : List[str]
        column names in the metadata to group the data

    Returns
    -------
    pd.DataFrame
        mutation statistics 
    """
    groupby_list = groupby + [groupby]
    result = []
    for col in groupby_list:
        for gp,sub_df in metadata.groupby(col):
            if sub_df.shape[0] == 1:
                print(sub_df)
            sub_mut = mutation_data.loc[mutation_data['sample'].isin(sub_df.index),:]
            n_gene = sub_mut[['sample','category']].drop_duplicates()['category'].value_counts()
            res_df = (n_gene.rename('Freq_mutation') / sub_df.shape[0]).to_frame()
            res_df['N_sample'] = sub_df.shape[0]
            res_df['N_Mut'] = n_gene
            res_df.index.name='Gene'
            if isinstance(col,str):
                res_df[col] = gp
            else:
                for i,g in enumerate(col):
                    res_df[g] = gp[i]
            result.append(res_df.reset_index())
    result = pd.concat(result,axis=0)
    ## Fill all categories for NA values, because it has not been counted in the calculation
    for g in groupby:
        fill_value = ','.join(metadata[g].unique().tolist())
        result[g].fillna(fill_value,inplace=True)

    return result
    
        
    
    
@click.command()
@click.argument('data_processed_folder', type=click.Path(exists=True))
def main(data_processed_folder):
    output_prefix = f"{data_processed_folder}/mutation"
    dataset = bulkDataset(data_processed_folder)
    mutation_data = dataset.load('somatic_mutation')
    metadata = dataset.metadata
    metadata = metadata.loc[(metadata.WES_combined_TiN<=0.3)&
                            (metadata.WES_combined_fracContam<=0.04)&
                            (metadata.WES_Profile==True),:]
    mutation_statistics = calculate_mutation_frequency(mutation_data=mutation_data,
                                 metadata=metadata,
                                 groupby=['Treatment_Arm','Timepoint'])
    print(metadata.shape)
    mutation_statistics.to_csv(f"{output_prefix}_frequency_stats.csv",index=False)
    genes_of_interests = mutation_statistics.Gene[(mutation_statistics.Freq_mutation > .1)&
                                                  (~mutation_statistics.Timepoint.str.contains(','))].unique()
    
    analyzer = bulkAnalysis(
        feature_mtrx =(pd.crosstab(mutation_data['sample'], mutation_data['category']) > 0).astype(int),#.apply(lambda x: pd.Categorical(x,categories=[True,False]))[genes_of_interests],
        metadata=metadata,
        timepoints_of_interested = ['Baseline'],
        confounders=['WES_absolute_purity']
    )
    
    logger.info(
        '--------------------Calculate clinical benefits of genes with mutations--------------------'
    )
    analyzer.analyze_arm_response_difference(output_prefix=output_prefix)

    analyzer = bulkAnalysis(
        feature_mtrx = (pd.crosstab(mutation_data['sample'], mutation_data['category']) > 0).apply(lambda x: pd.Categorical(x,categories=[True,False]))[genes_of_interests],
        metadata=metadata,
        confounders=None,
    )
    
