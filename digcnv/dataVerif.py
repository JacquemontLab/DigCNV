from typing import Optional
from digcnv.digCNV_logger import logger as dc_logger
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def computeNaPercentage(cnvs: pd.DataFrame, dimensions: list, remove_na_data=True) -> Optional[tuple[pd.DataFrame, pd.DataFrame]]:  
    """Will Check for each column given as input if their is any NA value and remove them if the option is set.
    Machine learning can't work with missing values so they must be removed before using the model.

    :param cnvs: list of CNVs with their scores
    :type cnvs: pd.DataFrame
    :param dimensions: list of specific columns used next for DigCNV model. 
    :type dimensions: list
    :param remove_na_data: Boolean indicating if the function removes NA value in addition of showing NA percentage by column, defaults to True
    :type remove_na_data: bool, optional
    :return: If remove na option is `True` will return a tuple of Two Dataframes, a first one containing CNVs with all 
    :rtype: Optional[tuple[pd.DataFrame, pd.DataFrame]]
    """    
    removed_cnvs = pd.DataFrame()
    if remove_na_data:
        cnvs_clean = cnvs.copy()
    split_cnvs = cnvs.loc[:, dimensions]
    dc_logger.info("Check if no data is missing:")
    i = 0
    while i < len(split_cnvs.columns):
        nb_line = split_cnvs.shape[0]-1
        count = split_cnvs.iloc[:, i].count()-1
        if count/nb_line == 1:
            dc_logger.info("{}: 100%".format(split_cnvs.columns[i]))
        else:
            dc_logger.info("{}: {:.1f}%".format(
                split_cnvs.columns[i], count/nb_line*100))
            if remove_na_data:
                removed_index = split_cnvs[split_cnvs[split_cnvs.columns[i]].isnull(
                )].index.tolist()
                dc_logger.info(
                    "{} CNV removed due to Null Data".format(len(removed_index)))
                split_cnvs = split_cnvs.drop(index=removed_index)
                removed_cnvs = pd.concat(
                    [removed_cnvs, cnvs_clean.iloc[removed_index]])
                cnvs_clean = cnvs_clean.drop(index=removed_index)
                cnvs_clean.reset_index(drop=True, inplace=True)
                split_cnvs.reset_index(drop=True, inplace=True)
        i += 1
    if remove_na_data:
        return cnvs_clean, removed_cnvs


def plotCorrelationHeatMap(cnvs:pd.DataFrame, list_dim:list, output_path=None, plot_fig = True):
    """Plot the correlation heatmap for the given list of features of the CNV list

    :param cnvs: list of CNVs with their scores
    :type cnvs: pd.DataFrame
    :param list_dim: list of specific columns used next for DigCNV model. 
    :type list_dim: list
    :param output_path: Pathway of the image could be a `PNG` or a `PDF` format, defaults to None
    :type output_path: str, optional
    :param plot_fig: Option to disable plotting figure in interactive session
    :type plot_fig: bool, optional
    """    
    data_clean = cnvs.loc[:, list_dim]
    cor_data = data_clean.corr()
    fig, ax = plt.subplots()
    # Set the width and height
    fig.set_figwidth(20)
    fig.set_figheight(20)
    sns.heatmap(cor_data.abs(), linewidth=0.5, vmin=0.0,
                vmax=1.0, cmap="OrRd", annot=round(cor_data, 2))

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")
    ax.set_title("Correlation table of CNV features")
    if output_path != None:
        plt.savefig(output_path)
    
    if(plot_fig):
        plt.show()



def main():
    cnvs = pd.read_csv('../data/UKBB_clean_for_DigCNV.tsv',
                       sep='\t', index_col=False)
    print('ok')


if __name__ == '__main__':
    main()
