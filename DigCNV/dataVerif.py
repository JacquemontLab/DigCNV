from DigCNV import digCNV_logger
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def computeNaPercentage(cnvs, dimensions):
    split_cnvs = cnvs.loc[:,dimensions]
    digCNV_logger.logger.info("Check if no data is missing:")
    i = 0
    nb_line = split_cnvs.shape[0]-1
    while i < len(split_cnvs.columns):
        count = split_cnvs.iloc[:, i].count()-1
        if count/nb_line == 1:
            digCNV_logger.logger.info("{}: 100%".format(split_cnvs.columns[i]))
        else:
            digCNV_logger.logger.info("{}: {:.1f}%".format(split_cnvs.columns[i], count/nb_line*100))
        i+=1

def plotCorrelationHeatMap(cnvs, list_dim, output_path=None):
    data_clean = cnvs.loc[:, list_dim]
    cor_data = data_clean.corr()
    fig, ax =  plt.subplots()
    # Set the width and height
    fig.set_figwidth(20)
    fig.set_figheight(20)
    sns.heatmap(cor_data.abs(), linewidth=0.5, vmin=0.0, vmax=1.0, cmap="OrRd", annot=round(cor_data,2))

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")
    ax.set_title("Correlation table of CNV features")
    if output_path != None:
        plt.savefig(output_path)
    plt.show()

def checkIfMandatoryColumnsExist(cnvs: pd.DataFrame, post_data_preparation=True):
    if post_data_preparation:
        mandatory_columns = ['START', 'STOP', 'CHR', 'SNP', 'SCORE', 'WF', 'TwoAlgs']
    else:
        mandatory_columns = ['START', 'STOP', 'CHR', 'SNP', 'SCORE', 'WF', 'TwoAlgs']
    if len(set(cnvs.columns.tolist()) & set(mandatory_columns)) != len(mandatory_columns):
        missing_col = list(set(mandatory_columns).difference(cnvs.columns.tolist()))
        raise Exception("\nSome columns are mandatory: {}\n{} are missing".format(mandatory_columns, missing_col))
    else:
        digCNV_logger.logger.info("All mandatory columns exist in the given dataframe")

def checkColumnsformats(cnvs:pd.DataFrame, post_data_preparation=True):
    # TODO
    if post_data_preparation == False:
        if cnvs[['START', 'STOP']].dtypes().tolist().unique() != int:
            print("pas OK")
    
    
def main():
    cnvs = pd.read_csv('../data/UKBB_clean_for_DigCNV.tsv', sep='\t', index_col=False)
    checkIfMandatoryColumnsExist(cnvs, False)
    print('ok')

if __name__ == '__main__':
    main()
    #     #validate format of CHR column
    # chr23 = "chr23" %in% cnv$CHR
    # if (chr23){
    # indexCHR23 = cnv$CHR %in% "chr23"
    # cnv$CHR[indexCHR23]= "chrX"
    # warnings("We converted the expression chr23 by chrX in the CHR column")
    # }