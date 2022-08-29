# from DigCNVlib import digCNV_logger
import dataVerif, dataPreparation, DigCnvPreProcessing
import pandas as pd


def main():
    # python setup.py bdist_wheel
    cnvs = pd.read_csv('../data/UKBB_clean_for_DigCNV.tsv', sep='\t', index_col=False)
    cnvs.rename(columns={'LociStart':'START',
                        'LociStop':'STOP',
                        'Chr':'CHR'}, inplace=True)
    print(cnvs.columns.tolist())
    cnvs_clean = cnvs[['SampleID','START', 'STOP', 'CHR', 'SCORE', 'SNP', 'TwoAlgs', 'overlapCNV_SegDup', 'overlapCNV_Centromere', 'LRR_mean', 'LRR_SD', 'BAF_mean', 'WF', 'GCWF', 'SnipPeep_Ok']]
    cnvs_clean.to_csv('../data/example_cnvs.tsv', sep='\t', index=False)
    callrates = cnvs[['SampleID', 'CallRate']]
    callrates.to_csv('../data/callrates.tsv', sep='\t', index=False)
    
    cnvs = pd.read_csv('../data/example_cnvs.tsv', sep='\t')
    print(cnvs.columns.tolist())

    # List of mandatories columns for the model.
    model_dimensions = ['SCORE', 'Score_SNP']

    dataVerif.checkIfMandatoryColumnsExist(cnvs, post_data_preparation=False)
    print('---------')
    cnvs = dataPreparation.addDerivedFeatures(cnvs)
    print('---------')
    cnvs = dataPreparation.addCallRateToDataset(cnvs, call_rate_path='../data/callrates.tsv', callrate_colname='CallRate', individual_colname='SampleID')
    print('---------')
    cnvs = dataPreparation.addNbProbeByTech(cnvs, nb_prob_tech=7777)
    cnvs = dataPreparation.addNbProbeByTech(cnvs, pfb_file_path='../data/UKBB_PFB.pfb')
    print('---------')
    cnvs = dataPreparation.transformTwoAlgsFeatures(cnvs)
    print('---------')
    cnvs = dataPreparation.addChromosomicAnnotation(cnvs)
    print('---------')
    dataVerif.checkIfMandatoryColumnsExist(cnvs, post_data_preparation=True)
    print('---------')
    dataVerif.checkColumnsformats(cnvs)
    print('---------')
    dataVerif.computeNaPercentage(cnvs, dimensions=model_dimensions)
    print('---------')
    dataVerif.plotCorrelationHeatMap(cnvs, list_dim=model_dimensions, output_path='../outputs/correlation.png')
    print('---------')

    cnvs, removed = DigCnvPreProcessing.removeLinesWithNA(cnvs, dimensions=model_dimensions)
    print('---------')
    # X_train, y_train, X_test, y_test = DigCnvPreProcessing.createTrainingTestingDatasets(cnvs, X_dimension='SnipPeep_Ok')
    print('---------')
    # DigCnvPreProcessing.uniformizeClassesSizes()

    print('ok')

if __name__ == '__main__':
    main()