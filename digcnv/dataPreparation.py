from digcnv.digCNV_logger import logger as dc_logger
import pandas as pd
import numpy as np
from os.path import exists, split, join

def addMicroArrayQualityData(cnvs, data_path:str):
    data = pd.read_csv(data_path, sep="\t")
    dc_logger.info("Micro-array quality data opened")
    dc_logger.info("Add {} columns".format(data.columns.tolist()))
    data.columns = ["SampleID", "LRR_mean", "LRR_median", "LRR_SD", "BAF_mean", "BAF_median", "BAF_SD", "BAF_DRIFT", "WF", "GCWF"]
    cnvs_qc = pd.merge(cnvs, data, how="left", on='SampleID')
    dc_logger.info("Micro-array quality data merged to CNV quality data")
    return cnvs_qc


def addDerivedFeatures(cnvs: pd.DataFrame) -> pd.DataFrame:
    cnvs["SIZE"] = cnvs["STOP"] - cnvs["START"] + 1 
    cnvs["DENSITY"] =  cnvs["SNP"] / cnvs["SIZE"]
    cnvs["Score_SNP"] = cnvs["SCORE"] / cnvs["SNP"]
    dc_logger.info("Derived features (DENSITY and Score_SNP columns) created")
    return cnvs

def addCallRateToDataset(cnvs: pd.DataFrame, call_rate_path: str, callrate_colname="callRate", individual_colname="FID") -> pd.DataFrame:
    if exists(call_rate_path):
        callrates = pd.read_csv(call_rate_path, sep='\t')
    else:
        raise Exception ("Given path doesn't exist")
    callrates.rename(columns = {callrate_colname:"CallRate"}, inplace=True)
    if len(callrates[individual_colname].unique()) != callrates.shape[0]:
        raise Exception("CallRates file must contains only unique individuals")
    if "CallRate" in cnvs.columns.tolist():
        dc_logger.info("Clean existing CallRate column")
        cnvs.drop(columns={"CallRate"}, inplace=True)

    cnvs_with_callrate = pd.merge(cnvs, callrates.loc[:,[individual_colname, "CallRate"]], how="left", left_on=individual_colname, right_on=individual_colname)

    if individual_colname != "SampleID":
        cnvs_with_callrate.drop(columns=[individual_colname], inplace = True)
    dc_logger.info("CallRate added to dataset")
    return cnvs_with_callrate

def addNbProbeByTech(cnvs: pd.DataFrame, nb_prob_tech=None, pfb_file_path=None) -> pd.DataFrame:
    if nb_prob_tech != None:
        cnvs["Nb_Probe_tech"] = nb_prob_tech
        dc_logger.info("Number of probes in technology added thanks to the given number")
    elif pfb_file_path != None:
        lines = 0
        with open(pfb_file_path) as f:
            for line in f:
                lines = lines + 1
        cnvs["Nb_Probe_tech"] = lines - 1
        dc_logger.info("Number of probes in technology added after counting number of lines in pfb file")
    else:
        dc_logger.info("You have to give at least one of these two parameters (nb_prob_tech or pfb_file_path)")

    dc_logger.info("Number of probes in technology added")
    return cnvs

def addChromosomicAnnotation(cnvs: pd.DataFrame, centromere_list_path=None, segdup_list_path=None) -> pd.DataFrame:
    this_dir, this_filename = split(__file__)
    if centromere_list_path == None:
        centromere_list_path = join(this_dir, "data", "Region_centromere_hg19.dat")
    if segdup_list_path == None:
        segdup_list_path = join(this_dir, "data", "SegDup_filtres_Ok_Oct.map")

    cnvs["overlapCNV_Centromere"] = getCentromereOverlap(cnvs, centromere_list_path)
    dc_logger.info("Centromere overlap added to CNVs")
    cnvs = getSegDupOverlap(cnvs, segdup_list_path)
    dc_logger.info("Both chromosomic annotation finished")
    return cnvs

def getCentromereOverlap(cnvs: pd.DataFrame, centromeres_list_path:str) -> list:
    if exists(centromeres_list_path):
        centromeres = pd.read_csv(centromeres_list_path, sep='\t')
    else:
        raise Exception ("File {} note found".format(centromeres_list_path))

    #validate the necessary titles of the regions of interest file
    nescessaryRegionTitles = ["CHR", "START", "STOP"]
    if len(centromeres.columns.intersection(nescessaryRegionTitles)) < 3:
        raise Exception("The input file for the regions of interest must contain the following columns: CHR, START and STOP")
    
    centromeres.rename(columns={"CHR":"CHR_centro", "START":"START_centro", "STOP":"STOP_centro"}, inplace=True)
    temp = pd.merge(cnvs, centromeres, left_on="CHR", right_on="CHR_centro")
    overlap = (temp[["STOP_centro", "STOP"]].min(axis=1) - temp[["START_centro", "START"]].max(axis=1) + 1)/(temp.STOP - temp.START + 1)
    overlap = np.where(overlap > 0, overlap, 0)
    dc_logger.info("Centromere overlap created")
    return overlap

def getSegDupOverlap(cnvs:pd.DataFrame, segdup_list_path:str):
    if exists(segdup_list_path):
        segdups = pd.read_csv(segdup_list_path, sep='\t')
        dc_logger.info("Segmental duplication file opened")
    else:
        raise Exception ("File note found")
    complete_cnvs = pd.DataFrame()
    for chr in cnvs.CHR.unique():
        cnvs_chr = cnvs[cnvs.CHR == chr].copy()
        segdups_chr = segdups[segdups.CHR == chr].copy()
        segdups_chr.reset_index(inplace=True)
        vals = segdups_chr.apply(lambda x: computeOneOverlap(cnvs_chr, x.START, x.STOP), axis=1)
        overlaps = pd.DataFrame(vals.tolist())
        cnvs_chr.loc[:,"overlapCNV_SegDup"] = overlaps.sum().tolist()
        complete_cnvs = pd.concat([complete_cnvs, cnvs_chr])
    dc_logger.info("Segmental duplication overlap computed and CNVs annotated")
    return complete_cnvs

def computeOneOverlap(cnvs, START, STOP):
    overlap = (np.where(cnvs['STOP']< STOP, cnvs['STOP'], STOP) - np.where(cnvs['START'] < START, START, cnvs['START']) + 1)/(cnvs['STOP'] - cnvs['START'] + 1)
    overlap = np.where(overlap > 0, overlap, 0)
    return overlap

def transformTwoAlgsFeatures(cnvs: pd.DataFrame) -> pd.DataFrame:
    if cnvs.TwoAlgs.dtype == object:
        cnvs.TwoAlgs = cnvs.TwoAlgs.str[:-1]
        cnvs.TwoAlgs = cnvs.TwoAlgs.astype(int)
    if cnvs.TwoAlgs.describe()[7] <= 1.0:
        cnvs.TwoAlgs = cnvs.TwoAlgs *100
        dc_logger.info("Transform TwoAlgs function into percentage format")
    elif cnvs.TwoAlgs.describe()[7] > 1.0 <= 100.0:
        dc_logger.info("Keep TwoAlgs function into percentage format")
    else: 
        dc_logger.info("Error in TwoAlgs format must be a percentage or a rate")
        quit()
    return cnvs

def main():
    cnvs = pd.read_csv('../data/UKBB_clean_for_DigCNV.tsv', sep='\t', index_col=False)
    cnvs = getSegDupOverlap(cnvs, '../data/SegDup_filtres_Ok_Oct.map')
    print('ok')

if __name__ == '__main__':
    main()