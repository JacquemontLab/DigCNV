from digcnv.digCNV_logger import logger as dc_logger
import pandas as pd
import numpy as np
from os.path import exists, split, join


def addMicroArrayQualityData(cnvs: pd.DataFrame, data_path: str) -> pd.DataFrame:
    """Add Microarray quality data to the list of merged CNVs

    :param cnvs: list of CNVs and their quality scores
    :type cnvs: pd.DataFrame
    :param data_path: pathway containing PennCNV quality chip output
    :type data_path: str
    :return: list of CNVs with sample quality aggregated
    :rtype: pd.DataFrame
    """
    data = pd.read_csv(data_path, sep="\t")
    dc_logger.info("Micro-array quality data opened")
    dc_logger.info("Add {} columns".format(data.columns.tolist()))
    data.columns = ["SampleID", "LRR_mean", "LRR_median", "LRR_SD",
                    "BAF_mean", "BAF_median", "BAF_SD", "BAF_DRIFT", "WF", "GCWF"]
    cnvs_qc = pd.merge(cnvs, data, how="left", on='SampleID')
    dc_logger.info("Micro-array quality data merged to CNV quality data")
    return cnvs_qc


def addDerivedFeatures(cnvs: pd.DataFrame) -> pd.DataFrame:
    """Compute and add 3 derived CNV scores 

    :param cnvs: list of CNVs with at least 4 mandatory columns, (`START`, `STOP`, `SNP`, `SCORE`) 
    :type cnvs: pd.DataFrame
    :return: same list of CNVs with 3 new column scores, (`SIZE`, `DENSITY`, `Score_SNP`)
    :rtype: pd.DataFrame
    """
    cnvs["SIZE"] = cnvs["STOP"] - cnvs["START"] + 1
    cnvs["DENSITY"] = cnvs["SNP"] / cnvs["SIZE"]
    cnvs["Score_SNP"] = cnvs["SCORE"] / cnvs["SNP"]
    dc_logger.info("Derived features (DENSITY and Score_SNP columns) created")
    return cnvs


def addCallRateToDataset(cnvs: pd.DataFrame, call_rate_path: str, callrate_colname="callRate", individual_colname="SampleID") -> pd.DataFrame:
    """Add CallRate information to the given list of CNVs

    :param cnvs: list of CNVs with their scores
    :type cnvs: pd.DataFrame
    :param call_rate_path: Pathway to the file listing CallRate data of each microarray. The file mut be a **tsv**.
    :type call_rate_path: str
    :param callrate_colname: CallRate column name, defaults to "callRate"
    :type callrate_colname: str, optional
    :param individual_colname: Individual column name, defaults to "FID"
    :type individual_colname: str, optional
    :raises Exception: If the given CallRate pathway doesn't exist
    :raises Exception: If the given CallRate file contains more than one value by sampleID
    :return: list of CNVs with CallRates aggregated
    :rtype: pd.DataFrame
    """
    if exists(call_rate_path):
        callrates = pd.read_csv(call_rate_path, sep='\t')
    else:
        raise Exception("Given path doesn't exist")
    callrates.rename(columns={callrate_colname: "CallRate"}, inplace=True)
    if len(callrates[individual_colname].unique()) != callrates.shape[0]:
        raise Exception("CallRates file must contains only unique individuals")
    if "CallRate" in cnvs.columns.tolist():
        dc_logger.info("Clean existing CallRate column")
        cnvs.drop(columns={"CallRate"}, inplace=True)

    cnvs_with_callrate = pd.merge(cnvs, callrates.loc[:, [
                                  individual_colname, "CallRate"]], how="left", left_on=individual_colname, right_on=individual_colname)

    if individual_colname != "SampleID":
        cnvs_with_callrate.drop(columns=[individual_colname], inplace=True)
    dc_logger.info("CallRate added to dataset")
    return cnvs_with_callrate


def addNbProbeByTech(cnvs: pd.DataFrame, nb_prob_tech=None, pfb_file_path=None) -> pd.DataFrame:
    """Add the number of Probes used to genotyped individuals. **Warning**: if your data comes from multiple genotyping technologies
    Please split your dataframe to get only one technology by dataframe before using this function.

    :param cnvs: list of CNVs with their scores
    :type cnvs: pd.DataFrame
    :param nb_prob_tech: Number of probe in the genotyping technology used, defaults to None
    :type nb_prob_tech: int, optional
    :param pfb_file_path: Pathway to the PFB file to compute get the number of probe used by the technology, defaults to None
    :type pfb_file_path: str, optional
    raises Exception: If none of the optional value is given. You have to set either nb_probe_tech or pfb_file_path parameter.
    :return: list of CNVs with Nb probe used in technology column aggregated
    :rtype: pd.DataFrame
    """
    if nb_prob_tech != None:
        cnvs["Nb_Probe_tech"] = nb_prob_tech
        dc_logger.info(
            "Number of probes in technology added thanks to the given number")
    elif pfb_file_path != None:
        lines = 0
        with open(pfb_file_path) as f:
            for line in f:
                lines = lines + 1
        cnvs["Nb_Probe_tech"] = lines - 1
        dc_logger.info(
            "Number of probes in technology added after counting number of lines in pfb file")
    else:
        raise Exception(
            "You have to give at least one of these two parameters (nb_prob_tech or pfb_file_path)")

    dc_logger.info("Number of probes in technology added")
    return cnvs


def addChromosomicAnnotation(cnvs: pd.DataFrame, centromere_list_path=None, segdup_list_path=None) -> pd.DataFrame:
    """Compute percentage of overlap for each CNV with specific chromosomic regions: centromere and Segmental Duplication. By default the overlap is map on Hg19 Human genome.

    :param cnvs: list of CNVs with their scores
    :type cnvs: pd.DataFrame
    :param centromere_list_path: Pathway to a tsv file containing a list of centromere regions coordinates. By default will use the list human centromere mapped on Hg19 genome, defaults to None
    :type centromere_list_path: str, optional
    :param segdup_list_path: Pathway to a tsv file containing a list of Segmental duplication regions coordinates. By default will use the list human segmental duplication mapped on Hg19 genome, defaults to None
    :type segdup_list_path: str, optional
    :return: list of CNVs with chromosomic features columns aggregated
    :rtype: pd.DataFrame
    """
    this_dir, this_filename = split(__file__)
    if centromere_list_path == None:
        centromere_list_path = join(
            this_dir, "data", "Region_centromere_hg19.dat")
    if segdup_list_path == None:
        segdup_list_path = join(this_dir, "data", "SegDup_filtres_Ok_Oct.map")

    cnvs = getCentromereOverlap(cnvs, centromere_list_path)
    dc_logger.info("Centromere overlap added to CNVs")
    cnvs = getSegDupOverlap(cnvs, segdup_list_path)
    dc_logger.info("Both chromosomic annotation finished")
    return cnvs


def getCentromereOverlap(cnvs: pd.DataFrame, centromeres_list_path: str) -> pd.DataFrame:
    """Compute percentage of overlap for each CNV with centromere.

    :param cnvs: list of CNVs with their scores
    :type cnvs: pd.DataFrame
    :param centromeres_list_path: Pathway to a tsv file containing a list of centromere regions coordinates
    :type centromeres_list_path: str
    :raises Exception: If the given centromere list file pathway doesn't exist
    :raises Exception: If the given file hasn't mandatory columns (`CHR`, `START`, `STOP`)
    :return: list of CNVs with centromere overlap feature aggregated
    :rtype: pd.DataFrame
    """
    if exists(centromeres_list_path):
        centromeres = pd.read_csv(centromeres_list_path, sep='\t')
    else:
        raise Exception("File {} note found".format(centromeres_list_path))

    # validate the necessary titles of the regions of interest file
    nescessaryRegionTitles = ["CHR", "START", "STOP"]
    if len(centromeres.columns.intersection(nescessaryRegionTitles)) < 3:
        raise Exception(
            "The input file for the regions of interest must contain the following columns: CHR, START and STOP")

    centromeres.rename(columns={
                       "CHR": "CHR_centro", "START": "START_centro", "STOP": "STOP_centro"}, inplace=True)
    temp = pd.merge(cnvs, centromeres, left_on="CHR", right_on="CHR_centro", how="left")
    overlap = (temp[["STOP_centro", "STOP"]].min(
        axis=1) - temp[["START_centro", "START"]].max(axis=1) + 1)/(temp.STOP - temp.START + 1)
    overlap = np.where(overlap > 0, overlap, 0)
    dc_logger.info("Centromere overlap created")
    cnvs["overlapCNV_Centromere"] = overlap
    return cnvs


def getSegDupOverlap(cnvs: pd.DataFrame, segdup_list_path: str) -> pd.DataFrame:
    """Compute percentage of overlap for each CNV with Segmental Duplication regions.

    :param cnvs: list of CNVs with their scores
    :type cnvs: pd.DataFrame
    :param segdup_list_path: Pathway to a tsv file containing a list of centromere regions coordinates
    :type segdup_list_path: str
    :raises Exception:  If the given sgmental duplication list file pathway doesn't exist
    :return: list of CNVs with sgmental duplication overlap feature aggregated
    :rtype: pd.DataFrame
    """
    if exists(segdup_list_path):
        segdups = pd.read_csv(segdup_list_path, sep='\t')
        dc_logger.info("Segmental duplication file opened")
    else:
        raise Exception("File note found")
    complete_cnvs = pd.DataFrame()

    for chr in cnvs.CHR.unique():
        if chr in segdups.CHR.unique():
            cnvs_chr = cnvs[cnvs.CHR == chr].copy()
            segdups_chr = segdups[segdups.CHR == chr].copy()
            segdups_chr.reset_index(inplace=True)
            vals = segdups_chr.apply(lambda x: computeOneOverlap(
                cnvs_chr, x.START, x.STOP), axis=1)
            overlaps = pd.DataFrame(vals.tolist())
            cnvs_chr.loc[:, "overlapCNV_SegDup"] = overlaps.sum().tolist()
            complete_cnvs = pd.concat([complete_cnvs, cnvs_chr])
    dc_logger.info("Segmental duplication overlap computed and CNVs annotated")
    return complete_cnvs


def computeOneOverlap(cnvs:pd.DataFrame, START:int, STOP:int):
    """Compute the percentage of overlap between a list of CNVs and START -- STOP coordinates of a given Segmental Duplication region. 
    Warning: CNVs must be on the same chromsome than the tested segmental duplication

    :param cnvs: list of CNVs coming from a unique chromosome similar to the segmental duplication analyzed 
    :type cnvs: pd.DataFrame
    :param START: Start coordinate of the segmental duplication region
    :type START: int
    :param STOP: Stop coordinate of the segmental duplication region
    :type STOP: int
    :return: An array listing the percentage of overlap for all given CNVs with the segmental duplication region
    :rtype: np.ArrayLike
    """
    overlap = (np.where(cnvs['STOP'] < STOP, cnvs['STOP'], STOP) - np.where(
        cnvs['START'] < START, START, cnvs['START']) + 1)/(cnvs['STOP'] - cnvs['START'] + 1)
    overlap = np.where(overlap > 0, overlap, 0)
    return overlap


def transformTwoAlgsFeatures(cnvs: pd.DataFrame) -> pd.DataFrame:
    """CHeck if the TwoAlgs column has the right format and correct it if not.

    :param cnvs: list of CNVs with their scores
    :type cnvs: pd.DataFrame
    :return: list of CNVs with the TwoAlgs column corrected if necessary
    :rtype: pd.DataFrame
    """
    if cnvs.TwoAlgs.dtype == object:
        cnvs.TwoAlgs = cnvs.TwoAlgs.str[:-1]
        cnvs.TwoAlgs = cnvs.TwoAlgs.astype(int)
    if cnvs.TwoAlgs.describe()[7] <= 1.0:
        cnvs.TwoAlgs = cnvs.TwoAlgs * 100
        dc_logger.info("Transform TwoAlgs function into percentage format")
    elif cnvs.TwoAlgs.describe()[7] > 1.0 <= 100.0:
        dc_logger.info("Keep TwoAlgs function into percentage format")
    else:
        dc_logger.info(
            "Error in TwoAlgs format must be a percentage or a rate")
        quit()
    return cnvs


def main():
    """ """
    cnvs = pd.read_csv('../data/UKBB_clean_for_DigCNV.tsv',
                       sep='\t', index_col=False)
    cnvs = getSegDupOverlap(cnvs, '../data/SegDup_filtres_Ok_Oct.map')
    print('ok')


if __name__ == '__main__':
    main()
