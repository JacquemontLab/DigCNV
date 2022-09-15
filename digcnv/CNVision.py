from digcnv.digCNV_logger import logger as dc_logger

import tempfile
from os.path import split, join, exists
from subprocess import Popen, PIPE
import pandas as pd


def formatPennCNVforCNVision(pennCNVfile_path: str, output_path: str):
    """Format PennCNV file to be readable by CNVision algorithm and write formatted CNVs into the given output path

    :param pennCNVfile_path: Pathway of the PennCNV output file listing all CNVs with their score. For File format example read README.md page on DigCNV github page.
    :type pennCNVfile_path: str
    :param output_path: Pathway of the PennCNV file formatted to use in CNVision function.
    :type output_path: str

    """
    dc_logger.info("format and filter PennCNV file to CNVision requirements")
    this_dir = split(__file__)[0]
    perl_script = join(this_dir, 'data', "CNVision_format_merge_ukbb.pl")
    proc1 = Popen(["perl", perl_script, "--PNformat", pennCNVfile_path,
                  output_path, 'tmp'], stdout=PIPE,  bufsize=1, text=True)
    while proc1.poll() is None:
        line = proc1.stdout.readline()
        dc_logger.info(line.rstrip())
    proc1.wait()
    dc_logger.info("PennCNV file formatted for CNVision")


def formatQuantiSNPforCNVision(quantiSNP_file_path: str, output_path: str):
    """Format QuantiSNP file to be readable by CNVision algorithm and write formatted CNVs into the given output path

    :param quantiSNP_file_path: Pathway of the QuantiSNP output file listing all CNVs with their score. For File format example read README.md page on DigCNV github page.
    :type quantiSNP_file_path: str
    :param output_path: Pathway of the QuantiSNP file formatted to use in CNVision function.
    :type output_path: str
    """
    dc_logger.info("format and filter QuantiSNP file to CNVision requirements")
    this_dir = split(__file__)[0]
    perl_script = join(this_dir, 'data', "CNVision_format_merge_ukbb.pl")
    proc1 = Popen(["perl", perl_script, "--QTformat", quantiSNP_file_path,
                  output_path, 'tmp'], stdout=PIPE, bufsize=1, text=True)
    while proc1.poll() is None:
        line = proc1.stdout.readline()
        dc_logger.info(line.rstrip())
    proc1.wait()
    dc_logger.info("PennCNV file formatted for CNVision")


def runMergingScript(formated_data_paths: list, output_path: str):
    """Create a merged CNV list file at the given output path. CNV files inputed must be in CNVision ready format.

    :param formated_data_paths: list of pathways to called CNVs and formatted for CNVision.
    :type formated_data_paths: list
    :param output_path: Pathway where the merged list of CNV will be written with the new algorithm overlap information added.
    :type output_path: str
    """
    dc_logger.info("format and filter QuantiSNP file to CNVision requirements")
    this_dir = split(__file__)[0]
    perl_script = join(this_dir, 'data', "CNVision_format_merge_ukbb.pl")
    proc1 = Popen(["perl", perl_script, "--merge", formated_data_paths[0],
                  formated_data_paths[1], output_path, 'tmp'], stdout=PIPE, bufsize=1, text=True)
    while proc1.poll() is None:
        line = proc1.stdout.readline()
        dc_logger.info(line.rstrip())
    proc1.wait()
    dc_logger.info("PennCNV file formatted for CNVision")


def mergeMultipleCNVCallingOutputs(list_calling_outputs_path: list, list_calling_softwares: list) -> pd.DataFrame:
    """Create a dataframe containing merged CNVs comming from the given algorithm outputs. Presently working only with PennCNV and QuantiSNP only.
    Please indicate the list of pathway and the list of algorithm names in the same order.

    :param list_calling_outputs_path: list of calling output pathways
    :type list_calling_outputs_path: list_calling_outputs_path: list
    :param list_calling_softwares: {`PennCNV`, `QuantiSNP`} list of calling algorithm names, indicated in the same order as the list of pathways 
    :type list_calling_softwares: list_calling_softwares: list
    :returns: Dataframe containing all CNVs coming from all caller outputs merged and with their overlap annotated
    :rtype: pd.DataFrame
    :raises Exception: If two given list haven't same sizes
    :raises Exception: If the calling algorithm name isn't supported

    """
    if len(list_calling_outputs_path) != len(list_calling_softwares):
        raise Exception("Both list must have same sizes")
    list_tmp_paths = []
    with tempfile.TemporaryDirectory() as tmp_dir:
        for i, calling_soft in enumerate(list_calling_softwares):
            if calling_soft == 'PennCNV':
                formatPennCNVforCNVision(list_calling_outputs_path[i], tmp_dir)
                if exists(join(tmp_dir, 'PC_tmp_CNVisionFormated.txt')):
                    list_tmp_paths.append(
                        join(tmp_dir, 'PC_tmp_CNVisionFormated.txt'))
            elif calling_soft == 'QuantiSNP':
                formatQuantiSNPforCNVision(
                    list_calling_outputs_path[i], tmp_dir)
                if exists(join(tmp_dir, 'QS_tmp_CNVisionFormated.txt')):
                    list_tmp_paths.append(
                        join(tmp_dir, 'QS_tmp_CNVisionFormated.txt'))
            else:
                raise Exception(
                    "Given calling software: {} in't support by the software please contact us if you want to add this software".format(calling_soft))
        runMergingScript(list_tmp_paths, tmp_dir)
        CNVs_list = pd.read_csv(
            join(tmp_dir, 'Sum_CNVisionMerged_PC_QS_tmp.txt'), sep='\t')

    return CNVs_list


def main():
    print(mergeMultipleCNVCallingOutputs(['/home/thomas/Documents/temp_data/PC_allCNV.txt',
          '/home/thomas/Documents/temp_data/QS_allCNV.txt'], ['PennCNV', 'QuantiSNP']).head())


if __name__ == "__main__":
    main()
