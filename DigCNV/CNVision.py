# from DigCNV import digCNV_logger

import tempfile
from os.path import split, join, exists
import subprocess
import pandas as pd

def formatPennCNVforCNVision(pennCNVfile_path:str, output_path:str):
    # digCNV_logger.logger.info("format and filter PennCNV file to CNVision requirements")
    this_dir = split(__file__)[0]
    perl_script = join(this_dir, 'scripts', "CNVision_format_merge_ukbb.pl")
    proc1=subprocess.Popen(["perl", perl_script, "--PNformat", pennCNVfile_path, output_path, 'tmp'])
    proc1.wait()
    # digCNV_logger.logger.info("PennCNV file formatted for CNVision")

    return 'ok'

def formatQuantiSNPforCNVision(quantiSNP_file_path:str, output_path:str):
    # digCNV_logger.logger.info("format and filter QuantiSNP file to CNVision requirements")
    this_dir = split(__file__)[0]
    perl_script = join(this_dir, 'scripts', "CNVision_format_merge_ukbb.pl")
    proc1=subprocess.Popen(["perl", perl_script, "--QTformat", quantiSNP_file_path, output_path, 'tmp'])
    proc1.wait()
    # digCNV_logger.logger.info("PennCNV file formatted for CNVision")

    return 'ok'

def runMergingScript(formated_data_paths:str, output_path:str):
    # digCNV_logger.logger.info("format and filter QuantiSNP file to CNVision requirements")
    this_dir = split(__file__)[0]
    perl_script = join(this_dir, 'scripts', "CNVision_format_merge_ukbb.pl")
    print(' '.join(formated_data_paths))
    proc1=subprocess.Popen(["perl", perl_script, "--merge", formated_data_paths[0], formated_data_paths[1], output_path, 'tmp'])
    proc1.wait()
    # digCNV_logger.logger.info("PennCNV file formatted for CNVision")

    return 'ok'

def mergeMultipleCNVCallingOutputs(list_calling_outputs_path=[], list_calling_softwares=[])-> pd.DataFrame:
    if len(list_calling_outputs_path) != len(list_calling_softwares):
        raise Exception("Both list must have same sizes")
    list_tmp_paths = []
    with tempfile.TemporaryDirectory() as tmp_dir:
        print(tmp_dir)
        for i, calling_soft in enumerate(list_calling_softwares):
            print(i)
            if calling_soft == 'PennCNV':
                formatPennCNVforCNVision(list_calling_outputs_path[i], tmp_dir)
                if exists(join(tmp_dir, 'PC_tmp_CNVisionFormated.txt')):
                    list_tmp_paths.append(join(tmp_dir, 'PC_tmp_CNVisionFormated.txt'))
            elif calling_soft == 'QuantiSNP':
                formatQuantiSNPforCNVision(list_calling_outputs_path[i], tmp_dir)
                if exists(join(tmp_dir, 'QS_tmp_CNVisionFormated.txt')):
                    list_tmp_paths.append(join(tmp_dir, 'QS_tmp_CNVisionFormated.txt'))
            else:
                raise Exception("Given calling software: {} in't support by the software please contact us if you want to add this software".format(calling_soft))
        runMergingScript(list_tmp_paths, tmp_dir)
        CNVs_list = pd.read_csv(join(tmp_dir,'Sum_CNVisionMerged_PC_QS_tmp.txt'), sep = '\t')
    return CNVs_list


def main():
    print(mergeMultipleCNVCallingOutputs(['/home/thomas/Documents/temp_data/PC_allCNV.txt', '/home/thomas/Documents/temp_data/QS_allCNV.txt'], ['PennCNV', 'QuantiSNP']).head())

if __name__ == "__main__":
    main()