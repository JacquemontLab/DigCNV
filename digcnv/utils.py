from digcnv import digCNV_logger

import configparser


def getConfigFileExample(exmpl_conf_output: str):
    """Write an example config file used for next steps

    Args:
        exmpl_conf_output (str): pathway where the example config file will be written
    """    
    config_file = configparser.ConfigParser()
    config_file["Inputs"] = {
        "PC_output_path": "<Path to PennCNV file>",
        "PC_QC_path" : "<Path to PennCNV microarray quality file>",
        "QS_output_path": "<Path to QuantiSNP file>"
    }
    # config_file["CallRates"] = {
    #     "CallRate_path": "<Path to CallRates file>",
    #     "individual_colname": "<Column name of individials>",
    #     "callrates_colname": "<Column name of Callrate values>"
    # }
    # config_file["PFBs"] = {
    #     "PFB_path": "<Path to PFB file>"
    # }
    config_file["Output"] = {
        "Save_to_file": "<True or False if you want to save the dataset to a file>",
        "Output_path": "<Path of the future output file>"
    }

    config_file["DigCNV"] = {
        "model_path": "<Path of the downloaded model. Available at : https://murena.io/s/xEsyae6gxfMEnWJ>"
    }

    with open(exmpl_conf_output, 'w') as configfileObj:
        config_file.write(configfileObj)
        configfileObj.flush()
        configfileObj.close()
    digCNV_logger.logger.info("Example Config file written to {}".format(exmpl_conf_output))


def readDigCNVConfFile(conf_file_path: str) -> dict:
    """Reads the given config file and return a dictionnary of named parameters

    Args:
        conf_file_path (str): pathway of the config file

    Returns:
        dict: Dictionnary containing all parameters values
    """ 
    parameters = {}
    config_file = configparser.ConfigParser()
    config_file.read(conf_file_path)
    digCNV_logger.logger.info("Config file {} opened".format(conf_file_path))
    parameters["PC"] = config_file.get('Inputs', "PC_output_path")
    parameters["QS"] = config_file.get('Inputs', "QS_output_path")
    parameters["QC"] = config_file.get('Inputs', "PC_QC_path")

    parameters["CallRate"] = config_file.get('Inputs', "callrate_path")
    # parameters["ind_name"] = config_file.get('CallRates', "individual_colname")
    # parameters["CR_name"] = config_file.get('CallRates', "callrates_colname")

    # parameters["PFB"] = config_file.get('PFBs', "PFB_path")
    parameters['DigCnvModel'] = config_file.get('DigCNV', 'model_path')
    parameters["save"] = config_file.get('Output', 'Save_to_file')
    parameters["output"] = config_file.get('Output', 'Output_path')
    digCNV_logger.logger.info("Set of parameters created for DigCNV")
    return parameters
