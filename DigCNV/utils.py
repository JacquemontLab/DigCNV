import configparser


def getConfigFileExample(exmpl_conf_output: str):
    """Write an example config file used for next steps

    Args:
        exmpl_conf_output (str): pathway where the example config file will be written
    """    
    config_file = configparser.ConfigParser()
    config_file["Inputs"] = {
        "PC_output_path": "<Path to PennCNV file>",
        "QS_output_path": "<Path to QuantiSNP file>"
    }
    config_file["CallRates"] = {
        "CallRate_path": "<Path to CallRates file>",
        "individual_colname": "<Column name of individials>",
        "callrates_colname": "<Column name of Callrate values>"
    }
    config_file["PFBs"] = {
        "PFB_path": "<Path to PFB file>"
    }

    with open(exmpl_conf_output, 'w') as configfileObj:
        config_file.write(configfileObj)
        configfileObj.flush()
        configfileObj.close()


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
    parameters["PC"] = config_file.get('Pathways', "PC_output_path")
    parameters["QS"] = config_file.get('Pathways', "QS_output_path")

    parameters["CR_path"] = config_file.get('CallRates', "CallRate_path")
    parameters["ind_name"] = config_file.get('CallRates', "individual_colname")
    parameters["CR_name"] = config_file.get('CallRates', "callrates_colname")

    parameters["PFB"] = config_file.get('PFBs', "PFB_path")
    return parameters
