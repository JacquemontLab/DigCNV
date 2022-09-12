import logging
import warnings
import pandas as pd

warnings.simplefilter(action='ignore', category=pd.errors.DtypeWarning)

logging.basicConfig(format='%(asctime)s :: %(funcName)s:\t %(message)s', level=logging.INFO)
logger = logging.getLogger()
logger.setLevel(logging.INFO) 

def changeLoggingLevel(verbose = False):
    if verbose == False:
        logger.setLevel(logging.WARNING)
    else:
        logger.setLevel(logging.INFO) 