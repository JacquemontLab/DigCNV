from DigCNV import digCNV_logger
import pandas as pd
from sklearn.model_selection import train_test_split
from imblearn.over_sampling import SMOTE
from imblearn.under_sampling import RandomUnderSampler
from imblearn.pipeline import Pipeline

def removeLinesWithNA(cnvs, dimensions):
    removed_cnvs = pd.DataFrame()
    cnvs_clean = cnvs.copy()
    split_cnvs = cnvs.loc[:,dimensions]
    nb_line = split_cnvs.shape[0]-1
    i = 0
    while i < len(split_cnvs.columns):
        if split_cnvs.iloc[:, i].count()-1 != nb_line:
            removed_index = split_cnvs[split_cnvs[split_cnvs.columns[i]].isnull()].index.tolist()
            digCNV_logger.logger.info("{} CNV removed due to Null Data".format(len(removed_index)))
            split_cnvs = split_cnvs.drop(index=removed_index)
            removed_cnvs = pd.concat([removed_cnvs, cnvs_clean.iloc[removed_index]])
            cnvs_clean = cnvs_clean.drop(index=removed_index)
            cnvs_clean.reset_index(drop=True, inplace = True)
            split_cnvs.reset_index(drop=True, inplace = True)
    return cnvs_clean, removed_cnvs

def createTrainingTestingDatasets(data, X_dimension = "SnipPeep_Ok", all_data_set=False):
    if all_data_set:
        X_train = data.drop(columns=[X_dimension])
        X_test = data.drop(columns=[X_dimension])
        y_train = data.SnipPeep_Ok.tolist()
        y_test = data.SnipPeep_Ok.tolist()
    else:
        data_rf = data.drop(columns=[X_dimension])
        labels = data.SnipPeep_Ok
        X_train, X_test, y_train, y_test = train_test_split(data_rf, labels, random_state=42, test_size = 0.33, shuffle=True)
    digCNV_logger.logger.info("Training dataset")
    digCNV_logger.logger.info(X_train.shape)
    digCNV_logger.logger.info("Testing dataset")
    digCNV_logger.logger.info(X_test.shape)
    return X_train, X_test, y_train, y_test

def uniformizeClassesSizes(X_train, y_train, k_neighbors, over_sampling, under_sampling):
    digCNV_logger.logger.info("Training dataset classes")
    digCNV_logger.logger.info(pd.Series(y_train).value_counts())
    over = SMOTE(sampling_strategy=over_sampling, k_neighbors = k_neighbors, n_jobs=-1, random_state=42)
    under = RandomUnderSampler(sampling_strategy=under_sampling, random_state=42)
    steps = [('o', over), ('u', under)]
    pipeline = Pipeline(steps=steps)
    X_train, y_train  = pipeline.fit_resample(X_train, y_train)
    digCNV_logger.logger.info("Training dataset after uniformizing classes")
    digCNV_logger.logger.info(pd.Series(y_train).value_counts())
    return X_train, y_train