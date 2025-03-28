from digcnv import digCNV_logger
import pandas as pd
from sklearn.model_selection import train_test_split
from imblearn.over_sampling import SMOTE
from imblearn.under_sampling import RandomUnderSampler
from imblearn.pipeline import Pipeline


def createTrainingTestingDatasets(cnvs:pd.DataFrame, dimensions:list, X_dimension:str, all_data_set=False) -> tuple[pd.DataFrame, pd.Series, pd.DataFrame, pd.Series]:
    """Split CNV list into two groups randomly. A first one used as training dataset and a second used as testing dataset.
    Created dataset will only contains useful features for the model and the column indicating the true classification will be put aside into two separate Series.

    :param cnvs: List of CNVs with all their features
    :type cnvs: pd.DataFrame
    :param dimensions: List of dimension to use in DigCNV model, must contains only integer or float columns.
    :type dimensions: list
    :param X_dimension: Column name of the true classification
    :type X_dimension: str
    :param all_data_set: Indicating how to split dataset, if `True` training and testing datasets will contains all CNVs while `False` will split dataset into 2/3 for training and 1/3 for testing, defaults to False
    :type all_data_set: bool, optional
    :return: A tuple containing 4 objects, the training dataframe and classes in first part and the testing dataframe and classes in second part
    :rtype: tuple[pd.DataFrame, pd.Series, pd.DataFrame, pd.Series]
    """
    data = cnvs.loc[:,dimensions + [X_dimension]]
    if all_data_set:
        X_train = data.drop(columns=[X_dimension])
        X_test = data.drop(columns=[X_dimension])
        y_train = data[X_dimension].tolist()
        y_test = data[X_dimension].tolist()
    else:
        data_rf = data.drop(columns=[X_dimension])
        labels = data[X_dimension]
        X_train, X_test, y_train, y_test = train_test_split(data_rf, labels, random_state=42, test_size = 0.33, shuffle=True)
    digCNV_logger.logger.info("Training dataset")
    digCNV_logger.logger.info(X_train.shape)
    digCNV_logger.logger.info("Testing dataset")
    digCNV_logger.logger.info(X_test.shape)
    return X_train, y_train, X_test, y_test

def uniformizeClassesSizes(X_train:pd.DataFrame, y_train:pd.Series, k_neighbors:int, over_sampling:float, under_sampling:float) -> tuple[pd.DataFrame, pd.Series]:
    """Wil uniformize classes for unbalanced datasets

    :param X_train: Training dataframe with only useful features.
    :type X_train: pd.DataFrame
    :param y_train: Training Series containing all true classes
    :type y_train: pd.Series
    :param k_neighbors: Number of K nearest neighbors used to create new Random CNV of the minority class
    :type k_neighbors: int
    :param over_sampling: Ratio of over-sampling (must be lower than under_sampling parameter)
    :type over_sampling: float
    :param under_sampling: Ratio of under sampling (percentage of CNVs to remove from the majority class)
    :type under_sampling: float
    :return: a tuple containing a dataframe and a series with more balanced classes
    :rtype: tuple[pd.DataFrame, pd.Series]
    """    
    digCNV_logger.logger.info("\nTraining dataset classes\n{}".format(y_train.value_counts()))
    over = SMOTE(sampling_strategy=over_sampling, k_neighbors = k_neighbors, random_state=42)
    under = RandomUnderSampler(sampling_strategy=under_sampling, random_state=42)
    steps = [('o', over), ('u', under)]
    pipeline = Pipeline(steps=steps)
    X_train, y_train  = pipeline.fit_resample(X_train, y_train)
    digCNV_logger.logger.info("\nTraining dataset after uniformizing classes\n{}".format(y_train.value_counts()))
    return X_train, y_train