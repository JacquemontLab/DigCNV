import matplotlib.pyplot as plt
import time
import pandas as pd
from digcnv import digCNV_logger

from sklearn.base import BaseEstimator
from sklearn.experimental import enable_halving_search_cv
from sklearn.model_selection import GridSearchCV, HalvingGridSearchCV

def evaluate_param(model:BaseEstimator, parameter:str, X_train:pd.DataFrame, y_train:pd.Series, num_range:list, index:int, validation_score:str) -> tuple[plt.plot, pd.DataFrame]:
    """Plot estimator performance when variating the hyperparameter value

    :param model: Scikit-learn model classifier to tune
    :type model: BaseEstimator
    :param parameter: hyperparameter name to test on the model
    :type parameter: str
    :param X_train: Training dataframe with only useful features.
    :type X_train: pd.DataFrame
    :param y_train: Training Series containing all true classes
    :type y_train: pd.Series
    :param num_range: hyperparameter value list to test and plot results
    :type num_range: list
    :param index: subplot index for the global plot
    :type index: int
    :param validation_score: Validation score used to evaluate model performance. must be a validation score available in scikit-learn package (`https://scikit-learn.org/stable/modules/model_evaluation.html#scoring-parameter`)
    :type validation_score: str
    :return: return two object, first a plot representing model performace based on hyperparamater value, second a dataframe representing these data.
    :rtype: tuple[plt.plot, pd.DataFrame]
    """    
    grid_search = GridSearchCV(model, param_grid = {parameter: num_range}, n_jobs = -1, scoring = validation_score)
    grid_search.fit(X_train, y_train)
    df = {}
    for i, score in enumerate(grid_search.cv_results_["mean_test_score"]):
        df[grid_search.cv_results_["params"][i][parameter]] = score
#         df[score[0][parameter]] = score[1]
    df = pd.DataFrame.from_dict(df, orient='index')
    df.reset_index(level=0, inplace=True)
    df = df.sort_values(by='index')

    plt.subplot(3,2,index)
    plot = plt.plot(df['index'], df[0])
    plt.title(parameter)
    return plot, df

def evaluateHyperparameterIndividually(model:BaseEstimator, param_dict:dict, X_train:pd.DataFrame, y_train:pd.Series, validation_score:str):
    """Plot multiple performance chart for each hyperparameter individually

    :param model: Scikit-learn model classifier to tune
    :type model: BaseEstimator
    :param param_dict: Dictionnary containing hyperparameters and their value range or list
    :type param_dict: dict
    :param X_train: Training dataframe with only useful features.
    :type X_train: pd.DataFrame
    :param y_train: Training Series containing all true classes
    :type y_train: pd.Series
    :param validation_score: Validation score used to evaluate model performance. must be a validation score available in scikit-learn package (`https://scikit-learn.org/stable/modules/model_evaluation.html#scoring-parameter`)
    :type validation_score: str
    """    
    start_time = time.time()
    index = 1
    plt.figure(figsize=(16,12))
    for parameter, param_range in param_dict.items():
        evaluate_param(model, parameter, X_train, y_train, param_range, index, validation_score)
        index += 1
    plt.show()
    digCNV_logger.logger.info("--- %s seconds ---" % (time.time() - start_time))
    
def evaluateHyperparameterGrid(X_train:pd.DataFrame, y_train:pd.Series, model:BaseEstimator, param_grid:dict, validation_score:str) -> dict:
    """_summary_

    :param X_train: Training dataframe with only useful features.
    :type X_train: pd.DataFrame
    :param y_train: Training Series containing all true classes
    :type y_train: pd.Series
    :param model: Scikit-learn model classifier to tune
    :type model: BaseEstimator
    :param param_grid: Dictionnary containing hyperparameters and their value range or list
    :type param_grid: dict
    :param validation_score: Validation score used to evaluate model performance. must be a validation score available in scikit-learn package (`https://scikit-learn.org/stable/modules/model_evaluation.html#scoring-parameter`)
    :type validation_score: str
    :return: Dictionnary containing each hyperparameters with and their best values 
    :rtype: dict
    """    
    start_time = time.time()
    grid_search = HalvingGridSearchCV(estimator = model, param_grid = param_grid, min_resources=100,
                              cv = 5, n_jobs = -1, verbose = 0, scoring = validation_score, random_state=42)
    grid_search.fit(X_train, y_train)
    digCNV_logger.logger.info(grid_search.best_params_)
    params = grid_search.best_params_
    digCNV_logger.logger.info("--- %s seconds ---" % (time.time() - start_time))
    return params
