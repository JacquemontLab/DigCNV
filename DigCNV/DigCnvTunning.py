import matplotlib.pyplot as plt
import time
import pandas as pd

from DigCNV import digCNV_logger

from sklearn.experimental import enable_halving_search_cv
from sklearn.model_selection import GridSearchCV, HalvingGridSearchCV

def evaluate_param(model, parameter, X_train, y_train, num_range, index, validation_score):
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

def evaluateHyperparameterIndividually(model, param_dict, X_train, y_train, validation_score):
    start_time = time.time()
    index = 1
    plt.figure(figsize=(16,12))
    for parameter, param_range in param_dict.items():
        evaluate_param(model, parameter, X_train, y_train, param_range, index, validation_score)
        index += 1
    plt.show()
    digCNV_logger.logger.info("--- %s seconds ---" % (time.time() - start_time))
    
def evaluateHyperparameterGrid(X_train, y_train, model, param_grid, validation_score):
    start_time = time.time()
    grid_search = HalvingGridSearchCV(estimator = model, param_grid = param_grid, min_resources=100,
                              cv = 5, n_jobs = -1, verbose = 0, scoring = validation_score, random_state=42)
    grid_search.fit(X_train, y_train)
    digCNV_logger.logger.info(grid_search.best_params_)
    params = grid_search.best_params_
    digCNV_logger.logger.info("--- %s seconds ---" % (time.time() - start_time))
    return params
