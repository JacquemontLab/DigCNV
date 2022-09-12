from sklearn.ensemble import RandomForestClassifier, VotingClassifier, GradientBoostingClassifier, AdaBoostClassifier
from sklearn.metrics import roc_auc_score, accuracy_score, f1_score, RocCurveDisplay
import joblib
import pandas as pd
import matplotlib.pyplot as plt
import warnings


from digcnv import digCNV_logger

class DigCnvModel:
    
    def __init__(self):
        self._rf_params = {"n_estimators": 142,
                  "max_depth": 24,
                  "min_samples_split": 2,
                  "min_samples_leaf": 1,
                  "max_leaf_nodes": 151,
                  "min_weight_fraction_leaf":0.1}
        self._gtb_params = {"n_estimators": 209,
                  "max_depth": 12,
                  "min_samples_split": 12,
                  "min_samples_leaf": 31,
                  "learning_rate" : 0.5}
        self._ada_params = {"n_estimators": 117,
                  "learning_rate" : 0.1}

        self._model = None
        self._dimensions = []
        digCNV_logger.logger.info("Empty DigCNV model created")

    @property
    def rf_params(self) -> dict:
         return self._rf_params
       
    # a setter function
    @rf_params.setter
    def rf_params(self, params: dict):
        mandatory_params = {"n_estimators", "max_depth", "min_samples_split", "min_samples_leaf", "max_leaf_nodes", "min_weight_fraction_leaf"}
        intersect = params.keys() & mandatory_params
        if intersect < 6:
            raise ValueError("Sorry at least of the mandatory hyperparameter is missing: {}".format(mandatory_params))
        self._rf_params = params
        digCNV_logger.logger.info("Random Forest hyperparameters set")

    @property
    def gtb_params(self) -> dict:
         return self._gtb_params
       
    # a setter function
    @gtb_params.setter
    def gtb_params(self, params: dict):
        mandatory_params = {"n_estimators", "max_depth", "min_samples_split", "min_samples_leaf", "learning_rate"}
        intersect = params.keys() & mandatory_params
        if intersect < 6:
            raise ValueError("Sorry at least of the mandatory hyperparameter is missing: {}".format(mandatory_params))
        self._gtb_params = params
        digCNV_logger.logger.info("Gradient Tree boosting hyperparameters set")

    @property
    def ada_params(self) -> dict:
         return self._ada_params
       
    # a setter function
    @ada_params.setter
    def ada_params(self, params: dict):
        mandatory_params = {"n_estimators", "learning_rate"}
        intersect = params.keys() & mandatory_params
        if intersect < 6:
            raise ValueError("Sorry at least of the mandatory hyperparameter is missing: {}".format(mandatory_params))
        self._ada_params = params
        digCNV_logger.logger.info("Ada boosting hyperparameters set")

    def createDigCnvClassifier(self, rf_params={}, gtb_params={}, ada_params={}) -> VotingClassifier:
        if rf_params is None:
            rf_params = self.rf_params

        rf_clf = RandomForestClassifier(n_estimators = rf_params["n_estimators"], 
                                 max_depth = rf_params["max_depth"], 
                                 min_samples_split = rf_params["min_samples_split"], 
                                 min_samples_leaf = rf_params["min_samples_leaf"], 
                                 max_leaf_nodes = rf_params["max_leaf_nodes"],
                                 min_weight_fraction_leaf = 0.0,
                                random_state=42)
        if gtb_params is None:
            gtb_params = self.gtb_params
        
        gtb_clf = GradientBoostingClassifier(n_estimators = gtb_params["n_estimators"], 
                                     max_depth = gtb_params["max_depth"], 
                                     min_samples_split = gtb_params["min_samples_split"], 
                                     min_samples_leaf = gtb_params["min_samples_leaf"], 
                                     learning_rate = gtb_params["learning_rate"],
                                     min_weight_fraction_leaf = 0.0,
                                    random_state=42)
        if ada_params is None:
            ada_params = self.ada_params
            
        ada_clf = AdaBoostClassifier(n_estimators=ada_params["n_estimators"],
                                     learning_rate=ada_params["learning_rate"],
                                    random_state=42)


        voting_clf = VotingClassifier(estimators=[('Random Forest', rf_clf), ("Gradient Tree Boosting", gtb_clf), ("Ada Boosting", ada_clf)], voting='soft')
        self._model = voting_clf
        return voting_clf
    
    def openPreTrainedDigCnvModel(self, model_path:str):
        has_warning = False
        with warnings.catch_warnings(record=True) as w:
            dimensions, model = joblib.load(model_path)
            digCNV_logger.logger.info("Pre trained model loaded from {}".format(model_path))
            if len(w) > 0:
                digCNV_logger.logger.info("Version Warning: {}".format(w[0]))
            self._dimensions = dimensions
            digCNV_logger.logger.info("Pre trained model will use {} as predictors".format(dimensions))
            self._model = model


    def checkIfDigCnvFitted(self) -> bool: 
        return hasattr(self._model, "classes_")
    
    def trainDigCnvModel(self, training_data: pd.DataFrame, training_cat: pd.Series):
        self._dimensions = training_data.columns.tolist()
        self._model.fit(training_data, training_cat)
    
    def saveDigCnvModelToPkl(self, output_path:str):
        if self.checkIfDigCnvFitted():
            joblib.dump([self._dimensions, self._model], output_path)
        else:
            raise Exception("DigCNV model not defined!\nSaving the model impossible")
    
    def predictCnvClasses(self, cnvs: pd.DataFrame, use_percentage=False) -> pd.DataFrame:
        split_cnvs = cnvs.loc[:,self._dimensions]
        if self.checkIfDigCnvFitted():
            predictions = self._model.predict(split_cnvs)
            digCNV_logger.logger.info("CNVs classes are now predicted by the model")
            cnvs["DigCNVpred"] = predictions
            if use_percentage:
                predict_proba = self._model.predict_proba(split_cnvs)
                digCNV_logger.logger.info("Classes probabilities added to CNV resutls")
                digCNV_logger.logger.info(predict_proba)
        else:
            raise Exception("DigCNV model not defined!")
        return cnvs
    
    def evaluateCnvClassification(self, testing_df: pd.DataFrame, expected_values: pd.Series, images_dir_path=""):
        split_cnvs = testing_df.loc[:,self._dimensions]
        if self.checkIfDigCnvFitted():
            predictions = self._model.predict(split_cnvs)
        else:
            raise Exception("DigCNV model isn't trained so you can't perform classifications")
        results = pd.DataFrame({"predict": predictions, "bon": expected_values})
        groups_count = results.groupby(['predict', 'bon']).size()
    
        digCNV_logger.logger.info(f"Sensibility : {groups_count[1][1]/groups_count[1].sum():.3f}")
        digCNV_logger.logger.info(f"Specificity : {groups_count[0][0]/groups_count[0].sum():.3f}")
        digCNV_logger.logger.info(f"AUC : {roc_auc_score(expected_values, predictions):.2f}%")
        digCNV_logger.logger.info(f"Accuracy : {accuracy_score(expected_values, predictions):.3f}")
        digCNV_logger.logger.info(f"F1 Score : {f1_score(expected_values, predictions):.3f}")
    
        RocCurveDisplay.from_estimator(self._model, testing_df, expected_values)
        if images_dir_path != "":
            plt.savefig("{}/ROC_curve.pdf".format(images_dir_path))
        plt.show()
        plt.close()
        
        proba = self._model.predict_proba(split_cnvs)
        proba = pd.DataFrame(proba)
        proba["predict"] = predictions
        proba["true_class"] = expected_values.tolist()
        proba.columns = ["C0","C1", "predict","true_class"]
        split_cnvs = split_cnvs.reset_index()
        full_results = pd.concat([split_cnvs, proba], axis=1)
        ax = full_results.pivot(columns='true_class', values='C1').plot.hist(bins=100, figsize=(12,8), color=["red", "green"], label=["CNV", "No CNV"], stacked=True)
        plt.legend()
        if images_dir_path != "":
            plt.savefig("{}/proba_distribution.pdf".format(images_dir_path))
        plt.show()
        plt.close()
