from sklearn.ensemble import RandomForestClassifier, VotingClassifier, GradientBoostingClassifier, AdaBoostClassifier
import joblib
import pandas as pd

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
        dimensions, model = joblib.load(model_path)
        self._dimensions = dimensions
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
            cnvs["DigCNVpred"] = predictions
            if use_percentage:
                predict_proba = self._model.predict_proba(split_cnvs)
        else:
            raise Exception("DigCNV model not defined!")
        return cnvs
    
    def evaluateCnvClassification(self, testing_df: pd.DataFrame, expected_valeues: pd.Series):
        split_cnvs = testing_df.loc[:,self._dimensions]
        if self.checkIfDigCnvFitted():
            predictions = self._model.predict(split_cnvs)
