from sklearn.ensemble import RandomForestClassifier, VotingClassifier, BaggingClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn import preprocessing

from sklearn.metrics import roc_auc_score, accuracy_score, f1_score, RocCurveDisplay
import joblib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings


from digcnv import digCNV_logger


class DigCnvModel:
    """Class of the DigCNV model
    """

    def __init__(self):
        """Creating the DigCNV instance and setting potential hyperparameters But no model is created you must create the classifier or load a pretrained model.
        """
        self.rf_params  = {"n_estimators":790,
                    "max_depth": 352,
                    "min_samples_split": 14,
                    "min_samples_leaf": 1,
                    "max_leaf_nodes": 495,
                    "min_weight_fraction_leaf":0.0}

        self.bg_knn_params = {'n_estimators': 191,
                        'max_samples' :0.35,
                        'estimator__n_neighbors':1}

        self.svm_params = {'C':49.8,
                    'gamma' : 'scale',
                    'tol': 0.008858667904100823}

        self._model = None
        self._dimensions = []
        self._dimensions_scales = {}
        digCNV_logger.logger.info("Empty DigCNV model created")

    @property
    def rf_params(self) -> dict:
        """get the list of Random forest hyperparameters

        :return: dictionnary of Random forest hyperparameters
        :rtype: dict
        """
        return self._rf_params

    # a setter function
    @rf_params.setter
    def rf_params(self, params: dict):
        """set the list of Random forest hyperparameters

        :param params: dictionary containing Random forest hyperparameters to set
        :type params: dict
        :raises ValueError: If at least one of the dictonnary value is missing
        """
        mandatory_params = {"n_estimators", "max_depth", "min_samples_split",
                            "min_samples_leaf", "max_leaf_nodes", "min_weight_fraction_leaf"}
        intersect = params.keys() & mandatory_params
        if len(intersect) < 6:
            raise ValueError(
                "Sorry at least one of the mandatory hyperparameter is missing: {}".format(mandatory_params))
        self._rf_params = params
        digCNV_logger.logger.info("Random Forest hyperparameters set")

    @property
    def bg_knn_params(self) -> dict:
        """get the list of Gradient tree boosting hyperparameters

        :return: dictionnary of Gradient tree boosting hyperparameters
        :rtype: dict
        """
        return self._bg_knn_params

    # a setter function
    @bg_knn_params.setter
    def bg_knn_params(self, params: dict):
        """set the list of Gradient tree boosting hyperparameters

        :param params: dictionary containing Gradient tree boosting hyperparameters to set
        :type params: dict
        :raises ValueError: If at least one of the dictonnary value is missing
        """
        mandatory_params = {"n_estimators", "max_samples",
                            "estimator__n_neighbors"}
        intersect = params.keys() & mandatory_params
        if len(intersect) < 3:
            raise ValueError(
                "Sorry at least of one the mandatory hyperparameter is missing: {}".format(mandatory_params))
        self._bg_knn_params = params
        digCNV_logger.logger.info("Bagging of KNN hyperparameters set")

    @property
    def svm_params(self) -> dict:
        """get the list of SVC hyperparameters

        :return: dictionnary of SVC hyperparameters
        :rtype: dict
        """
        return self._svm_params

    # a setter function
    @svm_params.setter
    def svm_params(self, params: dict):
        """set the list of SVC hyperparameters

        :param params: dictionary containing SVC hyperparameters to set
        :type params: dict
        :raises ValueError: If at least one of the dictonnary value is missing
        """
        mandatory_params = {"C", "gamma", "tol"}
        intersect = params.keys() & mandatory_params
        if len(intersect) < 3:
            raise ValueError(
                "Sorry at least one of the mandatory hyperparameter is missing: {}".format(mandatory_params))
        self._svm_params = params
        digCNV_logger.logger.info("SVC hyperparameters set")

    def createDigCnvClassifier(self, rf_params=None, bg_knn_params=None, svm_params=None) -> VotingClassifier:
        """Create the DigCNV classifier based on three Classifiers, a Random forest, a bagging of KNN and SVC.
        You can set dictionnaries of hyperparameters for each machine learning model or used hyperparameters stored in object by letting arguments empty

        :param rf_params: dictionary containing Random forest hyperparameters to set, defaults to None
        :type rf_params: dict, optional
        :param bg_knn_params: dictionary containing Bagging of KNN hyperparameters to set, defaults to None
        :type bg_knn_params: dict, optional
        :param svm_params: dictionary containing SVC hyperparameters to set, defaults to None
        :type svm_params: dict, optional
        :return: The DigCNV model created and ready for training.
        :rtype: VotingClassifier
        """
        if rf_params is None:
            rf_params = self.rf_params
        rf_clf = RandomForestClassifier(n_estimators = rf_params["n_estimators"],
                                        max_depth = rf_params["max_depth"],
                                        min_samples_split = rf_params["min_samples_split"],
                                        min_samples_leaf = rf_params["min_samples_leaf"],
                                        max_leaf_nodes = rf_params["max_leaf_nodes"],
                                        min_weight_fraction_leaf = 0.0,
                                        random_state=42)

        if bg_knn_params is None:
            bg_knn_params = self.bg_knn_params

        knn_clf = BaggingClassifier(estimator = KNeighborsClassifier(weights ="distance",
                                                                        n_neighbors = bg_knn_params["estimator__n_neighbors"]),
                                    bootstrap = True,
                                    bootstrap_features = True,
                                    n_estimators = bg_knn_params["n_estimators"],
                                    max_samples = bg_knn_params["max_samples"],
                                    random_state=42)
        if svm_params is None:
            svm_params = self.svm_params

        svm_clf = SVC(gamma = svm_params['gamma'],
                            C = svm_params['C'],
                            tol = svm_params['tol'],
                            probability=True,
                            random_state=42)

        voting_clf = VotingClassifier(estimators=[('Random Forest', rf_clf), (
            "Bagging KNN", knn_clf), ("SVC", svm_clf)], voting='soft')
        self._model = voting_clf
        return voting_clf

    def openPreTrainedDigCnvModel(self, model_path: str):
        """Open a pre-trained DigCNV model and update the DigCNV object

        :param model_path: Pathway to a pkl object containing the list of dimensions used in the model and the trained model. 
        :type model_path: str
        """        
        has_warning = False
        with warnings.catch_warnings(record=True) as w:
            dimensions, model = joblib.load(model_path)
            digCNV_logger.logger.info(
                "Pre trained model loaded from {}".format(model_path))
            if len(w) > 0:
                digCNV_logger.logger.info("Version Warning: {}".format(w[0]))

            self._dimensions = list(dimensions.keys())
            self._dimensions_scales = dimensions
            digCNV_logger.logger.info(
                "Pre trained model will use {} as predictors".format(dimensions.keys()))
            self._model = model

    def checkIfDigCnvFitted(self) -> bool:
        """Check if the DigCNV model object has been trained or not.
        A not trained model cannot be used for classification.

        :return: Boolean indicating if the model trained, `True` if the model is trained
        :rtype: bool
        """        
        return hasattr(self._model, "classes_")

    def trainDigCnvModel(self, training_data: pd.DataFrame, training_cat: pd.Series):
        """ Train the DigCNV model object thnaks to the given training data

        :param training_data: A dataframe containing all important features to describe CNVs
        :type training_data: pd.DataFrame
        :param training_cat: A list of binary annotation for CNVs indicating if each CNV is a True CNV or an artefact, must int values.
        :type training_cat: pd.Series
        """   
        scaler = preprocessing.StandardScaler()
        cols = training_data.columns
        X_train_scale = pd.DataFrame(scaler.fit_transform(training_data, training_cat), columns=cols)
        digCNV_logger.logger.info(
                "Predictor scaled for training")
        # for col in cols:
            # self._dimensions_scales[col] = [training_data[col].min(), training_data[col].max()]         
        for col in cols:
            self._dimensions_scales[col] = [training_data[col].mean(), training_data[col].std()] 
        self._dimensions = X_train_scale.columns.tolist()
        self._model.fit(X_train_scale, training_cat)
        digCNV_logger.logger.info(
                f"Model trained on the {X_train_scale.shape[0]} CNVs with {X_train_scale.shape[1]} features")

    def saveDigCnvModelToPkl(self, output_path: str):
        """Save a trained DigCNV model to a pkl file to be used later

        :param output_path: Pathway of the pkl file must be finishing by `,pkl`
        :type output_path: str
        :raises Exception: if the model isn't trained 
        """        
        if self.checkIfDigCnvFitted():
            joblib.dump([self._dimensions_scales, self._model], output_path)
        else:
            raise Exception(
                "DigCNV model not defined!\nSaving the model impossible")

    def predictCnvClasses(self, cnvs: pd.DataFrame, use_percentage=False) -> pd.DataFrame:
        """Will predict the CNVs classification based on the dataframe of CNV features given. For pre-trained models, classes are `0` for False CNVs and `1` for True CNVs

        :param cnvs: DataFrame containing describing features
        :type cnvs: pd.DataFrame
        :param use_percentage: Indicate if results must be binary or probabilities, `True` if you want to have probabilities, `False` if you want only classification, defaults to False
        :type use_percentage: bool, optional
        :raises Exception: if model isn't trained
        :return: CNVs with their classification aggregated
        :rtype: pd.DataFrame
        """        
        split_cnvs = cnvs.loc[:, self._dimensions]
        if self.checkIfDigCnvFitted():
            
            # Scale the data based on the training data
            for col in split_cnvs.columns:
                split_cnvs[col] = (split_cnvs[col] - self._dimensions_scales[col][0]) / self._dimensions_scales[col][1]

            predict_proba = self._model.predict_proba(split_cnvs)
            digCNV_logger.logger.info(
                "CNVs classes are now predicted by the model")
            if use_percentage:
                cnvs["class_1"] = predict_proba[:, 1]
                cnvs["class_0"] = predict_proba[:, 0]
                digCNV_logger.logger.info(
                    "Classes probabilities added to CNV resutls")
                digCNV_logger.logger.info(predict_proba)
            predictions = np.where(predict_proba[:, 1] > 0.5, 1, 0)
            cnvs["DigCNVpred"] = predictions
        else:
            raise Exception("DigCNV model not defined!")
        return cnvs

    def evaluateCnvClassification(self, testing_df: pd.DataFrame, expected_values: pd.Series, images_dir_path=None):
        """Evaluate a trained model with a list CNVs with already none classification.

        :param testing_df: List of CNVs with known classification
        :type testing_df: pd.DataFrame
        :param expected_values: list of the Classification. Must be the same as classes used for training. For example, pre-trained models have `0` for False CNVs and `1` for True CNVs
        :type expected_values: pd.Series
        :param images_dir_path: Pathway to the directory where figures will be plotted, Recommended if not runned within a Notebook, defaults to None
        :type images_dir_path: str, optional
        :raises Exception: if model isn't trained
        """        
        split_cnvs = testing_df.loc[:, self._dimensions]
        if self.checkIfDigCnvFitted():
            predictions = self._model.predict(split_cnvs)
        else:
            raise Exception(
                "DigCNV model isn't trained so you can't perform classifications")
        results = pd.DataFrame(
            {"predict": predictions, "bon": expected_values})
        groups_count = results.groupby(['predict', 'bon']).size()

        digCNV_logger.logger.info(
            f"Sensibility : {groups_count[1][1]/groups_count[1].sum():.3f}")
        digCNV_logger.logger.info(
            f"Specificity : {groups_count[0][0]/groups_count[0].sum():.3f}")
        digCNV_logger.logger.info(
            f"AUC : {roc_auc_score(expected_values, predictions):.2f}%")
        digCNV_logger.logger.info(
            f"Accuracy : {accuracy_score(expected_values, predictions):.3f}")
        digCNV_logger.logger.info(
            f"F1 Score : {f1_score(expected_values, predictions):.3f}")

        RocCurveDisplay.from_estimator(
            self._model, testing_df, expected_values)
        if images_dir_path != "":
            plt.savefig("{}/ROC_curve.pdf".format(images_dir_path))
        plt.show()
        plt.close()

        proba = self._model.predict_proba(split_cnvs)
        proba = pd.DataFrame(proba)
        proba["predict"] = predictions
        proba["true_class"] = expected_values.tolist()
        proba.columns = ["C0", "C1", "predict", "true_class"]
        split_cnvs = split_cnvs.reset_index()
        full_results = pd.concat([split_cnvs, proba], axis=1)
        ax = full_results.pivot(columns='true_class', values='C1').plot.hist(
            bins=100, figsize=(12, 8), color=["red", "green"], label=["CNV", "No CNV"], stacked=True)
        plt.legend()
        if images_dir_path != "":
            plt.savefig("{}/proba_distribution.pdf".format(images_dir_path))
        plt.show()
        plt.close()
