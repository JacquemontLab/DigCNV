from os.path import exists, split, join

import utils
import CNVision
import dataPreparation
import dataVerif
import digCnvModel

def runWholeAnnotationProcess(conf_file_path:str):
    parameters = utils.readDigCNVConfFile(conf_file_path)
    cnvs = CNVision.mergeMultipleCNVCallingOutputs([parameters["PC"],parameters["QS"]], ["PennCNV", "QuantiSNP"])

    cnvs = dataPreparation.addDerivedFeatures(cnvs)
    cnvs = dataPreparation.addCallRateToDataset(cnvs, call_rate_path=parameters["CR_path"], callrate_colname=parameters["CR_name"], individual_colname=parameters["ind_name"])
    
    cnvs = dataPreparation.addNbProbeByTech(cnvs, pfb_file_path=parameters["PFB"])
    cnvs = dataPreparation.transformTwoAlgsFeatures(cnvs)

    model = digCnvModel.DigCnvModel()
    model_path = join(split(__file__)[0], 'data', 'DigCnvModel_Trained_Mega_Spark_Ukbb.pkl')
    model.openPreTrainedDigCnvModel(model_path)

    dataVerif.checkIfMandatoryColumnsExist(cnvs, post_data_preparation=True)
    dataVerif.checkColumnsformats(cnvs, post_data_preparation=True)
    dataVerif.computeNaPercentage(cnvs, dimensions=model._dimensions)
    predicted_cnvs = model.predictCnvClasses(cnvs)
    if parameters["save"]:
        predicted_cnvs.to_csv(parameters["output_path"])
    return predicted_cnvs



def main():
    utils.getConfigFileExample("/home/thomas/Documents/temp_data/test.conf")
    runWholeAnnotationProcess("/home/thomas/Documents/temp_data/DC.conf")

if __name__ == "__main__":
    main()