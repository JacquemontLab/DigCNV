# DigCNV
### A machine learning software discrimating true from false CNVs
CNVs called by algorithms such as *PennCNV* or *QuantiSNP* have high false positive rates.
False CNVs bring noise to analysis and could distort a diagnosis
Thus we developped a algorithm analysing multiples CNVs features to remove these false CNVs 

## Classify CNVs with one line

``` sh
# To run only once to install package
pip install digcnv

python3 -m digcnv [-v] <Path to config file>
```
#### Config file
```
[Inputs]
pc_output_path = Path to the PennCNV output file
pc_qc_path = Path to the PennCNV microarray quality file 
qs_output_path = Path to the QuantiSNP output file

[CallRates]
callrate_path = Path to callarate file (.tsv)
individual_colname = Column name of the individuals (Matching names in PennCNV and QuantiSNP outputs)
callrates_colname = Column name of the callrates

[PFBs]
pfb_path = Path to the PFB file 

[Output]
Save_to_file = True
Output_path = /home/thomas/Documents/scripts/DigCNV/temp_data/DigCNV_pred.tsv
```


## Run DigCNV with more options

#### Prepare data

```python
from digcnv import CNVision, dataPreparation

# Import and merge CNV coming from two CNV calling algorthims (only PennCNV and QuantiSNP in this version)
cnvs = CNVision.mergeMultipleCNVCallingOutputs("<list of PennCNV and QuantiSNP output pathways>", ["PennCNV", "QuantiSNP"])

# Add microarray quality data to the list of CNVs
cnvs = dataPreparation.addMicroArrayQualityData(cnvs, "<path to the PennCNV microarray quality file>")

# Compute derived features and add it to the CNV list
cnvs = dataPreparation.addDerivedFeatures(cnvs)

# Add CallRate data to CNV list
cnvs = dataPreparation.addCallRateToDataset(cnvs, call_rate_path="<Pathway to the callrate file>", callrate_colname="<CallRate column name>", individual_colname="<Individual column name>")

# Add Chromosomic information such as centromere and Segmental Duplications overlap
cnvs = dataPreparation.addChromosomicAnnotation(cnvs)

# Add number of probes used with the technology (Useful for model trained on multiple datasets and multiple technologies)   
cnvs = dataPreparation.addNbProbeByTech(cnvs, pfb_file_path="<Pathway to the PFB file>")
 

```
#### Run DigCNV from the pre-trained model
More information at *"article"*

```python
from digcnv import digCnvModel, dataVerif

# Create an empty DigCNV model
model = digCnvModel.DigCnvModel()

# Open pre-trained model and update object
model_path = join(split(__file__)[0], 'data', 'DigCNV_model_multiple_technos.pkl')
model.openPreTrainedDigCnvModel(model_path)

# Check if mandatory columns for the DigCNV model exist and have right formats
dataVerif.checkIfMandatoryColumnsExist(cnvs, post_data_preparation=True)
dataVerif.checkColumnsformats(cnvs, post_data_preparation=True)

# Optional plot a correlation heatmap between different predictors used in model
dataVerif.plotCorrelationHeatMap(cnvs, list_dim=model._dimensions, output_path="<Pathway where output plot (.pdf or .png)>")

# Check NaN data within mandatory columns and split data into two dataframes: first for CNVs with all information available
# and a second one with all CNVs with at least one missing data (can't be used for prediction)
cnvs, cnvs_with_na = dataVerif.computeNaPercentage(cnvs, dimensions=model._dimensions, remove_na_data=True)

# Discriminate true from false CNVs from CNVs with all data, then produce a list of classes
predicted_cnvs = model.predictCnvClasses(cnvs)
cnvs["DigCNVpred"] = predicted_cnvs

```

#### Train your own DigCNV model

```python
from digcnv import digCNVModel
# Uses CNVs created in Prepare data section :

# Add your own annotation to your CNVs (0 for false CNVs and 1 for true)
cnvs["visualized_class"] = ["<classes of each CNV visualized>"]

# Dimensions used to classify CNVs (Recommended dimensions you can use your own)
predictors = ["",""]



# Remove CNVs with at least one missing values in used predictors or in visualized column
cnvs, removed = DigCnvPreProcessing.removeLinesWithNA(cnvs, dimensions=predictors + ["visualized_class"])

# Split dataset into two groups a training dataset and a testing dataset (70% - 30%)  
X_train, y_train, X_test, y_test = DigCnvPreProcessing.createTrainingTestingDatasets(cnvs, X_dimension="visualized_class")

# If ratio between the two classes is too unbalanced uniformize classes by split majoritary class and adding new pseudo CNVs to minority class
X_train, y_train = DigCnvPreProcessing.uniformizeClassesSizes(X_train, y_train, 17, 0.4, 0.5)

# Create a DigCNV model
model = digCnvModel.DigCnvModel()
model.createDigCnvClassifier()

# Train the DigCNV model with the given training dataset
model.trainDigCnvModel(training_data=X_train, training_cat=y_train)

# Analyse classification accuracy for the trained model
predicted_cnvs = model.evaluateCnvClassification(testing_df=X_test, expected_values=y_test, images_dir_path="<Path to the image output directory>")
```

#### Tune DigCNV hyperparameters

```python
from digcnv import DigCNVTunning
```

## Files format example
#### PennCNV output
#### PennCNV quality output
#### QuantiSNP output
#### CallRate file
#### PFB file