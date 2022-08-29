# DigCNV
### A machine learning software discrimating true from false CNVs
CNVs called by algorithms such as *PennCNV* or *QuantiSNP* have high false positive rates.
False CNVs bring noise to analysis and could distort a diagnosis
Thus we developped a algorithm analysing multiples CNVs features to remove these false CNVs 

## Data format

| **START** | **STOP** | **CHR** | **SNP** | **SCORE** | **TwoAlgs** | **LRR_mean** | **WF** | **** | **** |
|-----------|----------|---------|---------|-----------|-------------|--------------|--------|------|------|
| int       | int      | int     | int     | float     | int         | float        | float  |      |      |
|           |          |         |         |           |             |              |        |      |      |
|           |          |         |         |           |             |              |        |      |      |
|           |          |         |         |           |             |              |        |      |      |
|           |          |         |         |           |             |              |        |      |      |
|           |          |         |         |           |             |              |        |      |      |
|           |          |         |         |           |             |              |        |      |      |
|           |          |         |         |           |             |              |        |      |      |
|           |          |         |         |           |             |              |        |      |      |

## Prepare dataset

#### Add useful columns

```python
from DigCNVlib import dataPreparation, dataVerif
import pandas as pd

# Open cnv file as pandas dataframe
cnvs = pd.read_csv("data_pathway.csv")

# Check if all mandatory columns exists and have the right format
dataVerif.checkIfMandatoryColumnsExist(cnvs, post_data_preparation=False)

# Create some features used in the model
cnvs = dataPreparation.addDerivedFeatures(cnvs)
cnvs = dataPreparation.addCallRateToDataset(cnvs, call_rate_path='../data/callrates.tsv', callrate_colname='CallRate', individual_colname='SampleID')
cnvs = dataPreparation.addNbProbeByTech(cnvs, pfb_file_path='../data/UKBB_PFB.pfb')
cnvs = dataPreparation.transformTwoAlgsFeatures(cnvs)
```
#### Check data quality
```python
from DigCNVlib import dataVerif

dataVerif.checkColumnsformats(cnvs)
dataVerif.checkIfMandatoryColumnsExist(cnvs, post_data_preparation=False)
dataVerif.computeNaPercentage(cnvs, dimensions=model_dimensions)
dataVerif.plotCorrelationHeatMap(cnvs, list_dim=model_dimensions, output_path='../output_path/correlation.png')

```