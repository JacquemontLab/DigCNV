from DigCNV import dataPreparation
import pandas as pd
import pytest 


@pytest.fixture
def global_data():
    cnvs = pd.read_csv('data/UKBB_clean_for_DigCNV.tsv', sep='\t')
    return cnvs


def test_addDerivedFeatures(global_data):
    cnvs = dataPreparation.addDerivedFeatures(global_data)
    assert cnvs["Score_SNP"].mean() == 1.449798389612388
    assert cnvs["DENSITY"].mean() == 0.0007772926241674732

def test_addCallRateToDataset(global_data):
    # TODO
    #cnvs = dataPreparation.addCallRateToDataset(global_data, "data/callrates")
    assert global_data['CallRate'].mean() == 0.9977302081136972

def test_addNbProbeByTech(global_data):
    # TODO
    assert global_data["Nb_Probe_tech"].count() == global_data.shape[0]
    assert global_data.at[0, "Nb_Probe_tech"] == 733256

def test_addChromosomicAnnotation(global_data):
    # TODO
    assert global_data.shape[0] == 7713

def test_transformTwoAlgsFeatures(global_data):
    cnvs = dataPreparation.transformTwoAlgsFeatures(global_data)
    assert cnvs.TwoAlgs.max() == 100
    assert cnvs.TwoAlgs.min() == 0
    