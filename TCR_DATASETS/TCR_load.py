import pandas as pd
from typing import Literal

Integrity: bool = False

if not Integrity:
    J_DEPENDENT = 'TCR_DATASETS/TCR_sequencing_hedimed_DATASET_J_dependent.csv'
    V_DEPENDENT = 'TCR_DATASETS/TCR_sequencing_hedimed_DATASET_V_dependent.csv'
    VJ_DEPENDENT = 'TCR_DATASETS/TCR_sequencing_hedimed_DATASET_VJ_dependent.csv'
    VJ_INDEPENDENT = 'TCR_DATASETS/TCR_sequencing_hedimed_DATASET_VJ_independent.csv'


if Integrity:
    J_DEPENDENT = 'TCR_sequencing_hedimed_DATASET_J_dependent.csv'
    V_DEPENDENT = 'TCR_sequencing_hedimed_DATASET_V_dependent.csv'
    VJ_DEPENDENT = 'TCR_sequencing_hedimed_DATASET_VJ_dependent.csv'
    VJ_INDEPENDENT = 'TCR_sequencing_hedimed_DATASET_VJ_inependent.csv'

Base_columns = {
    "patient_id": str,
    "clonotype_sequence": str,
    "specific_seq_count": int,
    "ratio": float,
    "category": str
}

J_DEPENDENT_COLUMNS = Base_columns.copy()
J_DEPENDENT_COLUMNS["J_region"] = str

V_DEPENDENT_COLUMNS = Base_columns.copy()
V_DEPENDENT_COLUMNS["V_region"] = str

VJ_DEPENDENT_COLUMNS = Base_columns.copy()
VJ_DEPENDENT_COLUMNS["V_region"] = str
VJ_DEPENDENT_COLUMNS["J_region"] = str

VJ_INDEPENDENT_COLUMNS = Base_columns.copy()


data_sources = {
        'J_DEPENDENT': (J_DEPENDENT, J_DEPENDENT_COLUMNS),
        'V_DEPENDENT': (V_DEPENDENT, V_DEPENDENT_COLUMNS),
        'VJ_DEPENDENT': (VJ_DEPENDENT, VJ_DEPENDENT_COLUMNS),
        'VJ_INDEPENDENT': (VJ_INDEPENDENT, VJ_INDEPENDENT_COLUMNS),
    }

def tcr_info() -> None:
    # Return TCR data information

    print("TCR_DATASETS/TCR_sequencing_hedimed_DATASET_ratio_J_dependent.csv - NAME: J_DEPENDENT")
    print("TCR_DATASETS/TCR_sequencing_hedimed_DATASET_ratio_V_dependent.csv - NAME: V_DEPENDENT")
    print("TCR_DATASETS/TCR_sequencing_hedimed_DATASET_ratio_VJ_dependent.csv - NAME: VJ_DEPENDENT")
    print("TCR_DATASETS/TCR_sequencing_hedimed_DATASET_ratio_VJ_independent.csv - NAME: VJ_INDEPENDENT")


# Load TCR data
def load_tcr_data(NAME: Literal['J_DEPENDENT', 'V_DEPENDENT', 'VJ_DEPENDENT', 'VJ_INDEPENDENT'], check_integrity = False) -> pd.DataFrame:
    if NAME in data_sources:
        file_path, columns = data_sources[NAME]
        tcr_data = pd.read_csv(file_path, dtype=columns, usecols=columns.keys())
    else:
        raise ValueError("Invalid TCR data name")
    
    return pd.DataFrame(tcr_data)
