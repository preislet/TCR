import pandas as pd
import numpy as np


def check_duplicates(df: pd.DataFrame) -> bool:
    # Check for duplicates in the dataframe
    return df.duplicated().any()

def check_missing_values(df: pd.DataFrame) -> bool:
    # Check for missing values in the dataframe
    return df.isnull().values.any()

def check_for_anomalies(df: pd.DataFrame) -> bool:
    # Check for anomalies in the dataframe
    if "Unknown" in df["category"].unique():
        return True
    return False

def check_V_region(df: pd.DataFrame) -> bool:
    # Check for anomalies in the V region
    if "V_region" not in df.columns:
        return False
    # Start with the V and second one is number
    return not df["V_region"].str.match(r"V\d+").all()

def check_J_region(df: pd.DataFrame) -> bool:
    # Check for anomalies in the J region
    if "J_region" not in df.columns:
        return False
    # Start with the J and second one is number
    return not df["J_region"].str.match(r"J\d+").all()

def check_integrity(df: pd.DataFrame) -> bool:
    print("Checking for duplicates..." if check_duplicates(df) else "No duplicates found.")
    print("Checking for missing values..." if check_missing_values(df) else "No missing values found.")
    print("Checking for anomalies..." if check_for_anomalies(df) else "No anomalies found.")
    print("Checking for anomalies in V region..." if check_V_region(df) else "No anomalies found in V region.")
    print("Checking for anomalies in J region..." if check_J_region(df) else "No anomalies found in J region.")
    print("-"*100)

def start_integrity_check() -> None:
    for name in TCR_load.data_sources:
        print(f"Loading {name} data...")
        df = TCR_load.load_tcr_data(name)
        check_integrity(df)


if __name__ == "__main__":
    import TCR_load
    TCR_load.Integrity = True
    start_integrity_check()

else:
    import TCR_DATASETS.TCR_load as TCR_load
    TCR_load.Integrity = True
