import pandas as pd
import numpy as np
import pyreadr as pr
import os
from copy import deepcopy

BACKUP = True

def load_data() -> pd.DataFrame:
    if not os.path.exists("backup/productive_backup.csv"):
        if not os.path.exists("Databaze/table_productive.csv"):
            return pr.read_r("Databaze/table_productive.RData")["table_productive"] # load data from RData file
        return pd.DataFrame(pd.read_csv("Databaze/table_productive.csv"))
    
    return pd.DataFrame(pd.read_csv("backup/productive_backup.csv"))

def make_dataset(df, group_by = ["patient_id", "clonotype_sequence"],drop = None) -> pd.DataFrame:
    working_df = deepcopy(df) # copy the dataframe to avoid changing the original one
    if drop is not None: 
        working_df.drop(columns=drop, inplace=True) # drop columns that are not needed for ratio calculation

    working_df["all_seq_count"] =  1 
    working_df["specific_seq_count"] = 1

    working_df["all_seq_count"] = working_df.groupby(["patient_id"]).all_seq_count.transform("size") # count all sequences for each patient
    working_df["specific_seq_count"] = working_df.groupby(group_by).specific_seq_count.transform("size") # count specific sequences for each group

    working_df["ratio"] = working_df["specific_seq_count"] / working_df["all_seq_count"] # calculate ratio of specific sequences to all sequences

    # filter out sequences that do not start with C and end with F
    working_df = working_df[working_df["clonotype_sequence"].str.startswith("C", na=False) & working_df["clonotype_sequence"].str.endswith("F", na=False)] 
    
    # drop duplicates - sequences that are the same. We keep specific_seq_count and for ratio calculation, so we do not need duplicates
    working_df.drop_duplicates(inplace=True)
    working_df.drop(columns=["all_seq_count"], inplace=True)

    working_df.sort_values(by=["patient_id", "specific_seq_count"], inplace=True)

    return working_df

def export_dataframe(df, name) -> None:
    path = "TCR_DATASETS/"
    if not os.path.exists(path):
        os.makedirs(path)
    name = "TCR_sequencing_hedimed_DATASET_" + name + ".csv"
    df.to_csv(path + name, index=False)
    print(name + " exported")

#_______________________________________________________________________________________________________________________
data_df: pd.DataFrame = pd.DataFrame(pd.read_csv("Databaze/table_productive.csv"))
print("Database loaded. Shape: " + str(data_df.shape))
in_indexes = [496130, 533406,14096288,14468978, 15161098, 15452605, 16920138,18676245, 26100297 ]
print(data_df.loc[in_indexes])

replacements: list[str] = ["VJ:Vb-(Db)-Jb", "VJ:Va-Ja", "VJ:VJ:Vh-(Dh)-Jh", "VJ:Vl-Jl", "VJ:Vk-Jk", "VJ:Va-Jd", "VJ:Vh-(Dh)-Jh", "VJ:Vg-Jg"]

for replacement in replacements: 
    data_df["clonotype"] = data_df["clonotype"].str.replace(replacement, "") # remove types of regions from clonotype
# before: VJ:Vb-(Db)-Jb  V10-1 -0/14/-5 J2-7  CASSESAGVTHEQYF
# after: V10-1 -0/14/-5 J2-7  CASSESAGVTHEQYF
print("Types of regions removed")


data_df["clonotype"] = data_df.apply(lambda x: x["clonotype"].replace(x["clonotype_sequence"], ""), axis=1) # remove clonotype_sequence from clonotype
# before: V10-1 -0/14/-5 J2-7  CASSESAGVTHEQYF
# after: V10-1 -0/14/-5 J2-7
print("Clonotype sequences removed")

data_df["clonotype"] = data_df["clonotype"].str.replace("~", "") # remove ~ from clonotype
print("~ removed")
data_df["clonotype"] = data_df["clonotype"].str.lstrip() # remove leading whitespaces
print("Leading whitespaces removed")
data_df["clonotype"] = data_df["clonotype"].str.rstrip() # remove trailing whitespaces
print("Trailing whitespaces removed")
print("Clonotypes cleaned")
print("-"*50)

"""
Fix invalid clonotypes. There are rows: 496130, 533406, 14096288, 14468978, 15161098, 15452605, 16920138, 18676245, 26100297 that have invalid clonoty_sequence:
                    sample                                      clonotype                   clonotype_sequence
496130         1A7_TRB-VJ_S7_R1_001            VJ:Vb-(Db)-Jb  V5-6 -18/3/-22 J2-3  V                  V
533406         1A8_TRB-VJ_S8_R1_001             VJ:Vb-(Db)-Jb  V27 -34/0/-16 J2-7  V                  V
14096288   3G3-TRB-VJ_S89_merged_R1  VJ:Vb-(Db)-Jb  V5-5=V5-6=V5-7 -40/3/-26 J1-5  V                  V
14468978   3G8-TRB-VJ_S94_merged_R1     VJ:Vb-(Db)-Jb  V11-1=V11-3 -18/0/-16 J2-7  V                  V
15161098  3H7-TRB-VJ_S102_merged_R1             VJ:Vb-(Db)-Jb  V28 -18/0/-16 J2-7  V                  V
15452605  3I1-TRB-VJ_S105_merged_R1  VJ:Vb-(Db)-Jb  V5-5=V5-6=V5-7 -40/3/-26 J1-5  V                  V
16920138  4B3-TRB-VJ_S124_merged_R1             VJ:Vb-(Db)-Jb  V18 -41/1/-20 J1-5  V                  V
18676245   4D6-TRB-VJ_S21_merged_R1            VJ:Vb-(Db)-Jb  V7-2 -23/1/-21 J1-4  V                  V
26100297  5E8-TRB-VJ_S111_merged_R1  VJ:Vb-(Db)-Jb  V5-5=V5-6=V5-7 -40/3/-26 J1-5  V                  V

V will be added to the beginning of the clonotype, because it was removed by cleaning clonotypes:
data_df["clonotype"] = data_df.apply(lambda x: x["clonotype"].replace(x["clonotype_sequence"], ""), axis=1)
before: V5-6 -18/3/-22 J2-3 V     V
after: 5-6 -18/3/-22 J2-3

Sample will be removed no metter what, because clonoype_sequence does not start with C and end with F, but will be count for ratio calculation
"""

invalid_indices = [] 
for index, row in enumerate(data_df["clonotype"]):
    if not row.startswith("V") or not row[1].isdigit():
        data_df.at[index, "clonotype"] = "V" + data_df.at[index, "clonotype"]

print("Invalid clonotypes fixed")


data_df["sample"] = data_df["sample"].str.replace("HEDIMED-", "") # remove HEDIMED- from sample
data_df["sample"] = data_df["sample"].str[:3] # take first 3 characters from sample
data_df["sample"] = data_df["sample"].astype(str) 
data_df.rename(columns={"sample": "patient_id"}, inplace=True) # rename sample to patient_id
# before: HEDIMED-1A7_TRB-VJ_S7_R1_001
# after: 1A7
print("Samples fixed")

# extract V and J regions from clonotype and add them to the dataframe as separate columns
data_df["V_region"] = data_df["clonotype"].astype(str).str.split().str[0]
data_df["J_region"] = data_df["clonotype"].astype(str).str.split().str[2]
data_df= data_df.drop(columns=["clonotype"])
print("Regions extracted")


# add category to the dataframe
vial_code_data = pd.DataFrame(pd.read_excel('Databaze/HEDIMED_kodovani.xlsx', usecols= ['category', 'Vial code'], dtype = {'category': "string", 'Vial code': "string"}))


vial_code_data.rename(columns={"Vial code": "patient_id"}, inplace=True) # rename Vial code to patient_id for merging with data_df on patient_id
working_df = pd.merge(data_df, vial_code_data, on="patient_id", how="left") # merge data_df with vial_code_data on patient_id
data_df["category"] = data_df["category"].fillna("Unknown") # fill NaN values in category with Unknown (Never occured in the dataset)
print("Vial codes added")

if BACKUP:
    data_df.to_csv("backup/productive_backup.csv", index=False) # backup created. Load_data() will load this backup if it exists

print("Backup created")
print("Datframe ready" + str(data_df.shape))

# dataset_info = [group_by, drop, name]
VJ_INDEPENDENT_INFO = [["patient_id", "clonotype_sequence"], ["V_region", "J_region"], "VJ_independent"]
VJ_DEPENDENT_INFO = [["patient_id", "clonotype_sequence", "V_region", "J_region"], None, "VJ_dependent"]
V_DEPENDENT_INFO = [["patient_id", "clonotype_sequence", "V_region"], ["J_region"], "V_dependent"]
J_DEPENDENT_INFO = [["patient_id", "clonotype_sequence", "J_region"], ["V_region"], "J_dependent"]

DATASETS_INFO = [VJ_INDEPENDENT_INFO, VJ_DEPENDENT_INFO, V_DEPENDENT_INFO, J_DEPENDENT_INFO]
for dataset in DATASETS_INFO:
    export_dataframe(make_dataset(data_df, group_by=dataset[0], drop=dataset[1]), dataset[2])