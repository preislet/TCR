import pandas as pd
import numpy as np      
import warnings

import plotly.express as px
import plotly.graph_objects as go

with warnings.catch_warnings():
    warnings.simplefilter("ignore")

from TCR_DATASETS.TCR_load import load_tcr_data

from typing import Literal
from itertools import product
import os
import tqdm
import pickle

from deepdiff import DeepDiff

def make_category_dataframes(region_column: str, data_df: pd.DataFrame, type: str) -> dict:
    """
    Create dataframes for each category in the dataset based on a specific region column and type of analysis.
    The resulting dataframes are either loaded from a pickle file or generated and saved to a pickle file.

    Args:
    - region_column (str): The column name for regions in the dataset.
    - data_df (pd.DataFrame): The input dataframe containing TCR data.
    - type (str): The type of analysis to perform ('frequency' or 'sum').

    Returns:
    - dict: A dictionary of dataframes for each category.
    """

    # Check if the category dataframes already exist in a pickle file
    if os.path.exists(f"backup/category_dfs_{region_column}_{type}.pkl"):
        with open(f"backup/category_dfs_{region_column}_{type}.pkl", "rb") as f:
            category_dfs: dict = pickle.load(f)
            
    else:
        # Generate the category dataframes if not found in the pickle file
        categories: np.ndarray = data_df["category"].unique() # Get all categories in the dataset
        category_dfs: dict = {}
        for category in categories:
            category_df: pd.DataFrame = data_df[data_df["category"] == category] # Get the dataframe for the current category
            category_regions: np.ndarray = category_df[region_column].unique() # Get all regions in the current category
            category_regions_freq: list = []
            for region in category_regions:
                region_df: pd.DataFrame = category_df[category_df[region_column] == region] # Get the dataframe for the current region

                if type == "frequency": region_frequence_for_all_sequences = np.sum(region_df["ratio"]) # Calculate the frequency of the region for all sequences
                elif type == "sum": region_frequence_for_all_sequences = np.sum(region_df["specific_seq_count"]) # Calculate the sum of the region for all sequences

                category_regions_freq.append(region_frequence_for_all_sequences)

            category_dfs[category] = pd.DataFrame({region_column: category_regions, "output": category_regions_freq})
        
        # Save the generated dataframes to a pickle file
        with open(f"backup/category_dfs_{region_column}_{type}.pkl", "wb") as f:
            pickle.dump(category_dfs, f)

    return category_dfs


def make_full_plot(show, save, frequecy_dfs, categories, region_column, top, type, add_annotations_to_unique_regions = True):
    """
    Create a full plot for all categories combined, displaying the top regions based on the specified type of analysis.

    Args:
    - show (bool): Whether to display the plot.
    - save (bool): Whether to save the plot.
    - frequecy_dfs (dict): Dictionary of frequency dataframes for each category.
    - categoryes (list): List of categories.
    - region_column (str): The column name for regions in the dataset.
    - top (int): Number of top regions to display.
    - type (str): The type of analysis to perform ('frequency' or 'sum').
    """
    fig = go.Figure()
    for category in categories:
        fig.add_trace(go.Bar(x=frequecy_dfs[category][region_column], y=frequecy_dfs[category]["output"],name=category))


    path = f"Analysis/{type}_TOP_{top}_{region_column}"
    fig.update_layout(title=f"{top} {region_column}s for each category (U = unique region) ", 
                    xaxis_title=f"{region_column}s", 
                    yaxis_title=f"{type} of regions",
                    legend_title="Category", height=600, width=1800)
    

    if add_annotations_to_unique_regions:
        for category in categories:
            for i, region in enumerate(frequecy_dfs[category][region_column]):
                if frequecy_dfs[category].iloc[i]["unique"]:
                    fig.add_annotation(x=region, y=frequecy_dfs[category].iloc[i]["output"], text="U", showarrow=True, arrowhead=1, arrowwidth=2,
                                    arrowcolor="black", ax=0, ay=-40, yshift=10, opacity=0.8, font=dict(color="red"),)



    if save: fig.write_image(f"{path}/{region_column}_frequency_ALL_GROUP_{top}.png")
    if show: fig.show()
    fig.update_layout(barmode='stack')
    if save: fig.write_image(f"{path}/{region_column}_frequency_ALL_STACK_{top}.png")
    if show: fig.show()

def make_plot_for_each_category(show, save, frequecy_dfs, categories, region_column, top, type):
    """
    Create individual plots for each category, displaying the top regions based on the specified type of analysis.

    Args:
    - show (bool): Whether to display the plots.
    - save (bool): Whether to save the plots.
    - frequecy_dfs (dict): Dictionary of frequency dataframes for each category.
    - categories (list): List of categories.
    - region_column (str): The column name for regions in the dataset.
    - top (int): Number of top regions to display.
    - type (str): The type of analysis to perform ('frequency' or 'sum').
    """
    path = f"Analysis/{type}_TOP_{top}_{region_column}"
    for category in categories:
        fig = go.Figure()
        fig.add_trace(go.Bar(x=frequecy_dfs[category][region_column], y=frequecy_dfs[category]["output"], name=category,
                        marker_color=frequecy_dfs[category]["unique"].apply(lambda x: "red" if x else "blue")))
        
        fig.update_layout(
            title=f"{top} {region_column}s for {category} (U = unique region)",
            xaxis_title=f"{region_column}s",
            yaxis_title=f"{type} of regions",
            legend_title="Category",
            height=600,
            width=1200
        )

        # add annotations to unique regions dont overlap with other arows
        for i, region in enumerate(frequecy_dfs[category][region_column]):
            if frequecy_dfs[category]["unique"].iloc[i]:
                fig.add_annotation( x=region, y=frequecy_dfs[category]["output"].iloc[i], text="U", showarrow=True, arrowhead=1, yshift=10,
                    font=dict(family="Courier New, monospace", size=16,color="red"))
                
        if show: fig.show()
        if save: fig.write_image(f"{path}/{region_column}_frequency_{category}_{top}.png")

def Make_analysis(data_df: pd.DataFrame, region_column, type_of_analysis, save_plot = False, show_plot = False, regions_count = 10):
    """
    Perform analysis on the provided TCR data, creating and saving plots for the specified region column and type of analysis.

    Args:
    - data_df (pd.DataFrame): The input dataframe containing TCR data.
    - region_column (str): The column name for regions in the dataset.
    - type_of_analysis (str): The type of analysis to perform ('frequency' or 'sum').
    - save_plot (bool): Whether to save the plots.
    - show_plot (bool): Whether to display the plots.
    - regions_count (int): Number of top regions to include in the analysis.
    """     
    categories: np.ndarray = data_df["category"].unique()

    # Count the number of patients in each category
    patients_count = [len(data_df[data_df["category"] == category]["patient_id"].unique()) for category in categories]

    # Count the number of sequences in each category can be used for normalization too
    # patients_count = [len(data_df[data_df["category"] == category]) for category in categories]

    patients_counts_dict = dict(zip(data_df["category"].unique(), patients_count))
    category_dfs = make_category_dataframes(region_column=region_column, data_df=data_df, type=type_of_analysis)
    """
    code: print(category_dfs[categories[0]].head(10))
    Example of category_dfs[categories[0]].head(10) output:
            VJ_region    output
        0  V10-1_J2-7  0.048718
        1  V10-1_J1-1  0.113048
        2  V10-1_J1-4  0.004329
        3  V10-1_J1-5  0.019606
        4  V10-1_J1-2  0.019730
        5  V10-1_J1-6  0.009741
        6  V10-1_J1-3  0.007991
        7  V10-1_J2-3  0.001031
        8  V10-1_J2-2  0.001592
        9  V10-2_J1-4  0.008143
    """
    
    # Normalize the output by the number of patients in each category
    for category in categories:
        category_dfs[category] = category_dfs[category].sort_values(by="output", ascending=False).head(regions_count)
        category_dfs[category]["output"] /= patients_counts_dict[category]

    # Identify names of the regions for each category
    regions = { category: category_dfs[category][region_column].unique() for category in categories }
    # Identify names of unique regions for each category
    unique_regions = {}
    for category in categories:
        # Get all regions from other categories
        other_regions = np.concatenate([regions[other_category] for other_category in categories if category != other_category])
        unique_regions[category] = np.setdiff1d(regions[category], other_regions) # Get unique regions for the current category



    # Mark unique regions in the dataframes
    for category in categories:
        category_dfs[category]["unique"] = category_dfs[category][region_column].apply(lambda x: x in unique_regions[category]) # Mark unique regions in the dataframes


    # Generate and save/display plots
    make_full_plot(show_plot, save_plot, category_dfs, categories, region_column, regions_count, type_of_analysis)
    make_plot_for_each_category(show_plot, save_plot, category_dfs, categories, region_column, regions_count, type_of_analysis)

def make_dirs():
    # make dir for each combination of analysis TYPE_OF_ANALYSIS, REGIONS_COUNT_LIST, REGIONS_LIST
    if not os.path.exists("Analysis"):
            os.makedirs("Analysis")
    
    for analysis, count, regions_name in product(TYPE_OF_ANALYSIS, REGIONS_COUNT_LIST, REGIONS_LIST_NAMES):
        # Construct the directory name based on the current combination
        dir_name = f"Analysis/{analysis}_TOP_{count}_{regions_name}"
        
        # Create the directory if it doesn't exist
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

    

TYPE_OF_ANALYSIS = ["frequency", "sum"]
REGIONS_COUNT_LIST = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
REGIONS_LIST_NAMES = ["VJ_region", "V_region", "J_region"]
SAVE_PLOTS = True
SHOW_PLOTS = False

TCR_DATASETS = {
    "J_region": "J_DEPENDENT",
    "VJ_region": "VJ_DEPENDENT",
    "V_region": "V_DEPENDENT",
}

if __name__ == "__main__":
    make_dirs()
    for analysis, count, regions_name in product(TYPE_OF_ANALYSIS, REGIONS_COUNT_LIST, REGIONS_LIST_NAMES):
        data_df = load_tcr_data(TCR_DATASETS[regions_name])
        if regions_name == "VJ_region": data_df["VJ_region"] = data_df["V_region"] + "_" + data_df["J_region"]
        Make_analysis(data_df, regions_name, analysis, save_plot = SAVE_PLOTS,show_plot = SHOW_PLOTS, regions_count = count)



