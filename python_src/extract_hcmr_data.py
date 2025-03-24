import os
import pandas as pd
import numpy as np


ROOT_FOLDER = "/usr/local/scratch/metaGOflow-COMPLETED-results/Batch1and2/HCMR-data/data"
SUBFOLDERS = ["220223_JARVIS_HWLTKDRXY",
              "220622_JARVIS_HMGW5DSX3",
              "220712_JARVIS_HMNJKDSX3",
]

OUT_PATH = "/usr/local/scratch/nf-metaGOflow/wf-test/nf-testing/reuslts-hcmr"

BATCH1_RUN_INFO_PATH = (
    "https://raw.githubusercontent.com/emo-bon/sequencing-data/main/shipment/"
    "batch-001/run-information-batch-001.csv"
)
BATCH2_RUN_INFO_PATH = (
    "https://raw.githubusercontent.com/emo-bon/sequencing-data/main/shipment/"
    "batch-002/run-information-batch-002.csv"
)
    

## write a function to extract wll absolute paths which contain directory /results
def extract_results_paths(base_dir):
    """
    Extract all absolute paths containing the directory '/results' 
    and include all its subdirectories.

    :param base_dir: The base directory to search in.
    :return: A list of absolute paths containing '/results' and their subdirectories.
    """
    results_paths = []
    for root, dirs, files in os.walk(base_dir):
        if '/results' in root:
            results_paths.append(os.path.abspath(root))
            # Add all subdirectories under the current '/results' directory
            # for subdir in dirs:
            #     results_paths.append(os.path.abspath(os.path.join(root, subdir)))
    return results_paths

# search_string = "DBB_AACDOSDA_4_HMGW5DSX3"
def find_sample(search_string):
    for i, batch in enumerate([BATCH1_RUN_INFO_PATH, BATCH2_RUN_INFO_PATH]):
        df = pd.read_csv(BATCH1_RUN_INFO_PATH)
        df1 = df[["reads_name", "ref_code", "source_mat_id"]]

        # replace nan with empty string
        df1 = df1.replace(np.nan, '', regex=True)
        ans = df1[df1['reads_name'].str.contains(search_string)]
        if not ans.empty:
            return ans
    return None


# Example usage
if __name__ == "__main__":

    for folder in SUBFOLDERS:
        base_directory = os.path.join(ROOT_FOLDER, folder)
        paths = extract_results_paths(base_directory)

        top_paths = []
        for path in paths:
            # match paths where there is no subfolder of results
            if path.split('/')[-1] == 'results':
                top_paths.append(path)

        # extract archive name from the parent folder
        for path in top_paths:
            reads_name = path.split("/")[-2].split("_")[-1]
            print('Reads name', reads_name, 'of', path)
            if '.' not in reads_name:
                print('Skippint invalid archive name')
                continue

            # search for the sample in the batch run information
            ans = find_sample(reads_name)
            print(ans)

            # create folder with the reads name
            out_folder = os.path.join(OUT_PATH, f"{ans['source_mat_id'].values[0]}-tables")
            os.makedirs(out_folder, exist_ok=True)
            print(out_folder)

            # move some files
            # move the file with the reads name
            os.system(f"cp {os.path.join(path, "functional-annotation", "DBB.merged.summary.go")} {out_folder}")