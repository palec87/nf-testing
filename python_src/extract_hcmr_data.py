import os
import pandas as pd
import numpy as np


ROOT_FOLDER = "/usr/local/scratch/metaGOflow-COMPLETED-results/Batch1and2/HCMR-data/data"
SUBFOLDERS = ["220223_JARVIS_HWLTKDRXY",
              "220622_JARVIS_HMGW5DSX3",
              "220712_JARVIS_HMNJKDSX3",
]

OUT_PATH = "/usr/local/scratch/nf-metaGOflow/wf-test/nf-testing/results-hcmr"

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
    missing_data = []
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
            # print(ans)

            # create folder with the reads name
            out_folder = os.path.join(OUT_PATH, f"{ans['source_mat_id'].values[0]}-tables")
            os.makedirs(out_folder, exist_ok=True)
            # print(out_folder)

            # move some files
            ret = os.system(f"cp {os.path.join(path, 'functional-annotation', 'DBB.merged.summary.go')} {out_folder}")
            if ret != 0:
                missing_data.append((reads_name, 'DBB.merged.summary.go'))

            ret = os.system(f"cp {os.path.join(path, 'functional-annotation', 'DBB.merged.summary.pfam')} {out_folder}")
            if ret != 0:
                missing_data.append((reads_name, 'DBB.merged.summary.pfam'))

            ret = os.system(f"cp {os.path.join(path, 'functional-annotation', 'DBB.merged.summary.go_slim')} {out_folder}")
            if ret != 0:
                missing_data.append((reads_name, 'DBB.merged.summary.go_slim'))

            ret = os.system(f"cp {os.path.join(path, 'functional-annotation', 'DBB.merged.summary.ips')} {out_folder}")
            if ret != 0:
                missing_data.append((reads_name, 'DBB.merged.summary.ips'))

            ret = os.system(f"cp {os.path.join(path, 'functional-annotation', 'DBB.merged.summary.ko')} {out_folder}")
            if ret != 0:
                missing_data.append((reads_name, 'DBB.merged.summary.ko'))

            # taxonomy
            ret = os.system(f"cp {os.path.join(path, 'taxonomy-summary', "LSU", "DBB.merged_LSU.fasta.mseq.tsv")} {out_folder}")
            if ret != 0:
                missing_data.append((reads_name, 'DBB.merged_LSU.fasta.mseq.tsv'))

            ret = os.system(f"cp {os.path.join(path, 'taxonomy-summary', "SSU", "DBB.merged_SSU.fasta.mseq.tsv")} {out_folder}")
            if ret != 0:
                missing_data.append((reads_name, 'DBB.merged_SSU.fasta.mseq.tsv'))



            print(missing_data)