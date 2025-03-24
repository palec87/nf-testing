import os


ROOT_FOLDER = "/usr/local/scratch/metaGOflow-COMPLETED-results/Batch1and2/HCMR-data/data"
SUBFOLDERS = ["220223_JARVIS_HWLTKDRXY",
              "220622_JARVIS_HMGW5DSX3",
              "220712_JARVIS_HMNJKDSX3",
]

# for subfolder in SUBFOLDERS:
#     subfolder_path = os.path.join(ROOT_FOLDER, subfolder)
#     # print(subfolder_path)
#     # list folders in subfolder
#     folders = os.listdir(subfolder_path)
#     print("####### 1. layer folders ######")
#     print('subfolder', subfolder)
#     print("###############################")

#     folders_to_keep = []
#     for folder in folders:
#         if "DBH" in str(folder) or "DBB" in str(folder):
#             folders_to_keep.append(folder)
#             continue
#         if "results" not in (folder):
#             continue

#         # this is already the archive

#         parent_folder = subfolder_path.split("/")[-1].split("_")[-1] 
#         print('parent folder', parent_folder, 'of', folder)

#     #next layer
#     print("####### 2. layer folders ######")
#     folders_to_keep2 = []
#     print(folders_to_keep)
#     for folder in folders_to_keep:
#         print(os.listdir(os.path.join(subfolder_path, folder)))
#         if "DBH" in str(folder) or "DBB" in str(folder) and ".zip" not in str(folder):
#             folders_to_keep2.append(folder)
#             continue
#         if "results" not in str(folder):
#             continue

#         # this is already the archive
#         parent_folder = folder.split("/")[-1].split("_")[-1]
#         print('parent folder', parent_folder, 'of', folder)

#     if len(folders_to_keep2) > 0:
#         print("Folders to keep: ", folders_to_keep2)
#         raise ValueError(f"There are still folders to keep")
    

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
            for subdir in dirs:
                results_paths.append(os.path.abspath(os.path.join(root, subdir)))
    return results_paths

# Example usage
if __name__ == "__main__":
    for folder in SUBFOLDERS:
        base_directory = os.path.join(ROOT_FOLDER, folder)
        paths = extract_results_paths(base_directory)
        for path in paths:
            print(path)