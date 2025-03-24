import os


ROOT_FOLDER = "/usr/local/scratch/metaGOflow-COMPLETED-results/Batch1and2/HCMR-data/data"
SUBFOLDERS = ["220223_JARVIS_HWLTKDRXY",
              "220622_JARVIS_HMGW5DSX3",
              "220712_JARVIS_HMNJKDSX3",
]

for subfolder in SUBFOLDERS:
    subfolder_path = os.path.join(ROOT_FOLDER, subfolder)
    # print(subfolder_path)
    # list folders in subfolder
    folders = os.listdir(subfolder_path)
    print("####### 1. layer folders ######")
    print(folders)

    folders_to_keep = []
    for folder in os.path.abspath(folders):
        if "DBH" or "DBB" in folder:
            folders_to_keep.append(folder)
            continue
        if "results" not in folder:
            continue

        # this is already the archive
        parent_folder = os.path.join(subfolder_path, folder).split("/")[-1].split("_")[-1] 
        print(parent_folder)

    #next layer
    print("####### 2. layer folders ######")
    print(folders_to_keep)
    for folder in os.path.abspath(folders_to_keep):
        if "DBH" or "DBB" in folder:
            folders_to_keep.append(folder)
            continue
        if "results" not in folder:
            continue

        # this is already the archive
        parent_folder = os.path.join(subfolder_path, folder).split("/")[-1].split("_")[-1] 
        print(parent_folder)
