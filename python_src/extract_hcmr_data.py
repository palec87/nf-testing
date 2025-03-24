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
    print('subfolders', folders)

    folders_to_keep = []
    for folder in folders:
        if "DBH" in str(folder) or "DBB" in str(folder):
            folders_to_keep.append(folder)
            continue
        if "results" not in (folder):
            continue

        # this is already the archive

        parent_folder = subfolder_path.split("/")[-1].split("_")[-1] 
        print('parent folder', parent_folder, 'of', folder)

    #next layer
    print("####### 2. layer folders ######")
    folders_to_keep2 = []
    print(folders_to_keep)
    for folder in folders_to_keep:
        if "DBH" in str(folder) or "DBB" in str(folder):
            folders_to_keep2.append(folder)
            continue
        if "results" not in str(folder):
            continue

        # this is already the archive
        parent_folder = os.path.join(subfolder_path, folder).split("/")[-1].split("_")[-1] 
        print('parent folder', parent_folder, 'of', folder)
