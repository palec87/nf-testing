import os


ROOT_FOLDER = "/usr/local/scratch/metaGOflow-COMPLETED-results/Batch1and2/HCMR-data/data"
SUBFOLDERS = ["220223_JARVIS_HWLTKDRXY",
              "220622_JARVIS_HMGW5DSX3",
              "220712_JARVIS_HMNJKDSX3",
]

for subfolder in SUBFOLDERS:
    subfolder_path = os.path.join(ROOT_FOLDER, subfolder)
    print(subfolder_path)
    # list folders in subfolder
    folders = os.listdir(subfolder_path)
    for folder in folders:
        if "results" not in folder:
            continue


        folder_path = os.path.join(subfolder_path, folder)
        print(folder_path)
