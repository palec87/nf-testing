import io
import math
import os
import sys
import argparse
import textwrap
import urllib.request
import logging
from pathlib import Path
import pandas as pd
# from pyspark.sql import SparkSession

# spark = SparkSession.builder.getOrCreate()

desc = """
Combine tables in the specific extraction directory into a single tables
"""

logger = logging.getLogger(name="CombineTables")
logging.basicConfig(format="\t%(levelname)s: %(message)s", level=logging.INFO)
# sys.path.append(str(PROJECT_DIR / "src"))

# from minio import S3Error
# from minio import Minio


# with open("credentials.json") as f:
#     _creds = json.load(f)


# client = Minio(
#     "10.4.1.4:9000",
#     secure=False,
#     access_key=_creds["accessKey"],
#     secret_key=_creds["secretKey"],
# )


# Batch 1 has 80 records sent to squencing lab
# Batch 2 has 108 records sent to squencing lab
# BATCH1AND2_TOTAL = 188

# changed because blanks are removed
BATCH1AND2_TOTAL = 181
BATCH1_RUN_INFO_PATH = (
    "https://raw.githubusercontent.com/emo-bon/sequencing-data/main/shipment/"
    "batch-001/run-information-batch-001.csv"
)
BATCH2_RUN_INFO_PATH = (
    "https://raw.githubusercontent.com/emo-bon/sequencing-data/main/shipment/"
    "batch-002/run-information-batch-002.csv"
)
# Expected number of metaGOflow analyses in version 1 of the data release
EXPECTED_ANALYSES = 54
TAXONOMY_RANK_KEYS = {
    "sk": "superkingdom",
    "k": "kingdom",  # Note that prokaryotes do not have a kingdom entry
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}

VERSION = 2  # 11 Sept 2024


# this just for the minio uploading parquet files, I should do the same I guess
# TABLES = [
#     "metagoflow_analyses.SSU",
#     "metagoflow_analyses.LSU",
#     # "metagoflow_analyses.go_slim",
#     # "metagoflow_analyses.go"
#     # "metagoflow_analyses.ips"
#     # "metagoflow_analyses.ko",
#     # "metagoflow_analyses.pfam"
# ]

# The destination bucket and filename on the MinIO server
def extract_keys(data):
    d = {}
    for _, row in data.iterrows():
        try:
            # Not all entries are sequenced
            reads_name = row["reads_name"].split("_")[-1]
        except AttributeError as err:
            if not math.isnan(row["reads_name"]):
                raise err
        code = row["ref_code"]
        if code in d:
            raise ValueError(f"Duplicate reads_name: {reads_name}")
        prefix = row["ref_code_seq"].split("_")[0]

        d[reads_name] = (code, prefix, row['source_mat_id'])
    # print(f"Extracted {len(d)} records from batch sheets")
    return d

## method to merge tracker and shipments data
def merge_tracker_metadata(df_tracker):
    df = pd.read_csv(BATCH1_RUN_INFO_PATH)
    df1 = pd.read_csv(BATCH2_RUN_INFO_PATH)
    df = pd.concat([df, df1])
    print(len(df))

    # select only batch 1 and 2 from tracker data
    df_tracker = df_tracker[df_tracker['batch'].isin(['001', '002'])]

    # merge with tracker data on ref_code
    df = pd.merge(df, df_tracker, on='ref_code', how='inner')

    # split the sequenced tables
    df_discarded = df[df['lib_reads_seqd'].isnull()]
    df1 = df[~df['lib_reads_seqd'].isnull()]

    # deliver only filters and sediment samples, no blanks
    df_to_deliver = df1[df1['sample_type'].isin(['filters', 'sediment'])].reset_index(drop=True)
    #drop columns
    df_to_deliver.drop(columns='reads_name_y', inplace=True)
    df_to_deliver.rename(columns={'reads_name_x': 'reads_name'}, inplace=True)
    return df_to_deliver, df_discarded

# Uses code_keys dictionary from above
def parse_taxonomy(taxonomy):
    """
    12928>--1.0>sk__Bacteria;k__;p__Firmicutes;c__Erysipelotrichia;o__Erysipelotrichales;
    f__Erysipelotrichaceae;g__Turicibacter;s__Turicibacter_sanguinis>--154288
    [ebi_otu_id, abundance, (super_kingdon, phylum, class, order, family,
    geneus, species), ncbi_tax_id]
    """
    classification = {}
    for level in taxonomy.split(";"):
        rank, name = level.split("__")
        classification[TAXONOMY_RANK_KEYS[rank]] = name
    return classification


def parse_local_inventory(inv: str, code_keys: dict[tuple[str, str]], folders: list[Path] = None):
    count = 0
    all_objs_data = []  # list of dicts, each a taxonomic entry

    for key, val_tuple in code_keys.items():
        all_sample_data = []
        # For each of the 54 LSU inventories
        prefix = val_tuple[1]
        fn = f"{prefix}.merged_{inv}.fasta.mseq.tsv"
        fps = [os.path.join(folder, f"{val_tuple[2]}-tables", fn)
              for folder in folders]
        # chekc if the file exists in any of the folders
        for f in fps:
            if os.path.exists(f):
                fp = f
                break
        else:
            print(f"{val_tuple[2]}-tables not found in any of the folders, look for {key, val_tuple}")
            continue

        # fp = os.path.join(folder, f"{val_tuple[2]}-tables", fn)

        ############ this will be eventually removed ################
        # lst = ["EMOBON_OSD74_Wa_21-tables", "EMOBON_VB_Wa_43-tables"]
        # if any(x in fp for x in lst):
        #     logger.info(f'FULL path, {fp}, {val_tuple[2]}')
        #############################################################
    
        try:
            csv_data = pd.read_csv(fp, sep="\t", skiprows=1)
        except FileNotFoundError as e:
            # print(e)
            continue

        for _, row in csv_data.iterrows():
            data = {}  # For one taxonomic entry
            # For each row in the inventory
            taxonomy = parse_taxonomy(row["taxonomy"])
            data["ref_code"] = val_tuple[0]
            data["ncbi_tax_id"] = row.get("taxid")  # NCBI taxid
            data["abundance"] = row.get(f"{inv}_rRNA")  # Abundance
            data["superkingdom"] = taxonomy.get("superkingdom")
            data["kingdom"] = taxonomy.get("kingdom", None)
            data["phylum"] = taxonomy.get("phylum", None)
            data["class"] = taxonomy.get("class", None)
            data["order"] = taxonomy.get("order", None)
            data["family"] = taxonomy.get("family", None)
            data["genus"] = taxonomy.get("genus", None)
            data["species"] = taxonomy.get("species", None)
            all_sample_data.append(data)

        all_objs_data.extend(all_sample_data)

        count += 1
    logger.info(f"Found {count} inventories from {inv}")
    # if count != EXPECTED_ANALYSES:
    #     raise ValueError(f"Could not find all {EXPECTED_ANALYSES} records in v1")
    return all_objs_data


def parse_other_tax_tables(inv, code_keys, folders: list[Path] = None):
    """
    Parse IPS/KEGG/pfam summary files
    """
    count = 0
    all_objs_data = []  # list of dicts, each a taxonomic entry
    for _, val_tuple in code_keys.items():
        all_sample_data = []

        prefix = val_tuple[1]           # like DBB etc
        fn = f"{prefix}.merged.summary.{inv}"
        fps = [os.path.join(folder, f"{val_tuple[2]}-tables", fn)
              for folder in folders]
        # chekc if the file exists in any of the folders
        for f in fps:
            if os.path.exists(f):
                fp = f
                break
        else:
            print(f"{val_tuple[2]}-tables not found in any of the folders")
            continue

        # fp = os.path.join(folder, f"{val_tuple[2]}-tables", fn)  # this is the
        try:
            csv_data = pd.read_csv(fp, header=None)
        except FileNotFoundError as e:
            continue

        if inv == "ips":
            keys = ["ref_code", "accession", "description", "abundance"]
        else:
            keys = ["ref_code", "entry", "name", "abundance"]

        for _, row in csv_data.iterrows():
            data = {}
            data[keys[0]] = val_tuple[0]
            data[keys[1]] = row[1]
            data[keys[2]] = row[2]
            data[keys[3]] = int(row[0])
            all_sample_data.append(data)
            #     (val_tuple[0],              # ref_code
            #      row[1],                    # accession for IPS, or KEGG/PFAM entry
            #      row[2],                    # description for IPS, "name" for KEGG/PFAM
            #      int(row[0])                # abundance
            #     )
            # )
        all_objs_data.extend(all_sample_data)
        count += 1
    print(f"Found {count} samples for {inv} summary")
    return all_objs_data


def go_tables(inv, code_keys, folders: list[Path] = None):
    """
    Parse GO summary files
    """
    count = 0
    all_objs_data = []  # list of dicts, each a taxonomic entry
    for _, val_tuple in code_keys.items():
        all_sample_data = []

        prefix = val_tuple[1]           # like DBB etc
        fn = f"{prefix}.merged.summary.{inv}"
        # fp = os.path.join(folder, f"{val_tuple[2]}-tables", fn)  # this is the final path
        fps = [os.path.join(folder, f"{val_tuple[2]}-tables", fn)
              for folder in folders]
        # chekc if the file exists in any of the folders
        for f in fps:
            if os.path.exists(f):
                fp = f
                break
        else:
            print(f"{val_tuple[2]}-tables not found in any of the folders")
            continue
        try:
            csv_data = pd.read_csv(fp, header=None)
        except FileNotFoundError as e:
            continue

        for _, row in csv_data.iterrows():
            data = {}
            data['ref_code'] = val_tuple[0]
            data['id'] = row[0]
            data['name'] = row[1]
            data['aspect'] = row[2]
            data['abundance'] = row[3]
            all_sample_data.append(data)
        all_objs_data.extend(all_sample_data)
        count += 1
    print(f"Found {count} samples for {inv}")
    return all_objs_data


def main(project_dirs, out_dir=None):
    # batch1 = "https://raw.githubusercontent.com/emo-bon/sequencing-data/main/shipment/batch-001/run-information-batch-001.csv"
    # batch2 = "https://raw.githubusercontent.com/emo-bon/sequencing-data/main/shipment/batch-002/run-information-batch-002.csv"
    url_tracker = "https://raw.githubusercontent.com/emo-bon/momics-demos/refs/heads/main/wf0_landing_page/emobon_sequencing_master.csv"
    df_tracker = pd.read_csv(url_tracker ,index_col=0)
    df_to_deliver, df_discarded = merge_tracker_metadata(df_tracker)
    if out_dir is None:
        OUT_PATH = os.path.join(os.getcwd(), "combined_tables")
    else:
        OUT_PATH = os.path.join(out_dir, "combined_tables")

    logger.info(f"OUT_PATH: {OUT_PATH}")
    # tables_folder = os.path.join(project_dir, "results-tables")

    logger.info("creating code_keys")
    code_keys = {}
    # for batch in (batch1, batch2):
    #     with urllib.request.urlopen(batch) as f:  # noqa: S310
    #         data = pd.read_csv(f)  # CSV
    #         code_keys.update(extract_keys(data))

    code_keys.update(extract_keys(df_to_deliver))

    if len(code_keys) != BATCH1AND2_TOTAL:
        raise ValueError(f"Expected {BATCH1AND2_TOTAL} keys, got {len(code_keys)}")
    else:
        logger.info(f"Extracted the expected {len(code_keys)} records from batch sheets")

    logger.info(f"I have following number of keys: {len(code_keys)}")
    all_data = {}
    LSU_data = parse_local_inventory("LSU", code_keys, folders=project_dirs)
    SSU_data = parse_local_inventory("SSU", code_keys, folders=project_dirs)
    go_data = go_tables("go", code_keys, folders=project_dirs)
    go_slim_data = go_tables("go_slim", code_keys, folders=project_dirs)
    ips_data = parse_other_tax_tables("ips", code_keys, folders=project_dirs)
    ko_data = parse_other_tax_tables("ko", code_keys, folders=project_dirs)
    pfam_data = parse_other_tax_tables("pfam", code_keys, folders=project_dirs)
    logger.info(f"Parsed {len(LSU_data)} rows from LSU data")
    logger.info(f"Parsed {len(SSU_data)} rows from SSU data")
    logger.info(f"Parsed {len(go_data)} rows from GO data")
    logger.info(f"Parsed {len(go_slim_data)} rows from GO slim data")
    logger.info(f"Parsed {len(ips_data)} rows from IPS data")
    logger.info(f"Parsed {len(ko_data)} rows from KEGG data")
    logger.info(f"Parsed {len(pfam_data)} rows from PFAM data")
    # lsu_df = pd.DataFrame.from_records(LSU_data)
    # ssu_df = pd.DataFrame.from_records(SSU_data)
    # go_df = pd.DataFrame.from_records(go_data)
    # go_slim_df = pd.DataFrame.from_records(go_slim_data)

    all_data["metagoflow_analyses.LSU"] = LSU_data
    all_data["metagoflow_analyses.SSU"] = SSU_data
    all_data["metagoflow_analyses.go"] = go_data
    all_data["metagoflow_analyses.go_slim"] = go_slim_data
    all_data["metagoflow_analyses.ips"] = ips_data
    all_data["metagoflow_analyses.ko"] = ko_data
    all_data["metagoflow_analyses.pfam"] = pfam_data
    # create the output directory if it does not exist
    if not os.path.exists(OUT_PATH):
        os.makedirs(OUT_PATH)
    logger.info("SAVING PATh for the combined table.")
    logger.info(f"OUT_PATH: {OUT_PATH}")
    # logger.info(f"tables_folder: {project_dir}")
    for path, table in all_data.items():
        # logger.info("SAVING PATh for the combined table.")
        # logger.info(os.path.join(OUT_PATH, f"{path}.csv"))

        table_df = pd.DataFrame.from_records(table)
        table_df.to_csv(os.path.join(OUT_PATH, f"{path}.csv"), index=False)

        # parquet files here
        table_df.to_parquet(
            os.path.join(OUT_PATH, f"{path}.parquet"),
            engine="pyarrow",
            compression="snappy",
            index=False,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(desc),
    )
    parser.add_argument(
        '-d','--dirs',
        nargs='+',
        help="Target directories where the all the extracted tables are organized per ref-code",
    )
    (
        parser.add_argument(
            '-o',
            "--out_dir",
            help="Saving table for the output",
        ),
    )
    args = parser.parse_args()
    print(args.dirs)
    print(args.out_dir)
    main(args.dirs, args.out_dir)
