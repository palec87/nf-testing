
import io
import math
import os
import sys
import json
import urllib.request
from pathlib import Path

import pandas as pd

PROJECT_DIR = Path.cwd()
sys.path.append(str(PROJECT_DIR))
# sys.path.append(str(PROJECT_DIR / "src"))

from minio import S3Error
from minio import Minio


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
BATCH1AND2_TOTAL = 188

# Expected number of metaGOflow analyses in version 1 of the data release
EXPECTED_ANALYSES = 54
OUT_PATH = PROJECT_DIR / "data_tables"
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

TABLES = [
    "metagoflow_analyses.SSU",
    "metagoflow_analyses.LSU",
    # "metagoflow_analyses.go_slim",
    # "metagoflow_analyses.go"
    # "metagoflow_analyses.ips"
    # "metagoflow_analyses.ko",
    # "metagoflow_analyses.pfam"
]

# The destination bucket and filename on the MinIO server
bucket_name = "emo-bon-data"
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


def parse_inventories(inv):
    """
    Parse the LSU and SSU inventories
    """
    objects = client.list_objects(bucket_name)
    count = 0
    all_objs_data = []  # list of dicts, each a taxonomic entry

    for obj in objects:
        all_sample_data = []
        # For each of the 54 LSU inventories
        prefix = code_keys[obj.object_name[:-1]][1]
        fn = f"{prefix}.merged_{inv}.fasta.mseq.tsv"
        fp = os.path.join(obj.object_name, fn)
        try:
            response = client.get_object(bucket_name, fp)
        except S3Error:
            raise
        finally:
            csv_data = pd.read_csv(response, sep="\t", skiprows=1)
            response.close()
            response.release_conn()
        for _, row in csv_data.iterrows():
            data = {}  # For one taxonomic entry
            # For each row in the inventory
            taxonomy = parse_taxonomy(row["taxonomy"])
            data["ref_code"] = code_keys[obj.object_name[:-1]][0]
            # ref_code from Github Batch Run Information sheet
            data["reads_name"] = obj.object_name[:-1]  # reads_name from same
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
    print(f"Found {count} inventories from {inv}")
    if count != EXPECTED_ANALYSES:
        raise ValueError(f"Could not find all {EXPECTED_ANALYSES} records in v1")
    return all_objs_data


def parse_local_inventory(inv: str, code_keys: dict[tuple[str, str]], folder: Path = None):
    count = 0
    all_objs_data = []  # list of dicts, each a taxonomic entry

    for _, val_tuple in code_keys.items():
        all_sample_data = []
        # For each of the 54 LSU inventories
        prefix = val_tuple[1]
        fn = f"{prefix}.merged_{inv}.fasta.mseq.tsv"
        fp = os.path.join(folder, f"{val_tuple[2]}-tables", fn)
        lst = ["EMOBON_OSD74_Wa_21-tables", "EMOBON_VB_Wa_43-tables"]
        if any(x in fp for x in lst):
            print('FULL path', fp, val_tuple[2])
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
            # ref_code from Github Batch Run Information sheet
            data["reads_name"] = val_tuple[0]  # reads_name from same
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
    print(f"Found {count} inventories from {inv}")
    # if count != EXPECTED_ANALYSES:
    #     raise ValueError(f"Could not find all {EXPECTED_ANALYSES} records in v1")
    return all_objs_data


def main():
    batch1 = "https://raw.githubusercontent.com/emo-bon/sequencing-data/main/shipment/batch-001/run-information-batch-001.csv"
    batch2 = "https://raw.githubusercontent.com/emo-bon/sequencing-data/main/shipment/batch-002/run-information-batch-002.csv"


    code_keys = {}
    for batch in (batch1, batch2):
        with urllib.request.urlopen(batch) as f:  # noqa: S310
            data = pd.read_csv(f)  # CSV
            code_keys.update(extract_keys(data))

    if len(code_keys) != BATCH1AND2_TOTAL:
        raise ValueError(f"Expected {BATCH1AND2_TOTAL} keys, got {len(code_keys)}")
    else:
        print(f"Extracted the expected {len(code_keys)} records from batch sheets")

    # print(code_keys)
    LSU_data = parse_local_inventory("LSU", code_keys, folder=PROJECT_DIR / "results-tables/")


    # LSU_data = parse_inventories("LSU")
    # SSU_data = parse_inventories("SSU")
    # print(f"Parsed {len(LSU_data)} rows from LSU data")
    # print(f"Parsed {len(SSU_data)} rows from SSU data")
    # lsu_df = pd.DataFrame.from_records(LSU_data)
    # ssu_df = pd.DataFrame.from_records(SSU_data)
    # lsu_df.info()
    # ssu_df.info()
    # lsu_outfile = "metagoflow_analyses.LSU"
    # ssu_outfile = "metagoflow_analyses.SSU"
    # lsu_df.to_csv(OUT_PATH.joinpath(lsu_outfile), index=False)
    # ssu_df.to_csv(OUT_PATH.joinpath(ssu_outfile), index=False)


    # The destination bucket and filename on the MinIO server
    # bucket_name = "emo-bon-tables"

    # for table in TABLES:
    #     dfp = pd.read_csv(OUT_PATH / f"{table}")
    #     # https://www.iana.org/assignments/media-types/application/vnd.apache.parquet
    #     mime_application_type = "vnd.apache.parquet"
    #     df_bytes = dfp.to_parquet(None, engine="pyarrow", compression="snappy")
    #     buffer = io.BytesIO(df_bytes)
    #     bucket_name = Path(f"v{str(VERSION)}")
    #     file_name = Path(f"{table}.parquet")
    #     dest = bucket_name / file_name
    #     client.put_object(
    #         bucket_name,
    #         dest,
    #         data=buffer,
    #         length=len(df_bytes),
    #         content_type=mime_application_type,
    #     )
    #     print(f"\t{dest} successfully uploaded as object to {bucket_name}")


if __name__ == "__main__":
    main()


