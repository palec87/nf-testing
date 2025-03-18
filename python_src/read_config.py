#! /usr/bin/env python3

import os
import math
import argparse
import textwrap
import sys
import yaml
import datetime
import logging as log
from pathlib import Path

import pandas as pd

desc = """
DP 25/03/13: this outputs the template metadata file for the ro-crate and additional file containing folder name for the renaming process


Build a MetaGOflow Data Products ro-crate from a YAML configuration.

Invoke
$ create-ro-crate.py <target_directory> <yaml_configuration>

where:
    target_directory is the toplevel output directory of MetaGOflow
        This must be the MGF run_id, e.g. HWLTKDRXY.UDI210
    yaml_configuration is a YAML file of metadata specific to this ro-crate
        a template is here:
        https://raw.githubusercontent.com/emo-bon/MetaGOflow-Data-Products-RO-Crate/main/ro-crate-config.yaml

e.g.

$ create-ro-crate.py HWLTKDRXY-UDI210 config.yml

This script expects to be pointed to directory of MetaGOflow output.

When invoked, the MetaGOflow run_wf.sh script writes all output to a directory specified by
the "-d" parameter:

    $ run_wf.sh -n green -d  HWLTKDRXY-UDI210 -f input_data/${DATA_FORWARD} -r input_data/${DATA_REVERSE}

    $ tree -1
    HWLTKDRXY-UDI210
    ├── prov
    ├── results
    ├── green.yml
    └── tmp

    3 directories, 1 file

"""

#########################################################################################
# These are the paths to the external data files that are used to build the RO-Crate
# run-information files for each batch - these return the ref_code for a given MGF run_id
BATCH1_RUN_INFO_PATH = (
    "https://raw.githubusercontent.com/emo-bon/sequencing-data/main/shipment/"
    "batch-001/run-information-batch-001.csv"
)
BATCH2_RUN_INFO_PATH = (
    "https://raw.githubusercontent.com/emo-bon/sequencing-data/main/shipment/"
    "batch-002/run-information-batch-002.csv"
)
# ENA ACCESSSION INFO for each batch
BATCH1_ENA_ACCESSION_INFO_PATH = (
    "https://raw.githubusercontent.com/emo-bon/sequencing-data/refs/heads/main/"
    "shipment/batch-001/ena-accession-numbers-batch-001.csv"
)
BATCH2_ENA_ACCESSION_INFO_PATH = (
    "https://raw.githubusercontent.com/emo-bon/sequencing-data/refs/heads/main/"
    "shipment/batch-002/ena-accession-numbers-batch-002.csv"
)
# The combined sampling event logsheets for batch 1 and 2
COMBINED_LOGSHEETS_PATH = (
    "https://raw.githubusercontent.com/emo-bon/emo-bon-data-validation/"
    "refs/heads/main/validated-data/Batch1and2_combined_logsheets_2024-11-12.csv"
)
# The ro-crate metadata template from Github
TEMPLATE_URL = (
    "https://raw.githubusercontent.com/emo-bon/MetaGOflow-Data-Products-RO-Crate"
    "/main/ro-crate-metadata.json-template"
)
# The MGF _Run_Track Google Sheets
FILTERS_MGF_PATH = (
    "https://docs.google.com/spreadsheets/d/"
    "1j9tRRsRCcyViDMTB1X7lx8POY1P5bV7UijxKKSebZAM/gviz/tq?tqx=out:csv&sheet=FILTERS"
)
SEDIMENTS_MGF_PATH = (
    "https://docs.google.com/spreadsheets/d/"
    "1j9tRRsRCcyViDMTB1X7lx8POY1P5bV7UijxKKSebZAM/gviz/tq?tqx=out:csv&sheet=SEDIMENTS"
)

# This is the workflow YAML file, the prefix is the "-n" parameter of the
# "run_wf.sh" script:
WORKFLOW_YAML_FILENAME = "./{run_parameter}.yml"
CONFIG_YAML_PARAMETERS = [
    "datePublished",
    "run_parameter",
    "missing_files",
]

# This is a submodule of the current repository
# https://github.com/emo-bon/metaGOflow-rocrates-dvc
RO_CRATE_REPO_PATH = "metaGOflow-rocrates-dvc"

MANDATORY_FILES = [
    "./fastp.html",
    "./RNA-counts",
    "./final.contigs.fa.bz2",
    "./config.yml",
    "./functional-annotation/stats/go.stats",
    "./functional-annotation/stats/interproscan.stats",
    "./functional-annotation/stats/ko.stats",
    "./functional-annotation/stats/orf.stats",
    "./functional-annotation/stats/pfam.stats",
    "./taxonomy-summary/LSU/krona.html",
    "./taxonomy-summary/SSU/krona.html",
    "./functional-annotation/{prefix}.merged_CDS.I5.tsv.gz",
    "./functional-annotation/{prefix}.merged.hmm.tsv.gz",
    "./functional-annotation/{prefix}.merged.summary.go",
    "./functional-annotation/{prefix}.merged.summary.go_slim",
    "./functional-annotation/{prefix}.merged.summary.ips",
    "./functional-annotation/{prefix}.merged.summary.ko",
    "./functional-annotation/{prefix}.merged.summary.pfam",
    "./functional-annotation/{prefix}.merged.emapper.summary.eggnog",
    "./taxonomy-summary/SSU/{prefix}.merged_SSU.fasta.mseq.gz",
    "./taxonomy-summary/SSU/{prefix}.merged_SSU.fasta.mseq_hdf5.biom",
    "./taxonomy-summary/SSU/{prefix}.merged_SSU.fasta.mseq_json.biom",
    "./taxonomy-summary/SSU/{prefix}.merged_SSU.fasta.mseq.tsv",
    "./taxonomy-summary/SSU/{prefix}.merged_SSU.fasta.mseq.txt",
    "./taxonomy-summary/LSU/{prefix}.merged_LSU.fasta.mseq.gz",
    "./taxonomy-summary/LSU/{prefix}.merged_LSU.fasta.mseq_hdf5.biom",
    "./taxonomy-summary/LSU/{prefix}.merged_LSU.fasta.mseq_json.biom",
    "./taxonomy-summary/LSU/{prefix}.merged_LSU.fasta.mseq.tsv",
    "./taxonomy-summary/LSU/{prefix}.merged_LSU.fasta.mseq.txt",
]


def get_ref_code_and_prefix(conf):
    """Get the reference code for a given run_id.
    run_id is the last part of the reads_name in the run information file.
    e.g. 'DBH_AAAAOSDA_1_HWLTKDRXY.UDI235' and is the name used to label the target_directory.
    """
    for i, batch in enumerate([BATCH1_RUN_INFO_PATH, BATCH2_RUN_INFO_PATH]):
        df = pd.read_csv(batch)
        for row in df[["reads_name", "ref_code", "source_mat_id"]].values.tolist():
            if isinstance(row[0], str):
                # print(row)
                if row[0].split("_")[-1] == conf["run_id"]:
                    conf["ref_code"] = row[1]
                    conf["prefix"] = row[0].split("_")[0]
                    conf["batch_number"] = i + 1
                    conf["source_mat_id"] = row[2]
                    log.info(f"EMO BON ref_code: {conf['ref_code']}")
                    log.info(f"Source mat ID: {conf['source_mat_id']}")
                    log.info(f"Prefix: {conf['prefix']}")
                    return conf
            elif math.isnan(row[0]):
                # Not all samples with an EMO BON code were sent to sequencing
                continue
    log.error("Cannot find the ref_code for run_id %s" % conf["run_id"])
    sys.exit()


def read_yaml(yaml_config):
    # Read the YAML configuration
    if not os.path.exists(yaml_config):
        log.error(f"YAML configuration file does not exist at {yaml_config}")
        sys.exit()
    with open(yaml_config, "r") as f:
        conf = yaml.safe_load(f)
    # Check yaml parameters are formated correctly, but not necessarily sane
    for param in CONFIG_YAML_PARAMETERS:
        log.debug(f"Config paramater {param}: {conf[param]}")
        if param == "datePublished":
            if conf[param] == "None":
                # No specified date, delete from conf
                # Its absence will trigger formatting
                # with today's date
                del conf[param]
                continue
            else:
                if not isinstance(conf[param], str):
                    log.error(
                        "'dataPublished' should either be a string or 'None'. Bailing..."
                    )
                    sys.exit()
                try:
                    datetime.datetime.fromisoformat(conf[param])
                except ValueError:
                    log.error(f"'datePublished' must conform to ISO 8601: {param}")
                    log.error("Bailing...")
                    sys.exit()
        elif param == "missing_files":
            if param not in conf:
                continue
            else:
                for filename in conf[param]:
                    if not isinstance(filename, str):
                        log.error(
                            f"Parameter '{filename}' in 'missing_files' list in YAML file must be a string."
                        )
                        log.error("Bailing...")
                        sys.exit()
        else:
            if not conf[param] or not isinstance(conf[param], str):
                log.error(f"Parameter '{param}' in YAML file must be a string.")
                log.error("Bailing...")
                sys.exit()
    log.info("YAML configuration looks good...")
    log.info(f"Configuration file: {conf}")
    return conf


def main(
    target_directory,
    yaml_config,
    debug,
):
    # Logging
    if debug:
        log_level = log.DEBUG
    else:
        log_level = log.INFO
    log.basicConfig(format="\t%(levelname)s: %(message)s", level=log_level)

    # Read the YAML configuration
    log.info("Reading YAML configuration...")
    conf = read_yaml(yaml_config)

    # Check the target_directory name
    if not os.path.exists(target_directory):
        log.error("Cannot find the target directory %s" % target_directory)
        sys.exit()
    log.debug("Found target directory %s" % target_directory)
    run_id = Path(target_directory).name
    log.debug("run_id = %s" % run_id)
    if "UDI" not in run_id.split(".")[-1]:
        log.error("Target directory name does NOT appear to be correct format")
        log.error("It needs to match the format HWLTKDRXY.UDI210")
        log.error("Exiting...")
        sys.exit()
    conf["run_id"] = run_id

    # Get the emo bon ref_code, batch number, and prefix
    conf = get_ref_code_and_prefix(conf)

    # Check that an archive with the same name does not already exist
    ro_crate_name = Path(RO_CRATE_REPO_PATH, conf["source_mat_id"] + "-ro-crate")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(desc),
    )
    parser.add_argument(
        "target_directory",
        help="Name of target directory containing MetaGOflow output",
    )
    parser.add_argument(
        "yaml_config", help="Name of YAML config file for building RO-Crate"
    )
    parser.add_argument("-d", "--debug", action="store_true", help="DEBUG logging")

    args = parser.parse_args()
    main(
        args.target_directory,
        args.yaml_config,
        args.debug,
    )