#! /usr/bin/env python3

import os
import math
import argparse
import textwrap
import sys
import yaml
import json
import datetime
import re
import requests
import shutil
import glob
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


def concatenate_ips_chunks(path):
    """Concatenate the I5 files for the MGF functional-annotation"""

    prefix = path.name.split(".")[0]
    log.debug(f"Prefix: {prefix}")

    # Get the I5 search path
    dirname = path.parents[0]
    search = dirname / f"{prefix}.merged_CDS.I5_0*"
    I5_paths = glob.glob(str(search))  # Glob needs a string not Path object

    outpath = Path(dirname, f"{prefix}.merged_CDS.I5.tsv.gz")
    log.debug(f"Outpath: {outpath}")
    # Concatenate the I5 files
    with open(outpath, "wb") as wfp:
        for I5_path in I5_paths:
            with open(I5_path, "rb") as rfp:
                log.debug(f"Concatenating {I5_path}...")
                shutil.copyfileobj(rfp, wfp)
    log.info("I5 chunks concatenated")


def remove_hmm_chunk_file(target_directory, conf):
    """Remove the single HMM chunk file"""
    p = Path(
        target_directory,
        "results",
        "functional-annotation",
        f"{conf['prefix']}.merged.hmm.tsv.chunks",
    )
    if p.exists():
        log.info("Removing HMM chunk file")
        p.unlink()


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


def sequence_categorisation_stanzas(target_directory, template, conf):
    """Glob the sequence_categorisation directory and build a stanza for each
    zipped data file; add to MANDATORY_FILES list

    Return updated template, and list of sequence category filenames
    """

    seq_cat_dir_path = Path(target_directory, "results", "sequence-categorisation")
    seq_cat_paths = seq_cat_dir_path.glob("*.gz")
    # Just the file names as @ids changed later
    seq_cat_files = [sq.name for sq in seq_cat_paths]
    log.debug(f"Seq_cat_files: {seq_cat_files}")
    # Add the sequence categorisation files to the list of mandatory files
    # So that they can be used to build the upload script later
    qualified_paths = [
        "/".join(["./sequence-categorisation", str(sq)]) for sq in seq_cat_files
    ]
    # global MANDATORY_FILES
    MANDATORY_FILES.extend(qualified_paths)
    log.debug(f"MANDATORY_FILES (after seq_categorisation)= {MANDATORY_FILES}")

    # Sequence-categorisation stanza
    for i, stanza in enumerate(template["@graph"]):
        if stanza["@id"] == "./sequence-categorisation/":
            # NB with the qualified paths
            stanza["hasPart"] = [dict([("@id", f"{fn}")]) for fn in qualified_paths]
            log.debug(f"Seq catagoriastaion stanza['hasPart'] == {stanza['hasPart']}")
            sq_index = i
            break

    qualified_paths.reverse()
    for fn in qualified_paths:
        bits = Path(fn).name.split(".")
        if bits[0] == "5_8S":
            name_string = f"RNA prediction for {bits[0]}"
        else:
            name_string = (
                f"RNA prediction for {bits[0]} - Rfam accesssion number: {bits[1]}"
            )

        d = dict(
            [
                ("@id", f"{fn}"),
                ("@type", "File"),
                ("name", name_string),
                ("downloadUrl", ""),
                ("encodingFormat", "application/zip"),
            ]
        )
        template["@graph"].insert(sq_index + 1, d)
    return template


def add_sequence_data_stanzas(target_directory, template, conf):
    """
    Add the sequence data stanzas to the template:

    DBH_AAACOSDA_1_1_HWLTKDRXY.UDI211_clean.fastq.trimmed.fasta.bz2
    DBH_AAACOSDA_1_1_HWLTKDRXY.UDI211_clean.fastq.trimmed.qc_summary
    DBH_AAACOSDA_1_2_HWLTKDRXY.UDI211_clean.fastq.trimmed.fasta.bz2
    DBH_AAACOSDA_1_2_HWLTKDRXY.UDI211_clean.fastq.trimmed.qc_summary
    DBH.merged_CDS.faa.bz2
    DBH.merged_CDS.ffn.bz2
    DBH.merged.cmsearch.all.tblout.deoverlapped.bz2
    DBH.merged.fasta.bz2
    DBH.merged.motus.tsv.bz2
    DBH.merged.qc_summary
    DBH.merged.unfiltered_fasta.bz2

    """
    data = {
        r"^{prefix}_[A-Za-z0-9]+_[1,2]_1_[A-Za-z0-9]+\.[A-Za-z0-9]+_clean\.fastq\.trimmed\.fasta\.bz2$": (
            "Trimmed forward reads",
            "All forward reads after trimming in fasta format",
            "application/x-bzip2",
        ),
        r"^{prefix}_[A-Za-z0-9]+_[1,2]_1_[A-Za-z0-9]+\.[A-Za-z0-9]+_clean\.fastq\.trimmed\.qc_summary$": (
            "Trimmed forward reads QC summary",
            "Quality control summary of trimmed forward reads",
            "text/plain",
        ),
        r"^{prefix}_[A-Za-z0-9]+_[1,2]_2_[A-Za-z0-9]+\.[A-Za-z0-9]+_clean\.fastq\.trimmed\.fasta\.bz2$": (
            "Trimmed reverse reads",
            "All reverse reads after trimming in fasta format",
            "application/x-bzip2",
        ),
        r"^{prefix}_[A-Za-z0-9]+_[1,2]_2_[A-Za-z0-9]+\.[A-Za-z0-9]+_clean\.fastq\.trimmed\.qc_summary$": (
            "Trimmed reverse reads QC summary",
            "Quality control summary of trimmed reverse reads",
            "text/plain",
        ),
        "{prefix}.merged_CDS.faa.bz2": (
            "Protein coding amino acid sequences",
            "Coding sequences of merged reads in amino acid format",
            "application/x-bzip2",
        ),
        "{prefix}.merged_CDS.ffn.bz2": (
            "Protein coding nucleotide sequences",
            "Coding sequences of merged reads in nucleotide format",
            "application/x-bzip2",
        ),
        "{prefix}.merged.cmsearch.all.tblout.deoverlapped.bz2": (
            "Overlapped coding sequences",
            "Overlapped coding sequences (intermediate file)",
            "application/x-bzip2",
        ),
        "{prefix}.merged.fasta.bz2": (
            "Merged reads",
            "Merged forward and reverse reads in fasta format",
            "application/x-bzip2",
        ),
        "{prefix}.merged.motus.tsv.bz2": (
            "MOTUs",
            "Metagenomic Operational Taxonomic Units (MOTUs) in tab-separated format",
            "application/x-bzip2",
        ),
        "{prefix}.merged.qc_summary": (
            "QC summary of merged reads",
            "Quality control analysis summary of merged reads",
            "text/plain",
        ),
        "{prefix}.merged.unfiltered_fasta.bz2": (
            "Unfiltered merged reads",
            "All merged reads before fileting in fasta format",
            "application/x-bzip2",
        ),
    }
    seq_data_paths = Path(target_directory, "results").glob(f"{conf['prefix']}*")
    # Just the file names as @ids changed later
    seq_data_files = [sq.name for sq in seq_data_paths]
    if len(seq_data_files) != 11:
        log.error("Expected 11 sequence data files, found %s" % len(seq_data_files))
        log.error(f"sequence data files: {seq_data_files}")
        sys.exit()
    log.debug(f"Seq_data_files: {seq_data_files}")

    # global MANDATORY_FILES
    MANDATORY_FILES.extend([f"./{fn}" for fn in seq_data_files])
    log.debug(f"MANDATORY_FILES (after seq_data_files) = {MANDATORY_FILES}")

    for stanza in template["@graph"]:
        if stanza["@id"] == "./":
            # Update the hasPart dict with the sequence data files
            ns = [dict([("@id", f"./{fn.format(**conf)}")]) for fn in seq_data_files]
            stanza["hasPart"].extend(ns)
            break

    # Add the sequence data stanzas after config.yml
    # Get the config.yml index
    for sq_index, stanza in enumerate(template["@graph"]):
        if stanza["@id"] == "./config.yml":
            break

    for seq_file in seq_data_files:
        log.debug(f"Looking for seq_file: {seq_file}")
        for patt, value in data.items():
            found = False
            filename = patt.format(**conf)
            if re.match(filename, seq_file):
                d = dict(
                    [
                        ("@id", f"./{seq_file}"),
                        ("@type", "File"),
                        ("name", value[0]),
                        ("description", value[1]),
                        ("downloadUrl", ""),
                        ("encodingFormat", value[2]),
                    ]
                )
                template["@graph"].insert(sq_index + 1, d)
                found = True
                log.debug(f"Added stanza for {seq_file}")
                break
        if not found:
            log.error(f"Cannot find pattern for {seq_file}")
            sys.exit()

    return template


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
    return conf


def check_and_format_data_file_paths(target_directory, conf, check_exists=True):
    """Check that all mandatory files are present in the target directory"""

    workflow_yaml_path = WORKFLOW_YAML_FILENAME.format(**conf)
    filepaths = [f.format(**conf) for f in MANDATORY_FILES]
    if check_exists:
        path = Path(target_directory, workflow_yaml_path)
        log.debug(f"Looking for worflow YAML file at: {path}")
        if not os.path.exists(path):
            log.error(f"Cannot find workflow YAML file at {path}")
            sys.exit()
        else:
            MANDATORY_FILES.append(workflow_yaml_path)
        # The fixed file paths
        for filepath in filepaths:
            log.debug(f"File path: {filepath}")

            # Deal with config.yml
            if filepath == "./config.yml":
                path = Path(target_directory, filepath)
            else:
                path = Path(target_directory, "results", filepath)

            # Originally MGF did not concatenate the I5 files so this is needed for some of V1.0 and development runs
            if (
                filepath
                == f"./functional-annotation/{conf['prefix']}.merged_CDS.I5.tsv.gz"
                and not path.exists()
            ):
                log.info(
                    "Functional annotation I5 files are in chunks... concatenating..."
                )
                concatenate_ips_chunks(path)

            # Deal with the emapper summary file
            if (
                filepath
                == f"./functional-annotation/{conf['prefix']}.merged.emapper.summary.eggnog"
                and not path.exists()
            ):
                log.info("Eggnog emapper summary file is missing")
                fn = "./functional-annotation/{prefix}.merged.emapper.summary.eggnog"
                log.debug(f"Removing '{fn} from MANDATORY_FILES")
                MANDATORY_FILES.remove(fn)
                continue

            # And the rest
            if not os.path.exists(path):
                for filename in conf["missing_files"]:
                    if os.path.split(filepath)[1] == filename.format(**conf):
                        # This file is known to be missing, ignoring
                        log.info(
                            "Ignoring specified missing file: %s"
                            % os.path.split(filepath)[1]
                        )
                        parts = list(Path(filepath).parts[:-1])
                        fn = os.path.join("./", str(Path(*parts, filename)))
                        log.debug(f"Removed {fn} from MANDATORY_FILES")
                        MANDATORY_FILES.remove(fn)
                        break
                else:
                    log.error(
                        "Could not find the mandatory file '%s' at the following path: %s"
                        % (filepath, path)
                    )
                    log.error(
                        "Consider adding it to the 'missing_files' list in the YAML configuration."
                    )
                    log.error("Cannot continue...")
                    sys.exit()
            else:
                log.debug("Found %s" % path)

    # There's a single HMM file that needs to be removed
    remove_hmm_chunk_file(target_directory, conf)
    log.info("Data look good...")
    return conf


def get_creator_and_mgf_version_information(conf, overide_error=False):
    """Get the creator and mgf version information."""

    # Read the relevant row in sample sheet
    try:
        df_samp = pd.read_csv(COMBINED_LOGSHEETS_PATH)
    except requests.exceptions.RequestException:
        log.error("Cannot find the combined logsheets at %s" % COMBINED_LOGSHEETS_PATH)
        sys.exit()
    row_samp = df_samp.loc[df_samp["ref_code"] == conf["ref_code"]].to_dict()
    log.debug("Row in sample sheet: %s" % row_samp)
    # Get the env_package either water_column or soft_sediments
    try:
        env_package = list(row_samp["env_package"].values())[0]
    except IndexError:
        log.error("Cannot find the env_package for ref_code %s" % conf["ref_code"])
        sys.exit()
    log.debug("env_package: %s" % env_package)
    assert env_package in ["water_column", "soft_sediment"], (
        "env_package must be either 'water_column' or 'soft_sediment'"
        "Found: %s" % env_package
    )

    # Add MGF analysis creator_person
    mgf_path = FILTERS_MGF_PATH if env_package == "water_column" else SEDIMENTS_MGF_PATH
    data = pd.read_csv(mgf_path).to_dict(orient="records")
    for row in data:
        if row["ref_code"] == conf["ref_code"]:
            log.debug("Row in %s: %s" % (mgf_path, row))
            if row["who"] == "CCMAR":
                conf["creator_person_name"] = "Cymon J. Cox"
                conf["creator_person_identifier"] = (
                    "https://orcid.org/0000-0002-4927-979X"
                )

            elif row["who"] == "HCMR":
                conf["creator_person_name"] = "Stelios Ninidakis"
                conf["creator_person_identifier"] = (
                    "https://orcid.org/0000-0003-3898-9451"
                )

            else:
                log.error("Unrecognised creater of MGF data: %s" % row["who"])
                sys.exit()

            # Metagoflow
            log.info("MetaGOflow run at: %s" % row["who"])
            log.info("MetaGOflow version: %s" % row["version"])
            conf["metagoflow_version_id"] = row["version"]
            # Hard coding the metagoflow version URL here:
            if row["version"] == "develop (3cf3a7d)":
                conf["metagoflow_version"] = (
                    "https://github.com/emo-bon/MetaGOflow/commit/3cf3a7d39fabc6e75a8cb2971a711c2d781c84d0"
                )
            elif row["version"] == "1.0":
                conf["metagoflow_version"] = (
                    "https://github.com/emo-bon/MetaGOflow/releases/tag/v1.0.0"
                )
            else:
                log.error("Unrecognised MetaGOflow version: %s" % row["version"])
                sys.exit()

            break
    else:
        if not overide_error:
            log.error(
                f"Cannot find the creator of MGF data for ref_code {conf['ref_code']}"
            )
            log.error(
                "You can overide this error by setting the 'overide_error' flag to True"
            )
            log.error("This will set the default to CCMAR")
            sys.exit()
        else:
            log.info("Creator person missing: defaulting to CCMAR")
            conf["creator_person_name"] = "Cymon J. Cox"
            conf["creator_person_identifier"] = "https://orcid.org/0000-0002-4927-979X"
            log.info("MetaGOflow version missing: defaulting to develop (3cf3a7d)")
            conf["metagoflow_version_id"] = "develop (3cf3a7d)"
            conf["metagoflow_version"] = (
                "https://github.com/emo-bon/MetaGOflow/commit/3cf3a7d39fabc6e75a8cb2971a711c2d781c84d0"
            )

    return conf


def get_ena_accession_data(conf):
    """Get the ENA accession data for a given ref_code."""
    # Read the relevant row in sample sheet
    if conf["batch_number"] == 1:
        df_ena = pd.read_csv(BATCH1_ENA_ACCESSION_INFO_PATH)
    elif conf["batch_number"] == 2:
        df_ena = pd.read_csv(BATCH2_ENA_ACCESSION_INFO_PATH)
    else:
        log.error(f"Batch number not recognised {conf['batch_number']}")
        sys.exit()
    row_ena = df_ena.loc[df_ena["ref_code"] == conf["ref_code"]].to_dict()
    # Get the ENA accession data
    conf["ena_accession_number"] = list(
        row_ena["ena_accession_number_run_metag"].values()
    )[0]
    conf["ena_accession_number_url"] = (
        f"https://www.ebi.ac.uk/ena/browser/view/{conf['ena_accession_number']}"
    )
    return conf


def add_sequence_data_links(conf, override_error=False):
    """Add the links to the raw sequence data files in ENA"""
    # ENA ACCESSION filereport for sample
    filereport_url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={ena_accession_number}"
        "&result=read_run&fields=submitted_ftp&format=json&download=true&limit=-1"
    )
    log.debug(f"ENA filereport URL: {filereport_url.format(**conf)}")
    filereport = requests.get(filereport_url.format(**conf))
    if filereport.status_code == requests.codes.ok:
        filereport_json = filereport.json()
        log.debug(f"ENA filereport: {filereport_json}")
        # Get the FTP links to the raw sequence data
        links = filereport_json[0]["submitted_ftp"].split(";")
        if len(links) != 2:
            log.error("Cannot find the 2 raw sequence data links in the ENA filereport")
            sys.exit()
        # Add the links to the conf dictionary
        conf["forward_reads_link"] = f"https://{links[0]}"
        conf["reverse_reads_link"] = f"https://{links[1]}"
        log.info("ENA raw sequence data links added to conf")
    else:
        if not override_error:
            log.error("Cannot get the ENA filereport")
            log.error("Raw sequence data are not available from ENA")
            log.error("Use override_error flag to bypass this error")
            sys.exit()
        else:
            log.info("Raw sequence data links missing: defaulting to local")
            conf["forward_reads_link"] = "http://localhost/forward_reads_link"
            conf["reverse_reads_link"] = "http://localhost/reverse_reads_link"

    return conf


def write_metadata_json(
    target_directory, conf, without_sequence_data=False, override_error=False
):
    metadata_json_template = "ro-crate-metadata.json-template"
    if os.path.exists(metadata_json_template):
        log.debug("Using local metadata.json template")
        with open(metadata_json_template, "r") as f:
            template = json.load(f)
    else:
        # Grab the template from Github
        log.debug("Downloading metadata.json template from Github")
        req = requests.get(TEMPLATE_URL)
        if req.status_code == requests.codes.ok:
            template = req.json()
        else:
            log.error("Unable to download the metadata.json file from Github")
            log.error(f"Check {TEMPLATE_URL}")
            log.error("Exiting...")
            sys.exit()

    # Build the conf dictionary
    conf = get_creator_and_mgf_version_information(conf, override_error)
    conf = get_ena_accession_data(conf)
    conf = add_sequence_data_links(conf, override_error)
    log.debug("Conf dict: %s" % conf)

    log.info("Writing ro-crate-metadata.json...")

    # Add strings first
    # Add "ref_code"'s to "name", "title", and "description" fields
    template["@graph"][1]["name"] = template["@graph"][1]["name"].format(**conf)
    template["@graph"][1]["description"] = template["@graph"][1]["description"].format(
        **conf
    )

    # Add date to "datePublished"
    if "datePublished" in conf:
        template["@graph"][1]["datePublished"] = template["@graph"][1][
            "datePublished"
        ].format(**conf)
    else:
        template["@graph"][1]["datePublished"] = datetime.datetime.now().strftime(
            "%Y-%m-%d"
        )

    # Deal with the hasMember stanzas
    for section in template["@graph"]:
        if section["@id"] == "./":
            for member in section["pcdm:hasMember"]:
                if member["@id"] == "{ena_accession_number_url}":
                    member["@id"] = member["@id"].format(**conf)
                elif member["@id"] == "{metagoflow_version}":
                    member["@id"] = member["@id"].format(**conf)
                else:
                    log.error("Cannot find the hasMember stanzas")
                    sys.exit()
    # Add metadGOflow version id
    for section in template["@graph"]:
        if section["@id"] == "{metagoflow_version}":
            section["@id"] = section["@id"].format(**conf)
            section["softwareVersion"] = section["softwareVersion"].format(**conf)
            section["downloadUrl"] = section["downloadUrl"].format(**conf)
            break
    else:
        log.error("Cannot find the MetaGOflow version stanza")
        sys.exit()

    # Add ena_accession_number to the field
    for stanza in template["@graph"]:
        # Not yet formatted
        if stanza["@id"] == "{ena_accession_number_url}":
            stanza["@id"] = stanza["@id"].format(**conf)
            stanza["name"] = stanza["name"].format(**conf)
            stanza["downloadUrl"] = stanza["downloadUrl"].format(**conf)
            continue
        # Add the raw sequence data links
        if stanza["@id"] in ["{forward_reads_link}", "{reverse_reads_link}"]:
            stanza["@id"] = stanza["@id"].format(**conf)
            stanza["description"] = stanza["description"].format(**conf)
            stanza["downloadUrl"] = stanza["downloadUrl"].format(**conf)

    # creator  - the MGF data creator and institution
    # "creator": {}
    template["@graph"][1]["creator"] = template["@graph"][1]["creator"] = dict(
        [("@id", f"{conf['creator_person_identifier']}")]
    )

    # Add creater person stanza
    person_stanza = dict(
        [
            ("@id", f"{conf['creator_person_identifier']}"),
            ("@type", "Person"),
            ("name", f"{conf['creator_person_name']}"),
        ]
    )
    template["@graph"].insert(5, person_stanza)

    # Add eggnog summary file if present
    if (
        "./functional-annotation/{prefix}.merged.emapper.summary.eggnog"
        in MANDATORY_FILES
    ):
        log.info("Adding Eggnog emapper summary stanza to graph")
        eggnog_summary = dict(
            [
                (
                    "@id",
                    f"./functional-annotation/{conf['prefix']}.merged.emapper.summary.eggnog",
                ),
                ("@type", "File"),
                ("name", "Eggnog emapper summary"),
                ("description", "Summary of Eggnog emapper analysis"),
                ("downloadUrl", ""),
                ("encodingFormat", "text/plain"),
            ]
        )
        for stanza in template["@graph"]:
            if stanza["@id"] == "./functional-annotation/":
                fn = f"./{conf['prefix']}.merged.emapper.summary.eggnog"
                stanza["hasPart"].append(dict([("@id", fn)]))
                break
        else:
            log.error("Cannot find the functional-annotation stanza")
            sys.exit()
        # Add the eggnog summary stanza to the graph before the sequence categorisation stanzas
        for i, stanza in enumerate(template["@graph"]):
            if stanza["@id"] == "./sequence-categorisation/":
                template["@graph"].insert(i, eggnog_summary)
                log.debug(f"Added eggnog summary stanza at index {i}")
                break
    # Add sequence_categorisation stanza separately as they can vary in number and identity
    template = sequence_categorisation_stanzas(target_directory, template, conf)
    # Add sequence data stanzas
    if not without_sequence_data:
        template = add_sequence_data_stanzas(target_directory, template, conf)

    log.info("Metadata (first part) JSON written...")
    return template


def main(
    target_directory,
    yaml_config,
    debug,
    without_sequence_data=False,
    override_error=False,
):
    """ """
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
    log.debug(f"ro_crate_name is: {ro_crate_name}")
    if os.path.exists(ro_crate_name):
        log.error(f"An archive with the name {ro_crate_name} already exists")
        log.error("Exiting...")
        sys.exit()

    # Check all files are present
    log.info("Checking data files...")
    # filepaths includes the seq_cat files
    check_and_format_data_file_paths(target_directory, conf, check_exists=True)

    # Create the metadata.json file but dont write yet, need to add links later
    log.info("Formatting metadata.json...")
    metadata_json = write_metadata_json(
        target_directory, conf, without_sequence_data, override_error
    )

    #save metadata_josn to a file
    log.info("Writing metadata_part1.json to file.")
    with open("metadata_part1.json", "w") as f:
        f.write(json.dumps(metadata_json, indent=4))

    # save the path to a file
    log.info(f"Writing new folder name to 'path1.csv' file.")
    with open("path.csv", "w") as f:
        f.write(str(target_directory + ', ' + conf["source_mat_id"]))


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
    parser.add_argument(
        "-w",
        "--without_sequence_data",
        action="store_true",
        default=False,
        help="Do not add sequence data files (default: False)",
    )
    parser.add_argument(
        "-o",
        "--override_error",
        action="store_true",
        default=False,
        help="Override creator person missing error (default: False)",
    )
    args = parser.parse_args()
    main(
        args.target_directory,
        args.yaml_config,
        args.debug,
        args.without_sequence_data,
        args.override_error,
    )