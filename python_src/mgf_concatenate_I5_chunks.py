#! /usr/bin/env python3

import argparse
import glob
import logging as log
import shutil
import textwrap
from pathlib import Path

desc = """
Concatenate the I5 chunks of a MGF functional-annotation run

e.g.
Concatenate the I5 files:
<prefix>.merged_CDS.I5_001.tsv.gz
<prefix>.merged_CDS.I5_002.tsv.gz
into:
<prefix>.merged_CDS.I5.tsv.gz

Remove the chunks and the chunk list file
"""


def main(target_dir, debug=False):
    """Concatenate the I5 files for the MGF functional-annotation"""
    # Logging
    log_level = log.DEBUG if debug else log.INFO
    log.basicConfig(format="\t%(levelname)s: %(message)s", level=log_level)
    log.info("Concatenating I5 chunks...")

    # Get the directory path
    dirpath = Path(target_dir, "results", "functional-annotation")

    # Glob the chunks file to get the prefix:
    chunk_file_path = dirpath / "*.merged_CDS.I5.tsv.chunks"
    result = glob.glob(str(chunk_file_path))  # Glob needs a string not Path object

    if not result:
        log.info(f"No chunk file found at {chunk_file_path}")
        log.info("Exiting...")
        return
    chunk_file = result[0]
    log.debug(f"Found chunk file: {result[0]}")
    prefix = chunk_file.split("/")[-1].split(".")[0]
    log.debug(f"Prefix: {prefix}")

    # Get the I5 search path
    search = dirpath / f"{prefix}.merged_CDS.I5_0*"
    # log.debug(f"I5 search path: {search}")
    I5_paths = glob.glob(str(search))  # Glob needs a string not Path object
    # log.debug(f"I5 paths: {I5_paths}")

    outpath = Path(dirpath, f"{prefix}.merged_CDS.I5.tsv.gz")
    log.debug(f"Outpath: {outpath}")
    # Concatenate the I5 files
    with open(outpath, "wb") as wfp:
        for I5_path in I5_paths:
            with open(I5_path, "rb") as rfp:
                log.debug(f"Concatenating {I5_path}...")
                shutil.copyfileobj(rfp, wfp)
    log.info("I5 chunks concatenated")

    log.info("Removing I5 chunks and chunk list file")
    Path(chunk_file).unlink()
    for chunk in I5_paths:
        Path(chunk).unlink()
    log.info("Done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(desc),
    )
    (
        parser.add_argument(
            "target_directory",
            help="Name of target directory containing MetaGOflow output",
        ),
    )
    parser.add_argument("-d", "--debug", action="store_true", help="DEBUG logging")
    args = parser.parse_args()
    main(args.target_directory, args.debug)