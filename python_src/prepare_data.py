#! /usr/bin/env python3

from pathlib import Path
import sys
import os
import logging as log
import argparse
import textwrap
import subprocess
import psutil

from utils import find_bzip2, open_archive

desc = """
Prepare the MGF data archives for the ro-crate building. All files in the target
directory will be opened and the sequence files compressed.

This script opens the MGF results archive and compresses the individual data
files from building the ro-crate.
"""

FILE_PATTERNS = [
    "*.fastq.trimmed.fasta",
    "*.merged_CDS.faa",
    "*.merged_CDS.ffn",
    "*.merged.cmsearch.all.tblout.deoverlapped",
    "*.merged.fasta",
    "*.merged.motus.tsv",
    "*.merged.unfiltered_fasta",
    "final.contigs.fa",
]

def find_archive(home_dir, subfolders, archive_name):
    for subfolder in subfolders:
        path = Path(home_dir, subfolder)
        log.debug(f"Checking {path}")
        if path.exists():
            tarball_files = list(path.glob("*.tar.bz2"))
            log.debug(f"Found {len(tarball_files)} tarball files")
            tarball_file_names = [f.name for f in tarball_files]
            log.debug(f"Looking for {archive_name}in {tarball_file_names}")
            if archive_name in tarball_file_names:
                return path
    return None


def main(
    target_directory,
    archive_name,
    out_dir=os.getcwd(),
    debug=False,
):
    log.basicConfig(
        format="\t%(levelname)s: %(message)s", level=log.DEBUG if debug else log.INFO
    )
    # check first which folder contains the archive
    home_dir = Path.cwd()
    subfolders = ["BLANKS", "FILTERS", "MOCKS", "SEDIMENTS"]
    target_directory = Path(home_dir, target_directory)
    log.debug(f"Looking for {archive_name} in {target_directory}")

    target_directory = find_archive(target_directory, subfolders, archive_name)
    if not target_directory:
        log.error(f"Cannot find the target directory for {archive_name} in {target_directory}")
        sys.exit()

    log.debug(f"Found target directory {target_directory}")
    os.chdir(target_directory)

    # Where the open archives will go
    outpath = "prepared_archives"

    # concatenate the out_dir to the outpath
    log.info(f"Creating prepared archives in {out_dir} / {outpath}")
    outpath = os.path.join(out_dir, outpath)
    if Path(outpath).exists():
        log.debug("'prepared_archives' directory already exists")
    else:
        log.debug("Creating 'prepared_archives' directory")
        Path(outpath).mkdir()

    bzip2_program = find_bzip2()

    run_id = Path(str(archive_name).rsplit(".", 2)[0])

    oca_dir = Path(outpath, f"{run_id}")
    if oca_dir.exists():
        log.info(f"An prepared archive already exists for {run_id}")
        sys.exit()

    if run_id.exists():
        log.debug("Found open archive")
    else:
        # Open the archive
        log.info(f"Opening archive {archive_name}")
        open_archive(os.path.abspath(archive_name), bzip2_program, outpath)

    # Compress the sequence archive files
    log.info(f"Compressing sequence files for {run_id}")
    os.chdir(Path(run_id, "results"))
    log.debug(f"CWD: {os.getcwd()}")
    for fp in FILE_PATTERNS:
        sf = Path("./").glob(fp)
        for f in sf:
            log.debug(f"Compressing {f}")
            # Can't use f{} style formatting in subprocess call
            # of the program name
            if bzip2_program == "lbzip2":
                threads = psutil.cpu_count() - 4
                log.debug(f"Using lbzip2 with {threads} threads")
                subprocess.check_call(
                    [
                        "lbzip2",
                        "-9",
                        f"-n {threads}",
                        f"./{f}",
                    ]
                )
            elif bzip2_program == "bzip2":
                log.debug(f"bzip2 -9 {f}")
                subprocess.check_call(
                    [
                        "bzip2",
                        "-9",
                        f"./{f}",
                    ]
                )

    log.debug(f"Moving to {target_directory}")
    os.chdir(target_directory)

    # if len(tarball_files) == 0:
    #     log.error(f"Cannot find any tarball files in {target_directory}")
    #     sys.exit()

    # # Check the target_directory name
    # target_directory = Path(target_directory)
    # if not target_directory.exists():
    #     log.error(f"Cannot find the target directory {target_directory}")
    #     sys.exit()
    # log.debug(f"Found target directory {target_directory}")

    # # CD to target directory
    # log.info(f"Changing directory to {target_directory}")
    # home_dir = Path.cwd()
    # target_directory = Path(home_dir, target_directory)
    # os.chdir(target_directory)

    # # Get list of tarball files
    # tarball_files = list(Path.cwd().glob("*.tar.bz2"))

    # log.debug(f"Found {len(tarball_files)} tarball files")
    # for tarball in tarball_files:
    #     log.debug(f"Tarball name: {tarball.name}")
    # tarball_files = [f.name for f in tarball_files]
    # if len(tarball_files) == 0:
    #     log.error(f"Cannot find any tarball files in {target_directory}")
    #     sys.exit()

    # # Where the open archives will go
    # outpath = "prepared_archives"

    # # concatenate the out_dir to the outpath
    # log.info(f"Creating prepared archives in {out_dir} / {outpath}")
    # outpath = os.path.join(out_dir, outpath)
    # if Path(outpath).exists():
    #     log.debug("'prepared_archives' directory already exists")
    # else:
    #     log.debug("Creating 'prepared_archives' directory")
    #     Path(outpath).mkdir()

    # bzip2_program = find_bzip2()

    # if archive_name:
    #     tarball_files = [archive_name]

    # # Loop through the tarball files
    # for tarball_file in tarball_files:
    #     run_id = Path(str(tarball_file).rsplit(".", 2)[0])

    #     oca_dir = Path(outpath, f"{run_id}")
    #     if oca_dir.exists():
    #         log.info(f"An prepared archive already exists for {run_id}")
    #         continue

    #     if run_id.exists():
    #         log.debug("Found open archive")
    #     else:
    #         # Open the archive
    #         log.info(f"Opening archive {tarball_file}")
    #         open_archive(os.path.abspath(tarball_file), bzip2_program, outpath)

    #     # Compress the sequence archive files
    #     log.info(f"Compressing sequence files for {run_id}")
    #     os.chdir(Path(run_id, "results"))
    #     log.debug(f"CWD: {os.getcwd()}")
    #     for fp in FILE_PATTERNS:
    #         sf = Path("./").glob(fp)
    #         for f in sf:
    #             log.debug(f"Compressing {f}")
    #             # Can't use f{} style formatting in subprocess call
    #             # of the program name
    #             if bzip2_program == "lbzip2":
    #                 threads = psutil.cpu_count() - 4
    #                 log.debug(f"Using lbzip2 with {threads} threads")
    #                 subprocess.check_call(
    #                     [
    #                         "lbzip2",
    #                         "-9",
    #                         f"-n {threads}",
    #                         f"./{f}",
    #                     ]
    #                 )
    #             elif bzip2_program == "bzip2":
    #                 log.debug(f"bzip2 -9 {f}")
    #                 subprocess.check_call(
    #                     [
    #                         "bzip2",
    #                         "-9",
    #                         f"./{f}",
    #                     ]
    #                 )

        # # Move the compressed files to the open-compressed-samples directory
        # log.debug(f"Moving to {target_directory}")
        # os.chdir(target_directory)
        # # log.info(f"Moving compressed files to '{outpath}'")
        # # log.debug(f"Moving {run_id} to {oca_dir}")
        # # Path(run_id).rename(oca_dir)

    # CD back to home directory
    # log.info(f"Changing directory to {home_dir}")
    # os.chdir(home_dir)

    # log.info("Done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(desc),
    )
    parser.add_argument(
        "target_directory",
        help=(
            "Name of target directory containing MetaGOflow output"
            " archives to prepare relative to the current working directory"
        ),
    )
    parser.add_argument("archive_name", help="process single archive only")
    parser.add_argument(
        "-o",
        "--out_dir",
        help="Name of directory to store the prepared archives",
        default=os.getcwd(),
    )
    parser.add_argument("-d", "--debug", action="store_true", help="DEBUG logging")
    args = parser.parse_args()
    main(
        args.target_directory,
        args.archive_name,
        args.out_dir,
        args.debug,
    )