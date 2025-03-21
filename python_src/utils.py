import os
import sys
import subprocess
import logging
from pathlib import Path
import shutil

log = logging.getLogger(__name__)


"""
Utilities for dealing with MGF data archives


"""


def find_bzip2():
    # Find best bzip2 programme
    if shutil.which("lbzip2"):
        # Get nunmber of threads/cpu cores
        bzip2_program = "lbzip2"
    elif shutil.which("bzip2"):
        log.info("Using bzip2")
        bzip2_program = "bzip2"
    else:
        log.error("Cannot find lbzip2 or bzip2")
        log.error("Exiting...")
        sys.exit()
    return bzip2_program


def open_archive(tarball_file, bzip2_program, outpath=None, hcmr=False):
    """
    You are expected to be in the dir with the tarball when calling this function
    """
    # Create temp dir in which to untar archives
    if outpath is not None:
        os.chdir(outpath)
    try:
        Path("temp").mkdir(exist_ok=False)
    except FileExistsError as e:
        log.error(
            f"A directory called 'temp' in the target directory already exists: {e}"
        )
        log.error("Exiting...")
        sys.exit()
    os.chdir("temp")
    
    # tar --use-compress-program lbunzip2 -xvf ../HMNJKDSX3.UDI200.tar.bz2
    if hcmr:
        log.debug(f"Opening archive {tarball_file} with unzip")
        subprocess.check_call(
            [
                "unzip",
                f"{tarball_file}",
            ]
        )
    else:
        log.debug(f"Opening archive {tarball_file} with {bzip2_program}")
        subprocess.check_call(
            [
                "tar",
                "--use-compress-program",
                f"{bzip2_program}",
                "-xf",
                f"{tarball_file}",
            ]
        )
    # Check archive
    run_id = Path(str(tarball_file).split('/')[-1].rsplit(".", 2)[0])
    log.debug(f"run_id = {run_id}")
    if run_id.exists():
        os.chdir("..")
        # Move archive up to target directory
        Path("temp", run_id).rename(run_id)
        shutil.rmtree("temp")
    else:
        # Archive with top level directory
        log.debug("Checking archive in ./temp")
        yml_files = list(Path.cwd().glob("*.yml"))
        log.debug(f"Found {yml_files} yml files")
        if Path("results").exists() and len(yml_files) == 2:
            log.debug("Found archive without top level directory")
            os.chdir("..")
            log.debug(f"Renaming ./temp to {run_id}")
            Path("temp").rename(run_id)
        else:
            # Deal with broken archive
            log.error(f"Archive looks completely broken at {tarball_file}")
            os.chdir("..")
            shutil.rmtree("temp")
            log.debug("Renaming broken archive")
            Path(tarball_file).rename(f"{tarball_file}-broken")


# Legacy code dont use; just open_archive
def fix_all_archives(target_directory, fix_archive, debug):
    """MGF archives should have a top level directory with the run_id
    e.g. HMNJKDSX3.UDI20 below which sits the results directory

    """

    log.basicConfig(
        format="\t%(levelname)s: %(message)s", level=log.DEBUG if debug else log.INFO
    )

    # Check the target_directory name
    target_directory = Path(target_directory)
    if not target_directory.exists():
        log.error(f"Cannot find the target directory {target_directory}")
        sys.exit()
    log.debug(f"Found target directory {target_directory}")

    # CD to target directory
    log.info(f"Changing directory to {target_directory}")
    home_dir = Path.cwd()
    os.chdir(target_directory)

    # Get list of tarball files
    tarball_files = list(Path.cwd().glob("*.tar.bz2"))

    log.debug(f"Found {len(tarball_files)} tarball files")
    for tarball in tarball_files:
        log.debug(f"Tarball: {tarball}")
        log.debug(f"Tarball name: {tarball.name}")
    tarball_files = [f.name for f in tarball_files]
    if len(tarball_files) == 0:
        log.error(f"Cannot find any tarball files in {target_directory}")
        sys.exit()
    else:
        # Where the open archives will go
        log.debug("Creating fixed-archives directory")
        Path("fixed-archives").mkdir(exist_ok=True)

    bzip2_program = find_bzip2()

    for tarball_file in tarball_files:
        run_id = Path(str(tarball_file).rsplit(".", 2)[0])
        log.debug(f"run_id = {run_id}")
        open_archive(tarball_file, bzip2_program)

    os.chdir(home_dir)