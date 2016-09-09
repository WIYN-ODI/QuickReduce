#!/usr/bin/env python

#
# Update this variable below when changing versions
#
pipeline_plver = "QuickReduce 2.0"
pipeline_name = "QuickReduce"
pipeline_version = "2.0.5"

odi_help_email = "odi-ppa-help@googlegroups.com"


#
#   |
#   |  Do NOT change the code below when updating version information
#   v
#

import subprocess
import os
import sys
from podi_definitions import add_fits_header_title


def get_branch_name(src_directory=None):

    if (src_directory is None):
        src_directory,_ = os.path.split(__file__)

    cmd = "git rev-parse --abbrev-ref HEAD"
    get_branch = subprocess.Popen(
        cmd.split(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=src_directory
    )
    out,err = get_branch.communicate()
    branch_name = out.strip() if (get_branch.returncode == 0) else "UNKNOWN"
    return branch_name

def get_last_commit(src_directory=None):

    if (src_directory is None):
        src_directory,_ = os.path.split(__file__)

    cmd = "git rev-parse --verify HEAD"
    get_commit = subprocess.Popen(
        cmd.split(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=src_directory
    )
    out, err = get_commit.communicate()
    commit_hash = out.strip() if (get_commit.returncode == 0) else "??????????????????"
    return commit_hash



def record_pipeline_versioning(header):

    header["PLVER"] = (pipeline_plver, "name and version")
    header["PIPELINE"] = (pipeline_name, "pipeline name")
    header["PLVERSIO"] = (pipeline_version, "pipeline version")
    header["PLAUTHOR"] = ("Ralf Kotulla", "pipeline author")
    header["PLEMAIL"] = ("kotulla@wisc.edu", "contact email")
    header['FOR_HELP'] = (odi_help_email, "contact ODI help with questions")

    src_directory,_ = os.path.split(__file__)

    # get some info from GIT - branch name and last commit
    branch_name = get_branch_name(src_directory)
    commit_hash = get_last_commit(src_directory)

    header['GITBRNCH'] = (branch_name, "git branch name")
    header['GITCOMIT'] = (commit_hash, "git commit hash")

    add_fits_header_title(header, "Pipeline information", 'PLVER')



if __name__ == "__main__":
    pass
