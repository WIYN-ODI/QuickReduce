#!/usr/bin/env python
from __future__ import print_function

#
# Update this variable below when changing versions
#
pipeline_plver = "QuickReduce 2.2rc1"
pipeline_name = "QuickReduce"
pipeline_version = "2.2rc1"

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


pl_author = "Ralf Kotulla"
pl_email = "kotulla@wisc.edu"
def record_pipeline_versioning(header):

    header["PLVER"] = (pipeline_plver, "name and version")
    header["PIPELINE"] = (pipeline_name, "pipeline name")
    header["PLVERSIO"] = (pipeline_version, "pipeline version")
    header["PLAUTHOR"] = (pl_author, "pipeline author")
    header["PLEMAIL"] = (pl_email, "contact email")
    header['FOR_HELP'] = (odi_help_email, "contact ODI help with questions")

    src_directory,_ = os.path.split(__file__)

    # get some info from GIT - branch name and last commit
    branch_name = get_branch_name(src_directory)
    commit_hash = get_last_commit(src_directory)

    header['GITBRNCH'] = (branch_name, "git branch name")
    header['GITCOMIT'] = (commit_hash, "git commit hash")

    add_fits_header_title(header, "Pipeline information", 'PLVER')


if __name__ == "__main__":

    git_branch = get_branch_name()
    git_commit = get_last_commit()

    if (len(sys.argv) > 1 and sys.argv[1] == "batch"):
        print(pipeline_version, git_branch, git_commit)
        sys.exit(0)

    print("This is %s version %s (GIT:: branch:%s commit:%s)" % (
        pipeline_name, pipeline_version, git_branch, git_commit))
    print("Developed by: %s (%s)" % (pl_author, pl_email))
    print("For questions, email %s" % (odi_help_email))
    pass
