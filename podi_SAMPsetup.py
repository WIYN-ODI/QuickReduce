
output_format = "/scratch/%OBSID.fits"

use_ssh = False #True
ssh_user = "rkotulla"
ssh_host = "localhost"
ssh_executable = "/work/podi_devel/podi_collectcells.py"

def translate_filename_local2remote(filename):

    # This is a template to re0write the filenames from local paths to the
    # path used on the remote system

    if (filename.startswith("/Volumes")):
        # This is on a Mac
        new_filename = filename.replace("/Volumes/", "/mnt/")
    else:
        new_filename = filename

    return new_filename
