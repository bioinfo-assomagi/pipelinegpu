import os
import config
import utils


def setup_directories(**kwargs):
    """ If principal_directory already exists, dont build it. """
    principal_directory = kwargs.pop("principal_directory", None)
    proj_name = kwargs.pop("project", None)

    if proj_name is None:
        proj_name = utils.createProjectName()

    principal_directory = os.path.join(path, config.RESULTS_DIRECTORY_NAME, "{}_{}".format(proj_name, panel))

    download_folder = create_downloadfolders(principal_directory, dest)
    
    config.update_cfg("DEFAULTS", "DOWNLOAD_FOLDER", download_folder)




def create_downloadfolders(principal_directory, server_id):
    download_folder = None
    if server_id == 'r':
        download_folder = '/home/magi/VIRTUAL/EUREGIO/DOWNLOADS/NGSINFO/BED_INFO/%s/%s/' % ('ROVERETO',principal_directory.split('/')[-1])
        if not os.path.exists(download_folder):
            os.makedirs(download_folder)
    elif server_id == 'b':
        download_folder = '/home/magi/VIRTUAL/EUREGIO/DOWNLOADS/NGSINFO/BED_INFO/%s/%s/' % ('BOLZANO',principal_directory.split('/')[-1])
        if not os.path.exists(download_folder):
            os.makedirs(download_folder)
    elif server_id == 's':
        download_folder = '/home/magi/VIRTUAL/SANFELICE/DOWNLOADS/NGSINFO/BED_INFO/%s/%s/' % ('SANFELICE',principal_directory.split('/')[-1])
        if not os.path.exists(download_folder):
            os.makedirs(download_folder)
    elif server_id == 'z':
        dest = 'RICERCA'
        download_folder = '/home/magi/VIRTUAL/RICERCA/DOWNLOADS/NGSINFO/BED_INFO/%s/' % (principal_directory.split('/')[-1])
        if not os.path.exists(download_folder):
            os.makedirs(download_folder)
