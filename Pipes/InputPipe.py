from Pipes.Pipe import Pipe
import os
import config
import sys
import utils
from DBContext import DBContext

class Setup(Pipe):

    def __init__(self):
        pass

    def process(self, **kwargs):
        path = kwargs.pop('path', None)
        if path is None:
            sys.exit("No path to create the principal directory!")

        panel = kwargs.pop('panel')
        panel = panel if panel == 'trusightone' else panel.upper()
        genome = kwargs.pop('genome', 'geno38')
        over = kwargs.pop('over', True)
        dest = kwargs.pop('dest')
        proj_name = kwargs.pop('project')
        if proj_name is None:
            proj_name = utils.createProjectName()

        principal_directory = os.path.join(path, config.RESULTS_DIRECTORY_NAME, "{}_{}".format(proj_name, panel))
        #TODO: add logging
        download_folder = self.create_downloadfolders(principal_directory, dest)

        config.update_cfg("DEFAULTS", "PRINCIPAL_DIRECTORY", principal_directory)
        config.update_cfg("DEFAULTS", "DOWNLOAD_FOLDER", download_folder)
        #self.db_setup(dest)

        kwargs.update({"proj_name": proj_name, "panel": panel, "path": path, "genome": genome, "over": over, "name_folder": principal_directory, "dest": dest})
        #print(kwargs)
        return kwargs

    def create_downloadfolders(self, principal_directory, server_id):
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

        

        return download_folder

