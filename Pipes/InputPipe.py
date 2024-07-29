from Pipes.Pipe import Pipe
import os
import config
import sys
import utils
import shutil
from DBContext import DBContext
import dir_tree


""" Initiation of the project organization structure. The name of the results directory will be constructed, and be searched for its existence. Additionally,
the config will be updated with additional information, that will be only useful for the current project execution. """

class Setup(Pipe):

    def __init__(self):
        pass

    def process(self, **kwargs):
        upd_config_params = []
        path = kwargs.pop('path', None)
        panel = kwargs.pop('panel')
        panel = panel if panel == 'trusightone' else panel.upper()
        genome = kwargs.pop('genome', 'geno38')
        over = kwargs.pop('over', True)
        dest = kwargs.pop('dest')
        proj_name = kwargs.pop('project', None)
    

        #TODO: add logging 

        principal_directory = self.create_principal_directory(path, proj_name, over, panel) # TODO: maybe not a good idea to keep principal_directory (as well as dir_tree) global.
        download_folder = self.create_downloadfolders(principal_directory, dest)

        upd_config_params += [("DEFAULTS", "PRINCIPAL_DIRECTORY", principal_directory), ("DEFAULTS", "DOWNLOAD_FOLDER", download_folder)]
        config.update(upd_config_params)

        dir_tree.build(principal_directory)
        print(dir_tree.principal_directory)
        kwargs.update({"proj_name": proj_name, "panel": panel, "path": path, "genome": genome, "over": over, "principal_directory": principal_directory, "dest": dest})
        #print(kwargs)
        return kwargs



    def create_principal_directory(self, path, proj_name, over, panel):
        if path is None:
            sys.exit("No path to create the principal directory!")

        if proj_name is None:
            proj_name = utils.createProjectName()

        principal_directory = os.path.join(path, config.RESULTS_DIRECTORY_NAME, "{}_{}".format(proj_name, panel))

        if os.path.exists(principal_directory):
            if not over:
                # sys.exit("Folder exist just - Choose another project name, please!!! You can run with <-proj> option")
                return principal_directory
            else:
                shutil.rmtree(principal_directory)

        return principal_directory


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
