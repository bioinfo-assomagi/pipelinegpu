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

"""
Module Description:
-------------------

This module initiates the project organization structure. It constructs the name of the results directory
and checks for its existence. Additionally, it updates the configuration with information that is only
useful for the current project execution.

Classes:
- `Setup`: Initializes the project structure and updates configurations.

Functions:
- None

**Usage Example:**

.. code-block:: python

    setup_pipe = Setup()
    kwargs = {
        'path': '/path/to/results',
        'panel': 'CANCER',
        'genome': 'hg38',
        'over': False,
        'dest': 'r',
        'project': 'ProjectName'
    }
    result = setup_pipe.process(**kwargs)
"""

class Setup(Pipe):

    """
    The `Setup` class is a subclass of `Pipe` that initializes the project organization structure.
    It constructs the results directory, updates configurations, and prepares the environment
    for the pipeline execution.

    **Methods:**

    - `process(**kwargs)`: Main method to process the setup with given arguments.
    - `create_principal_directory(path, proj_name, over, panel)`: Creates the principal directory for results.
    - `create_downloadfolders(principal_directory, server_id)`: Creates the download folders based on server ID.
    """

    def __init__(self):
        pass

    def process(self, **kwargs):

        """
        Processes the setup by creating the principal directory, updating configurations,
        and building the directory tree.

        :param kwargs: Keyword arguments containing setup parameters.
            - `path` (str): The base path where the results directory will be created.
            - `panel` (str): The panel name (e.g., 'trusightone', 'CANCER').
            - `genome` (str, optional): The genome version (default is 'geno38').
            - `over` (bool, optional): Whether to overwrite existing directories (default is True).
            - `dest` (str): Destination server ID ('r', 'b', 's', 'z').
            - `project` (str, optional): The project name.

        :return: Updated keyword arguments with additional setup information.
        :rtype: dict

        **Example:**

        .. code-block:: python

            setup_pipe = Setup()
            kwargs = {
                'path': '/path/to/results',
                'panel': 'CANCER',
                'genome': 'hg38',
                'over': False,
                'dest': 'r',
                'project': 'ProjectName'
            }
            result = setup_pipe.process(**kwargs)
        """

        upd_config_params = []
        path = kwargs.pop('path', None)
        panel = kwargs.pop('panel')
        panel = panel if panel == 'trusightone' else panel.upper()
        genome = kwargs.pop('genome', 'geno38')
        over = kwargs.pop('over', True)
        dest = kwargs.pop('dest')
        proj_name = kwargs.pop('project', None)
        # TODO: take fastq names so you can build the phenotype file
    

        #TODO: add logging 

        principal_directory = self.create_principal_directory(path, proj_name, over, panel) # TODO: maybe not a good idea to keep principal_directory (as well as dir_tree) global.
        download_folder = self.create_downloadfolders(principal_directory, dest)

        upd_config_params += [("DEFAULTS", "PRINCIPAL_DIRECTORY", principal_directory), ("DEFAULTS", "DOWNLOAD_FOLDER", download_folder)]
        config.update(upd_config_params)

        dir_tree.build(principal_directory)
        print(dir_tree.principal_directory)

        #TODO: create the phenotype file here
        kwargs.update({"proj_name": proj_name, "panel": panel, "path": path, "genome": genome, "over": over, "principal_directory": principal_directory, "dest": dest})
        #print(kwargs)
        return kwargs



    def create_principal_directory(self, path, proj_name, over, panel):
        """
        Creates the principal directory for storing results.

        :param path: The base path where the results directory will be created.
        :type path: str
        :param proj_name: The project name.
        :type proj_name: str
        :param over: Whether to overwrite the existing directory if it exists.
        :type over: bool
        :param panel: The panel name (e.g., 'CANCER').
        :type panel: str
        :return: The path to the principal directory.
        :rtype: str

        :raises SystemExit: If the `path` is not provided.

        **Example:**

        .. code-block:: python

            principal_dir = self.create_principal_directory('/path/to/results', 'ProjectName', False, 'CANCER')
        """

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
        """
        Creates the download folders based on the server ID.

        :param principal_directory: The principal directory path.
        :type principal_directory: str
        :param server_id: The server identifier ('r', 'b', 's', 'z').
        :type server_id: str
        :return: The path to the download folder.
        :rtype: str

        **Example:**

        .. code-block:: python

            download_folder = self.create_downloadfolders('/path/to/principal_directory', 'r')
        """
         
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
