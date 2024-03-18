import os
import logging
import re
import sys
import configparser

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIGNAME = os.path.join(ROOT_DIR, 'config.cfg')  # not really used at all


def _get_cfg(cfg):
    if os.path.exists(cfg):
        # 1) is it already an accessible file?
        pass
    elif os.path.exists(os.path.join(ROOT_DIR, cfg)):
        # 2) search GENEMAGI
        return os.path.join(ROOT_DIR, cfg)
    else:
        raise ValueError("Failed to locate config file {}".format(cfg))

    return os.path.realpath(cfg)
    


class CustomConfigParser():
    """No config is created here, only read from files. Default config to be added. TODO: add a default config, that creates a .cfg file if none exists. 
    If a file exists, defaults should already be there, but will be overwritten by the file configs. """

    cfg_template = "config.cfg"

    def __init__(self, *args, **kwargs):
        self.parser = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation()) # Composition over inheritance
        self.filename = kwargs.pop("filename", CONFIGNAME)
        default_cfg = _get_cfg(self.cfg_template)
        self.parser.read_file(open(default_cfg))
    
    def __getattr__(self, attr_name: str):
        return getattr(self.parser, attr_name)
    
    def write(self):
        with open(self.cfg_template, 'w') as configfile:
            self.parser.write(configfile)

    @property
    def configuration(self):
        """Serves only to represent the configurations we are interested in as a dictionary of variables that will be
        made available as globals in the module."""
        configuration = {
            'configfilename': self.filename,
            'HOME_DIR': self.get('DEFAULTS', 'HOME_DIR'),
            'WORKING_PATH': self.get('DEFAULTS', 'WORKING_PATH'),
            'PROJECT_PATH': self.get('DEFAULTS', 'PROJECT_PATH'),
            'RESULTS_DIRECTORY_NAME': self.get('DEFAULTS', 'RESULTS_DIRECTORY_NAME'),
            'GATK': self.get('TOOLS', 'GATK'),
            'VEPS': self.get('TOOLS', 'VEPS'),
            'FASTQX': self.get('TOOLS', 'FASTQX'),
            'FASTQC': self.get('TOOLS', 'FASTQC'),
            'BWA': self.get('TOOLS', 'BWA'),
            'SAMTOOLS': self.get('TOOLS', 'SAMTOOLS'),
            'DB_EUREGIO': self.get('DATABASES', 'DB_EUREGIO'),
            'DB_MAGI': self.get('DATABASES', 'DB_MAGI'),
            'DB_RICERCA': self.get('DATABASES', 'DB_RICERCA'),
            'DB_PATH': self.get('DATABASES', 'DB_PATH'),
            'REF': self.get('PARABRICKS', 'REF'),
            'PARABRICKS_VERSION': self.get('PARABRICKS', 'VERSION'),
            'PARABRICKS_VERSION_DEEPVARIANT': self.get('PARABRICKS', 'DEEPVARIANT_VERSION'),
            'DOCKER_WORKDIR': self.get('DOCKER', 'WORKDIR_VOLUME'),
            'DOCKER_OUTPUTDIR': self.get('DOCKER', 'OUTPUTDIR_VOLUME'),
            'DOCKER_REFDIR': self.get('DOCKER', 'REFDIR_VOLUME'),
            'SERVER_EUREGIO': self.get('SERVER', 'EUREGIO'),
            'SERVER_MAGIS': self.get('SERVER', 'MAGIS'),
            'SERVER_RICERCA': self.get('SERVER', 'RICERCA')
            # 'PARABRICKS_INPUTDIR': self.get('PARABRICKS', 'INPUTDIR_NAME')
        }

        return configuration

    def get_configuration(filename=CONFIGNAME):  # not used, globals updated when the file is imported
        """Reads and parses the configuration file."""
        cfg = CustomConfigParser(filename=filename)  # update module-level cfg
        globals().update(cfg.configuration)  # update configdir, templatesdir ...
        # configuration = cfg.configuration          # update module-level configuration
        return cfg


def update_cfg(section, parameter, value):
    cfg = CustomConfigParser()
    cfg.set(section, parameter, value)
    #cfg.write() # can be ignored if we don't want to save dynamic config options in the file
    globals().update({parameter: value})


globals().update(CustomConfigParser().configuration)