import os
import logging
import re
import sys
import configparser

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIGNAME = os.path.join(ROOT_DIR, 'config.cfg')  # not really used at all


def _get_cfg(cfg):
    if os.path.exists(cfg):
        # 1) is it already an accessible file? ::: if a config.cfg is already in the path from where you are executing the pipeline (many times 
        # bin/diagnosys_pipleline.py which in turn calls bin_jurgen/launcher2.py, will be executed) it will get that config.cfg. e.g.
        # since the cfg argument is specified inside the CustomConfigParser constructor and is set to be just "config.cfg"
        # if there is such named filed inside the calling directory (e.g. PROJECT/diagnosys/) it will point there. i.e. enter in this condition
        # and proceed reading the config.cfg already in PROJECT/diagnosys/
        pass
    elif os.path.exists(os.path.join(ROOT_DIR, cfg)):
        # 2) search in the execution directory ::: it will always loook for bin_jurgen/config.cfg
        # it will enter this condition if there is no config.cfg file in the e.g. /PROJECT/diagnosys/ directory
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
            'SERVER_RICERCA': self.get('SERVER', 'RICERCA'),
            'REF_GENOME_NAME': self.get('PARABRICKS', 'REF_GENOME_NAME'),
            'BUCHIARTIFICIALI': self.get('FILES', 'BUCHIARTIFICIALI'),
            'GENO37': self.get('FILES', 'GENO37'),
            'GENO38': self.get('FILES', 'GENO38'),
            'EREDITA37': self.get('FILES', 'EREDITA37'),
            'EREDITA38': self.get('FILES', 'EREDITA38'),
            'HGMD37': self.get('FILES', 'HGMD37'),
            'HGMD38': self.get('FILES', 'HGMD38'),
            'APPRIS': self.get('FILES', 'APPRIS'),
            'ECCEZIONI': self.get('FILES', 'ECCEZIONI'),
            'dbNSFP38_gz': self.get('FILES', 'dbNSFP38_gz'),
            'dbscSNV11_gz': self.get('FILES', 'dbscSNV11_gz'),
            'OCULAR': self.get('VERTICAL', 'OCULAR'),
            'VASCULAR': self.get('VERTICAL', 'VASCULAR'),
            'NEUROLOGY': self.get('VERTICAL', 'NEUROLOGY'),
            'INTEGRACARDIOSTANCHEZZA': self.get('VERTICAL', 'INTEGRACARDIOSTANCHEZZA'),
            'MIXED': self.get('VERTICAL', 'MIXED'),
            'LYMPHOBESITY': self.get('VERTICAL', 'LYMPHOBESITY'),
            'INFERTILIT': self.get('VERTICAL', 'INFERTILIT'),
            'CHERATOCONO': self.get('VERTICAL', 'CHERATOCONO'),
            'PCDH19': self.get('VERTICAL', 'PCDH19'),
            'GENEOBNEW': self.get('VERTICAL', 'GENEOBNEW'),
            'CANCER': self.get('VERTICAL', 'CANCER'),
            'OVERWRITE': self.get('MISC', 'OVERWRITE'),
            'PIPELINETYPE': self.get('MISC', 'PIPELINETYPE'),
            'SANGER_MODEL': self.get('MODELS', 'SANGER_MODEL')

            # 'PARABRICKS_INPUTDIR': self.get('PARABRICKS', 'INPUTDIR_NAME')
        }

        return configuration

    # def get_configuration(filename=CONFIGNAME):  # not used, globals updated when the file is imported
    #     """Reads and parses the configuration file."""
    #     cfg = CustomConfigParser(filename=filename)  # update module-level cfg
    #     globals().update(cfg.configuration)  # update configdir, templatesdir ...
    #     #configuration = cfg.configuration          # update module-level configuration
    #     return cfg



def update(upd_config_params):
    for params in upd_config_params:
        update_cfg(*params)

# TODO: rename to private
def update_cfg(section, parameter, value):
    cfg = CustomConfigParser()
    cfg.set(section, parameter, value)
    cfg.write() # can be ignored if we don't want to save dynamic config options in the file
    globals().update({parameter: value})


globals().update(CustomConfigParser().configuration)