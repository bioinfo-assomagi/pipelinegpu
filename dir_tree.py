import json
import os
import config

""" Parse the directory tree of the project results. All the results of the project are stored in the RESULT directory, under the name of the principal_directory, which can be specified
by the user, or automatically by considering the current date and the name of the panel. """

#config.PRINCIPAL_DIRECTORY = "/home/magi/PROJECT/diagnosys/RESULTS_jurgen/27_Mar_2024_CANCER"

class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
    

class DirectoryTreeBuilder():

    tree = None
    list_dir = None

    def __init__(self, principal_directory):
        self.principal_directory = principal_directory
        self.directories = None

    def parse(self): #TODO: create l with a separate function
        principal_directory = self.principal_directory

        if principal_directory == None:
            raise Exception("Principal Directory not specified!")

        with open("/home/magi/PROJECT/diagnosys/bin_jurgen/directory_tree.json") as f:
            d_tree = json.load(f)

        d_tree['name'] = ""

        l = []
        d = dotdict({"principal_directory": dotdict({})})
        explore_directory_tree(d_tree, principal_directory, l, d.principal_directory)
        self.directories = l
        return d
    
    def run(self):
        d = self.parse()
        
        if self.directories is not None:
            for directory in self.directories:
                os.makedirs(directory, exist_ok = True)
        else:
            raise Exception("Error building directory tree!")

        return d

def explore_directory_tree(directory_tree, parent_path, l, d):
    """ Explores the directory stored in the json file. Save the paths of each directory or file into a list, which will be used by the program to create those paths.
    Saves the paths in a dict as well, to be able to allow programmatic access in the form of e.g. dir_tree.temp.to_count. """

    current = directory_tree['name']

    current_path = os.path.join(parent_path, current)
    l.append(current_path)
    d.path = current_path
 
    if directory_tree["type"] == "folder":
        if "children" in directory_tree.keys():
            for child in directory_tree["children"]:
                child_name = child["name"]
                d[child_name] = dotdict({})
                explore_directory_tree(child, current_path, l, d[child_name])

def build(principal_directory):
    d = DirectoryTreeBuilder(principal_directory).run()
    globals().update(d)

""" In the parse() mehtod of DirectoryTreeBuilder there is a tight coupling between the creation of directories and the construction of the directory tree dictionary. 
There will be no dictionary of directories unless those directories were created successfully. Therefore, checking if the dictionary is in the global variables, also checks if the directories
were created (NOTE: this will be a problem if we want to run different pipes of the pipelines separately. Actually, no, since the directories will be not overwritten by os.makedirs(directory, exist_ok = True)). 
The only issue that perists is regarding the config.PRINCIPAL_DIRECTORY, a variable which will be created only if the InputPipe() is executed. However, if the user wants to run, say, CoveragePipe()
only, he may edit the config file manually.
"""

# if "principal_directory" not in globals():
#     build()

if __name__ == "__main__":
    
    d = DirectoryTreeBuilder("/home/magi/PROJECT/diagnosys/RESULTS_jurgen/test_build_tree").run()
    print(d)