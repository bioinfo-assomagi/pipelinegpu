import json
import os
import config

""" Parse the directory tree of the project results. All the results of the project are stored in the RESULT directory, under the name of the principal_directory, which can be specified
by the user, or automatically by considering the current date and the name of the panel. """

"""
Module Description:
-------------------

This module parses the directory tree of the project results. All the results of the project are stored
in the RESULT directory, under the name of the `principal_directory`, which can be specified
by the user or automatically determined by considering the current date and the name of the panel.

It provides the following classes and functions:

- `dotdict`: A dictionary subclass that allows dot notation access to dictionary attributes.
- `DirectoryTreeBuilder`: A class that builds the directory tree and creates directories accordingly.
- `explore_directory_tree`: A recursive function that explores the directory structure defined in a JSON file.
- `build`: A function that instantiates `DirectoryTreeBuilder` and updates global variables with the directory paths.

**Usage Example:**

Assuming you have a JSON file defining the directory structure, you can build and create the directory tree as follows:

.. code-block:: python

    d = DirectoryTreeBuilder("/path/to/principal_directory").run()
    print(d)

"""

#config.PRINCIPAL_DIRECTORY = "/home/magi/PROJECT/diagnosys/RESULTS_jurgen/27_Mar_2024_CANCER"

class dotdict(dict):
    """
    dot.notation access to dictionary attributes
    
    **Usage Example:**

    .. code-block:: python

        d = dotdict()
        d.key = 'value'
        print(d['key'])  # Output: 'value'

    """
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
    

class DirectoryTreeBuilder():

    """
    The `DirectoryTreeBuilder` class parses a directory tree defined in a JSON file and creates the directory structure
    starting from a specified principal directory.

    **Attributes:**

    - `principal_directory` (str): The root directory where the tree will be built.
    - `directories` (list): A list of directories to be created.

    **Methods:**

    - `parse()`: Parses the directory tree JSON file and builds a dictionary representing the directory structure.
    - `run()`: Parses the directory tree and creates the directories.

    **Usage Example:**

    .. code-block:: python

        builder = DirectoryTreeBuilder("/path/to/principal_directory")
        d = builder.run()
        # Now, the directories are created, and 'd' contains the directory paths.
    """

    tree = None
    list_dir = None

    def __init__(self, principal_directory):
        self.principal_directory = principal_directory
        self.directories = None

    def parse(self): #TODO: create l with a separate function
        """
        Parses the directory tree JSON file and builds a dictionary representing the directory structure.

        :return: A `dotdict` object representing the directory structure with paths.
        :rtype: dotdict
        :raises Exception: If the principal directory is not specified or if the JSON file cannot be loaded.
        """
        
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
        """
        Parses the directory tree and creates the directories.

        :return: A `dotdict` object representing the directory structure with paths.
        :rtype: dotdict
        :raises Exception: If there is an error building the directory tree.
        """
         
        d = self.parse()
        
        if self.directories is not None:
            for directory in self.directories:
                os.makedirs(directory, exist_ok = True)
        else:
            raise Exception("Error building directory tree!")

        return d

def explore_directory_tree(directory_tree, parent_path, l, d):
    """ 
    Recursively explores the directory structure defined in the JSON file. Saves the paths of each directory
    into a list, which will be used to create those directories. Also saves the paths in a dictionary to allow
    programmatic access in the form of, e.g., `dir_tree.temp.to_count`.

    :param directory_tree: The current directory tree node.
    :type directory_tree: dict
    :param parent_path: The path of the parent directory.
    :type parent_path: str
    :param directories_list: A list to collect directory paths.
    :type directories_list: list
    :param directory_dict: A dictionary to collect directory paths for programmatic access.
    :type directory_dict: dotdict
    """

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
    """
    Builds the directory tree starting from the specified principal directory and updates global variables
    with the directory paths.

    :param principal_directory: The root directory where the tree will be built.
    :type principal_directory: str
    """
     
    d = DirectoryTreeBuilder(principal_directory).run()
    globals().update(d)

"""
Note:
-----

In the `parse()` method of `DirectoryTreeBuilder`, there is a tight coupling between the creation of directories
and the construction of the directory tree dictionary. There will be no dictionary of directories unless those
directories were created successfully. Therefore, checking if the dictionary is in the global variables also
checks if the directories were created.

This design ensures that when running different pipes of the pipelines separately, the directories will not be
overwritten due to `os.makedirs(directory, exist_ok=True)`. However, the only issue that persists is regarding
the `config.PRINCIPAL_DIRECTORY`, a variable which will be created only if the `InputPipe()` is executed.
If the user wants to run, say, `CoveragePipe()` only, they may need to edit the config file manually.
"""

# if "principal_directory" not in globals():
#     build()

if __name__ == "__main__":
    
    d = DirectoryTreeBuilder("/home/magi/PROJECT/diagnosys/RESULTS_jurgen/test_build_tree").run()
    print(d)