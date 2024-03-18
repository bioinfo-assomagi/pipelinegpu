from Pipes.Pipe import Pipe

class QualityReadsPipe(Pipe):
     
    def __init__(self):
        print("QualityReadsPipe init!")

    def process(self, **kwargs):
        # reading some args
        name_folder = kwargs.pop("name_folder", None) # the usual name_folder
        dummy_data = kwargs.pop("dummy_var", 'default_dummy_var')

        if name_folder:
            print("QualityReadsPipe is reading the {} directory ...".format(name_folder))

        print("QualityReadsPipe has received the following dummy_data: {}".format(dummy_data))

        return {'name_folder': name_folder, "dummy_var" : 'dummy_var2'}