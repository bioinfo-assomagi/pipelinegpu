from Pipes.Pipe import Pipe
import os
import config
import sys
import utils

class EndPipe(Pipe):

    def __init__(self):
        pass

    def process(self, **kwargs):
        print("End Pipe reached. The pipeline will now be terminated ...")

        return kwargs
        
