from Pipes.Pipe import Pipe
import os
import config
import sys
import utils

class EndPipe(Pipe):

    def __init__(self):
        pass

    def process(self, **kwargs):
        self.thread_id = kwargs.pop("thread_id", None)
        print("End Pipe reached. The process with id:{} will now be terminated ...".format(self.thread_id))
        
        return kwargs
        
