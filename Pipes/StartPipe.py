
from Pipes.Pipe import Pipe
import os
import config
import sys
import utils



""" Start and End pipe are used as tail wrappers of the pipelines to log messages, or execute pre/post functionalities. """
class StartPipe(Pipe):

    def __init__(self):
        pass

    def process(self, **kwargs):
        self.thread_id = kwargs.pop("thread_id", None)
        if self.thread_id is not None:
            print("Start Pipe started. The process with id:{} will now start execution ...".format(self.thread_id))
        
        kwargs.update({"thread_id": self.thread_id})
        return kwargs
        
