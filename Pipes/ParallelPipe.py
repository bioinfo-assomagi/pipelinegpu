

from Pipes.Pipe import Pipe
import os
import config
import sys
import utils

class ParallelPipe(Pipe):

    def __init__(self,):
        super().__init__()
        
    def process(self, **kwargs):
        return kwargs

    def thread_print(self, msg, name=__name__):
        utils.thread_print(self.thread_id, msg, name)

    def thread_filewrite(self, msg, file_path, mode, lock):
        # The lock is provided by the child
        if lock is None:
            raise ("No lock present!")
        
        lock.acquire()
        file = open(file_path, mode)
        file.write(msg)
        file.close()
        lock.release()

