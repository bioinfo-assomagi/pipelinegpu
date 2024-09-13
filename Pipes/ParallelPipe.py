

from Pipes.Pipe import Pipe
import os
import config
import sys
import utils


class ParallelPipe(Pipe):

    def __init__(self):
        super().__init__()

    def process(self, **kwargs):
        return kwargs

    def thread_print(self, msg):
        utils.thread_print(self.thread_id, msg)
