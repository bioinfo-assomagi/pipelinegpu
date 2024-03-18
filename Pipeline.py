
from Pipes.Pipe import Pipe
        
""" The pipeline class. Command and decorator pattern. An instance of the Pipeline class will be constitued by many Pipes assembled together sequentially. 
Each pipe will use the output of the previous pipe as input. """

class Pipeline():

    lastPipe = None

    def __init__(self, pipe : Pipe):
        self.lastPipe = pipe

    def assemblePipe(self, newPipe : Pipe):
        pipe = Pipe()
        pipe.process = lambda **kwargs : newPipe.process(**self.lastPipe.process(**kwargs))
        
        # now pipe will be the lastPipe for the returned Pipeline. i.e. we are returning an extended Pipeline, with pipe as the last assembled pipeline containing the new call stack, with newPipe.process() on top
        return Pipeline(pipe) 

    def start(self, **kwargs):
        return self.lastPipe.process(**kwargs)
