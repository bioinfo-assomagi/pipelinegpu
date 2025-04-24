
from Pipes.Pipe import Pipe
        

class Pipeline():

    """
    The pipeline class. Command (Pipe class) and Decorator pattern. An instance of the Pipeline class will be constitued by many Pipes assembled together sequentially. 
    Each pipe will use the output of the previous pipe as input. 

    Command Pattern 
    Encapsulation of Actions:

    Each Pipe subclass encapsulates a specific action or processing step by implementing the process method. This method represents a command that can be executed.
    Decoupling Invoker and Receiver:

    The Pipeline class acts as the invoker that calls the process method on Pipe objects (the receivers) without needing to know the details of their implementations.
    Parameterization of Requests:

    You can create different Pipe subclasses with various process implementations and assemble them in different orders within the Pipeline, effectively parameterizing the sequence of commands to execute.
    Flexible and Extensible Processing Chain:

    The Pipeline allows you to add, remove, or reorder Pipe objects, providing a flexible way to manage a sequence of commands.

    Decorator Pattern 
    Dynamic Addition of Responsibilities:

    Each Pipe wraps the previous one by defining a new process method that calls the previous Pipe's process method and then applies additional processing. This adds new responsibilities dynamically.
    Wrapping Objects:

    In the assemblePipe method, a new Pipe is created where its process method wraps the process methods of the existing Pipe objects, effectively decorating them.
    Transparency to Client Code:

    The client interacts with the Pipeline in the same way, regardless of how many Pipe objects have been assembled. This is characteristic of the Decorator Pattern, where the decoration is transparent to the client.
    How Both Patterns Work Together in Your Code
    Command Aspect: Each Pipe acts as a command by encapsulating a processing operation. The Pipeline invokes these commands in sequence.
    Decorator Aspect: The Pipeline assembles Pipe objects in a way that each new Pipe decorates the previous one, enhancing or modifying its behavior without altering the existing Pipe implementations.

    **Examples:**

    .. code-block:: python

        class Pipe:
            def __init__(self):
                pass

            def process(self, **kwargs):
                raise NotImplementedError

        class PipeA(Pipe):
            def process(self, **kwargs):
                # Do something with kwargs
                kwargs['a'] = 'Processed by PipeA'
                return kwargs

        class PipeB(Pipe):
            def process(self, **kwargs):
                # Do something else with kwargs
                kwargs['b'] = 'Processed by PipeB'
                return kwargs

        # Assemble the pipeline
        pipeline = Pipeline(PipeA()).assemblePipe(PipeB())
        result = pipeline.start(data='Initial data')
        print(result)
        # Output: {'data': 'Initial data', 'a': 'Processed by PipeA', 'b': 'Processed by PipeB'}

    In this example:

    - **Command Pattern:** `PipeA` and `PipeB` encapsulate actions as commands.
    - **Decorator Pattern:** `PipeB` decorates `PipeA` by wrapping its `process` method.
    """

    lastPipe = None

    def __init__(self, pipe : Pipe):
        """
        Initializes the Pipeline with a given Pipe.

        :param pipe: The initial Pipe to start the Pipeline.
        :type pipe: Pipe
        """
        self.lastPipe = pipe

    def assemblePipe(self, newPipe : Pipe):
        """
        Assembles a new Pipe onto the Pipeline.

        :param newPipe: The new Pipe to be added to the Pipeline.
        :type newPipe: Pipe
        :return: A new Pipeline with the new Pipe assembled.
        :rtype: Pipeline
        """
          
        pipe = Pipe()
        pipe.process = lambda **kwargs : newPipe.process(**self.lastPipe.process(**kwargs))
        
        # now pipe will be the lastPipe for the returned Pipeline. i.e. we are returning an extended Pipeline, with pipe as the last assembled pipeline containing the new call stack, with newPipe.process() on top
        return Pipeline(pipe) 

    def start(self, **kwargs):
        """
        Starts the Pipeline processing.

        :param kwargs: Keyword arguments to be passed to the first Pipe.
        :return: The result of processing through the Pipeline.
        """

        return self.lastPipe.process(**kwargs)
