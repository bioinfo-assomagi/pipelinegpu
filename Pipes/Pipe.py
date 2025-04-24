

class Pipe():

    """
    The `Pipe` class serves as an abstract base class for pipeline components.
    It defines the interface for processing data and is intended to be subclassed
    by concrete implementations. 
    
    A `Pipe` object can be executed individually, but assembling it into a `Pipeline` (even of one `Pipe`) is preferred.

    **Usage:**

    Subclass the `Pipe` class and implement the `process` method to define
    specific processing behavior.

    **Example:**

    .. code-block:: python

        class Pipe:
            def __init__(self):
                pass

            def process(self, **kwargs):
                raise NotImplementedError

        class PipeA(Pipe):
            def process(self, **kwargs):
                # Implement processing logic
                kwargs['a'] = 'Processed by PipeA'
                return kwargs

        class PipeB(Pipe):
            def process(self, **kwargs):
                # Implement additional processing logic
                kwargs['b'] = 'Processed by PipeB'
                return kwargs

        # Create instances of the pipes
        pipe_a = PipeA()
        pipe_b = PipeB()

        # Process data through PipeA
        result_a = pipe_a.process(data='Initial data')
        print(result_a)
        # Output: {'data': 'Initial data', 'a': 'Processed by PipeA'}

        # Process data through PipeB using the output of PipeA
        result_b = pipe_b.process(**result_a)
        print(result_b)
        # Output: {'data': 'Initial data', 'a': 'Processed by PipeA', 'b': 'Processed by PipeB'}

    In this example:

    - **Abstract Base Class:** `Pipe` defines the `process` method as an interface.
    - **Subclassing:** `PipeA` and `PipeB` implement the `process` method with specific logic.
    - **Chaining:** The output of `PipeA` is used as the input for `PipeB`, demonstrating how pipes can be chained together.
    """

    def __init__(self) -> None:
        pass

    def process(self, **kwargs):
        """
        Processes the input data and returns the result.

        :param kwargs: Keyword arguments containing the data to process.
        :return: The processed data.
        :raises NotImplementedError: If the method is not implemented in a subclass.
        """
        raise NotImplementedError("Subclasses must implement the 'process' method.")