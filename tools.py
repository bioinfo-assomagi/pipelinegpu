import sys
import subprocess
import six

""" Wrapper classes for tools that will be used by the pipeline. """

class CustomCommand():

    command_name = None

    def __init__(self, *args, **kwargs) -> None:
        self.args = args
        self.kwargs = kwargs

    def _combine_arglist(self, args, kwargs): # in case the execution will contain additional arguments from initation. Not used.
        #print(f'Combining arguments ... [{self.args}] and [{args}], as well as [{self.kwargs} and {kwargs}]')

        _args = self.args + args
        if sys.version_info < (3, 9, 0):
            _kwargs = {**self.kwargs, **kwargs}
        else:
            _kwargs = self.kwargs | kwargs
        
        _kwargs.update(kwargs)

        return _args, _kwargs
    
    def _run_command(self, *args, **kwargs): # TODO: capture the output of commands), also to chain it to the next command
        
        stderr = kwargs.pop("stderr", None)
        stdout = kwargs.pop("stdout", None)
        capture_output = kwargs.pop("capture_output", False)

        if capture_output is True:
            stderr = subprocess.PIPE
            stdout = subprocess.PIPE

        stdin = kwargs.pop("stdin", None)
        input = kwargs.pop("input", None)

        if input: # we make the p.Communicate input argument a "property" of Popen.
            stdin = subprocess.PIPE
            if isinstance(input, six.string_types) and not input.endswith("\n"):
                # make sure that input is a simple string with \n line endings
                input = six.text_type(input) + "\n"
            # else:
            #     try:
            #         # make sure that input is a simple string with \n line endings
            #         input = "\n".join(map(six.text_type, input)) + "\n"
            #     except TypeError:
            #         # so maybe we are a file or something ... and hope for the best
            #         pass


        #print(input)
        cmd = self._commandline(*args, **kwargs)

        try:
            # print(cmd)
            p = subprocess.Popen(cmd, stdin=stdin, stderr=stderr, stdout=stdout, universal_newlines=True)
            out, err = p.communicate(input=input) 
        except:
            raise

        rc = p.returncode
        return (rc, out, err), p

    def buildProcess(self, *args, **kwargs):
        pass



    def _commandline(self, *args, **kwargs):
        return [self.command_name] + self.format_args(*args, **kwargs)
    
    def format_args(self, *args, **kwargs): # put the commandline arguments in a GNU or POSIX style
        options = []
        for option, value in kwargs.items():
            if not option.startswith("-"):
                option = "-{}".format(option)
            if value is True:
                options.append(option)
                continue
            elif value is False:
                raise ValueError(
                    "A False value is ambiguous for option {0!r}".format(option)
                )

            options.extend((option, str(value)))

        return options + list(args)


    def expose(self):
        return self._commandline(*self.args, **self.kwargs)
    
    def __call__(self, *args, **kwargs):
        _args, _kwargs = self._combine_arglist(args, kwargs)
        results, p = self._run_command(*_args, **_kwargs)
        return results
    

# TODO: add input handling from user interaction
""" Not used, not needed, to be deleted. """
class PopenWithInput(subprocess.Popen): 

    def __init__(self, *args, **kwargs):
        self.input = kwargs.pop('input', None)

        self.command = args[0]

        super(PopenWithInput, self).__init__(*args, **kwargs)

    def communicate(self, use_input=True) -> tuple:
        return super().communicate(self.input) if use_input else super().communicate()
    

""" I have to define this class here, since it deviates from the default CustomCommand, the fastq file path is defined without an option argument. 
While this special case can be handled in the CustomCommand class, better to keep things seperate, who knows what other tool would use a fastq file. """
class FastQC(CustomCommand):

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)

    def format_args(self, *args, **kwargs): # put the commandline arguments in a GNU or POSIX style, for FastQC, when it is a fastq file, it should treat is as an option
        options = []
        for option, value in kwargs.items():
            if option.endswith(".fastq"):
                options.append(option)
            else:
                if not option.startswith("-"):
                    option = "-{}".format(option)
                if value is True:
                    options.append(option)
                    continue
                elif value is False:
                    raise ValueError(
                        "A False value is ambiguous for option {0!r}".format(option)
                    )

                options.extend((option, str(value)))

        return options + list(args)

class ToolFactory():

    def __init__(self, tool_name, command_name, executable=None) -> None:
        self.tool_name = tool_name
        self.command_name = command_name
        self.executable = executable

    def get_tool(self):
        """ Here's the whole idea of the CommandWrapper; to wrap up all the argument preparation code and make it reusable
        for every tool. """
        tool_dict = {
            'command_name': self.command_name,
            'executable': self.executable
        }
        tool = type(self.tool_name, (CustomCommand,), tool_dict) # returns a Class. e.g. if u have a class as following: class FastQC: ..., tool will hold a reference to that
        return tool() # here, we instaitate that class, the same as saying return FastQC(). Then, we run the FastQC command, by calling its run method. It should look like FastQC()(), but since here we return an instatiated object FastQC(), to actually call and run the command we type only FastQC() in the code where we execute the command.


# Return FASTQ and FASTX tools as CustomCommands
def get_tools():
    import config
    import os

    tools = {}
    tools['FastQC'] = type("FastQC", (FastQC,), {'command_name': os.path.join(config.FASTQC, "fastqc"), 'executable': None})()
    tools['FastQualityFilter'] = ToolFactory('FastQualityFilter', os.path.join(config.FASTQX, "fastq_quality_filter"), None).get_tool()
    tools['FastXTrimmer'] = ToolFactory('FastXTrimmer', os.path.join(config.FASTQX, "fastx_trimmer"), None).get_tool()
    tools['BWA'] = ToolFactory('BWA', os.path.join(config.BWA, 'bwa'), None).get_tool()
    tools['ls'] = ToolFactory('Ls', 'ls', None).get_tool()
    tools['grep'] = ToolFactory("Grep", "grep", None).get_tool()
    tools['echo'] = ToolFactory("Echo", "echo", None).get_tool()

    return tools
    
globals().update(get_tools())

# class FastQC(CustomCommand):
     
#     def __init__(self, *args, **kwargs) -> None:
#         super().__init__(*args, **kwargs)

# class FastQualityFilter(CustomCommand):

#     def __init__(self, *args, **kwargs) -> None:
#         super().__init__(*args, **kwargs)

# class FastXTrimmer(CustomCommand):

#     def __init__(self, *args, **kwargs) -> None:
#         super().__init__(*args, **kwargs)