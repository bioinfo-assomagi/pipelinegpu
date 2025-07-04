import logging
import sys
from colorama import Fore, Style, init

# Init colorama on Windows (no‐op elsewhere)
init()

class ColorFormatter(logging.Formatter):
    # Map log levels to colors
    LEVEL_COLORS = {
        'DEBUG': Fore.CYAN,
        'INFO': Fore.GREEN,
        'WARNING': Fore.YELLOW,
        'ERROR': Fore.RED,
        'CRITICAL': Fore.MAGENTA,
    }

    THREAD_COLOR = Fore.CYAN
    MODULE_COLOR = Fore.MAGENTA
    PROC_COLOR   = Fore.CYAN
    #MODULE_COLOR = Fore.MAGENTA

    # how wide each column should be on screen
    COL_WIDTH = {
        "module": 65,     # 55 chars for “[module.path:line]”
        "sample": 25,     # 15 chars for “[Sample 123.456]”
    }

    def format(self, record):
        # Base timestamp + level part
        ts = self.formatTime(record, "%H:%M:%S")
        thread_name = record.threadName
        levelname = record.levelname.ljust(8)
        modulename = f"[{record.name}:{record.lineno}]".ljust(self.COL_WIDTH["module"])
        color = self.LEVEL_COLORS.get(record.levelname, Fore.WHITE)
        funcName = f"[{record.funcName}]".ljust(30) 

        # grab process name & pid instead of threadName
        pname  = record.processName      # e.g. "ForkPoolWorker-1"
        pid    = record.process          # e.g.  12345

        # Build colored prefix
        prefix = (
            f"{Fore.BLUE}{ts}{Style.RESET_ALL} "
            f"{color}{levelname}{Style.RESET_ALL} "
            f"{self.PROC_COLOR}[{pname}:{pid}]{Style.RESET_ALL} "
            f"{self.MODULE_COLOR}{modulename}{Style.RESET_ALL} "
            f"{color}{funcName}{Style.RESET_ALL}"
        )

        # The actual message stays uncolored (or you can choose another color)
        message = record.getMessage()
        return f"{prefix} {message}"

# Logger helper: just for log clarity, helps debugging in the future - for now hardcoded at the top of each script; TODO: generalize it so that all scripts
# theat report a sample name can use it
class SampleAlignedFormatter(ColorFormatter):
    """
    Same colours as CLI.ColorFormatter, but gives Sample a fixed-width column.
    """
    # def format(self, record):
    #     # If 'sample' is present, pad it to 12 chars
    #     if hasattr(record, "sample_id"):
    #         record.sample_id = f"{record.sample_id.ljust(self.COL_WIDTH['sample'])}"
    #     else:
    #         record.sample_id = ""          # keep format() happy
    #     return super().format(record)
    pass