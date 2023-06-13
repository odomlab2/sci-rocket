import logging
from rich.console import Console
from rich.logging import RichHandler

def init_logger(verbose: bool = False):
    """
    Initial a Logger with rich handler.

    Returns:
        logging.Logger: logger
    """
    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)

    console = Console(force_terminal=True)
    ch = RichHandler(show_path=False, console=console, show_time=True)
    formatter = logging.Formatter("snakemake-sciseq: %(message)s")
    ch.setFormatter(formatter)
    log.addHandler(ch)
    log.propagate = False

    if verbose:
        log.setLevel(logging.DEBUG)

    return log