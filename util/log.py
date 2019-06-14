import inspect
import logging
import os
import select
import subprocess

FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

logger = logging.getLogger()
logging.basicConfig(level=logging.INFO, format=FORMAT)

def error(error):
    logger.error(error)

def warning(warning):
    logger.warning(warning)

def info(info):
    logger.info(info)

def debug(debug):
    logger.debug(debug)

def logging_call(popenargs, output_stream, error_function):

    """
    Call a subprocess and log the output to stderr to a logfile.

    Args:
        popenargs: The arguments to run in the subprocess
        output_stream: The stream to pipe stdout to

    Returns:
        A path to the SAM output file created.
    """
    process = subprocess.Popen(" ".join(popenargs), 
                               stdout=output_stream, 
                               stderr=subprocess.PIPE,
                               shell=True,
                               env=os.environ)

    # Get the caller of this function, 
    # so we can report where we're logging from.
    current_frame = inspect.currentframe()
    calling_frame = inspect.getouterframes(current_frame, 2)

    def check_io():
           while True:
                err = process.stderr.readline().decode()
                if err:
                    error_function("From " + calling_frame[1][3] + " - " + err.rstrip())
                else:
                    break

    # keep checking stdout/stderr until the child exits
    while process.poll() is None:
        check_io()