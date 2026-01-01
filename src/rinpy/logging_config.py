# -*- coding: utf-8 -*-
"""
    __author__ = 'Zehra Sarica'
    __email__ = ['sarica16@itu.edu.tr','zehraacar559@gmail.com']
"""

import logging
import pathlib
import shutil

RINPY_DEFAULT_FOLDER = ".rinpy"
DEBUG_LOG = "debug.log"
PROCESS_LOG = "process.log"
LOG_TXT = "log.txt"
UTF_8_ENCODING = 'utf-8'


class InfoOnlyFilter(logging.Filter):
    """Filter only INFO messages to be shown in the log files."""

    def filter(self, record):
        return record.levelno == logging.INFO


def get_rinpy_log_dir():
    """Get the rinpy log directory path."""
    home_dir = pathlib.Path.home()
    rinpy_dir = home_dir / RINPY_DEFAULT_FOLDER
    rinpy_dir.mkdir(parents=True, exist_ok=True)
    return rinpy_dir


log_dir = get_rinpy_log_dir()

info_log_file = log_dir / PROCESS_LOG
all_log_file = log_dir / DEBUG_LOG

formatter = logging.Formatter('%(asctime)s.%(msecs)03d %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

info_handler = logging.FileHandler(info_log_file, encoding=UTF_8_ENCODING)
info_handler.setLevel(logging.INFO)
info_handler.addFilter(InfoOnlyFilter())
info_handler.setFormatter(formatter)

all_handler = logging.FileHandler(all_log_file, encoding=UTF_8_ENCODING)
all_handler.setLevel(logging.DEBUG)
all_handler.setFormatter(formatter)

stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.DEBUG)
stream_handler.setFormatter(formatter)

logging.basicConfig(
    level=logging.DEBUG,
    handlers=[info_handler, all_handler, stream_handler]
)

logging.getLogger("matplotlib").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)


def clear_logs():
    """Clear the content of log files for each run."""
    for file in [info_log_file, all_log_file]:
        try:
            if file.exists():
                with open(file, "r+", encoding=UTF_8_ENCODING) as f:
                    f.truncate()
        except Exception as e:
            logging.warning(f"Log file cannot be cleared: {file} ({e})")


def save_logs(output_dir: str):
    """ Save .rinpy/process.log to the given output directory as log.txt.
    If .rinpy/process.log does not exist, log a warning.
    """
    process_log_file = log_dir / PROCESS_LOG  # .rin/process.log path

    try:
        if not process_log_file.exists():
            logging.warning(f"Process log not found at: {process_log_file}")
            return

        destination = pathlib.Path(output_dir)
        destination.mkdir(parents=True, exist_ok=True)

        output_file = destination / LOG_TXT
        shutil.copyfile(process_log_file, output_file)

        logging.info(f"Process log saved to: {output_file}")

    except Exception as e:
        logging.error(f"Failed to save process log to {output_dir}: {e}")
