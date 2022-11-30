from errno import ENOENT as errno_ENOENT
from os import strerror as os_strerror
import pathlib


def non_empty_existing_path(file_or_dir):
    input_path = pathlib.Path(file_or_dir)
    if not input_path.is_file() and not input_path.is_dir():
        raise FileNotFoundError(errno_ENOENT, os_strerror(errno_ENOENT), file_or_dir)
    elif not input_path.is_dir() and not input_path.stat().st_size:
        raise OSError("File is empty.")
    return file_or_dir
