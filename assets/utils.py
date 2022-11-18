from errno import ENOENT as errno_ENOENT
from os import strerror as os_strerror
import pathlib


def non_empty_existing_file(file):
    path_file = pathlib.Path(file)
    if not path_file.is_file() and not path_file.is_dir():
        raise FileNotFoundError(errno_ENOENT, os_strerror(errno_ENOENT), file)
    elif not path_file.is_dir() and not path_file.stat().st_size:
        raise OSError("File is empty.")
    return file
