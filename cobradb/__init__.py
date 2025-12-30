from typing import Iterable
from .cobradb_settings import settings
import os

__version__ = "3.0.0"

_ROOT = os.path.abspath(os.path.dirname(__file__))


def get_data(*args):
    return os.path.join(_ROOT, "data", *args)
