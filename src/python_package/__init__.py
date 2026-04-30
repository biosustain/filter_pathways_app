# The __init__.py file is loaded when the package is loaded.
# It is used to indicate that the directory in which it resides is a Python package
from importlib import metadata

__version__ = metadata.version("python_package")

from .mockup import hello_world, saved_world
from .filter_pathways import query_uniprot, filter_annotations, export_to_csv
from .uniprot_fields import UNIPROT_FIELDS, FIELD_DISPLAY_NAMES, fields_to_api_string

# The __all__ variable is a list of variables which are imported
# when a user does "from example import *"
__all__ = [
    "hello_world",
    "saved_world",
    "query_uniprot",
    "filter_annotations",
    "export_to_csv",
    "UNIPROT_FIELDS",
    "FIELD_DISPLAY_NAMES",
    "fields_to_api_string",
]
