# The __init__.py file is loaded when the package is loaded.
# It is used to indicate that the directory in which it resides is a Python package
from importlib import metadata

from .filter_pathways import export_to_csv, filter_annotations, query_uniprot
from .uniprot_fields import FIELD_DISPLAY_NAMES, UNIPROT_FIELDS, fields_to_api_string

__version__ = metadata.version("filter_pathways_app")


# The __all__ variable is a list of variables which are imported
# when a user does "from example import *"
__all__ = [
    "query_uniprot",
    "filter_annotations",
    "export_to_csv",
    "UNIPROT_FIELDS",
    "FIELD_DISPLAY_NAMES",
    "fields_to_api_string",
]
