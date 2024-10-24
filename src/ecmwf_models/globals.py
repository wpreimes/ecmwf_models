import os
from pathlib import Path


IMG_FNAME_TEMPLATE = "{product}_{type}_{datetime}.{ext}"
IMG_FNAME_DATETIME_FORMAT = "%Y%m%d_%H%M"

# ERA5 products supported by the reader.
SUPPORTED_PRODUCTS = ['era5', 'era5-land']

CDS_API_URL = "https://cds.climate.copernicus.eu/api"

# CDSAPI_RC variable must be set or we use home dir
DOTRC = os.environ.get('CDSAPI_RC', os.path.join(Path.home(), '.cdsapirc'))

EXPVER = {'0005': 'T', '0001': ''}


class CdoNotFoundError(ModuleNotFoundError):
    def __init__(self, msg=None):
        _default_msg = ("cdo and/or python-cdo not installed. "
                        "Pleas run `conda install -c conda-forge cdo` and also "
                        "`pip install cdo`.")
        self.msg = _default_msg if msg is None else msg