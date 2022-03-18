import logging
import cdsapi
import typing as tp

import numpy as np


class CDSStatusTracker:
    """
    Track the status of the CDS download by using the CDS callback functions
    """

    statuscode_ok = 0
    statuscode_error = -1
    statuscode_unavailable = 10

    def __init__(self):
        self.download_statuscode = self.statuscode_ok

    def handle_error_function(self, *args, **kwargs):
        message_prefix = args[0]
        message_body = args[1]
        if self.download_statuscode != self.statuscode_unavailable:
            if (message_prefix.startswith("Reason:") and
                    message_body == "Request returned no data"):
                self.download_statuscode = self.statuscode_unavailable
            else:
                self.download_statuscode = self.statuscode_error
        logging.error(*args, **kwargs)


class CDSDataRequest:
    def __init__(self,
                 variables: tp.Sequence[str],
                 year: tp.Sequence[int],
                 month: tp.Sequence[int],
                 days: tp.Sequence[int],
                 Client: cdsapi.Client = None,
                 hours: tp.Optional[tp.Tuple[int]] = range(1, 13),
                 product: tp.Optional[tp.Literal['era5', 'era5_land']]='era5',
                 dry_run: bool=False
        ):
        """
        Wrapper around request expected by CDS API.

        Parameters
        ----------
        Client: cdsapi.Client or None
            Client (with authentification) to handle requests
            Passing None will not raise any errors, but download wont work.
        variables: Lisŧ[str, ...]
            Variable names to download.
        year: int or List[int,...]
            Year in which `month` is.
        month: int or List[int,...]
            Month in which `days` are
        days: List[int, ...]
            Days for which `hours` are downloaded
        hours: List[int, ...]
            Hours for which images are downloaded
        product: str, optional (default: era5)
            Product to request data for
            era5 or era5_land
        """
    request = {
        "format": "grib",
        "variable": np.atleast_1d(variables).tolist(),
        "year": [str(y) for year in np.atleast_1d(year)],
        "month": [str(m).zfill(2) for m in np.atleast_1d(month)],
        "day": [str(d).zfill(2) for d in days],
        "time": [time(h, 0).strftime("%H:%M") for h in h_steps],
    }
    if bbox is not None:
        if (bbox[0] > bbox[2]) or (bbox[1] > bbox[3]):
            raise ValueError(f"Invalid bounding box passed: {bbox} "
                             f"Expected order: "
                             f"min_lat, min_lon, max_lat, max_lon")
        request["bbox"] = bbox
    if grid_size is not None:
        request["grid"] = [grid_size, grid_size]

    request.update(cds_kwds)
    result = {'success': False}
    if product == "era5":
        request["product_type"] = "reanalysis"
        if not dry_run:
            api_result = c.retrieve("reanalysis-era5-single-levels", request,
                                    target)
            result['success'] = True
    elif product == "era5-land":
        if not dry_run:
            api_result = c.retrieve("reanalysis-era5-land", request, target)
            result['success'] = True
    else:
        raise ValueError(
            product, "Unknown product, choose either 'era5' or 'era5-land'")

    if dry_run:
        result['success'] = True