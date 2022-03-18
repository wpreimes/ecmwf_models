
class CdoNotFoundError(ModuleNotFoundError):
    def __init__(self, msg=None):
        _default_msg = "cdo and/or python-cdo not installed. " \
                       "Use conda to install it them under Linux."
        self.msg = _default_msg if msg is None else msg

_expver_lut = {1: '', 5: 'T'}
