__version__ = '0.0.3'
__version_info__ = tuple([int(num) for num in __version__.split('.')])

from liana import method as mt, plotting as pl, resource as rs

# done after everything has been imported (adapted from scanpy)
import sys

sys.modules.update({f'{__name__}.{m}': globals()[m] for m in ['mt', 'pl', 'rs']})

del sys
