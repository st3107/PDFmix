import copy
import time

import tqdm
from pkg_resources import resource_filename

import pdfmix.core as core

NI_CIF = resource_filename("pdfmix", "data/Ni.cif")


def test_nCr():
    assert 1 is core.nCr(3, 3)
    assert 10 is core.nCr(5, 2)
