from pathlib import Path

import numpy as np
import xarray as xr
from pkg_resources import resource_filename

import pdfmix.core as core

CIF_DIR = Path(resource_filename("pdfmix", "data/"))


def test_nCr():
    assert 1 == core.nCr(3, 3)
    assert 10 == core.nCr(5, 2)


def test_show_default_config():
    print()
    core.show_default_config()


def test_write_default_config(tmp_path: Path):
    filename = "test_default_config.yaml"
    filepath = tmp_path.joinpath(filename)
    core.write_default_config(str(filepath))
    print()
    print(filepath.read_text())


def test_create_mixture_pdf_files_from_cif_directory(tmp_path: Path):
    # create temp output directory
    temp_dir = tmp_path.joinpath("temp_output_ncs")
    temp_dir.mkdir()
    temp_config_file = tmp_path.joinpath("temp_config.yaml")
    # create config file
    config = core.PDFMixConfigParser()
    qmaxs = [20.0, 30.0]
    fracs = [[0., 1.], [1., 0.]]
    config.read_dict({"qmax": qmaxs, "fracs": fracs, "rmax": [10.0], "verbose": 0})
    config.write(str(temp_config_file))
    # run functions
    core.create_mixture_pdf_files_from_cif_directory(
        str(temp_dir), str(CIF_DIR), config_file=str(temp_config_file)
    )
    # check output files
    assert len(list(temp_dir.glob("*.nc"))) == len(qmaxs) * len(fracs)


def test_r_range_creation(tmp_path: Path):
    # create temp output directory
    temp_dir = tmp_path.joinpath("temp_output_ncs")
    temp_dir.mkdir()
    temp_config_file = tmp_path.joinpath("temp_config.yaml")
    # create config file
    config = core.PDFMixConfigParser()
    rmax, rmin, rstep = 30.0, 1.5, 0.1
    config.read_dict(
        {
            "fracs": [[0.5, 0.5]],
            "rmax": [rmax],
            "rstep": [rstep],
            "rmin": [rmin],
            "qmin": [0.5],
            "qmax": [23.0, 25.0],
            "qdamp": [0.04],
            "qbroad": [0.01],
            "iso_adp": [0.008]
        }
    )
    config.write(str(temp_config_file))
    # run functions
    core.create_mixture_pdf_files_from_cif_directory(
        str(temp_dir), str(CIF_DIR), config_file=str(temp_config_file)
    )
    # check output files
    fs = list(temp_dir.glob("*.nc"))
    assert len(fs) > 0
    r = np.arange(rmin, rmax, rstep)
    ds0 = xr.load_dataset(str(fs[0]))
    assert ds0["G"].shape == r.shape
    for f in fs[1:]:
        ds = xr.load_dataset(str(f))
        assert ds["G"].shape == r.shape
        assert not ds0["G"].equals(ds["G"])
