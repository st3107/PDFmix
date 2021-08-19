import time
import typing
import math
import itertools as it
import copy
import collections
import tqdm
from pprint import pformat
from pathlib import PurePath, Path
from configparser import ConfigParser
from pkg_resources import resource_filename

import fire
import yaml
import xarray as xr
import numpy as np
from pyobjcryst.crystal import Crystal, CreateCrystalFromCIF
from pyobjcryst.scatteringpower import ScatteringPower
from diffpy.srreal.pdfcalculator import PDFCalculator
from diffpy.srfit.pdf import PDFGenerator
from diffpy.srfit.structure.basestructureparset import BaseStructureParSet
from diffpy.srfit.fitbase.parameterset import ParameterSet

VarDict = typing.Dict[str, typing.Tuple[typing.List[str], typing.Any, typing.Dict[str, typing.Any]]]
ConfigDict = typing.Dict[
    str, typing.Union[str, int, float, list, typing.List[str], typing.List[int], typing.List[float], typing.List[
        list]]
]
NAME2ATTRS = {}


class PDFMixError(Exception):
    pass


class CalculatorSetting(dict):

    def to_var_dict(self) -> VarDict:
        dct = {}
        for k, v in self.items():
            dct[k] = ([k], [v], NAME2ATTRS.get(k, {}))
        return dct


class StructureSetting(dict):

    def to_var_dict(self) -> VarDict:
        dct = {}
        for k, v in self.items():
            dct[k] = ([k], [v], NAME2ATTRS.get(k, {}))
        return dct


def get_CalculatorSettings(
        keys: typing.Iterable[str],
        dct: typing.Union[dict, collections.ChainMap]
) -> typing.List:
    grids = [dct[key] for key in keys]
    coords = np.column_stack([c.flatten() for c in np.meshgrid(*grids)])
    return [CalculatorSetting(zip(keys, c)) for c in coords]


def get_StructureSettings(
        keys: typing.Iterable[str],
        dct: typing.Union[dict, collections.ChainMap]
) -> typing.List:
    grids = [dct[key] for key in keys]
    coords = np.column_stack([c.flatten() for c in np.meshgrid(*grids)])
    return [StructureSetting(zip(keys, c)) for c in coords]


def load_yaml(filename: str) -> dict:
    with Path(filename).open("r") as f:
        dct = yaml.safe_load(f)
    return dct


def dump_yaml(filename: str, dct: dict) -> None:
    with Path(filename).open("w") as f:
        yaml.safe_dump(dct, f, width=80, indent=2, default_flow_style=False)
    return


class PDFMixConfigParser:

    def __init__(self, dct: dict = None):
        if dct is None:
            dct = {}
        # default configuration
        self._calc_config: ConfigDict = {
            "rmin": [0.0],
            "rmax": [20.0],
            "rstep": [0.01],
            "qmin": [0.0],
            "qmax": [23.0],
            "qdamp": [0.0],
            "qbroad": [0.0],
            "delta1": [0.0],
            "delta2": [0.0]
        }
        self._stru_config: ConfigDict = {
            "iso_adp": [0.02],
            "lat_scale": [1.0]
        }
        self._fracs_config: ConfigDict = {
            "fracs": [[1.0]]
        }
        self._other_config: ConfigDict = {
            "verbose": 1
        }
        # attributes to hold the data
        self._dct = dict()
        self.calc_settings: typing.List[CalculatorSetting] = []
        self.stru_settings: typing.List[StructureSetting] = []
        self.frac_combs: typing.List[typing.List[float]] = []
        self.verbose: int = 0
        self.n_phase: int = 0
        # populate the attributes
        self.read_dict(dct)

    def read_dict(self, dct: dict) -> None:
        dct = dict(
            collections.ChainMap(
                dct,
                self._calc_config,
                self._stru_config,
                self._fracs_config,
                self._other_config
            )
        )
        self._dct = dct
        # calculator settings
        self.calc_settings = get_CalculatorSettings(self._calc_config.keys(), dct)
        # structure settings
        self.stru_settings = get_StructureSettings(self._stru_config.keys(), dct)
        # phases and fractions
        self.frac_combs = dct["fracs"]
        if len(self.frac_combs) == 0:
            raise PDFMixError("There is no fraction list. Please add fraction config 'fracs'.")
        self.n_phase = len(self.frac_combs[0])
        # set verbose
        self.verbose = dct["verbose"]

    def read(self, filename: str = None, **kwargs) -> None:
        dct = load_yaml(filename) if filename else {}
        dct.update(kwargs)
        self.read_dict(dct)

    def show(self) -> None:
        print(pformat(self._dct))
        return

    def write(self, filename: str) -> None:
        dump_yaml(filename, self._dct)
        return


def load_config(config_file: str, **kwargs) -> PDFMixConfigParser:
    config = PDFMixConfigParser()
    config.read(config_file, **kwargs)
    return config


def find_all_files(
        directory: str,
        pattern: str
) -> typing.List[str]:
    _directory = Path(directory)
    fs = []
    for f in _directory.rglob(pattern):
        if f.is_file():
            fs.append(str(f))
    return fs


def load_crystal(filename: str) -> Crystal:
    with Path(filename).open("rb") as fb:
        c = CreateCrystalFromCIF(fb)
        # record the input cif in the cif_str
        c.cif_str = fb.read()
    return c


def serialize_crystals(
        crystals: typing.List[Crystal]
) -> typing.List[str]:
    return [c.cif_str for c in crystals]


def nCr(n: int, r: int) -> int:
    f = math.factorial
    return f(n) // f(r) // f(n - r)


def gen_file_combs_from_directory(
        config: PDFMixConfigParser,
        files: typing.List[str]
) -> typing.Generator[typing.List[str], None, None]:
    n = config.n_phase
    for fc in it.combinations(files, n):
        yield list(fc)


def create_crystals(
        file_comb: typing.List[str]
) -> typing.List[Crystal]:
    return [load_crystal(f) for f in file_comb]


def store_in_dataset(
        r: np.ndarray,
        g: np.ndarray,
        crystals: typing.List[Crystal],
        fracs: typing.List[float],
        stru_setting: StructureSetting,
        calc_setting: CalculatorSetting
) -> xr.Dataset:
    n = len(fracs)
    m = len(crystals)
    if n != m:
        raise PDFMixError("Number of phases doesn't match the number of fractions: {} != {}.".format(n, m))
    phase = list(range(n))
    s_crystals = serialize_crystals(crystals)
    data = {
        "G": (["r"], g, NAME2ATTRS.get("G", {})),
        "fraction": (["phase"], fracs, NAME2ATTRS.get("fraction", {})),
        "structure": (["phase"], s_crystals, NAME2ATTRS.get("structure", {}))
    }
    coords = {
        "r": (["r"], r, NAME2ATTRS.get("r", {})),
        "phase": (["phase"], phase, NAME2ATTRS.get("phase", {}))
    }
    data.update(
        stru_setting.to_var_dict()
    )
    data.update(
        calc_setting.to_var_dict()
    )
    ds = xr.Dataset(data, coords=coords)
    return ds


def create_pg(
        crystal: Crystal,
        frac: float,
        stru_setting: StructureSetting,
        calc_setting: CalculatorSetting
) -> PDFGenerator:
    pg = PDFGenerator()
    pg.setStructure(crystal)
    pg.scale.setValue(frac)
    pg.delta1.setValue(calc_setting["delta1"])
    pg.delta2.setValue(calc_setting["delta2"])
    pg.setQmax(calc_setting["qmax"])
    pg.setQmin(calc_setting["qmin"])
    pg.qdamp.setValue(calc_setting["qdamp"])
    pg.qbroad.setValue(calc_setting["qbroad"])
    phase: BaseStructureParSet = pg.phase
    lat: ParameterSet = phase.getLattice()
    lat.a.setValue(lat.a.getValue() * stru_setting["lat_scale"])
    lat.b.setValue(lat.b.getValue() * stru_setting["lat_scale"])
    lat.c.setValue(lat.c.getValue() * stru_setting["lat_scale"])
    atoms = phase.getScatterers()
    for atom in atoms:
        atom.Biso.setValue(stru_setting["iso_adp"])
    return pg


def create_r(
        calc_setting: CalculatorSetting
) -> np.ndarray:
    return np.arange(
        calc_setting["rmin"],
        calc_setting["rmax"],
        calc_setting["rstep"]
    )


def calc_mixed_pdf(
        crystals: typing.List[Crystal],
        fracs: typing.List[float],
        stru_setting: StructureSetting,
        calc_setting: CalculatorSetting
) -> typing.Tuple[np.ndarray, np.ndarray]:
    r = create_r(calc_setting)
    g = np.zeros_like(r)
    for crystal, frac in zip(crystals, fracs):
        pg = create_pg(crystal, frac, stru_setting, calc_setting)
        g += pg(r)
    return r, g


def create_mixture_pdf(
        file_comb: typing.List[str],
        fracs: typing.List[float],
        stru_setting: StructureSetting,
        calc_setting: CalculatorSetting
) -> xr.Dataset:
    crystals = create_crystals(file_comb)
    r, g = calc_mixed_pdf(crystals, fracs, stru_setting, calc_setting)
    ds = store_in_dataset(r, g, crystals, fracs, stru_setting, calc_setting)
    return ds


def save_mixture_pdf(
        mixture_pdf: xr.Dataset,
        directory: str,
        pattern: str,
        count: int
) -> None:
    filename = pattern.format(count)
    filepath = PurePath(directory).joinpath(filename)
    mixture_pdf.to_netcdf(filepath)
    return


def create_progress_bar(files: typing.List[str], config: PDFMixConfigParser) -> tqdm.tqdm:
    nfs = len(files)
    nph = config.n_phase
    if nfs < nph:
        raise PDFMixError("Number of cif files is smaller than number of phases: {} < {}".format(nfs, nph))
    ncomb = nCr(nfs, nph)
    nf = len(config.frac_combs)
    ns = len(config.stru_settings)
    nc = len(config.calc_settings)
    counts = ncomb * nf * ns * nc
    verbose = config.verbose
    return tqdm.tqdm(total=counts, disable=(verbose == 0))


def create_mixture_pdf_files_from_cif_directory(
        output_directory: str,
        input_directory: str = r"./",
        config_file: str = None,
        output_pattern: str = r"{:016d}.nc",
        input_pattern: str = r"[!.]*.cif",
        **kwargs
) -> None:
    _output_directory = Path(output_directory).expanduser()
    _input_directory = Path(input_directory).expanduser()
    config = load_config(config_file, **kwargs)
    files = find_all_files(str(_input_directory), input_pattern)
    file_combs = gen_file_combs_from_directory(config, files)
    frac_combs = config.frac_combs
    stru_settings = config.stru_settings
    calc_settings = config.calc_settings
    _output_directory.mkdir(parents=True, exist_ok=True)
    pb = create_progress_bar(files, config)
    count = 0
    for file_comb in file_combs:
        for frac_comb in frac_combs:
            for stru_setting in stru_settings:
                for calc_setting in calc_settings:
                    pb.update()
                    mixture_pdf = create_mixture_pdf(file_comb, frac_comb, stru_setting, calc_setting)
                    save_mixture_pdf(mixture_pdf, str(_output_directory), output_pattern, count)
                    count += 1
    pb.close()
    return


def show_default_config() -> None:
    config = PDFMixConfigParser()
    config.show()
    return


def write_default_config(filename: str) -> None:
    config = PDFMixConfigParser()
    config.write(filename)
    return


def cli():
    fire.Fire(
        {
            "create": create_mixture_pdf_files_from_cif_directory,
            "show": show_default_config,
            "write": write_default_config
        }
    )
    return
