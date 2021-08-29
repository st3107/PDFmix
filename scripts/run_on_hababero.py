from pdfmix.core import create_mixture_pdf_files_from_cif_directory
import math


create_mixture_pdf_files_from_cif_directory(
    "../notebooks/data",
    "../notebooks/cifs",
    output_pattern="{:08d}.nc",
    fracs=[
        [0.0, 1.0],
        [1.0, 0.0],
        [0.05, 0.95],
        [0.95, 0.05],
        [0.2, 0.8],
        [0.8, 0.2],
        [0.3, 0.7],
        [0.7, 0.3],
        [0.4, 0.6],
        [0.6, 0.4],
        [0.5, 0.5]
    ],
    rmin=[1.5],
    rmax=[30.0],
    rstep=[math.pi / 23.0],
    iso_adp=[0.008, 0.032],
    ncpu=24
)
