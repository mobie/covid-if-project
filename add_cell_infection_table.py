import os
from glob import glob

import numpy as np
import h5py
import pandas as pd

ROOT = "./data"
DEFAULT_PLATE = "/g/kreshuk/data/covid/data-processed/20200406_164555_328"


def parse_im_name(path):
    fname = os.path.split(path)[1]
    parts = fname.split('_')
    well = parts[0]
    im = parts[2]
    return f"{well}_{im}"


def decode_col(col):
    return np.array([int(v.decode()[:-2]) for v in col])


def add_infection_cols():
    files = glob(os.path.join(DEFAULT_PLATE, "*.h5"))
    for ff in files:
        name = parse_im_name(ff)
        table_path = os.path.join(f"./data/20200406_164555_328/tables/cell_segmentation_{name}/default.tsv")
        assert os.path.exists(table_path), table_path
        table = pd.read_csv(table_path, sep="\t")
        with h5py.File(ff, "r") as f:
            values = f["/tables/cell_classification/cell_segmentation/marker_tophat/cells"][:]
            label_id = decode_col(values[:, 0])
            infected = decode_col(values[:, 1])
            control = decode_col(values[:, 2])

        if label_id[0] == 0:
            label_id = label_id[1:]
            infected = infected[1:]
            control = control[1:]
        assert np.allclose(label_id, table.label_id.values.astype("int"))
        table["cell"] = ["infected" if infect else "control" for infect in infected]
        table.to_csv(table_path, sep="\t", index=False)


add_infection_cols()
