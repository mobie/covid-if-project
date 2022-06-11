import json
import os
from glob import glob
from subprocess import run

import mobie
import pandas as pd
from tqdm import tqdm

DS_FOLDER = "./data/20200406_164555_328"
SCALE = 2 * [0.65]
UNIT = "micrometer"


def write_scale(path):
    attrs_path = os.path.join(path, ".zattrs")
    with open(attrs_path, "r") as f:
        multiscales = json.load(f)["multiscales"][0]

    # write scale info
    datasets = multiscales["datasets"]
    new_datasets = []
    for ii, ds in enumerate(datasets):
        scale = [sc * 2**ii for sc in SCALE]
        ds["coordinateTransformations"] = [{"scale": scale, "type": "scale"}]
        new_datasets.append(ds)
    multiscales["datasets"] = new_datasets

    # write units
    new_axes = [
        {"name": "y", "type": "space", "unit": UNIT},
        {"name": "x", "type": "space", "unit": UNIT}
    ]
    multiscales["axes"] = new_axes

    # update version
    multiscales["version"] = "0.4"

    attrs = {"multiscales": [multiscales]}
    with open(attrs_path, "w") as f:
        json.dump(attrs, f)


def write_scale_info_local():
    metadata = mobie.metadata.read_dataset_metadata(DS_FOLDER)
    sources = metadata["sources"]
    for name, src in tqdm(sources.items(), total=len(sources)):
        src_type = list(src.keys())[0]
        local_path = src[src_type]["imageData"]["ome.zarr"]["relativePath"]
        path = os.path.join(DS_FOLDER, local_path)
        assert os.path.exists(path), path
        write_scale(path)


def update_attrs_s3():
    metadata = mobie.metadata.read_dataset_metadata(DS_FOLDER)
    sources = metadata["sources"]
    for name, src in tqdm(sources.items(), total=len(sources)):
        src_type = list(src.keys())[0]
        local_path = src[src_type]["imageData"]["ome.zarr"]["relativePath"]
        attrs_path = os.path.join(DS_FOLDER, local_path, ".zattrs")
        assert os.path.exists(attrs_path), attrs_path

        s3_addr = src[src_type]["imageData"]["ome.zarr.s3"]["s3Address"]
        s3_path = "/".join(s3_addr.split("/")[3:])
        s3_path = f"embl/{s3_path}/.zattrs"
        cmd = ["mc", "cp", attrs_path, s3_path]
        run(cmd)


def update_tables():
    update_cols = ["anchor_y", "anchor_x", "bb_min_y", "bb_min_x", "bb_max_y", "bb_max_x"]

    def _update(pattern):
        folders = glob(pattern)
        for folder in tqdm(folders, desc=f"Update tables with pattern {pattern}"):
            table_path = os.path.join(folder, "default.tsv")
            table = pd.read_csv(table_path, sep="\t")
            for col in update_cols:
                table[col] = table[col] * SCALE[0]
            table.to_csv(table_path, sep="\t", index=False)

    _update("./data/20200406_164555_328/tables/cell_segmentation_*")
    _update("./data/20200406_164555_328/tables/nucleus_segmentation_*")


def main():
    # write_scale_info_local()
    # update_attrs_s3()
    update_tables()


if __name__ == "__main__":
    main()
