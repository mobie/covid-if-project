import os
from glob import glob
from mobie.htm.data_import import _add_tables


def fix_tables():
    ds_folder = "./data/20200406_164555_328"
    seg_paths = glob(os.path.join(ds_folder, "images", "ome-zarr", "cell_segmentation*"))
    seg_paths.sort()
    seg_names = [
        os.path.split(path)[1].rstrip(".ome.zarr") for path in seg_paths
    ]
    res = (1.0, 1.0)
    _add_tables(
        file_format="ome.zarr",
        paths=seg_paths,
        segmentation_names=seg_names,
        resolution=res,
        ds_folder=ds_folder,
        tmp_folder="./tmp_fix_tables",
        target="local",
        max_jobs=24
    )


if __name__ == "__main__":
    fix_tables()
