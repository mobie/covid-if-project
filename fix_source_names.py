import os
import mobie
import z5py


def fix_source_names(ds_folder):
    metadata = mobie.metadata.read_dataset_metadata(ds_folder)
    sources = metadata["sources"]
    for name, source in sources.items():
        storage = source["image"]["imageData"]["ome.zarr"]["relativePath"]
        path = os.path.join(ds_folder, storage)
        assert os.path.exists(path)
        with z5py.File(path, mode="a") as f:
            mscales= f.attrs["multiscales"]
            mscales[0]["name"] = name
            f.attrs["multiscales"] = mscales
        with z5py.File(path, mode="r") as f:
            assert name == f.attrs["multiscales"][0]["name"]


if __name__ == '__main__':
    fix_source_names('./data/20200406_164555_328')
