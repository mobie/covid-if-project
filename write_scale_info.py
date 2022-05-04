import os
import mobie
import zarr

DS_FOLDER = "./data/20200406_164555_328"
SCALE = 2 * [0.65]
UNIT = "micrometer"


def write_scale(path):
    with zarr.open(path, "a") as f:
        attrs = f.attrs


def write_scale_info_local():
    metadata = mobie.metadata.read_dataset_metadata(DS_FOLDER)
    sources = metadata["sources"]
    for name, src in sources.items():
        src_type = list(src.keys())[0]
        local_path = src[src_type]["imageData"]["ome.zarr"]["relativePath"]
        path = os.path.join(DS_FOLDER, local_path)
        assert os.path.exists(path), path
        write_scale(path)


def main():
    write_scale_info_local()


if __name__ == "__main__":
    main()
