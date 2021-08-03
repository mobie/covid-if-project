import os
from shutil import rmtree
from subprocess import run
import mobie


def clear_seg(ds_folder, source):
    source = source["segmentation"]
    data_path = os.path.join(ds_folder, source["imageData"]["ome.zarr"]["relativePath"])
    table_folder = os.path.join(ds_folder, source["tableData"]["tsv"]["relativePath"])
    rmtree(data_path)
    run(["git", "rm", "-r", table_folder])


def clear_segs():
    ds_folder = "./data/20200406_164555_328"
    metadata = mobie.metadata.read_dataset_metadata(ds_folder)
    sources = metadata["sources"]
    new_sources = {}

    for name, source in sources.items():
        if list(source.keys())[0] == "segmentation":
            print("Clearing", name)
            clear_seg(ds_folder, source)
            continue
        new_sources[name] = source

    metadata["sources"] = new_sources
    # also clear the view
    metadata["views"] = {}
    mobie.metadata.write_dataset_metadata(ds_folder, metadata)


if __name__ == '__main__':
    clear_segs()
