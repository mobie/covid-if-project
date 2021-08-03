import os
import mobie
import napari
from elf.io import open_file


def check_pos_mobie(ds_folder, pos_sources, sources):
    with napari.gui_qt():
        v = napari.Viewer()

        for source_name in pos_sources:
            print("Loading source", source_name)
            source_type = list(sources[source_name].keys())[0]
            path = sources[source_name][source_type]["imageData"]["ome.zarr"]["relativePath"]
            path = os.path.join(ds_folder, path)
            assert os.path.exists(path)

            with open_file(path, "r") as f:
                data = f["s0"][:]

            if source_type == "image":
                v.add_image(data, name=source_name)
            else:
                v.add_labels(data, name=source_name)


def check_pos_original(well, site):
    pass


def check_site_mobie(well_name, pos):
    ds_folder = "./data/20200406_164555_328"
    metadata = mobie.metadata.read_dataset_metadata(ds_folder)

    sources = metadata["sources"]
    view = metadata["views"]["default"]
    grid_trafo = view["sourceTransforms"][0]["grid"]
    well_sources = grid_trafo["sources"]
    pos_sources = well_sources[f"{well_name}_{pos}"]

    check_pos_mobie(ds_folder, pos_sources, sources)


def check_site_original(well_name, pos):
    root = "/g/kreshuk/data/covid/data-processed"
    name = f"20200406_164555_328/Well{well_name}_Point{well_name}_{pos}_ChannelDAPI,WF_GFP,TRITC_Seq0477.h5"
    path = os.path.join(root, name)

    with open_file(path, "r") as f:
        serum = f["serum_IgG/s0"][:]
        nuclei = f["nuclei/s0"][:]
        cell_seg = f["cell_segmentation/s0"][:]
        nuc_seg = f["nucleus_segmentation/s0"][:]

    with napari.gui_qt():
        v = napari.Viewer()
        v.add_image(serum)
        v.add_image(nuclei)
        v.add_labels(cell_seg)
        v.add_labels(nuc_seg)


if __name__ == "__main__":
    # check_site_mobie("E06", "0000")
    check_site_original("E06", "0000")
