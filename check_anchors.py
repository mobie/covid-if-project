from glob import glob

import numpy as np
import pandas as pd
import z5py


def get_anchor(anchors, seg_id):
    anchor = (
        int(np.round(anchors["anchor_y"].loc[anchors["label_id"] == seg_id].values[0])),
        int(np.round(anchors["anchor_x"].loc[anchors["label_id"] == seg_id].values[0])),
    )
    return anchor


def check_anchors_visually(im, seg, anchors, n_cells=5):
    import napari
    seg_ids = np.unique(seg)[1:]
    seg_ids = np.random.choice(seg_ids, n_cells, replace=False)
    for seg_id in seg_ids:
        anchor = get_anchor(anchors, seg_id)
        v = napari.Viewer()
        v.add_image(im)
        v.add_labels(seg == seg_id)
        v.add_points([anchor])
        napari.run()


def validate_anchors(seg, anchors):
    seg_ids = np.unique(seg)[1:]
    assert len(anchors) == len(seg_ids)
    assert np.allclose(seg_ids, anchors["label_id"].values)
    for seg_id in seg_ids:
        anchor = get_anchor(anchors, seg_id)
        anchor_id = seg[anchor]
        assert anchor_id == seg_id, f"{seg_id}, {anchor_id}, {anchor}"
    print("Passed")


def check_anchors(name):
    im_path = f"./data/20200406_164555_328/images/ome-zarr/nuclei_{name}.ome.zarr"
    seg_path = f"./data/20200406_164555_328/images/ome-zarr/cell_segmentation_{name}.ome.zarr"
    table_path = f"./data/20200406_164555_328/tables/cell_segmentation_{name}/default.tsv"

    with z5py.File(seg_path, "r") as f:
        seg = f["s0"][:]
    with z5py.File(im_path, "r") as f:
        im = f["s0"][:]
    anchors = pd.read_csv(table_path, sep="\t")[["label_id", "anchor_y", "anchor_x"]]
    check_anchors_visually(im, seg, anchors)
    validate_anchors(seg, anchors)


def check_all_anchors():
    paths = glob("./data/20200406_164555_328/images/ome-zarr/cell_segmentation_*.ome.zarr")
    table_paths = glob("./data/20200406_164555_328/tables/cell_segmentation_*/default.tsv")
    assert len(paths) == len(table_paths)
    paths.sort()
    table_paths.sort()
    print("Checking anchors for", len(paths), "segmentations")
    for path, tab in zip(paths, table_paths):
        print("Check for", path)
        with z5py.File(path, "r") as f:
            seg = f["s0"][:]
        anchors = pd.read_csv(tab, sep="\t")[["label_id", "anchor_y", "anchor_x"]]
        validate_anchors(seg, anchors)
    print("All anchors are valid")


if __name__ == "__main__":
    # name = "WellE06_0002"
    # check_anchors(name)
    check_all_anchors()
