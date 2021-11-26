import argparse
import napari
import zarr


def view_segment(seg_path, label_id):
    nuc_path = seg_path.replace("cell_segmentation", "nuclei")
    with zarr.open(nuc_path, "r") as f:
        nuc = f["s0"][:]
    with zarr.open(seg_path, "r") as f:
        seg = f["s0"][:]
    v = napari.Viewer()
    v.add_image(nuc, name="nuclei")
    v.add_labels(seg == label_id, name=f"label-mask-{label_id}")
    v.add_labels(seg, name="cell-segmentation")
    napari.run()


def main():
    default_well = "/g/kreshuk/pape/Work/mobie/covid-if-project/data/20200406_164555_328/images/ome-zarr/cell_segmentation_WellE06_0000.ome.zarr"
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--path", default=default_well)
    parser.add_argument("-l", "--label", default=132)
    args = parser.parse_args()
    view_segment(args.path, args.label)


if __name__ == "__main__":
    main()
