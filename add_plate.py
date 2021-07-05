import argparse
import os
from glob import glob
import mobie

ROOT = "./data"
DEFAULT_PLATE = "/g/kreshuk/data/covid/data-processed/20200406_164555_328"
RESOLUTION = [1., 1.]  # TODO what's the resolution in microns?
SCALE_FACTORS = 3 * [[2, 2]]
CHUNKS = (512, 512)


def parse_im_name(path):
    fname = os.path.split(path)[1]
    parts = fname.split('_')
    well = parts[0]
    im = parts[2]
    return f"{well}_{im}"


def parse_plate_name(plate_folder):
    return os.path.split(plate_folder)[1]


# TODO support ome.zarr data format
def add_image_data(input_files, plate_name):
    tmp_root = f'./tmp_{plate_name}'
    os.makedirs(tmp_root, exist_ok=True)

    file_format = "ome.zarr"
    for in_file in input_files:
        name = parse_im_name(in_file)

        # add nucleus channel
        nuc_key = "nuclei/s0"
        im_name = f"nuclei_{name}"
        mobie.add_image(in_file, nuc_key, ROOT, plate_name, im_name, RESOLUTION, SCALE_FACTORS, CHUNKS,
                        file_format=file_format, menu_name="images", target="local", max_jobs=8,
                        tmp_folder=os.path.join(tmp_root, f'tmp_{im_name}'))

        # also handle IgA and IgG channels?
        # add serum channel
        serum_key = "serum_IgG/s0"
        im_name = f"serumIgG_{name}"
        mobie.add_image(in_file, serum_key, ROOT, plate_name, im_name, RESOLUTION, SCALE_FACTORS, CHUNKS,
                        file_format=file_format, menu_name="images", target="local", max_jobs=8,
                        tmp_folder=os.path.join(tmp_root, f'tmp_{im_name}'))

        # add marker channel
        serum_key = "marker/s0"
        im_name = f"marker_tophat_{name}"
        mobie.add_image(in_file, serum_key, ROOT, plate_name, im_name, RESOLUTION, SCALE_FACTORS, CHUNKS,
                        file_format=file_format, menu_name="images", target="local", max_jobs=8,
                        tmp_folder=os.path.join(tmp_root, f'tmp_{im_name}'))


def add_plate(plate_folder):
    plate_name = parse_plate_name(plate_folder)
    print("Adding plate", plate_name)
    input_files = glob(os.path.join(plate_folder, "*.h5"))
    input_files.sort()
    print("with", len(input_files), "images")
    add_image_data(input_files, plate_name)

    # TODO
    # add segmentations
    # add tables
    # add the full grid view


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # add default argument for debugging
    parser.add_argument('--input', '-i', default=DEFAULT_PLATE)
    args = parser.parse_args()
    add_plate(args.input)
