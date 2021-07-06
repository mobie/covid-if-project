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


def add_image_data(input_files, plate_name):
    tmp_root = f'./tmp_{plate_name}'
    os.makedirs(tmp_root, exist_ok=True)

    ds_folder = os.path.join('./data', plate_name)
    ds_meta = mobie.metadata.read_dataset_metadata(ds_folder)
    sources = ds_meta.get('sources', {})

    file_format = "ome.zarr"
    for in_file in input_files:
        name = parse_im_name(in_file)

        # add nucleus channel
        nuc_key = "nuclei/s0"
        im_name = f"nuclei_{name}"
        if im_name not in sources:
            mobie.add_image(in_file, nuc_key, ROOT, plate_name, im_name, RESOLUTION, SCALE_FACTORS, CHUNKS,
                            file_format=file_format, menu_name="images", target="local", max_jobs=8,
                            tmp_folder=os.path.join(tmp_root, f'tmp_{im_name}'))

        # also handle IgA and IgG channels?
        # add serum channel
        serum_key = "serum_IgG/s0"
        im_name = f"serumIgG_{name}"
        if im_name not in sources:
            mobie.add_image(in_file, serum_key, ROOT, plate_name, im_name, RESOLUTION, SCALE_FACTORS, CHUNKS,
                            file_format=file_format, menu_name="images", target="local", max_jobs=8,
                            tmp_folder=os.path.join(tmp_root, f'tmp_{im_name}'))

        # add marker channel
        serum_key = "marker/s0"
        im_name = f"marker_tophat_{name}"
        if im_name not in sources:
            mobie.add_image(in_file, serum_key, ROOT, plate_name, im_name, RESOLUTION, SCALE_FACTORS, CHUNKS,
                            file_format=file_format, menu_name="images", target="local", max_jobs=8,
                            tmp_folder=os.path.join(tmp_root, f'tmp_{im_name}'))


def make_2d(plate_name):
    ds_folder = os.path.join('./data', plate_name)
    ds_meta = mobie.metadata.read_dataset_metadata(ds_folder)
    ds_meta['is2D'] = True
    mobie.metadata.write_dataset_metadata(ds_folder, ds_meta)


def remove_single_source_views(ds_meta, remove_prefixes):
    new_views = {k: v for k, v in ds_meta['views'].items() if not any(k.startswith(pref) for pref in remove_prefixes)}
    ds_meta['views'] = new_views
    return ds_meta


def create_raw_views(plate_name):
    ds_folder = os.path.join('./data', plate_name)
    ds_meta = mobie.metadata.read_dataset_metadata(ds_folder)

    # remove single source views
    ds_meta = remove_single_source_views(ds_meta, remove_prefixes=['marker', 'serumIgG', 'nuclei'])

    # TODO add grid views

    # TODO replace the default view

    mobie.metadata.write_dataset_metadata(ds_folder, ds_meta)


def add_plate(plate_folder):
    plate_name = parse_plate_name(plate_folder)
    print("Adding plate", plate_name)
    input_files = glob(os.path.join(plate_folder, "*.h5"))
    input_files.sort()
    print("with", len(input_files), "images")
    add_image_data(input_files, plate_name)

    make_2d(plate_name)

    create_raw_views(plate_name)

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
