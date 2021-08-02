import argparse
import os
from glob import glob
import h5py

DEFAULT_PLATE = "/g/kreshuk/data/covid/data-processed/20200406_164555_328"


def check_plate(plate_folder):
    input_files = glob(os.path.join(plate_folder, "*.h5"))
    input_files.sort()
    ff = input_files[0]

    print("Check input plate @", plate_folder)
    with h5py.File(ff, "r") as f:
        keys = list(f.keys())
    print("Available image data:")
    print(keys)

    # TODO check the tables


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", default=DEFAULT_PLATE)
    args = parser.parse_args()
    check_plate(args.input)
