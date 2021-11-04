import argparse
import json
import os
import string
from concurrent.futures import ThreadPoolExecutor
from functools import partial
from glob import glob

import mobie
import mobie.htm as htm

import numpy as np
import pandas as pd
from elf.io import open_file
from tqdm import tqdm


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
    file_format = "ome.zarr"

    print("Add dapi / nuclei")
    tmp_folder = os.path.join(tmp_root, "ims_nuclei")
    image_names = [f"nuclei_{parse_im_name(path)}" for path in input_files]
    htm.add_images(input_files, ROOT, plate_name, image_names,
                   key="nuclei/s0", file_format=file_format,
                   resolution=RESOLUTION, scale_factors=SCALE_FACTORS, chunks=CHUNKS,
                   tmp_folder=tmp_folder, target="local", max_jobs=24)

    print("Add serum-IgG")
    tmp_folder = os.path.join(tmp_root, "ims_igg")
    image_names = [f"serumIgG_{parse_im_name(path)}" for path in input_files]
    htm.add_images(input_files, ROOT, plate_name, image_names,
                   key="serum_IgG/s0", file_format=file_format,
                   resolution=RESOLUTION, scale_factors=SCALE_FACTORS, chunks=CHUNKS,
                   tmp_folder=tmp_folder, target="local", max_jobs=24)

    print("Add marker")
    tmp_folder = os.path.join(tmp_root, "ims_marker")
    image_names = [f"marker_tophat_{parse_im_name(path)}" for path in input_files]
    htm.add_images(input_files, ROOT, plate_name, image_names,
                   key="marker/s0", file_format=file_format,
                   resolution=RESOLUTION, scale_factors=SCALE_FACTORS, chunks=CHUNKS,
                   tmp_folder=tmp_folder, target="local", max_jobs=24)


def add_segmentations(input_files, plate_name):
    tmp_root = f'./tmp_{plate_name}'
    os.makedirs(tmp_root, exist_ok=True)
    file_format = "ome.zarr"

    print("Add cell segmentation")
    tmp_folder = os.path.join(tmp_root, "segs_cell")
    seg_names = [f"cell_segmentation_{parse_im_name(path)}" for path in input_files]
    htm.add_segmentations(input_files, ROOT, plate_name, seg_names,
                          key="cell_segmentation/s0", file_format=file_format,
                          resolution=RESOLUTION, scale_factors=SCALE_FACTORS, chunks=CHUNKS,
                          tmp_folder=tmp_folder, target="local", max_jobs=24)

    print("Add nucleus segmentation")
    tmp_folder = os.path.join(tmp_root, "segs_nuc")
    seg_names = [f"nucleus_segmentation_{parse_im_name(path)}" for path in input_files]
    htm.add_segmentations(input_files, ROOT, plate_name, seg_names,
                          key="nucleus_segmentation/s0", file_format=file_format,
                          resolution=RESOLUTION, scale_factors=SCALE_FACTORS, chunks=CHUNKS,
                          tmp_folder=tmp_folder, target="local", max_jobs=24)


def compute_clims(source, ds_folder, tmp_folder):
    tmp_path = os.path.join(tmp_folder, f'clims_{source}.json')
    if os.path.exists(tmp_path):
        with open(tmp_path) as f:
            return json.load(f)

    sources = mobie.metadata.read_dataset_metadata(ds_folder)["sources"]

    def compute_clim_im(source_name):
        path = os.path.join(
            ds_folder,
            sources[source_name]['image']['imageData']['ome.zarr']['relativePath']
        )
        with open_file(path, 'r') as f:
            data = f['s0'][:]
            cmin = np.percentile(data, 2)
            cmax = np.percentile(data, 98)
        return cmin, cmax

    source_names = [name for name in sources.keys() if name.startswith(source)]
    with ThreadPoolExecutor(32) as tp:
        results = list(tqdm(
            tp.map(compute_clim_im, source_names),
            total=len(source_names),
            desc=f"Compute contrast limits for {source}"
        ))

    cmin = np.median([res[0] for res in results])
    cmax = np.median([res[1] for res in results])
    clim = [float(cmin), float(cmax)]
    with open(tmp_path, 'w') as f:
        json.dump(clim, f)
    return clim


def format_tab(data):
    formatted = [[rr.decode("utf-8") for rr in row] for row in data]
    return np.array(formatted)


def create_site_table(ds_folder, table_file):
    table_path = os.path.join(ds_folder, "tables", "sites")
    os.makedirs(table_path, exist_ok=True)
    table_path = os.path.join(table_path, "default.tsv")

    rel_table_path = "tables/sites"
    if os.path.exists(table_path):
        return rel_table_path

    col_names_in = ["score", "IgG_is_outlier", "n_outlier_cells",
                    "n_infected", "n_control",
                    "background_IgG_median", "background_IgG_mad"]
    col_names_out = ["score", "is_outlier", "n_outlier_cells",
                     "n_infected_cells", "n_non_infected_cells",
                     "serum_background_intensity", "serum_background_variance"]
    assert len(col_names_in) == len(col_names_out)
    with open_file(table_file, "r") as f:
        g = f["tables/images/default"]
        cols = [col.decode("utf-8") for col in g["columns"][:]]
        data = g["cells"][:]

    # get the site names from the table and order by site names (order is not preserved in the table!)
    site_names = [name.decode("utf-8") for name in data[:, 1]]
    data = data[np.argsort(site_names)]
    site_names = [name.decode("utf-8").replace("-", "_") for name in data[:, 1]]

    col_ids = np.array([cols.index(col_name) for col_name in col_names_in])
    tab = data[:, col_ids]
    assert tab.shape == (len(data), len(col_ids))

    wells = np.array([site_name.split("_")[0] for site_name in site_names])
    col_names_out = ["annotation_id", "well"] + col_names_out
    tab = np.concatenate([
        np.array(site_names)[:, None],
        wells[:, None],
        format_tab(tab)
    ], axis=1)
    tab = pd.DataFrame(tab, columns=col_names_out)
    tab.to_csv(table_path, sep="\t", index=False, na_rep="nan")

    return rel_table_path


def create_well_table(ds_folder, table_file, well_names):
    table_path = os.path.join(ds_folder, "tables", "well")
    os.makedirs(table_path, exist_ok=True)
    table_path = os.path.join(table_path, "default.tsv")

    rel_table_path = "tables/well"
    if os.path.exists(table_path):
        return rel_table_path

    col_names_in = ["score", "IgG_is_outlier", "n_outlier_cells",
                    "n_infected", "n_control",
                    "background_IgG_median", "background_IgG_mad"]
    col_names_out = ["score", "is_outlier", "n_outlier_cells",
                     "n_infected_cells", "n_non_infected_cells",
                     "serum_background_intensity", "serum_background_variance"]
    assert len(col_names_in) == len(col_names_out)
    with open_file(table_file, "r") as f:
        g = f["tables/wells/default"]
        cols = [col.decode("utf-8") for col in g["columns"][:]]
        data = g["cells"][:]
    assert well_names == [name.decode("utf-8") for name in data[:, 0]]

    col_ids = np.array([cols.index(col_name) for col_name in col_names_in])
    tab = data[:, col_ids]
    assert tab.shape == (len(well_names), len(col_ids))

    col_names_out = ["annotation_id"] + col_names_out
    tab = np.concatenate([
        np.array(well_names)[:, None],
        format_tab(tab)
    ], axis=1)
    tab = pd.DataFrame(tab, columns=col_names_out)
    tab.to_csv(table_path, sep="\t", index=False, na_rep="nan")

    return rel_table_path


def get_all_wells(ds_meta):
    source_names = list(ds_meta["sources"].keys())
    wells = [name[name.find("Well"):].split("_")[0] for name in source_names]
    all_wells = list(set(wells))
    all_wells.sort()
    all_wells = [name[4:] for name in all_wells]
    return all_wells


def filter_function(source, well_names):
    return any(wname in source for wname in well_names)


def create_plate_view(view_name, plate_name, site_table, well_table,
                      wells=None, use_transform_grid=False,
                      add_source_tables=True, add_segmentations=True):
    ds_folder = f"./data/{plate_name}"
    tmp_folder = f"./tmp_{plate_name}"
    os.makedirs(tmp_folder, exist_ok=True)

    source_prefixes = ["nuclei", "serumIgG", "marker_tophat"]
    source_types = ["image", "image", "image"]
    source_settings = [
        {"color": "blue", "contrastLimits": compute_clims("nuclei", ds_folder, tmp_folder), "visible": True},
        {"color": "green", "contrastLimits": compute_clims("serumIgG", ds_folder, tmp_folder), "visible": False},
        {"color": "red", "contrastLimits": compute_clims("marker_tophat", ds_folder, tmp_folder), "visible": False},
    ]
    if add_segmentations:
        source_prefixes += ["cell_segmentation", "nucleus_segmentation"]
        source_types += ["segmentation", "segmentation"]
        source_settings += [
            {"lut": "glasbey", "tables": ["default.tsv"], "visible": False, "showTable": False},
            {"lut": "glasbey", "tables": ["default.tsv"], "visible": False, "showTable": False}
        ]

    def to_site_name(source_name, prefix):
        return source_name[len(f"{prefix}_Well"):]

    def to_well_name(site_name):
        return site_name.split("_")[0]

    def to_position(well_name):
        r, c = well_name[0], well_name[1:]
        r = string.ascii_uppercase.index(r)
        c = int(c) - 1
        return [r, c]

    name_filter = None if wells is None else partial(filter_function, well_names=wells)
    mobie.htm.add_plate_grid_view(ds_folder, view_name, menu_name="plate",
                                  source_prefixes=source_prefixes, source_types=source_types,
                                  source_settings=source_settings,
                                  source_name_to_site_name=to_site_name,
                                  site_name_to_well_name=to_well_name,
                                  well_to_position=to_position,
                                  site_table=site_table,
                                  well_table=well_table,
                                  name_filter=name_filter,
                                  sites_visible=False,
                                  add_annotation_displays=add_source_tables)


# add all view combinations
def add_all_views(plate_name, site_table, well_table):
    selected_wells = ["E06", "E07"]

    def to_name(wells, use_transform_grid, add_source_tables, add_segmentations):
        name = "full" if wells is None else "default"
        if use_transform_grid:
            name += "_transform_grid"
        if add_source_tables:
            name += "_with_tables"
        if add_segmentations:
            name += "_with_segmentations"
        return name

    for wells in (None, selected_wells):
        for use_transform_grid in (False, True):
            for add_source_tables in (False, True):
                for add_segmentations in (False, True):
                    view_name = to_name(wells, use_transform_grid, add_source_tables, add_segmentations)
                    create_plate_view(view_name, plate_name, site_table, well_table,
                                      wells=wells, use_transform_grid=use_transform_grid,
                                      add_source_tables=add_source_tables, add_segmentations=add_segmentations)


# this is the view we want to have in the end
def add_default_view(plate_name, site_table, well_table):
    selected_wells = ["E06", "E07"]
    create_plate_view("default", plate_name, site_table, well_table, wells=selected_wells)
    # TODO make this the default view
    create_plate_view("full_grid", plate_name, site_table, well_table)


def add_plate(plate_folder, all_views):
    plate_name = parse_plate_name(plate_folder)
    ds_folder = f"./data/{plate_name}"

    print("Adding plate", plate_name)
    input_files = glob(os.path.join(plate_folder, "*.h5"))
    input_files.sort()
    print("with", len(input_files), "images")

    # add images and segmentations
    add_image_data(input_files, plate_name)
    add_segmentations(input_files, plate_name)
    mobie.metadata.set_is2d(ds_folder, True)

    # create site and well tables
    table_file = os.path.join(plate_folder, f"{plate_name}_table.hdf5")
    assert os.path.exists(table_file)
    all_wells = get_all_wells(mobie.metadata.read_dataset_metadata(ds_folder))
    site_table = create_site_table(ds_folder, table_file)
    well_table = create_well_table(ds_folder, table_file, all_wells)

    # adding all grids for test purposes
    if all_views:
        add_all_views(plate_name, site_table, well_table)
    # add only the default grids
    else:
        add_default_view(plate_name, site_table, well_table)

    # validate the project
    print("Validating the project ...")
    mobie.validation.validate_project('./data')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", default=DEFAULT_PLATE)
    parser.add_argument("--all_views", "-a", default=0, type=int)
    args = parser.parse_args()
    add_plate(args.input, args.all_views)
