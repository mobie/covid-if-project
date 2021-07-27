import argparse
import json
import os
import string
from concurrent.futures import ThreadPoolExecutor
from glob import glob

import mobie
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


def compute_clims(source, this_source_names, sources, tmp_folder, ds_folder):
    tmp_path = os.path.join(tmp_folder, f'clims_{source}.json')
    if os.path.exists(tmp_path):
        with open(tmp_path) as f:
            return json.load(f)

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

    with ThreadPoolExecutor(32) as tp:
        results = list(tqdm(
            tp.map(compute_clim_im, this_source_names),
            total=len(this_source_names),
            desc=f"Compute contrast limits for {source}"
        ))

    cmin = np.median([res[0] for res in results])
    cmax = np.median([res[1] for res in results])
    clim = [float(cmin), float(cmax)]
    with open(tmp_path, 'w') as f:
        json.dump(clim, f)
    return clim


def format_tab(data):
    formatted = [[rr.decode('utf-8') for rr in row] for row in data]
    return np.array(formatted)


def create_well_table(ds_folder, table_file):
    table_path = os.path.join(ds_folder, 'tables', 'wells')
    os.makedirs(table_path, exist_ok=True)
    table_path = os.path.join(table_path, 'default.tsv')

    rel_table_path = 'tables/wells'
    if os.path.exists(table_path):
        return rel_table_path

    col_names_in = ['score', 'IgG_is_outlier', 'n_outlier_cells',
                    'n_infected', 'n_control',
                    'background_IgG_median', 'background_IgG_mad']
    col_names_out = ['score', 'is_outlier', 'n_outlier_cells',
                     'n_infected_cells', 'n_non_infected_cells',
                     'serum_background_intensity', 'serum_background_variance']
    assert len(col_names_in) == len(col_names_out)
    with open_file(table_file, 'r') as f:
        g = f['tables/images/default']
        cols = [col.decode('utf-8') for col in g['columns'][:]]
        data = g['cells'][:]

    # get the site names from the table and order by site names (order is not preserved in the table!)
    site_names = [name.decode('utf-8') for name in data[:, 1]]
    data = data[np.argsort(site_names)]
    site_names = [name.decode('utf-8') for name in data[:, 1]]

    col_ids = np.array([cols.index(col_name) for col_name in col_names_in])
    tab = data[:, col_ids]
    assert tab.shape == (len(data), len(col_ids))

    wells = np.array([site_name.split('-')[0] for site_name in site_names])
    col_names_out = ['annotation_id', 'well'] + col_names_out
    tab = np.concatenate([
        np.arange(tab.shape[0])[:, None],
        wells[:, None],
        format_tab(tab)
    ], axis=1)
    tab = pd.DataFrame(tab, columns=col_names_out)
    tab.to_csv(table_path, sep='\t', index=False, na_rep="nan")

    return rel_table_path


def create_plate_table(ds_folder, table_file, well_names):
    table_path = os.path.join(ds_folder, 'tables', 'plate')
    os.makedirs(table_path, exist_ok=True)
    table_path = os.path.join(table_path, 'default.tsv')

    rel_table_path = 'tables/plate'
    if os.path.exists(table_path):
        return rel_table_path

    col_names_in = ['score', 'IgG_is_outlier', 'n_outlier_cells',
                    'n_infected', 'n_control',
                    'background_IgG_median', 'background_IgG_mad']
    col_names_out = ['score', 'is_outlier', 'n_outlier_cells',
                     'n_infected_cells', 'n_non_infected_cells',
                     'serum_background_intensity', 'serum_background_variance']
    assert len(col_names_in) == len(col_names_out)
    with open_file(table_file, 'r') as f:
        g = f['tables/wells/default']
        cols = [col.decode('utf-8') for col in g['columns'][:]]
        data = g['cells'][:]
    assert well_names == [name.decode('utf-8') for name in data[:, 0]]

    col_ids = np.array([cols.index(col_name) for col_name in col_names_in])
    tab = data[:, col_ids]
    assert tab.shape == (len(well_names), len(col_ids))

    col_names_out = ['annotation_id', 'well'] + col_names_out
    tab = np.concatenate([
        np.arange(tab.shape[0])[:, None],
        np.array(well_names)[:, None],
        format_tab(tab)
    ], axis=1)
    tab = pd.DataFrame(tab, columns=col_names_out)
    tab.to_csv(table_path, sep='\t', index=False, na_rep="nan")

    return rel_table_path


def filter_wells(sources, well_names):
    filtered_sources = {}
    for name, this_sources in sources.items():
        filtered_sources[name] = [source for source in this_sources if any(wname in source for wname in well_names)]
    return filtered_sources


def get_all_wells(ds_meta):
    source_names = list(ds_meta["sources"].keys())
    wells = [name[name.find("Well"):].split("_")[0] for name in source_names]
    all_wells = list(set(wells))
    all_wells.sort()
    all_wells = [name[4:] for name in all_wells]
    return all_wells


def add_plate_view(ds_meta, well_names, plate_table, well_table,
                   image_names, image_types, image_colors,
                   menu_name, view_name, exclusive,
                   ds_folder, tmp_folder):
    assert len(image_names) == len(image_types) == len(image_colors)

    all_sources = ds_meta['sources']
    this_sources = {
        im_name: [name for name in all_sources if name.startswith(im_name)]
        for im_name in image_names
    }
    this_sources = filter_wells(this_sources, well_names)

    n_positions = len(this_sources[image_names[0]])
    assert all(len(sources) == n_positions for sources in this_sources.values()), f"{this_sources} -:- {n_positions}"

    all_source_positions = {
        im_name: [source[len(f'{im_name}_Well'):] for source in this_sources[im_name]]
        for im_name in image_names
    }
    source_positions = all_source_positions[image_names[0]]
    assert all(source_pos == source_positions for source_pos in all_source_positions.values())

    wells = [pos.split('_')[0] for pos in source_positions]
    source_displays = []
    for im_name, im_type, color in zip(image_names, image_types, image_colors):

        if im_type == 'image':
            im_sources = this_sources[im_name]
            clims = compute_clims(im_name, im_sources, all_sources, tmp_folder, ds_folder)
            source_display = {"imageDisplay": {
                "color": color,
                "contrastLimits": clims,
                "opacity": 1.0,
                "name": im_name,
                "sources": im_sources}
            }
        else:
            assert im_type == 'segmentation'
            seg_sources = this_sources[im_name]
            source_display = {"segmentationDisplay": {
                "lut": "glasbey",
                "opacity": 0.5,
                "name": im_name,
                "sources": seg_sources}
            }

        source_displays.append(source_display)

    # well grid transforms
    source_transforms = []
    sources_per_well = {}
    wells = np.array(wells)
    for well in well_names:
        source_ids = np.where(wells == well)[0]
        well_sources = {
            ii: [sources[sid] for sources in this_sources.values()] for ii, sid in enumerate(source_ids)
        }
        well_trafo = {
            "grid": {
                "name": well,
                "sources": well_sources
            }
        }
        well_display = {
            "sourceAnnotationDisplay": {
                "lut": "glasbey",
                "name": well,
                "sources": well_sources,
                "opacity": 0.5,
                "tableData": {"tsv": {"relativePath": well_table}},
                "tables": ["default.tsv"]
            }
        }
        source_transforms.append(well_trafo)
        source_displays.append(well_display)
        sources_per_well[well] = [source for sources in well_sources.values() for source in sources]

    def _to_pos(well_name):
        r, c = well_name[0], well_name[1:]
        r = string.ascii_uppercase.index(r)
        c = int(c) - 1
        return [r, c]

    # plate: grid transform for wells
    plate_sources = {ii: sources_per_well[well] for ii, well in enumerate(well_names)}
    plate_trafo = {
        "grid": {
            "name": f"plate-{view_name}",
            "sources": plate_sources,
            "positions": {
                ii: _to_pos(well) for ii, well in enumerate(well_names)
            }
        }
    }
    plate_display = {
        "sourceAnnotationDisplay": {
            "name": f"plate-{view_name}",
            "sources": plate_sources,
            "opacity": 0.5,
            "lut": "glasbey",
            "tableData": {"tsv": {"relativePath": plate_table}},
            "tables": ["default.tsv"]
        }
    }
    source_transforms.append(plate_trafo)
    source_displays.append(plate_display)

    grid_view = {
        "isExclusive": exclusive,
        "sourceDisplays": source_displays,
        "sourceTransforms": source_transforms,
        "uiSelectionGroup": menu_name
    }

    views = ds_meta['views']
    views[view_name] = grid_view
    ds_meta['views'] = views
    return ds_meta


def create_raw_views(plate_name, plate_table, well_table):
    ds_folder = os.path.join('./data', plate_name)
    ds_meta = mobie.metadata.read_dataset_metadata(ds_folder)

    # remove single source views
    ds_meta = remove_single_source_views(ds_meta,
                                         remove_prefixes=['marker', 'serumIgG', 'nuclei'])

    all_wells = get_all_wells(ds_meta)
    # add the individual grid views
    tmp_folder = f'./tmp_{plate_name}'
    add_plate_view(ds_meta, all_wells, plate_table, well_table, ['nuclei'], ['image'], ['blue'],
                   menu_name="images", view_name="nuclei", exclusive=False,
                   ds_folder=ds_folder, tmp_folder=tmp_folder)
    add_plate_view(ds_meta, all_wells, plate_table, well_table, ['serumIgG'], ['image'], ['green'],
                   menu_name="images", view_name="serumIgG", exclusive=False,
                   ds_folder=ds_folder, tmp_folder=tmp_folder)
    add_plate_view(ds_meta, all_wells, plate_table, well_table, ['marker_tophat'], ['image'], ['red'],
                   menu_name="images", view_name="marker", exclusive=False,
                   ds_folder=ds_folder, tmp_folder=tmp_folder)

    # TODO replace the default view with a composite view of marker + serumIgG + nuclei

    mobie.metadata.write_dataset_metadata(ds_folder, ds_meta)


def create_test_views(plate_name, plate_table, well_table):
    ds_folder = os.path.join('./data', plate_name)
    ds_meta = mobie.metadata.read_dataset_metadata(ds_folder)

    well_names = ["E06", "E07"]
    tmp_folder = f'./tmp_{plate_name}'

    # add nucleus view
    add_plate_view(ds_meta, well_names, plate_table, well_table, ['nuclei'], ['image'], ['blue'],
                   menu_name="images", view_name="test-nuclei", exclusive=True,
                   ds_folder=ds_folder, tmp_folder=tmp_folder)

    # add full view
    add_plate_view(ds_meta, well_names, plate_table, well_table,
                   ['nuclei', 'serumIgG', 'marker_tophat'],
                   ['image', 'image', 'image'],
                   ['blue', 'green', 'red'],
                   menu_name="images", view_name="test-full", exclusive=True,
                   ds_folder=ds_folder, tmp_folder=tmp_folder)

    mobie.metadata.write_dataset_metadata(ds_folder, ds_meta)


def add_plate(plate_folder):
    plate_name = parse_plate_name(plate_folder)
    ds_folder = f"./data/{plate_name}"
    print("Adding plate", plate_name)

    input_files = glob(os.path.join(plate_folder, "*.h5"))
    input_files.sort()
    print("with", len(input_files), "images")
    # add_image_data(input_files, plate_name)
    # make_2d(plate_name)

    # create the tables
    table_file = os.path.join(plate_folder, f'{plate_name}_table.hdf5')
    assert os.path.exists(table_file)
    all_wells = get_all_wells(mobie.metadata.read_dataset_metadata(ds_folder))
    plate_table = create_plate_table(ds_folder, table_file, all_wells)
    well_table = create_well_table(ds_folder, table_file)

    # create_raw_views(plate_name, plate_table, well_table)
    create_test_views(plate_name, plate_table, well_table)

    # TODO
    # add segmentations
    # add image / well tables
    # add the full grid view

    # validate the project
    mobie.validation.validate_project('./data')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # add default argument for debugging
    parser.add_argument('--input', '-i', default=DEFAULT_PLATE)
    args = parser.parse_args()
    add_plate(args.input)
