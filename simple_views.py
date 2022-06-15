import os
import mobie
import numpy as np
import pandas as pd


def single_well_view(well, source_prefixes):
    ds_folder = "./data/20200406_164555_328"
    ds_meta = mobie.metadata.read_dataset_metadata(ds_folder)
    all_sources = ds_meta["sources"]
    sources = [
        [source for source in all_sources
         if (f"{well}_000{ii}" in source and any(pref in source for pref in source_prefixes))]
        for ii in range(9)
    ]

    display_groups = {}
    for source_list in sources:
        display_groups.update(**{source: [pref for pref in source_prefixes if pref in source][0]
                                 for source in source_list})

    reference_settings = ds_meta["views"]["Figure3a"]["sourceDisplays"]
    display_group_settings = {}
    for ref in reference_settings:
        ref_settings = next(iter(ref.values()))
        ref_name = ref_settings.pop("name")
        if ref_name in source_prefixes:
            source_index = [i for i, pref in enumerate(source_prefixes) if ref_name == pref][0]
            ref_settings.pop("sources")
            ref_settings["sources"] = [f"single_well-{source_index}"]
            ref_settings["visible"] = source_index == 1
            display_group_settings[ref_name] = ref_settings

    mobie.create_grid_view(
        ds_folder, "single_well", sources,
        menu_name="bookmarks", display_groups=display_groups,
        display_group_settings=display_group_settings,
        use_transformed_grid=False, overwrite=True
    )


def single_well_view_with_extra_table(well, source_prefixes):
    ds_folder = "./data/20200406_164555_328"
    ds_meta = mobie.metadata.read_dataset_metadata(ds_folder)
    all_sources = ds_meta["sources"]
    sources = [
        [source for source in all_sources
         if (f"{well}_000{ii}" in source and any(pref in source for pref in source_prefixes))]
        for ii in range(9)
    ]

    display_groups = {}
    for source_list in sources:
        display_groups.update(**{source: [pref for pref in source_prefixes if pref in source][0]
                                 for source in source_list})

    reference_settings = ds_meta["views"]["Figure3a"]["sourceDisplays"]
    display_group_settings = {}
    for ref in reference_settings:
        ref_settings = next(iter(ref.values()))
        ref_name = ref_settings.pop("name")
        if ref_name in source_prefixes:
            source_index = [i for i, pref in enumerate(source_prefixes) if ref_name == pref][0]
            ref_settings.pop("sources")
            ref_settings["sources"] = [f"single_well-{source_index}"]
            ref_settings["visible"] = source_index == 1
            display_group_settings[ref_name] = ref_settings

    view = mobie.create_grid_view(
        ds_folder, "single_well", sources,
        menu_name="bookmarks", display_groups=display_groups,
        display_group_settings=display_group_settings,
        use_transformed_grid=False, return_view=True,
    )

    table_folder = "./data/20200406_164555_328/tables/single_well"
    tab = pd.read_csv(os.path.join(table_folder, "default.tsv"), sep="\t")
    tab = tab.drop(columns="source")
    tab["extra-col1"] = np.random.rand(len(tab))
    tab["extra-col2"] = np.arange(3, 3 + len(tab))
    tab["extra-col3"] = np.random.rand(len(tab))
    tab["lorem"] = "ipsum dolor sit amet, consetetur sadipscing elitier, sed diam nonumy eirdom temporume indiviunt ut labore et dolore magna".split()[:len(tab)]
    tab.to_csv(os.path.join(table_folder, "extra.tsv"), sep="\t")

    # create the extra table
    view["sourceDisplays"][-1]["regionDisplay"]["tables"].append("extra.tsv")

    ds_meta = mobie.metadata.read_dataset_metadata("./data/20200406_164555_328")
    ds_meta["views"]["single-well-extra-table"] = view
    mobie.metadata.write_dataset_metadata("./data/20200406_164555_328", ds_meta)


if __name__ == "__main__":
    # single_well_view(well="F09", source_prefixes=["marker_tophat", "nuclei", "serumIgG"])
    single_well_view_with_extra_table(well="F09", source_prefixes=["marker_tophat", "nuclei", "serumIgG"])
