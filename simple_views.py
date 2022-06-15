import mobie


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
        if ref_settings["name"] in source_prefixes:
            name = ref_settings.pop("name")
            ref_settings.pop("sources")
            display_group_settings[name] = ref_settings

    mobie.create_grid_view(
        ds_folder, "single_well", sources,
        menu_name="simple_views", display_groups=display_groups,
        display_group_settings=display_group_settings,
        use_transformed_grid=False
    )


if __name__ == "__main__":
    single_well_view(well="F09", source_prefixes=["nuclei", "serumIgG", "marker_tophat"])
