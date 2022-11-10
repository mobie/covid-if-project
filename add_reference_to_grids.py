import mobie


def add_reference_to_grids_in_view(view):
    if "sourceTransforms" not in view:
        return view
    source_transforms = view["sourceTransforms"]
    new_source_transforms = []
    references = {}
    for trafo in source_transforms:
        trafo_type, trafo = next(iter(trafo.items()))
        if trafo_type == "mergedGrid":
            sources = trafo["sources"]
            source_prefix = sources[0].split("_")[0]
            if source_prefix.startswith("A01"):
                trafo["metadataSource"] = sources[0]
            else:
                if source_prefix not in references:
                    references[source_prefix] = sources[0]
                trafo["metadataSource"] = references[source_prefix]
        trafo = {trafo_type: trafo}
        new_source_transforms.append(trafo)
    view["sourceTransforms"] = new_source_transforms
    return view


def add_reference_all_views(ds_folder):
    metadata = mobie.metadata.read_dataset_metadata(ds_folder)
    views = metadata["views"]
    for name, view in views.items():
        updated_view = add_reference_to_grids_in_view(view)
        views[name] = updated_view
    mobie.metadata.write_dataset_metadata(ds_folder, metadata)


def update_source_names_in_view(view):
    if "sourceTransforms" not in view:
        return view

    source_transforms = view["sourceTransforms"]
    source_displays = view["sourceDisplays"]

    updated_source_names = {}
    for trafo in source_transforms:
        trafo_type, trafo = next(iter(trafo.items()))
        if trafo_type != "mergedGrid":
            continue
        sources = trafo["sources"]
        grid_name = trafo["mergedGridSourceName"]
        updated_source_names.update({name: f"{name}_{grid_name}" for name in sources})

    # this happens for the view with a transformed grid
    if not updated_source_names:
        return view

    new_source_displays = []
    for display in source_displays:
        display_type, display = next(iter(display.items()))
        if display_type == "regionDisplay":
            sources = display["sources"]
            try:
                sources = {pos: [updated_source_names[source] for source in pos_sources]
                           for pos, pos_sources in sources.items()}
            except KeyError as e:
                breakpoint()
                raise e
            display["sources"] = sources
        display = {display_type: display}
        new_source_displays.append(display)
    view["sourceDisplays"] = new_source_displays

    return view


def update_source_names_all_views(ds_folder):
    metadata = mobie.metadata.read_dataset_metadata(ds_folder)
    views = metadata["views"]
    for name, view in views.items():
        updated_view = update_source_names_in_view(view)
        views[name] = updated_view
    mobie.metadata.write_dataset_metadata(ds_folder, metadata)


def main():
    ds_folder = "./data/20200406_164555_328"
    # add_reference_all_views(ds_folder)
    update_source_names_all_views(ds_folder)


if __name__ == "__main__":
    main()
