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


def main():
    ds_folder = "./data/20200406_164555_328"
    metadata = mobie.metadata.read_dataset_metadata(ds_folder)
    views = metadata["views"]
    for name, view in views.items():
        updated_view = add_reference_to_grids_in_view(view)
        views[name] = updated_view
    mobie.metadata.write_dataset_metadata(ds_folder, metadata)


if __name__ == "__main__":
    main()
