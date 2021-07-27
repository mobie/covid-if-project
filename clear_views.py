import mobie


def clear_views():
    ds_folder = "./data/20200406_164555_328"
    metadata = mobie.metadata.read_dataset_metadata(ds_folder)

    views = metadata["views"]
    names = list(views.keys())
    for name in names:
        if name != "default":
            views.pop(name)

    metadata["views"] = views
    mobie.metadata.write_dataset_metadata(ds_folder, metadata)


clear_views()
