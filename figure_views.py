import os
import json
import mobie
import pandas as pd


def save_view(name, view, save_to_project):
    view["uiSelectionGroup"] = "paper"
    view["isExclusive"] = True

    if save_to_project:
        metadata = mobie.metadata.read_dataset_metadata("./data/20200406_164555_328")
        metadata["views"][name] = view
        mobie.metadata.write_dataset_metadata("./data/20200406_164555_328", metadata)

    else:
        view_root = "./data/20200406_164555_328/misc/bookmarks"
        os.makedirs(view_root, exist_ok=True)
        view_file = os.path.join(view_root, "mobie-paper.json")
        if os.path.exists(view_file):
            with open(view_file) as f:
                views = json.load(f)
        else:
            views = {"views": {}}
        views["views"][name] = view

        with open(view_file, "w") as f:
            json.dump(views, f, indent=2, sort_keys=True)


def panel_a():
    metadata = mobie.metadata.read_dataset_metadata("./data/20200406_164555_328")
    views = metadata["views"]
    grid = views["full_grid"]

    # make green channel visible
    serum = grid["sourceDisplays"][1]["imageDisplay"]
    serum["visible"] = True
    grid["sourceDisplays"][1]["imageDisplay"] = serum

    # select well F9
    wells = grid["sourceDisplays"][-1]["sourceAnnotationDisplay"]
    table = pd.read_csv("./data/20200406_164555_328/tables/well/default.tsv", sep="\t")
    well_names = table.annotation_id.values.tolist()
    well_id = well_names.index("F09")
    wells["selectedAnnotationIds"] = [f"0;{well_id}"]
    wells["showAsBoundaries"] = True
    wells["boundaryThickness"] = 100
    grid["sourceDisplays"][-1]["sourceAnnotationDisplay"] = wells

    # save the new view
    save_view("Figure 3 a)", grid, save_to_project=True)


def panel_b():
    pass


def panel_c():
    pass


def panel_d():
    pass


if __name__ == "__main__":
    panel_a()
    panel_b()
    panel_c()
    panel_d()
