import mobie

ds_folder = "./data/20200406_164555_328"
view = mobie.metadata.get_view(["nuclei", "cell-seg"], ["image", "segmentation"],
                               [["nuclei_WellE06_0000"], ["cell_segmentation_WellE06_0000"]],
                               [{"color": "blue"}, {"tables": ["default.tsv"], "lut": "glasbey"}],
                               is_exclusive=True, menu_name="bookmarks")
mobie.metadata.add_view_to_dataset(ds_folder, "single_image", view)
