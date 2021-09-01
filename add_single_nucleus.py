import mobie


in_path = "./data/20200406_164555_328/images/ome-zarr/nuclei_WellE05_0001.ome.zarr"
in_key = "s0"

scale_factors = [[2, 2], [2, 2], [2, 2], [2, 2]]
chunks = (256, 256)

mobie.add_image(in_path, in_key,
                root="./data", dataset_name="20200406_164555_328",
                image_name="single-nucleus", resolution=[1, 1], scale_factors=scale_factors,
                chunks=chunks, file_format="ome.zarr", menu_name="test")
