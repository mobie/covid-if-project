

def add_metadata():
    from mobie.metadata import add_remote_project_metadata

    bucket_name = "covid-if-project"
    endpoint = "https://s3.embl.de"
    add_remote_project_metadata("./data", bucket_name, endpoint)


def upload():
    from subprocess import run

    cmd = ["mc", "cp", "-r", "data/", "embl/covid-if-project/"]
    # print(cmd)
    run(cmd)


def update():
    cmd = ["mc", "cp", "data/20200406_164555_328/dataset.json", "embl/covid-if-project/20200406_164555_328/dataset.json"]


if __name__ == '__main__':
    # add_metadata()
    upload()
    # update()
