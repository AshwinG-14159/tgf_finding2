import subprocess as subp
import argparse
import glob
import numpy as np
import pandas as pd
from astropy.time import Time
import pexpect
import os
import requests


def getOrbitinfo(work_dir, download=True):
    """
    Download and get information for the orbits
    The function returns obsid_orbit and the start and end times of all the orbits
    """
    subp.call(f"cd {work_dir}", shell=True)
    if download:
        print("Downloading the orbitinfo.csv file")
        orbitinfo_url = "https://web.iucaa.in/~astrosat/czti_dqr/orbitinfo.csv"
        r = requests.get(orbitinfo_url)
        with open(f"{work_dir}/orbitinfo.csv".format(work_dir=work_dir), "wb") as f:
            f.write(r.content)
    orbitinfo = pd.read_csv(f"{work_dir}/orbitinfo.csv")
    print("Orbit info file read")

    def fix_times(time_array):
        """
        Some of the earlier orbits have times that look like this 2015-12-03T12-03-10 instead of 2015-12-03T12:03:10
        This function returns the correct format for isot (astropy.time)
        """
        time_array_colons = np.char.replace(time_array, "-", ":")
        time_array_isot = np.char.replace(
            time_array_colons, ":", "-", 2
        )  # only replace the first 2 ':' with '-' for the date part
        return time_array_isot

    tstart = fix_times(np.array(orbitinfo["tstart"], "str"))
    tend = fix_times(np.array(orbitinfo["tend"], "str"))

    obsnames = np.array(orbitinfo["obsname"], "str")
    ind_not_combined = np.logical_not(np.char.endswith(obsnames, "level2"))

    return obsnames[ind_not_combined], tstart[ind_not_combined], tend[ind_not_combined]


def find_orbits(triggerTime_UTC, work_dir, download=True, num_neighbours=5):
    """
    Find the event and neighbouring orbit indices of a given transient
    """
    triggertime_UTC = Time(triggerTime_UTC).isot
    obsnames, tstart, tend = getOrbitinfo(work_dir, download)

    ind_tend_greater = np.where(np.char.greater(tend, triggertime_UTC))[0]
    ind_tstart_lesser = np.where(np.char.less(tstart, triggertime_UTC))[0]

    matching_orbit_ind = np.intersect1d(ind_tend_greater, ind_tstart_lesser)

    if matching_orbit_ind.size == 0:
        print("Orbit not found! Please check http://www.iucaa.in/~astrosat/czti_dqr")
        raise Exception("Orbit Not Found")
    elif matching_orbit_ind.size > 1:
        print("Multiple Matching Orbits!")
    else:
        print("Orbit found")

    neighbouring_orbit_ind = np.append(
        np.arange(matching_orbit_ind.min(), matching_orbit_ind.min(), 1),
        np.arange(
            matching_orbit_ind.max() + 1,
            matching_orbit_ind.max() + num_neighbours + 1,
            1,
        ),
    )

    return obsnames[matching_orbit_ind], obsnames[neighbouring_orbit_ind]


def get_path_remote(root, obsname):
    """
    Get the path on nas given the obsname
    """
    date_obsid, orbit = obsname.rsplit("_", maxsplit=1)
    base_path = f"{root}/{date_obsid}/czti/orbit/{orbit}"
    obsid = date_obsid.split("_", maxsplit=1)[1].rsplit("_", maxsplit=1)[0]

    mkf_file_path = f"{base_path}/AS1{obsid}_{orbit}czt_level2.mkf"
    bc_file_path = f"{base_path}/modeM0/AS1{obsid}_{orbit}cztM0_level2_bc.evt"
    livetime_file_path = (
        f"{base_path}/modeM0/AS1{obsid}_{orbit}cztM0_level2_bc_livetime.fits"
    )
    return mkf_file_path, bc_file_path, livetime_file_path


def get_path_local(root, obsname):
    date_obsid, orbit = obsname.rsplit("_", maxsplit=1)
    obsid = date_obsid.split("_", maxsplit=1)[1].rsplit("_", maxsplit=1)[0]
    mkf_file_path = f"{root}/AS1{obsid}_{orbit}czt_level2.mkf"
    bc_file_path = f"{root}/AS1{obsid}_{orbit}cztM0_level2_bc.evt"
    livetime_file_path = f"{root}/AS1{obsid}_{orbit}cztM0_level2_bc_livetime.fits"

    return mkf_file_path, bc_file_path, livetime_file_path


def rsync_command(remote, local):
    # proxy_command = '"ProxyCommand ssh -A czti@10.112.5.88 -W %h:%p"'
    # command = f"rsync -vzrutP -e 'ssh -o {proxy_command}' cztipoc@192.168.11.37:{remote} {local}"
    command = f"rsync -vzrutP cztipoc@192.168.11.37:{remote} {local}"
    print(command)
    return command


def download(command):
    returncode = -1
    # while returncode!=0:
    # ret = subp.run(command, shell=True)
    # returncode = ret.returncode
    # czti_token = os.environ["czti_vb_token"]
    arun_token = os.environ["arun_token"]
    child = pexpect.spawn(command)
    # child.expect("czti@10.112.5.88's password:")
    # child.sendline(czti_token)
    # child.expect("cztipoc@192.168.11.37's password:")
    # child.sendline(arun_token)
    child.interact()
    print("File Downloaded")


def download_command(obsname, trans_name, event, work_dir, checknas, only_mkf=False):
    path_poc = "/data2/czti/level2"
    if event:
        path_proc = f"{work_dir}/{trans_name}/Event_files"
    else:
        path_proc = f"{work_dir}/{trans_name}"
    # path_proc = '.'
    mkf_l, bc_l, livet_l = get_path_local(path_proc, obsname)

    # Checking if the folder exists or not
    if not os.path.exists(work_dir):
        # print('Hola - ', work_dir, trans_name)
        subp.call(f"mkdir -p {work_dir}/{trans_name}/Event_files", shell=True)

    # Checking if the data exists in NAS2
    localroot = "/home/cift/grb_search/data/level2"
    date_obsid, orbit = obsname.rsplit("_", maxsplit=1)
    base_path = f"{localroot}/{date_obsid}/czti/orbit/{orbit}_V1.0"
    filelist = glob.glob(f"{base_path}/modeM0/*_bc.evt")
    bc_file = bc_l.split("/")[-1]
    all_bc_files = [path.split("/")[-1] for path in filelist]
    if bc_file in all_bc_files:
        print("Found in grb_search folder, Rsyncing")
        print(f"rsync -vzrutP {base_path}/modeM0/*_bc.evt {bc_l}")
        subp.call(f"rsync -vzrutP {base_path}/modeM0/*_bc.evt {bc_l}", shell=True)
        subp.call(
            f"rsync -vzrutP {base_path}/modeM0/*bc_livetime.fits {livet_l}", shell=True
        )
        subp.call(f"rsync -vzrutP {base_path}/*.mkf {mkf_l}", shell=True)
    # if not then download
    else:
        if checknas == 'yes':
            print("Checking in NAS2")
            nasroot = "/mnt/nas2_czti/czti/level2"
            date_obsid, orbit = obsname.rsplit("_", maxsplit=1)
            base_path = f"{nasroot}/{date_obsid}/czti/orbit/{orbit}_V1.0"
            filelist = glob.glob(f"{base_path}/modeM0/*_bc.evt")
            livetime_list = glob.glob(f"{base_path}/modeM0/*_bc_livetime.fits")
            if len(filelist) + len(livetime_list) > 1:
                print("Downloading from NAS2")
                bc_file = bc_l.split("/")[-1]
                all_bc_files = [path.split("/")[-1] for path in filelist]
                if bc_file in all_bc_files:
                    print("Found in NAS2, Rsyncing")
                    print(f"rsync -vzrutP {base_path}/modeM0/*_bc.evt {bc_l}")
                    subp.call(f"rsync -vzrutP {base_path}/modeM0/*_bc.evt {bc_l}", shell=True)
                    subp.call(
                        f"rsync -vzrutP {base_path}/modeM0/*bc_livetime.fits* {livet_l}",
                        shell=True,
                    )
                    subp.call(f"rsync -vzrutP {base_path}/*.mkf {mkf_l}", shell=True)
            else:
                print("Not all files found in NAS2, Downloading from Arun")
                mkf_r, bc_r, livet_r = get_path_remote(path_poc, obsname)
                print(mkf_r, bc_r, livet_r)

                command = rsync_command(mkf_r, mkf_l)
                download(command)

                if only_mkf == False:
                    command = rsync_command(bc_r, bc_l)
                    download(command)
                    command = rsync_command(livet_r, livet_l)
                    download(command)
        else:
            print("Downloading from Arun")
            mkf_r, bc_r, livet_r = get_path_remote(path_poc, obsname)
            print(mkf_r, bc_r, livet_r)

            command = rsync_command(mkf_r, mkf_l)
            download(command)

            if only_mkf == False:
                command = rsync_command(bc_r, bc_l)
                download(command)
                command = rsync_command(livet_r, livet_l)
                download(command)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("triggertime", help="Trigger time in UTC", type=str)
    parser.add_argument("transname", help="Transient Name", type=str)
    parser.add_argument(
        "--only_mkf",
        dest="only_mkf",
        help="Download only mkf if given",
        action="store_true",
    )
    parser.add_argument(
        "--event", dest="event", help="Download Event Files", action="store_true"
    )
    parser.add_argument(
        "--neighbour",
        dest="neighbour",
        help="Download Neighbour Files",
        action="store_true",
    )
    parser.add_argument(
        "--download",
        dest="download",
        help="Download the orbitinfo.csv file",
        action="store_true",
    )
    parser.add_argument(
        "--work_dir",
        dest="work_dir",
        help="Full path to the working or the download destination  directory",
        type=str,
        default="/home/czti/transient_search/",
    )
    parser.add_argument(
        "--code_dir",
        dest="code_dir",
        help="Full path to the code directory of fluxlimits codes",
        type=str,
        default="/home/czti/CZTI-fast-transients/fluxlimit codes/",
    )
    parser.add_argument(
        "--checknas",
        dest="checknas",
        help="Say yes if you want to do a NAS check",
        type=str,
        default="yes",
    )
    args = parser.parse_args()

    event_obs, neighbour_obs = find_orbits(
        args.triggertime, args.work_dir, args.download
    )
    print("event obs::", event_obs)

    for ev in event_obs:
        if args.event:
            # download_command(ev, args.transname, event=True, only_mkf=args.only_mkf)
            download_command(
                ev, args.transname, True, args.work_dir, args.checknas, only_mkf=args.only_mkf
            )
            # print(ev)
        else:
            mkf, bc, lv = get_path_remote("/data2/czti/level2/", ev)
            print(mkf)
            # print(bc)
            # print(lv)
    if args.neighbour:
        print("neighbour obs::", neighbour_obs)
        for ne in neighbour_obs:
            print("ne", ne)
            download_command(
                ne, args.transname, False, args.work_dir, args.checknas, only_mkf=args.only_mkf
            )