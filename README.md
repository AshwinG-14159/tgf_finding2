First ideas:

python3 get orbits.py --only_mkf --work_dir data_download 2023-06-16T07:13:37 name
This command will open the orbitinfo.csv and inform which mkf file contains data on the timestamp given
data_download should be the directory containing orbitinfo.csv for this guy to read.


Feed the mkf file to give_ra_dec.py and it will tell the pointing angle of czti at that point of time.


python3 angles.py /mnt/nas2_czti/czti/level2/20230615_A12_077T12_9000005692_level2/czti/orbit/41718_V1.0 123d 10d  023-06-16T07:13:37

This command gives the separation angle according to the provided mkf file between the location given and the earth.

Should combine these three so that given timestamp, one can get angle between czti pointing and earth.



Making this combined file:
1)

python3 time_to_sep.py  --only_mkf --work_dir ../data 2023-06-16T07:13:37 name

output:
Orbit info file read
Orbit found
event obs:: ['20230615_A12_077T12_9000005692_level2_41718']
/mnt/nas2_czti/czti/level2/20230615_A12_077T12_9000005692_level2/czti/orbit/41718_V1.0/AS1A12_077T12_9000005692_41718czt_level2.mkf
211.3632d
25.91804d

Now to add the angles functionality:

python3 time_to_sep.py --work_dir ../data 2023-06-16T07:13:37

output:
Orbit info file read
Orbit found
/mnt/nas2_czti/czti/level2/20230615_A12_077T12_9000005692_level2/czti/orbit/41718_V1.0/AS1A12_077T12_9000005692_41718czt_level2.mkf
211.3632d
25.91804d
Theta of earth is: 107.26292215755154
Phi of Earth is: 74.09527422250595
Earth Transient Angle:  107.26622613827726



















