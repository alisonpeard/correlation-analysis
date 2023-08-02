#!/Users/alison/mambaforge/envs/snakemake/bin/python
"""
Requires:
---------
./config.json
datadir/.../wrz.shp
datadir/.../wrz_code.xlsx

Output:
-------
tempdir/wrz.gpkg
tempdir/wrz_buffer.gpkg
"""
import os
os.environ['USE_PYGEOS'] = '0'
import utils
import pandas as pd
import geopandas as gpd

print("\nProcessing WRZ files.")

def main(config):
    datadir = config['paths']['datadir']
    outdir = config['paths']["tempdir"]
    uk_crs = config['config']['uk_crs']
    buffers = config['config']['buffers']

    # load and merge WRZ dataframes
    wrz = gpd.read_file(os.path.join(datadir, 'WRZ', 'WRZ.shp'), crs=uk_crs).to_crs('EPSG:4326')
    wrz_code = pd.read_excel(os.path.join(datadir, 'WRZ', 'wrz_code.xlsx'))
    wrz_code['RZ ID'] = pd.to_numeric(wrz_code['RZ ID'])
    wrz = wrz.merge(wrz_code, left_on='RZ_ID', right_on='RZ ID')

    # extract centroids of wrzs
    wrz_centroid = wrz.copy().to_crs(uk_crs)
    wrz_centroid['centroid_x'] = wrz_centroid.geometry.centroid.x
    wrz_centroid['centroid_y'] = wrz_centroid.geometry.centroid.y
    wrz_centroid = gpd.GeoDataFrame(wrz_centroid, geometry=gpd.points_from_xy(wrz_centroid.centroid_x, wrz_centroid.centroid_y), crs=uk_crs)

    # create buffer dataframes
    wrz_buffers = []
    for buffer in buffers:
        wrz_buffer = wrz_centroid.copy()
        wrz_buffer['geometry'] = wrz_buffer['geometry'].buffer(buffer).to_crs(4326)
        wrz_buffer['buffer'] = [buffer / 1000] * len(wrz_buffer)
        wrz_buffer = wrz_buffer[['RZ_ID', 'geometry', 'buffer']].drop_duplicates()
        wrz_buffers.append(wrz_buffer)
    wrz_buffer = pd.concat(wrz_buffers).sort_values(by=['RZ_ID', 'buffer'], ascending=True)[['RZ_ID', 'buffer', 'geometry']]

    wrz.to_file(os.path.join(outdir, 'wrz.gpkg'), driver="GPKG")
    wrz_buffer.to_file(os.path.join(outdir, 'wrz_buffer.gpkg'), driver="GPKG")


if __name__ == "__main__":
    config = utils.load_config()
    main(config)




