###############################################################################
#
# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.
#
###############################################################################
__version__ = "0.0.1"
from copy import deepcopy
from datetime import datetime
import json
import logging
from pathlib import Path
import numpy as np
from pyproj import Geod
from typing import Union

from bufr2geojson import transform as as_geojson
from wis2box.data.base import BaseAbstractData
from wis2box.api import upsert_collection_item

LOGGER = logging.getLogger(__name__)

g = Geod(ellps="WGS84")

bearingToName = {
    "0.0-90.0": "RAD1",
    "90.0-180.0": "RAD2",
    "180.0-270.0": "RAD3",
    "270.0-0.0": "RAD4"
}

class BUFR2JTWC(BaseAbstractData):
    """Model data"""
    def transform(self, input_data: Union[Path, bytes], filename: str = '') -> bool:

        LOGGER.debug('Procesing BUFR data')
        input_bytes = self.as_bytes(input_data)

        LOGGER.debug('Generating GeoJSON features')
        results = as_geojson(input_bytes, serialize=False)

        LOGGER.debug('Processing GeoJSON features')
        count = 0
        for collection in results:
            # results is an iterator, for each iteration we have:
            # - dict['id']
            # - dict['id']['_meta']
            # - dict['id']
            for id, item in collection.items():
                count += 1
                # extract identification
                # check if we have ensemble member number, if not use subset
                subset = item['geojson']['properties']['subset']
                for metadata in item['geojson']['properties']['metadata']:
                    if metadata['name'] == "ensemble_member_number":
                        subset = metadata['value']
                stormIdentifier = item['geojson']['properties']['wigos_station_identifier']
                stormName, stormNumber = stormIdentifier.split("-")
                # now warning time and forecast (tau)
                warningTime = item['geojson']['properties']['phenomenonTime']
                forecastTime = None
                tau = None
                if "/" in warningTime:
                    forecastTime, tau = warningTime.split("/")
                else:
                    forecastTime = warningTime
                    tau = item['geojson']['properties']['resultTime']
                try:
                    forecastTime = datetime.strptime(forecastTime,"%Y-%m-%dT%H:%M:%SZ")
                    tau = datetime.strptime(tau,"%Y-%m-%dT%H:%M:%SZ")
                except Exception as e:
                    LOGGER.debug(f"Error setting time in geojson processing, error: {e}")
                tau = int((tau - forecastTime).total_seconds() / 3600)
                # now construct key for completing data.frame
                #LOGGER.debug("extracting data")
                key1 = f"{stormIdentifier}-{subset}-{forecastTime}-{tau}"
                geojson_out = deepcopy(item['geojson'])
                if geojson_out['properties']['name'] ==  "pressure_reduced_to_mean_sea_level":
                    geojson_out = self.extract_MSLP(geojson_out)
                    key2 = "MSLP"
                elif geojson_out['properties']['name'] ==  "wind_speed_at10m":
                    geojson_out = self.extract_vmax(geojson_out)
                    key2 = "Vmax"
                elif geojson_out['properties']['name'] ==  "effective_radius_with_respect_to_wind_speeds_above_threshold":
                    bearing = None
                    for metadata in geojson_out['properties']['metadata']:
                        if metadata['name'] == "bearing_or_azimuth":
                            bearing = metadata['value']
                            bearing = f"{bearing[0]}-{bearing[1]}"
                    assert bearing is not None
                    geojson_out = self.extract_wind_polygon(geojson_out)
                    key2 = bearingToName[bearing]
                else:
                    assert False

                #LOGGER.debug("publishing")
                data_date = forecastTime

                self.output_data[f"{key1}-{key2}"] = {
                    'geojson': geojson_out, # json.dumps(geojson_out),
                    '_meta': {
                        'data_date': data_date.strftime('%Y-%m-%d %H:%M'),
                        'relative_filepath': self.get_local_filepath(data_date)
                    }
                }

        LOGGER.debug('Successfully finished transforming BUFR data')
        LOGGER.debug(f"{count} features processed")
        return True

    def get_local_filepath(self, date_):
        yyyymmdd = date_.strftime('%Y-%m-%d')
        return Path(yyyymmdd) / 'wis' / self.topic_hierarchy.dirpath

    def extract_vmax(self, feature):
        forecastTime = feature['properties']['phenomenonTime']
        if "/" in forecastTime:
            t1,t2 = forecastTime.split("/")
        else:
            t1 = forecastTime
            t2 = feature['properties']['resultTime']
        feature['properties']['resultTime'] = t1
        feature['properties']['phenomenonTime'] = t2
        return feature

    def extract_MSLP(self, feature):
        #LOGGER.debug("Extracting MSLP as GeoJSON")
        forecastTime = feature['properties']['phenomenonTime']
        if "/" in forecastTime:
            t1,t2 = forecastTime.split("/")
        else:
            t1 = forecastTime
            t2 = feature['properties']['resultTime']
        feature['properties']['resultTime'] = t1
        feature['properties']['phenomenonTime'] = t2
        return feature

    def extract_wind_polygon(self, feature):
        keep = ("centre", "generating_application", "storm_identifier", "long_storm_name",
                "technique_for_making_up_initial_perturbations", "ensemble_member_number", "ensemble_forecast_type",
                "meteorological_attribute_significance")

        radius = feature['properties']['value']
        parameters = feature['properties']['metadata']
        lon = feature['geometry']['coordinates'][0]
        lat = feature['geometry']['coordinates'][1]

        # get bearing
        bearing = None
        for parameter in parameters:
            if parameter["name"] == "bearing_or_azimuth":
                bearing = parameter["value"]
        assert (bearing is not None)
        if bearing[1] == 0:
            bearing[1] = 360
        # get wind speed
        wind_speed = None
        units = None
        for parameter in parameters:
            if parameter["name"] == "wind_speed_threshold":
                wind_speed = parameter["value"]
                units = parameter["units"]
        assert (wind_speed is not None)

        # drop unwanted / used parameters
        feature['properties']['parameters'] = list()
        for parameter in parameters:
            if parameter['name'] in keep:
                feature['properties']['parameters'].append(parameter)
        del feature['properties']['metadata']

        x = list(map(lambda b: g.fwd(lon, lat, b, radius)[0:2], np.arange(bearing[0], bearing[1] + 2.5, 2.5)))
        x.insert(0, (lon, lat))
        x.append((lon, lat))
        if radius > 0:
            feature['geometry']['type'] = "Polygon"
            feature['geometry']['coordinates'] = list()
            feature['geometry']['coordinates'].insert(0, x)
        feature['properties']['name'] = "wind_speed_threshold"
        feature['properties']['value'] = wind_speed
        feature['properties']['units'] = units
        forecastTime = feature['properties']['phenomenonTime']
        if "/" in forecastTime:
            t1,t2 = forecastTime.split("/")
        else:
            t1 = forecastTime
            t2 = feature['properties']['resultTime']
        feature['properties']['resultTime'] = t1
        feature['properties']['phenomenonTime'] = t2
        return feature

    def list_test(self, items):
        LOGGER.debug(len(items))
        assert len(items) > 1

    def publish(self) -> bool:
        LOGGER.info('Publishing output data')
        upsert_list = []
        for identifier, item in self.output_data.items():
            # now iterate over formats
            for format_, the_data in item.items():
                if format_ == '_meta':  # not data, skip
                    continue

                #LOGGER.debug(f'Processing format: {format_}')
                # check that we actually have data
                if the_data is None:
                    msg = f'Empty data for {identifier}-{format_}; not publishing'  # noqa
                    LOGGER.warning(msg)
                    continue
                upsert_list.append(deepcopy(the_data))
        LOGGER.debug('Publishing data to API')
        LOGGER.debug(f"{len(upsert_list)} items to publish")

        upsert_collection_item(self.topic_hierarchy.dotpath, upsert_list)

        return True