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
import logging
from pathlib import Path
from typing import Union

from bufr2geojson import transform as as_geojson
from wis2box.data.base import BaseAbstractData

LOGGER = logging.getLogger(__name__)


bearingToName = {
    "0-90": "RAD1",
    "90-180": "RAD2",
    "180-270": "RAD3",
    "270-0": "RAD4"
}

class BUFR2JTWC(BaseAbstractData):
    """Model data"""
    def transform(self, input_data: Union[Path, bytes], filename: str = '') -> bool:

        LOGGER.debug('Procesing BUFR data')
        input_bytes = self.as_bytes(input_data)

        LOGGER.debug('Generating GeoJSON features')
        results = as_geojson(input_bytes, serialize=False)

        LOGGER.debug('Processing GeoJSON features')
        for collection in results:
            # results is an iterator, for each iteration we have:
            # - dict['id']
            # - dict['id']['_meta']
            # - dict['id']
            JTWC = {}
            datarow = {
                "basin": "",
                "cyclone_number": None,
                "yyyymmddhh": None,
                "tech_number": None,
                "tech": "ECMF",
                "lat": None,
                "lon": None,
                "Vmax": None,
                "MSLP": None,
                "TY": "XX",
                "RAD": {}
            }
            rads = {
                "RAD": None,
                "wind_code": "NEQ",
                "RAD1": None,
                "RAD2": None,
                "RAD3": None,
                "RAD4": None
            }

            for id, item in collection.items():
                # extract identification
                stormIdentifier = item['geojson']['properties']['wigos_station_identifier']
                stormName, stormNumber = stormIdentifier.split("-")
                # now warning time and forecast (tau)
                warningTime = item['geojson']['properties']['phenomenonTime']
                forecastTime = None
                tau = None
                if "/" in warningTime:
                    forecastTime, tau = warningTime.split("/")
                else:
                    tau = item['geojson']['properties']['resultTime']
                forecastTime = datetime.strptime(forecastTime,"%Y-%m-%dT%H:%M:%SZ")
                tau = datetime.strptime(tau,"%Y-%m-%dT%H:%M:%SZ")
                tau = int((tau - forecastTime).total_seconds() / 3600)
                # now construct key for completing data.frame
                key1 = f"{stormIdentifier}-{forecastTime}-{tau}"

                if key1 not in JTWC:
                    JTWC[key1] = deepcopy(datarow)

                if item['geojson']['properties']['name'] == "pressure_reduced_to_mean_sea_level":
                    JTWC[key1]["MSLP"] = item['geojson']['properties']['value']
                    # this element also gives location
                    JTWC[key1]['lat'] = item['geojson']['geometry']['coordinates'][1]
                    JTWC[key1]['lon'] = item['geojson']['geometry']['coordinates'][2]
                    JTWC[key1]['cyclone_number'] = stormNumber
                    JTWC[key1]['yyyymmdd'] = forecastTime
                    JTWC[key1]['tau'] = tau
                elif item['geojson']['properties']['name'] == "wind_speed_at10m":
                    JTWC[key1]["Vmax"] = item['geojson']['properties']['value']
                elif item['geojson']['properties']['name'] == "effective_radius_with_respect_to_wind_speeds_above_threshold":  # noqc
                    bearing = None
                    threshold = None
                    for metadata in item['geojson']['properties']['metadata']:
                        if metadata['name'] == "bearing_or_azimuth":
                            bearing = metadata['value']
                            bearing = f"{bearing[0]}-{bearing[1]}"
                        elif metadata['name'] == "wind_speed_threshold":
                            threshold = metadata['value']
                    assert bearing is not None
                    assert threshold is not None
                    colname = bearingToName[bearing]
                    key2 = threshold

                    if key2 in JTWC[key1]['RAD']:
                        JTWC[key1]['RAD'][key2][colname] = item['geojson']['properties']['value']
                    else:
                        JTWC[key1]['RAD'][key2] = deepcopy(rads)
                        JTWC[key1]['RAD'][key2][colname] = item['geojson']['properties']['value']
            # end of items in collection
            # each item int JTCW contains 3 rows of dat file output and 3 objects

            for key1, data in JTWC.items():
                output_string = ''
                head = (data['basin'], data['cyclone_number'], data['yyyymmddhh'], data['tech_number'],
                        data['tech'], data['lat'], data['lon'], data['Vmax'], data['MSLP'], data['TY'])
                for key2, rad in data['RAD']:
                    quadrants = (rad['RAD'], rad["wind_code"], rad["RAD1"], rad["RAD2"], rad["RAD3"], rad["RAD4"])
                    # add line to output file
                    output_string += ', '.join(head + quadrants) + "\n"

                    # write geojson(s) for line (TODO)
                    # 1: MSLP, point
                    # 2: Vmax, point
                    # 3: quadrant, polygons
                self.output[key1] = {
                    'dat': output_string,
                    '_meta': {
                        'data_date': data['yyyymmddhh'],
                        'relative_filepath': self.get_local_filepath(data['yyyymmddhh'])
                    }
                }
        # end of collections

        LOGGER.debug('Successfully finished transforming BUFR data')
        return True

    def get_local_filepath(self, date_):
        yyyymmdd = date_[0:10]  # date_.strftime('%Y-%m-%d')
        return Path(yyyymmdd) / 'wis' / self.topic_hierarchy.dirpath