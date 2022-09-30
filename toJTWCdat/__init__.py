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
        for collection in results:
            # results is an iterator, for each iteration we have:
            # - dict['id']
            # - dict['id']['_meta']
            # - dict['id']
            JTWC = {}
            datarow = {
                "basin": "",
                "cyclone_number": "",
                "yyyymmddhh": "",
                "tech_number": "00",
                "tech": "ECMF",
                "clat": "",
                "clon": "",
                "flat": "",
                "flon": "",
                "Vmax": "",
                "MSLP": "",
                "TY": "XX",
                "RAD": {},
                "MSLP_metadata": None,
                "Vmax_metadata": None
            }
            rads = {
                "RAD": "",
                "wind_code": "NEQ",
                "RAD1": "",
                "RAD2": "",
                "RAD3": "",
                "RAD4": ""
            }
            geojsons_out = {}
            for id, item in collection.items():
                #LOGGER.debug(item)
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
                #LOGGER.debug(warningTime)
                forecastTime = None
                tau = None
                if "/" in warningTime:
                    forecastTime, tau = warningTime.split("/")
                else:
                    forecastTime = warningTime
                    tau = item['geojson']['properties']['resultTime']
                forecastTime = datetime.strptime(forecastTime,"%Y-%m-%dT%H:%M:%SZ")
                #LOGGER.debug("Calculating tau")
                tau = datetime.strptime(tau,"%Y-%m-%dT%H:%M:%SZ")
                tau = int((tau - forecastTime).total_seconds() / 3600)
                # now construct key for completing data.frame
                key1 = f"{stormIdentifier}-{subset}-{forecastTime}-{tau}"

                if key1 not in JTWC:
                    JTWC[key1] = deepcopy(datarow)
                if key1 not in geojsons_out:
                    geojsons_out[key1] = {}

                #LOGGER.debug(f"Processing {item['geojson']['properties']['name']}")
                if item['geojson']['properties']['name'] == "pressure_reduced_to_mean_sea_level":
                    geojsons_out[key1]["MSLP"] = self.extract_MSLP(item['geojson'])
                    JTWC[key1]['MSLP_metadata'] = deepcopy(item['geojson']['properties']['metadata'])
                    JTWC[key1]['MSLP_geometry'] = deepcopy(item['geojson']["geometry"])
                    JTWC[key1]["MSLP"] = int( item['geojson']['properties']['value'] )
                    # this element also gives location
                    lat = item['geojson']['geometry']['coordinates'][1]
                    if lat >=0:
                        sign = "N"
                    else:
                        sign = "S"
                    clat = int(abs(lat * 10))
                    JTWC[key1]['clat'] = f"{clat}{sign}"
                    JTWC[key1]['lat'] = lat
                    lon = item['geojson']['geometry']['coordinates'][0]
                    if lon >=0:
                        sign = "E"
                    else:
                        sign = "W"
                    clon = int(abs(lon * 10))
                    JTWC[key1]['clon'] = f"{clon}{sign}"
                    JTWC[key1]['lon'] = lon
                    JTWC[key1]['basin'] = stormNumber[2]
                    JTWC[key1]['cyclone_number'] = stormNumber[0:2]  # characters 1 - 2, 3rd is letter
                    JTWC[key1]['yyyymmddhh'] = forecastTime.strftime("%Y%m%d%H")
                    JTWC[key1]['tau'] = tau
                    #LOGGER.debug(JTWC[key1])
                elif item['geojson']['properties']['name'] == "wind_speed_at10m":
                    geojsons_out[key1]["vmax"] = self.extract_Vmax(item['geojson'])
                    JTWC[key1]["Vmax"] = int(item['geojson']['properties']['value']/0.51444444)
                    JTWC[key1]['Vmax_metadata'] = deepcopy(item['geojson']['properties']['metadata'])
                    JTWC[key1]['Vmax_geometry'] = deepcopy(item['geojson']["geometry"])
                elif item['geojson']['properties']['name'] == "effective_radius_with_respect_to_wind_speeds_above_threshold":  # noqc
                    bearing = None
                    threshold = None
                    #LOGGER.debug("Iterating over metadata")
                    #LOGGER.debug(item['geojson']['properties']['metadata'])
                    for metadata in item['geojson']['properties']['metadata']:
                        if metadata['name'] == "bearing_or_azimuth":
                            bearing = metadata['value']
                            bearing = f"{bearing[0]}-{bearing[1]}"
                        elif metadata['name'] == "wind_speed_threshold":
                            threshold = metadata['value']
                    assert bearing is not None
                    assert threshold is not None
                    colname = bearingToName[bearing]
                    key2 = int(threshold/0.51444444)
                    val = int(item['geojson']['properties']['value'] * 0.000539957)
                    if key2 in JTWC[key1]['RAD']:
                        JTWC[key1]['RAD'][key2][colname] = val
                    else:
                        JTWC[key1]['RAD'][key2] = deepcopy(rads)
                        JTWC[key1]['RAD'][key2][colname] = val
                    geojsons_out[key1][f"{key2}-{colname}"] = self.extract_wind_polygon(item['geojson'])

            # end of items in collection
            # each item int JTCW contains 3 rows of dat file output and 3 objects
            #LOGGER.debug("Building output string")
            for key1, data in JTWC.items():
                output_string = ''
                head = (f"{data['basin']:>2}",
                        f"{data['cyclone_number']:>2}",
                        f"{data['yyyymmddhh']:>10}",
                        f"{data['tech_number']:>2}",
                        f"{data['tech']:>4}",
                        f"{data['tau']:>3}",
                        f"{data['lat']:>4}",
                        f"{data['lon']:>5}",
                        f"{data['Vmax']:>3}",
                        f"{data['MSLP']:>4}",
                        f"{data['TY']:>2}")
                #LOGGER.debug(f"{head}")
                for key2, rad in data['RAD'].items():
                    #LOGGER.debug(key2)
                    #LOGGER.debug(rad)
                    quadrants = (f"{key2:>3}",
                                 f"{rad['wind_code']:>3}",
                                 f"{rad['RAD1']:>4}",
                                 f"{rad['RAD2']:>4}",
                                 f"{rad['RAD3']:>4}",
                                 f"{rad['RAD4']:>4}",
                                 "")
                    # add line to output file
                    output_string += ', '.join(head + quadrants) + "\n"

                    # write geojson(s) for line (TODO)
                    # 1: MSLP, point
                    # 2: Vmax, point
                    # 3: quadrant, polygons

                # add .dats to output
                self.output_data[key1] = {
                    'dat': output_string,
                    '_meta': {
                        'data_date': data['yyyymmddhh'],
                        'relative_filepath': self.get_local_filepath(data['yyyymmddhh'])
                    }
                }

            # now geojsons
            for key1, obj1 in geojsons_out.items():
                for key2, obj2 in obj1.items():
                    self.output_data[f"{key1}-{key2}"] = {
                        'geojson': json.dumps(obj2),
                        '_meta': {
                            'data_date': data['yyyymmddhh'],
                            'relative_filepath': self.get_local_filepath(data['yyyymmddhh'])
                        }
                    }
                    LOGGER.debug(self.output_data[f"{key1}-{key2}"])
                    LOGGER.debug(self.output_data[key1])

        LOGGER.debug('Successfully finished transforming BUFR data')
        return True

    def get_local_filepath(self, date_):
        yyyymmdd = date_[0:10]  # date_.strftime('%Y-%m-%d')
        return Path(yyyymmdd) / 'wis' / self.topic_hierarchy.dirpath


    def extract_Vmax(self, feature):
        #LOGGER.debug("Extracting Vmax as GeoJSON")
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
        #LOGGER.debug(feature)
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
        #LOGGER.debug("Extracting polygons")
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

