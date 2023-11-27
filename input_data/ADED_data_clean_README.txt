time: Time of observation in UTC time zone. Format can be anything recognized by the "as_datetime" function from the "lubridate" package in R. Note that the simulation file currently only supports data collected every minute. Time values need to be rounded to the minute (no seconds).

wind.direction: Wind direction measurements. The angle is defined as the direction the wind is coming from, with zero degrees meaning wind coming from the north and moving clockwise. The code can handle NA values.

wind.speed: Wind speed measurements in meters per second. The code can handle NA values.

methane: Methane concentration measurements in parts per million. The code can handle NA values.

name: Uniquely defined sensor names. Needs to be the same as the names in the "sensor_locations.csv" file.