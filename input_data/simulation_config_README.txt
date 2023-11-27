num.cores.to.use: Number of cores to parallelize over. Should be less than or equal to the number of cores on your machine.

dt: Simulation frequency. There will be dt seconds between subsequent puffs. dt = 1 is recommended. NOTE: 60/dt must be an even number, so viable options are dt=1,2,3,4,5,6,10,12,15,20,30,60 or fractions.

cutoff.t: Time after which puffs are ignored. A good rule of thumb is to take the 5th percentile of wind speeds during your experiment, and divide the diameter of your study area by this wind speed value. In minutes.

run.mode: Should the simulation be run in constant mode or dynamic mode. Input options are: "constant" or "dynamic". Constant mode assumes that puff characteristics (dispersion parameters, direction of travel, speed of travel) are fixed based on atmospheric conditions when the puff was created. Dynamic mode allows for the puff characteristics to be updated as the puff moves. Constant mode is used in the paper as it can be derived from the advection-diffusion equation. Dynamic mode cannot be rigorously derived from the governing PDE. Experimentation has shown that there is little difference between the two.

emission.rate: Assumed emission rate in kilograms per second. Default is to simulate at an arbitrary emission rate of 1 gram per second, as the event detection, localization, and quantification framework is designed to run on any arbitrary but constant emission rate.

start.time: Start time of simulation. Should be in format: YYYY-MM-DDTHH:MM:00, where Y is year, M is month, D is day, H is hour, M is minute, and T is not a variable but rather stands for "time". Note that simulation must be started on an even minute (no seconds).

end.time: End time of simulation. Should be in format: YYYY-MM-DDTHH:MM:00, where Y is year, M is month, D is day, H is hour, M is minute, and T is not a variable but rather stands for "time". Note that simulation must be started on an even minute (no seconds).

tz: Local time zone for simulation. Use the "area/location" naming convention. See https://en.wikipedia.org/wiki/List_of_tz_database_time_zones

raw.sensor.observations.path: Location (including file name) of continuous monitoring sensor data

source.locations.path: Location (including file name) of potential emission sources. Multiple sources can be listed in the csv file, and a separate simulation will be run for each.

sensor.locations.path: Location (including file name) of continuous monitoring sensors.

output.file.path: Location (including file name) of output data file that contains the simulation predictions.

helper.distance.conversions.path: Location (including file name) of the HELPER file that contains the distance conversion functions. Default name is "HELPER_distance_conversions.R"

helper.gpuff.function.path: Location (including file name) of the HELPER file that contains the Gaussian puff functions. Default name is "HELPER_gpuff_function.R"