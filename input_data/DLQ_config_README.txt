gap.time: Distinct events that are separated by gap.time minutes or less will be combined into one event.

length.threshold: Events shorter than or equal to length.threshold minutes will be discarded.

do.event.detection: Set to T if you want to estimate emission start and end times, and then perform localization and quantification on the resulting events. Set to F if you want to perform localization and quantification on non-overlapping 30-minute intervals.

first.sim.ind: Index (starting at 1) of the first simulation file in the "data" object.

n.samples: Number of times to sample from predictions and observations.

forward.model.path: Location (including file name) of simulation output (output of MAIN_1_simulate.R)

output.file.path: Location (including file name) of output data file that will contain the event detection, localization, and quantification results.

helper.spike.detection.alg.path: Location (including file name) of the HELPER file that contains the spike detection algorithm. Default name is "HELPER_spike_detection_algorithm.R"
