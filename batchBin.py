#! /usr/bin/python

from  runBinningPrograms import merge_and_run_binning_programs

#merge_and_run_binning_programs takes a single parameter samples which is a list of the form ['MPH', 'MPC', ... , 'OPC'] and returns a single process, the one which takes longest to run

samples_to_run = [
	'MPH MPC MSH MSC MBH MBC OPH OPC OSH OSC OBH OBC', 
	'MPH MPC MSH MSC MBH MBC',
	'OPH OPC OSH OSC OBH OBC',
	'MP MS MB OP OS OB',
	'PH SH BH PO SO BO',
	'PH SH BH',
	'PO SO BO',
	'MP MS MB',
	'OP OS OB',
	'MH OH',
	'M O'
	]

list_of_slow_processes = []

def run_processes(list_of_samples):
	for samples in list_of_samples:
		list_of_slow_processes.append(merge_and_run_binning_programs(samples))

if __name__ == '__main__':
	run_processes(samples_to_run)