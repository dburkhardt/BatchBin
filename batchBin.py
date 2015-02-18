#! /usr/bin/python
import argparse,os
from  runBinningPrograms import run_binning_pipeline

#run_binning_pipeline(argparser, samples) takes a single parameter 
#which is a and argparse object and returns a single process, 
#the one which takes longest to run

argparser = None

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
cwd = os.getcwd()
def run_processes(list_of_samples):
    try:
        for samples in list_of_samples:
            os.chdir(cwd)
            list_of_slow_processes.append(
                run_binning_pipeline(argparser, samples.split())
                )
            if len(list_of_slow_processes) > 0:
                print "Waiting..."
                list_of_slow_processes[0].wait()
    except KeyboardInterrupt:
        for p in list_of_slow_processes: p.kill()
        sys.exit()

def initializeArgparse():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--barcode_file', default = '~/barcode_table.tsv',
        help ='the barcode file'
        )
    parser.add_argument(
        '-b',
        '--bamdir', dest = 'bamdir', default = '~/bamfiles/',
        help = 'directory containing bamfiles, default is ~/bamfiles'
        )
    parser.add_argument(
        '-a','--assembly_file',dest='assembly', 
        default = '~/binning_files/1018256.scaffolds.fasta', 
        help = "path to assembly.fa, default is "
        "~/binning_files/1018256.scaffolds.fasta')"
    return parser.parse_args()


if __name__ == '__main__':
    argparser = initializeArgparse()
    run_processes(samples_to_run)
