#! /usr/bin/python
import re,sys,argparse,subprocess,os,time

parser=None
assembly=''
bamdir= ''
rundir= '' #Will eventually look like MPH_MPC_MSH

def initializeArgparse():
        global parser
        parser = argparse.ArgumentParser()
        parser.add_argument('--barcode_file',default = '~/barcode_table.tsv',help='the barcode file')
        parser.add_argument('samples', nargs='+',help='a list of codes of the form PHO (Prospect hill Heated)')
        parser.add_argument('-b','--bamdir',dest='bamdir', default= '~/bamfiles/',help= 'directory containing bamfiles, default is ~/bamfiles')
        parser.add_argument('-a','--assembly_file',dest='assembly', default='~/binning_files/1018256.scaffolds.fasta', help='path to assembly.fa, default is ~/binning_files/1018256.scaffolds.fasta')
        return parser.parse_args()

#parses the arguments and assigns values to the global variables
def initializeVariables(args):
        print "Binning: " + ' '.join(args.samples)
        barcode_table_asList = load_barcodeFile(args.barcode_file)
        global bamdir
        bamdir = args.bamdir
        global assembly
        assembly = args.assembly
        global rundir
        rundir = '_'.join(args.samples)
        os.system('mkdir -p %s' % rundir)
        os.chdir('./%s' % rundir)
        print "Placing files in " + os.getcwd()
        return args

#currently uses '/home/dan/barcode_table.tsv'
def load_barcodeFile(filename):
        barcode_file = open(os.path.expanduser(filename))
        barcode_file.seek(1)
        barcodes = barcode_file.read().splitlines()
        barcode_table_asList = []
        for line in barcodes:
                barcode_table_asList.append(line.split('\t'))

        del barcode_table_asList[0]
        for line in barcode_table_asList:
                del line[3]
        return barcode_table_asList

#take a string ID like 'HOB' or 'PO' and returns a list of the barcodes (e.g. ['AAATTT','GGGCCC'])
def get_barcodes( subset , barcode_table_asList ):
        subset_list = list(subset) #turns 'MPH' (subset) into ['M','P','H'] (subset_list)
        subset_of_barcode_table = barcode_table_asList 
        for code in subset_list:
                subset_of_barcode_table = [record for record in subset_of_barcode_table for entry in record if re.search(code,entry) and not re.search('[ACTG]{6}',entry)]
        return [record[0] for record in subset_of_barcode_table] #returns only the first column of the table

#returns a single Popen object running samtools merge on a set of bamfiles
def mergeBamFilesPopen(list_of_barcodes, subset):
        list_of_bamfiles = ['%s%s.bam' % (bamdir,barcode) for barcode in list_of_barcodes]
        if not os.path.exists('./mergedBamfiles/%s.tmp.bam' % subset): #check to see if this particular set of files has been created
                return subprocess.Popen('samtools merge ./mergedBamfiles/%s.tmp.bam ' % subset + ' '.join(list_of_bamfiles), shell=True)

def indexBamFilesPopen(list_of_barcodes, subset):
        list_of_bamfiles = ['%s/%s.bam' % (bamdir,barcode) for barcode in list_of_barcodes]
        if not os.path.exists('./mergedBamfiles/%s.tmp.bam.bai' % subset): #check to see if this particular set of files has been created
                return subprocess.Popen('samtools index ./mergedBamfiles/%s.tmp.bam' % subset, shell=True)

#spawns a bunch of samtools merge subprocesses and writes merged bamfiles to rundir/mergedBamfiles/ (BLOCKING)
def mergeBamfiles(samples, barcode_table_asList):
        print "Merging bamfiles"
        os.system('mkdir -p ./mergedBamfiles')
        list_of_bamfiles = ['./mergedBamfiles/'+sample+'.tmp.bam' for sample in samples]
        print list_of_bamfiles
        try:
                print [subset for subset in samples]
                print [get_barcodes(subset, barcode_table_asList) for subset in samples]
                processes = [mergeBamFilesPopen(get_barcodes(subset, barcode_table_asList),subset) for subset in samples]
                if processes[0] is None:
                        print 'Found merged bamfiles in ./%s/mergedBamfiles' % rundir
                        return list_of_bamfiles
                else:
                        print "Merging bamfiles..."
                        for p in processes:
                                p.wait()
                        for p in processes:
                                if p.returncode is not 0:
                                        print "Merging failed"
                                        sys.exit()
                                else: 
                                        print "Successfully merged bamfiles!"
                                        return list_of_bamfiles
        except KeyboardInterrupt:
                for p in processes: p.kill()
                sys.exit()

#spawns a bunch of samtools merge subprocesses and writes merged bamfiles to rundir/mergedBamfiles/ (BLOCKING)
def indexBamfiles(samples, barcode_table_asList):
        list_of_bamfiles = ['./mergedBamfiles/'+sample+'.tmp.bam' for sample in samples]
        try:
                processes = [indexBamFilesPopen(get_barcodes(subset, barcode_table_asList),subset) for subset in samples]
                if processes[0] is None:
                        print 'Found indexed bamfiles in ./%s/mergedBamfiles' % rundir
                        return list_of_bamfiles
                else:
                        print "Indexing bamfiles..."
                        for p in processes:
                                p.wait()
                                if p.returncode is not 0:
                                        print "Indexing failed"
                                        sys.exit()
#                                else:
                                        #print "Successfully indexed bamfiles!"
        except KeyboardInterrupt:
                for p in processes: p.kill()
                       sys.exit()

#returns a List of Popen objects running MetaBat using the Specific Sensitive and Specific_Paired settings
def runMetaBat(list_of_samples, list_of_merged_bamfiles):
        print 'Binning using MetaBat...' 
        makeDepthFile(list_of_merged_bamfiles) #will block until completed
        os.system('mkdir -p ./MetaBAT/Sensitive ./MetaBAT/Specific ./MetaBAT/SpecificPair')
        try:
                log_sensitive = open('./MetaBAT/Sensitive/MetaBat_sensititve.log','a')
                log_specific = open('./MetaBAT/Specific/MetaBat_specific.log','a')
                log_specific_paired = open('./MetaBAT/SpecificPair/MetaBat_SpecificPair.log','a')
                p = [subprocess.Popen('time nice -n 15 metabat -i %s -a MetaBat.depth.txt -o ./MetaBAT/Sensitive/bin --sensitive -v --saveTNF saved.tnf --saveDistance saved.gprob' % assembly,
                        stdout=log_sensitive, stderr=log_sensitive, shell=True)]
                p.append(subprocess.Popen('time nice -n 15 metabat -i %s -a MetaBat.depth.txt -o ./MetaBAT/Specific/bin --specific -v --saveTNF saved.tnf --saveDistance saved.gprob' % assembly,
                        stdout=log_specific, stderr=log_specific, shell=True))
                p.append(subprocess.Popen('time nice -n 15 metabat -i %s -a MetaBat.depth.txt -p MetaBat.paired.txt -o ./MetaBAT/SpecificPair/bin --specific -v --saveTNF saved.tnf --saveDistance saved.gprob' % assembly, 
                        stdout=log_specific_paired, stderr=log_specific_paired, shell=True, executable='/bin/bash'))
                return p
        except KeyboardInterrupt:
                for a in p: a.kill()
                sys.exit()

#returns a Popen object with a single subprocess running Concoct
def runConcoct():
        print 'Processing depth file for Concoct'
        os.system('mkdir -p ./Concoct ./Concoct/bins')
        log_concoct = open('./Concoct/concoct.log','a')
        if not os.path.exists('./depth_concoct.txt'):
                subprocess.check_call("awk \'NR > 1 {for(x=1;x<=NF;x++) if(x == 1 || (x >= 4 && x % 2 == 0)) printf \"%s\", $x (x == NF || x == (NF-1) ? \"\\n\":\"\\t\")}' ./MetaBat.depth.txt > ./depth_concoct.txt", shell=True)
        else: print "found Concoct depth file"
        call_to_source_concoct_env = ". ~/software/anaconda/envs/concoct_env/bin/activate concoct_env"
        call_to_bin_with_concoct = "concoct --composition_file %s --coverage_file ../depth_concoct.txt" % assembly
        extract_bins = "/home/dan/software/CONCOCT-0.4.0/scripts/extract_fasta_bins.py --output_path ./bins %s clustering_gt1000.csv" % assembly
        concoct = subprocess.Popen('cd Concoct ; %s ; time nice -n 15 %s ; time nice -n 15 %s' % (call_to_source_concoct_env, call_to_bin_with_concoct, extract_bins),
                stdout=log_concoct, stderr=log_concoct, shell=True, executable='/bin/bash')
        print 'Binning using CONCOCT'
        return concoct

#returns a Popen object with a single subprocess running GroopM pipeline
def runGroopM(list_of_bamfiles):
        database = "%s.gm" % rundir
        os.system('mkdir -p ./GroopM')
        log_groopM = open('./GroopM/groopm.log','a')
        parse = "groopm parse -t 16 %s %s " % (database, assembly) + ' '.join(list_of_bamfiles)
        core = "groopm core %s " % database
        extract = "groopm extract -t 32 --prefix ./GroopM/core_only/bin_groopm %s %s" % (database , assembly)
        recruit = "groopm recruit %s" % database
        extract_recruit = "groopm extract -t 32 --prefix ./GroopM/recruited/bin_groopm %s %s" % (database , assembly)
        print "Binning using GroopM..."
        return subprocess.Popen('time nice -n 15 %s ; time nice -n 15 %s ; time nice -n 15 %s ; time nice -n 15 %s ; time nice -n 15 %s' % (parse, core, extract, recruit, extract_recruit),
                stdout = log_groopM, stderr=log_groopM,shell = True , executable = '/bin/bash')

#creates a file called 'MetaBat.depth.txt' in the cwd (BLOCKING)
def makeDepthFile(list_of_merged_bamfiles):
        try:
                if os.path.exists('./MetaBat.depth.txt'): 
                        print 'Found MetaBat depth file. Skipping...'
                        return
                else: 
                        depth_log = open("./depth.log",'a')
                        p = subprocess.Popen([
                                'jgi_summarize_bam_contig_depths',
                                '--outputDepth',
                                'MetaBat.depth.txt',
                                '--pairedContigs',
                                'MetaBat.paired.txt'] + list_of_merged_bamfiles, 
                                stdout=depth_log, stderr=depth_log)
                        p.wait()
        except KeyboardInterrupt:
                       p.kill()
                       sys.exit()

def monitorProcesses_returnLast(list_of_processes):
        num_processes = len(list_of_processes)
        done_procs = 0
        try:
                while (num_processes - done_procs) > 1:
                        for p in list_of_processes:
                                if p.poll() is not None:
                                        list_of_processes.remove(p)
                                        done_procs+=1
                        time.sleep(5)
                return list_of_processes[0]
        except KeyboardInterrupt:
                for p in list_of_processes: p.kill
                sys.exit()

def merge_and_run_binning_programs(samples, barcode_table_asList):
        #bin using MetaBat
        list_of_bamfiles = mergeBamfiles(samples, barcode_table_asList)
        indexBamfiles(samples, barcode_table_asList)
        try:
                metabat_process_list = runMetaBat(samples,list_of_bamfiles)
                concoct_process = runConcoct() # will take the longest
                if len(list_of_bamfiles) >= 3:
                        groopm_process = runGroopM(list_of_bamfiles)
                else:
                        print "Too few bamfiles to bin using GroopM. Skipping..."
                        groopm_process = subprocess.Popen(['sleep', '0'])
                all_processes = []
                all_processes.extend(metabat_process_list)
                all_processes.append(concoct_process)
                all_processes.append(groopm_process)
                return monitorProcesses_returnLast(all_processes) #TODO put keyboard interrupt in here, handle end of function processes
        except KeyboardInterrupt:
                for p in all_processes: p.kill
                sys.exit()
        print ('done!')

def run_binning_pipeline(argparser, samples):
        argparser.samples = samples
        initializeVariables(argparser)
        barcode_table_asList = load_barcodeFile(argparser.barcode_file)
        return merge_and_run_binning_programs(samples, barcode_table_asList)

if __name__ == '__main__':
        argparser = initializeArgparse()
        initializeVariables(argparser)
        barcode_table_asList = load_barcodeFile(argparser.barcode_file)
        merge_and_run_binning_programs(argparser.samples, barcode_table_asList).wait()
        print "Done binning!"
