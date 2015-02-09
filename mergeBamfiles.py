#! /usr/bin/python
import re,sys,argparse,subprocess,os

parser = argparse.ArgumentParser()
parser.add_argument('barcode_file',help='the barcode file')
parser.add_argument('samples', nargs='+',help='a list of codes of the form PHO (Prospect hill Heated)')
parser.add_argument('-b','--bamdir',dest='bamdir', default= '~/bamfiles',help= 'directory containing bamfiles, default is ~/bamfiles')
parser.add_argument('-a','--assembly_file',dest='assembly', default='~/binning_files/1018256.scaffolds.fasta', help='path to assembly.fa, default is ~/binning_files/1018256.scaffolds.fasta')

assembly=''
bamdir = ''
rundir = ''

#currently uses '/home/dan/barcode_table.tsv'
def load_barcodeFile(filename):
        barcode_file = open(filename)
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
                return subprocess.Popen(['samtools', 'merge', './mergedBamfiles/%s.tmp.bam' % subset] + list_of_bamfiles)

#returns a Popen object which will run the whole MetaBat pipeline
def runMetaBat(list_of_samples, list_of_merged_bamfiles):
        makeDepthFile(list_of_merged_bamfiles) #will block until completed
        os.system('mkdir -p ./MetaBAT/Sensitive ./MetaBAT/Specific ./MetaBAT/SpecificPair')
        try:
        	log_sensitive = open('./MetaBAT/Sensitive/MetaBat_sensititve.log','a')
		log_specific = open('./MetaBAT/Specific/MetaBat_specific.log','a')
		log_specific_paired = open('./MetaBAT/SpecificPair/MetaBat_SpecificPair.log','a')
	        p = [subprocess.Popen('nice -n 15 metabat -i %s -a MetaBat.depth.txt -o ./MetaBAT/Sensitive/bin --sensitive -v --saveTNF saved.tnf --saveDistance saved.gprob' % assembly, stdout=log_sensitive, stderr=log_sensitive, shell=True)]
                p.append(subprocess.Popen('nice -n 15 metabat -i %s -a MetaBat.depth.txt -o ./MetaBAT/Specific/bin --specific -v --saveTNF saved.tnf --saveDistance saved.gprob' % assembly, stdout=log_specific, stderr=log_specific, shell=True))
                p.append(subprocess.Popen('nice -n 15 metabat -i %s -a MetaBat.depth.txt -p MetaBat.paired.txt -o ./MetaBAT/SpecificPair/bin --specific -v --saveTNF saved.tnf --saveDistance saved.gprob' % assembly, stdout=log_specific_paired, stderr=log_specific_paired, shell=True))
		print 'Binning using MetaBAT...'
		return p
        except KeyboardInterrupt:
                for a in p: a.kill()
                sys.exit()


def makeDepthFile(list_of_merged_bamfiles):
        try:
		if os.path.exists('./MetaBat.depth.txt'): 
			print 'Found MetaBat depth file. Skipping...'
			return
		p = subprocess.Popen([
			'jgi_summarize_bam_contig_depths',
			'--outputDepth',
			'MetaBat.depth.txt',
			'--pairedContigs',
			'MetaBat.paired.txt'] + list_of_merged_bamfiles)
                p.wait()
        except KeyboardInterrupt:
                p.kill()
                sys.exit()




def main():
        args = parser.parse_args()
        barcode_table_asList = load_barcodeFile(args.barcode_file)
        global bamdir
        bamdir = args.bamdir
        global assembly
        assembly = args.assembly
        global rundir
        rundir = '_'.join(args.samples)
        os.system('mkdir -p %s' % rundir)
        os.chdir('./%s' % rundir)
        print os.getcwd()
        print "Merging bamfiles..."
        os.system('mkdir -p ./mergedBamfiles')
        try:
                processes = [mergeBamFilesPopen(get_barcodes(subset, barcode_table_asList),subset) for subset in args.samples]
                if processes[0] is None: print 'Found merged bamfiles in ./%s/mergedBamfiles' % rundir
		else:
			for p in processes: p.wait()
			for p in processes:
				if p.returncode is not 0:
					print "Merging failed"
					sys.exit()
				else: 
					print "Successfully merged bamfiles!"
        except KeyboardInterrupt:
                for p in processes: p.kill()
                sys.exit()
        
        #bin using MetaBat
        list_of_bamfiles = ['./mergedBamfiles/'+sample+'.tmp.bam' for sample in args.samples]
        
	metabat_processes = runMetaBat(args.samples,list_of_bamfiles)
	print 'Binning using MetaBat...' # TODO put all the above in one large try catch block
	#bin using Concoct
	print 'Processing depth file for Concoct'
	subprocess.check_call("awk 'NR > 1 {for(x=1;x<=NF;x++) if(x == 1 || (x >= 4 && x % 2 == 0)) printf "%s", $x (x == NF || x == (NF-1) ? "\n":"\t")}' ./MetaBat.depth.txt > ./depth_concoct.txt", shell=True)
	call_to_source_concoct_env = ". ./software/anaconda/envs/concoct_env/bin/activate concoct_env"
	call_to_bin_with_concoct = "concoct --composition_file %s --coverage_file depth_concoct.txt" % assembly
	concoct = subprocess.Popen('%s ; %s' % (call_to_source_concoct_env, call_to_bin_with_concoct), shell=True, executable='/bin/bash')
	print 'Binning using CONCOCT'

        print ('done!')

if __name__ == '__main__':
        main()
