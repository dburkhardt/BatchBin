#! /usr/bin/python
import re,sys,argparse,subprocess,os

parser = argparse.ArgumentParser()
parser.add_argument('barcode_file',help='the barcode file')
parser.add_argument('samples', nargs='+',help='a list of codes of the form PHO (Prospect hill Heated)')
parser.add_argument('-b','--bamdir',dest='bamdir', default= '~/bamfiles',help= 'directory containing bamfiles, default is ~/bamfiles')
parser.add_argument('-a','--assembly_file',dest='assembly' default='~/binning_files/1018256.scaffolds.fasta', help='path to assembly.fa, default is ~/binning_files/1018256.scaffolds.fasta')

assembly=''
bamdir = ''
run = ''

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
def mergeBamFilesPopen(list_of_bamfiles, subset):
        if not os.path.exists('./mergedBamfiles/%s.tmp.bam'): #check to see if this particular set of files has been created
                return subprocess.Popen(['samtools', 'merge', './mergedBamfiles/%s.tmp.bam' % subset] + list_of_bamfiles) 

#returns a Popen object which will run the whole MetaBat pipeline
def runMetaBat(list_of_samples, list_of_merged_bamfiles,):
        makeDepthFile(list_of_merged_bamfiles) #will block until completed
        os.system('mkdir ./MetaBAT/Sensitive ./MetaBAT/Specific ./MetaBAT/SpecificPair')
        try:
                p = [subprocess.Popen('nice -n 15 metabat -i %s -a MetaBat.depth.txt -o ./MetaBAT/Sensitive/bin --sensitive -v --saveTNF saved.tnf --saveDistance saved.gprob' % assembly, Shell=True)]
                p.append(subprocess.Popen('nice -n 15 metabat -i %s -a MetaBat.depth.txt -o ./MetaBAT/Specific/bin --specific -v --saveTNF saved.tnf --saveDistance saved.gprob' % assembly, Shell=True))
                p.append(subprocess.Popen('nice -n 15 metabat -i %s -a MetaBat.depth.txt -p MetaBat.paired.txt -o ./MetaBAT/SpecificPair/bin --specific -v --saveTNF saved.tnf --saveDistance saved.gprob' % assembly, Shell=True))
        except KeyboardInterrupt:
                p.kill()
                sys.exit()


def makeDepthFile(list_of_merged_bamfiles):
        try:
                p = subprocess.Popen(['jgi_summarize_bam_contig_depths','--outputDepth','MetaBat.depth.txt','--pairedContigs','MetaBat.paired.txt'] + list_of_merged_bamfiles)
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
        global run
        rundir = '_'.join(args.samples)
        os.system('mkdir %s' % rundir)
        os.chdir('./%s' % rundir)
        os.getcwd()
        print "Merging bamfiles..."
        os.system('mkdir ./mergedBamfiles')
        try:
                processes = [mergeBamFilesPopen(get_barcodes(subset, barcode_table_asList),subset) for subset in args.samples]
                for p in processes: p.wait()
        except KeyboardInterrupt:
                for p in processes: p.kill()
                sys.exit()
        print "Successfully merged bamfiles!"
# at this point there exists in the working directory a set of files called "subset.tmp" for subset in args.samples
        
        #bin using MetaBat
        list_of_bamfiles = ['./mergedBamfiles/'+sample+'.tmp.bam' for sample in args.samples]
        process = runMetaBat(args.samples,list_of_bamfiles)


        print ('done!')

if __name__ == '__main__':
        main()
