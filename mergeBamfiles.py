#! /usr/bin/python
import re,sys,argparse,subprocess,os

parser = argparse.ArgumentParser()
parser.add_argument('barcode_file',help='the barcode file')
parser.add_argument('samples', nargs='+',help='a list of codes of the form PHO (Prospect hill Heated)')
parser.add_argument('-b','--bamdir',dest='bamdir', default= '~/bamfiles',help= 'directory containing bamfiles, default is ~/bamfiles')

bamdir = ''

#currently uses '/home/dan/barcode_table.tsv'
def load_barcodeFile(filename):
        barcode_file = open(filename)
        barcode_file.seek(1)
        barcodes = barcode_file.read().splitlines()
        barcode_list = []
        for line in barcodes:
                barcode_list.append(line.split('\t'))

        del barcode_list[0]
        for line in barcode_list:
                del line[3]
        return barcode_list

#take a string ID like 'HOB' or 'PO' and returns a list of the barcodes (e.g. ['AAATTT','GGGCCC'])
def get_barcodes( c , barcode_list ):
        subset_codes = list(c)
        subset_ids = barcode_list
        for code in subset_codes:
                subset_ids = [record for record in subset_ids for entry in record if re.search(code,entry) and not re.search('[ACTG]{6}',entry)]
        return [record[0] for record in subset_ids]

def mergeBamFilesPopen(barcodesToMerge, setName):
        print [bamdir + '/' + barcode + '.bam' for barcode in barcodesToMerge]
	list_of_bamfiles = [os.path.expanduser(bamdir + '/' + barcode + '.bam') for barcode in barcodesToMerge]
        if not os.path.exists('./%s.tmp'): #check to see if this particular set of files has been created
                return subprocess.Popen(['samtools', 'merge', '%s.tmp' % setName] + list_of_bamfiles) 


def main():
        args = parser.parse_args()
        barcode_list = load_barcodeFile(args.barcode_file)
        global bamdir 
	bamdir = args.bamdir
	print "Merging bamfiles..."
	try:
	        processes = [mergeBamFilesPopen(get_barcodes(subset, barcode_list),subset) for subset in args.samples]
        	for p in processes: p.wait()
	except KeyboardInterrupt:
		for p in processes: p.kill()
		sys.exit()
	print "Successfully merged bamfiles!"

# at this point there exists in the working directory a set of files called "subset.tmp" for subset in args.samples
        print ('done!')

if __name__ == '__main__':
	main()
