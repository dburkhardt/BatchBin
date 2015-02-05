#! /usr/bin/python
import re,sys,argparse

parser = argparse.ArgumentParser()
parser.add_argument("barcode_file",dest='filename',help="the barcode file")
parser.add_argument("samples", nargs='+',help="a list of codes of the form PHO (Prospect hill Heated Organic)")


#currently uses "/home/dan/barcode_table.tsv"
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

#take a string ID like "HOB" or "PO" and returns the barcodes for all the samples matching those
def get_barcodes( c ):
	subset_codes = list(c)
	subset_ids = barcode_list
	for code in subset_codes:
		subset_ids = [record for record in subset_ids for entry in record if re.search(code,entry)]
	return [record[0] for record in subset_ids]

def main():
	args = parser.parse_args()
	barcode_list = load_barcodeFile(args.filename)
	for subset in args.samples:
		print "Barodes for subset %s" % subset
		print get_barcodes(subset)

if __name__ == "__main__":
	main()