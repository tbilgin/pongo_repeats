import sys
import os
from os import listdir
from os.path import isfile, join

def get_file_paths(dir_path):
	return [join(dir_path,file_path) for file_path in listdir(dir_path) if isfile(join(dir_path, file_path))]

def get_motif_to_str_count(dir):
	files = get_file_paths(dir)
	motif_length_strs = {}
	for file in files:
		print(file)
		with open(file) as f:
			next(f)
			for line in f:
				entries = line.split(' ')
				if entries[4] != '.':
					coord = (entries[0], entries[1], entries[2])
					motif_length = len(entries[5].replace('\n',''))
					if motif_length not in motif_length_strs.keys():
						motif_length_strs[motif_length] = set([coord])
					else:
						motif_length_strs[motif_length].update([coord])
	for k,v in motif_length_strs.items():
		print(str(k) + ',' + str(len(v)))
def main():
	dir = '../genome-wide-txt-hq-cov20-processed'
	get_motif_to_str_count(dir)
main()
