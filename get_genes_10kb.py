import sys
from bs4 import BeautifulSoup
from urllib import urlopen
import os
from os import listdir
from os.path import isfile, join
import operator 
from collections import defaultdict

def get_file_paths(dir_path):
	return [join(dir_path,file_path) for file_path in listdir(dir_path) if isfile(join(dir_path, file_path))]


def get_variant_genes(file_path):
	file_content = open(file_path).readlines()[1:]
	genome_name = file_path.split("/")[-1].replace("_10kb.txt","")
	variant_genes = set()
	for entry in file_content:
		elements = entry.split()
		gene_name = elements[-2]
		variant = elements[4]
		if variant != '.' and gene_name != 'NA':
			variant_genes.add(gene_name)
	return genome_name, variant_genes
	
def get_island_genes(file_paths, island_genomes):
	island_genes = set()
	for file in file_paths:
		genes = set()
		genome = file.split('/')[-1].replace('_10kb.txt','')
		if genome in island_genomes:
			file_content = open(file).readlines()[1:]
			for line in file_content:
				entries = line.split(' ')
				if entries[4] != '.':
					chrom_nr = entries[0].replace('chr','')
					start = int(entries[1])
					end = int(entries[2])
					hugo_human_id = entries[9].replace('\n','')
					ensembl_human_id = entries[8]
#					if ensembl_human_id != 'NA':
#						genes.update([ensembl_human_id])
					dist = int(entries[5])
					if dist <= 10000 and ensembl_human_id != 'NA':
						genes.update([ensembl_human_id])
		island_genes.update(genes)
	return island_genes


def get_top_level_pathways_count_and_genes(file_path):
	pathway_id_to_genes = get_pathway_id_to_genes(file_path)
	top_level_pathways_to_genes = {}
	top_level_pathways_count = defaultdict(int)
	for pathway_id in pathway_id_to_genes.keys():
		top_level_pathway = get_top_level_pathway(pathway_id)
		top_level_pathways_count[top_level_pathway] += 1
		if top_level_pathway not in top_level_pathways_to_genes.keys():
			top_level_pathways_to_genes[top_level_pathway] = set(pathway_id_to_genes[pathway_id])
		else:
			top_level_pathways_to_genes[top_level_pathway].update(pathway_id_to_genes[pathway_id])
	return top_level_pathways_count, top_level_pathways_to_genes

def get_pathway_id_to_genes(file_path):
	file_content = open(file_path).readlines()[1:]
	pathway_id_to_genes = {}
	for line in file_content:
		entries = line.split(',\"')
		pathway_id = entries[0]
		genes = entries[2]
		pathway_id_to_genes[pathway_id] = genes.replace('\"','').split(';')
	return pathway_id_to_genes


def get_top_level_pathway(new_reactome_id):
#	new_reactome_id = "R-HSA-111448"
	query_id = new_reactome_id.split('-')[-1]
	url = "http://www.reactome.org/content/detail/" + query_id
	html = BeautifulSoup(urlopen(url), "html.parser")
	try:
		pathway_tag = html.find_all('span',class_="plus tree-root")[0]
		return pathway_tag.text.replace("\n","").strip().replace('(Homo sapiens)', '')
	except IndexError:
		for item in html.find_all('a'):
			href = item['href']
			if href.startswith('/PathwayBrowser/#/R-HSA'):
				return item.text.replace('(Homo sapiens)', '')
	
def main():

	dir = '10kb-gw'
	file_paths = get_file_paths(dir)
	
	borneo_genomes = ['PP_KB5404','PP_A940','PP_A941','PP_A943','PP_A944','PP_A938', 'PP_A983', 'PP_A939', 'PP_A942', 'PP_A946','PP_A987', 'PP_A988', 'PP_5062', 'PP_A989', 'PP_A984', 'PP_A985']
	sumatra_genomes = ['PA_A947', 'PA_A948', 'PA_A950', 'PA_A952', 'PA_A949','PA_B018', 'PA_B017','PA_B020', 'PA_A953','PA_A955','PA_A964','PA_B019']

	files = get_file_paths(dir)
	sumatra_genes = []
	borneo_genes = []
	for file in files:
		genome_name, variant_genes = get_variant_genes(file)
		print(genome_name)
		if len(variant_genes) > 0:
			if genome_name in sumatra_genomes:
				sumatra_genes.append(list(set(variant_genes)))  #set instead of list
			else:
				borneo_genes.append(list(set(variant_genes)))
				
	sumatra_genes = reduce(operator.add, sumatra_genes)
	borneo_genes = reduce(operator.add, borneo_genes)
	sumatra_dict = defaultdict(int)
	borneo_dict = defaultdict(int)
	for gene in sumatra_genes:
		sumatra_dict[gene] += 1
	for gene in borneo_genes:
		borneo_dict[gene] += 1

	sumatra_min2 = {k:v for k,v in sumatra_dict.items() if v >= 2}
	f1 = open('SumatraVariantGenes10kb_2_occurrences.txt', 'w')
	for k,v in sumatra_min2.items():
		f1.write(k + '\n')
	f1.close()
	borneo_min2 = {k:v for k,v in borneo_dict.items() if v >= 2}
	f2 = open('BorneoVariantGenes10kb_2_occurrences.txt', 'w')
	for k,v in borneo_min2.items():
		f2.write(k + '\n')
	f2.close()

#	borneo_genes = get_island_genes(file_paths, borneo_genomes)
#	b = open('BorneoUniqueGenesEnsembl10kb.txt', 'w')
#	for gene in borneo_genes:
#		b.write(gene + '\n')
#	b.close()
	
#	top_level_pathways_count, top_level_pathways_to_genes = get_top_level_pathways_count_and_genes('Borneo5kbPathways.csv')
#	print(top_level_pathways_count)
#	print(top_level_pathways_to_genes)

#	sumatra_genes = get_island_genes(file_paths, sumatra_genomes)
#	s = open('SumatraUniqueGenesEnsembl10kb.txt', 'w')
#	for gene in sumatra_genes:
#		s.write(gene + '\n')
#	s.close()

#	top_level_pathways_count, top_level_pathways_to_genes = get_top_level_pathways_count_and_genes('Sumatra5kbPathways.csv')
#	print(top_level_pathways_count)
#	print(top_level_pathways_to_genes)


main()
