def getLocusDict(filePath):
	fileContent = open(filePath).readlines()
	count = 0
	genomes = fileContent[0].split()[3:]
	selectedLoci = {}
	for line in fileContent[1:]:
		locusInfo = line.split()
		if locusInfo[0] != 'chrX':
			locusCoords = (locusInfo[0], locusInfo[1], locusInfo[2]) 
			validPairs = [elem for elem in zip(locusInfo[3:], genomes) if elem[0] != '-']
			diffFromRefs = [int(x) for elem in validPairs for x in elem[0].split(',')]
			if len([x for x in diffFromRefs if x!=0]) >= 2 and len(diffFromRefs) >= 32:
				selectedLoci[locusCoords] = validPairs
				count += 1
	print(count)
	return selectedLoci
	
def getLocusRefLength(filePath):
	fileContent = open(filePath).readlines()[1:]
	locusRefLength = {}
	for line in fileContent:
		elements = line.split()
		chrom = elements[0]
		start = elements[1]
		end = elements[2]
		refLength = int(elements[3])
		locusRefLength[(chrom, start, end)] = refLength
	
	return locusRefLength	

def main():
	filePath = "STR_Genotypes_Genome_Wide.txt"
	selectedLoci = getLocusDict(filePath)

	sumatraGenomes = ['PA_A947', 'PA_A948', 'PA_A950', 'PA_A952', 'PA_A949','PA_B018', 'PA_B017','PA_B020', 'PA_A953','PA_A955','PA_A964','PA_B019']
	borneoGenomes = ['PP_KB5404','PP_A940','PP_A941','PP_A943','PP_A944','PP_A938', 'PP_A983', 'PP_A939', 'PP_A942', 'PP_A946', 
		'PP_A987', 'PP_A988', 'PP_5062', 'PP_A989', 'PP_A984', 'PP_A985']
	allGenomes = sumatraGenomes + borneoGenomes
	
	locusRefLength = getLocusRefLength("LocusRefLengthGenomeWide.txt")
		
	out = open("StructureSpecies.txt",'w')

#	islandsMap = {}
#	for genome in sumatraGenomes:
#		islandsMap[genome] = "1"
#	for genome in borneoGenomes:
#		islandsMap[genome] = "2"

	
	pongoAbelii = ['PA_A947', 'PA_A948', 'PA_A950', 'PA_A952', 'PA_A949','PA_B018', 'PA_B017','PA_B020', 'PA_A953','PA_A955','PA_A964','PA_B019']
	pongoPygmaeusWurmbii = ['PP_KB5404','PP_A940','PP_A941','PP_A943','PP_A944','PP_A938', 'PP_A983']
	pongoPygmaeusPygmaeus = ['PP_A939', 'PP_A942', 'PP_A946']
	pongoPygmaeusMorio = ['PP_A987', 'PP_A988', 'PP_5062', 'PP_A989', 'PP_A984', 'PP_A985']
	
	allGenomes = pongoAbelii + pongoPygmaeusWurmbii + pongoPygmaeusPygmaeus + pongoPygmaeusMorio

	speciesMap = {}
	for genome in pongoAbelii:
		speciesMap[genome] = "1"
	for genome in pongoPygmaeusWurmbii:
		speciesMap[genome] = "2"
	for genome in pongoPygmaeusPygmaeus:
		speciesMap[genome] = "3"
	for genome in pongoPygmaeusMorio:
		speciesMap[genome] = "4"

	out.write("	")
	for i in range(1,len(selectedLoci.keys()) + 1):
		out.write("	Locus"+str(i))
	out.write("\n")

	i = 1
	for genome in allGenomes:
		firstRow = str(i) + "	" + speciesMap[genome]
		secondRow = str(i)  + "	" + speciesMap[genome]
		for locusCoord in selectedLoci.keys():		
			locusInfo = selectedLoci[locusCoord]
			refLength = locusRefLength[locusCoord]
			locusGenomes = [elem[1] for elem in locusInfo]
			if genome in locusGenomes:
				for genotype in locusInfo:
					alleles = genotype[0]
					a1, a2 = alleles.split(',')
					na1 = int(a1) + refLength
					na2 = int(a2) + refLength
					genomeName = genotype[1]
					if genome == genomeName:
						firstRow += "	" + str(na1)
						secondRow += "	" + str(na2)
			else:
				firstRow += "	-9"
				secondRow += "	-9"
		firstRow += "\n"
		secondRow += "\n"
		out.write(firstRow)
		out.write(secondRow)
		i += 1
	out.close()
	
main()
