#!/usr/bin/python
# Author: Gabriel dos Santos Goncalves
# Rio de Janeiro, Brasil - march 2016


# QuasiPycies

import re 
import math

class Clusters(object):

    def __init__(self, fasta):
        
        #Load the file and read lines
        self.fasta_file = open(fasta, "r").readlines()
        
        # Define the instances
        self.dict_fasta = {}
        self.names = []
        self.sequences = []
        self.percentages = []
        self.number_reads = []

        # Iterate over the lines of the file
        for number, line in enumerate(self.fasta_file):
            if line.startswith(">"):
                
                # Search for the number of reads in sequence name 
                n_search = re.search(";size=(\d+)", line)
                n = int(n_search.group(1))
                # Store the number of reads in a list
                self.number_reads.append(n)
                
                # Search for the percentage in sequence name
                percentage_search = re.search("(\d+.\d+)%", line)
                percentage = float(percentage_search.group(1))
                # Store the percentage in a list
                self.percentages.append(percentage)
                    
                # Define the name of the sequence without ">", number of reads and percentage
                name = line[1:n_search.start()]
                # Store the name in a list
                self.names.append(name)
                
                # Store the sequence in a list
                sequence = self.fasta_file[number+1]
                self.sequences.append(sequence)
                
                self.dict_fasta [name] = [sequence, percentage, n]
                
      
    
    def shannon_entropy(self, minimun_perc):
    	''' Shannon Entropy is a diversity index that measures the evenness of
    	the population. It uses the values of size of haplotypes 
    	(subpopulations) in the sample.  
    	'''
        
        number_reads = 0
        list_sizes = []

        for key, value in self.dict_fasta.items():
            
            if value[1] >=  minimun_perc:
                list_sizes.append(value[2])
                number_reads += value[2]
        
        #print(list_sizes)
        
        entropy = 0
        for v in list_sizes:
            subpop = (v/number_reads) * math.log(v/number_reads)
            entropy += subpop

        S = -(entropy/math.log(number_reads))
        #print(input1, "Shannon entropy :", round(S,3))
        return(str(round(S,3)))

    # 
    def simpsons_index(self):
    	''' Simpson's index measures the probability that two individuals 
    	(in this case reads) randomly selected from a sample will belong 
    	to the same species (in this case cluster).
    	'''
    	
        D = 0
	    for element in self.number_reads:
	        d= (element/sum(self.number_reads))**2
	        D += d
	        
	    return(D)

	def simpsons_index_of_diversity(self):
		''' Simpson's index of diversity measures the probability that two
		individuals (in this case reads) randomly selected from a sample will
		belong to different species (in this case cluster). 

		Simpson's index of diversity = 1 - Simpson's index
		'''

	    D = 0
	    for element in self.number_reads:
	        d= (element/sum(self.number_reads))**2
	        D += d
	        
	    return(1-D)


	def simpsons_Reciprocal_Index(self):
    ''' Simpson's Reciprocal Index 

    '''

	    D = 0
	    for element in self.number_reads:
	        d= (element/sum(self.number_reads))**2
	        D += d
	        
	    return(1/D)


    # This function writes to a fasta file the information of self.dict_fasta
    def clusters_object_to_fasta(self, fasta_name):
    	''' This function creates a fasta file from the dict_fasta object
    	'''
        
        # First we create a file to receive the data
        fasta_output = open(fasta_name, "w")
        
        # Then we iterate through the self.dict_fasta dictionary 
        for key, value in self.dict_fasta.items():
            fasta_output.write(">" + key + "_size=" + str(value[1]) + "_" + str(value[1]) + "%" + "\n")
            fasta_output.write(value[0])

    
    def filtering_clusters_by_percentage(self, minimun_perc, output1):
    	''' This function filters the clusters by percentage.
    	'''

        #file_in = open (input1, "r")
        #fasta_in = file_in.readlines()
        file_out = open (output1, "w")
        number_reads = 0
        list_sizes = []

        for key, value in self.dict_fasta.items():
            n = int(value[2])
            number_reads += n
        
        print(number_reads)

        for key, value in self.dict_fasta.items():
            #print(number_reads)
            #print(value[1])
            if value[2] >= (number_reads * float(minimun_perc)):
                print((number_reads * float(minimun_perc)))
                print(value[2])
                #print(">" + name + "_" + value[2]+ "_" + str(round(float((value[2]/number_reads)*100),2)) + "%" + "\n")
                
                #file_out.write(">" + name + "_" + value[2]+ "_" + str(round(float((value[2]/number_reads)*100),2)) + "%" + "\n")
                file_out.write(value[0])
                
    def number_haplotypes(self, value_type, cutoff):
    	''' This function counts the number of haplotypes greater than a 
    	cutoff. The cutoff can be "percentage" or "number_reads".

    	'''

        # Define the two variables to hold the number of haplotypes
        counter_haplotypes = 0
        counter_haplotypes_over_min_number = 0

        # If percentage is used to filter the haplotypes
        if value_type == "percentage":
            for key in self.dict_fasta.keys():
                counter_haplotypes += 1
                if self.dict_fasta [key][1] >= cutoff:
                    counter_haplotypes_over_min_number +=1
                    #print(self.dict_fasta [key][1])

        elif value_type == "number_reads":
            for key in self.dict_fasta.keys():
                counter_haplotypes += 1
                if self.dict_fasta [key][2] >= cutoff:
                    counter_haplotypes_over_min_number +=1
                    #print(self.dict_fasta [key][2])

        print ("From a total of " + str(counter_haplotypes) + " haplotypes only " + str(counter_haplotypes_over_min_number) + " are above the cutoff")
        return (counter_haplotypes_over_min_number)
    # This function filter the clusters base the percentage or number of reads cutoof
    # and returns the dict_fasta updated
    
    def filtering_Clusters(self, value_type, cutoff):

            self.filtered_dict_fasta = {}

            # Iterate over the keys of dict_fasta
            for key in self.dict_fasta.keys():

                if value_type == "percentage":
                    if self.dict_fasta[key][1] > cutoff:
                        print(self.dict_fasta[key][1])
                        self.filtered_dict_fasta[key] = self.dict_fasta[key]
                        #del self.dict_fasta[key]

                elif value_type == "number_reads":
                    if self.dict_fasta[key][2] > cutoff:
                        print(self.dict_fasta[key][2])
                        self.filtered_dict_fasta[key] = self.dict_fasta[key]
                        #del self.dict_fasta[key]


            return(self.filtered_dict_fasta)





'''		
# This class is in charge of dealing with pileup files, parsing it and calculating
# diversity indices like nucleotide diversity
class Pileup(object):

    def __init__(self, pileup_file, min_reads_per_position=0):
        
        # Create the dictionary to store each base data
        self.dict_pileup = {}

        #Load the file and read lines
        self.pileup_file = open(pileup_file, "r").readlines()
        
        # Iterate over the lines of the pileup file storing the information 
        # of each base position  
        for n,line in enumerate(self.pileup_file):
            line_elements =  line.split("\t")
            ref= line_elements[0]
            line_number = line_elements[1]
            ref_base = line_elements[2]
            pos_depth= line_elements[3]
            base_qualities= line_elements[5]
            pos_bases_original = line_elements[4]
            if n < min_reads_per_position:
                pos_bases_sub_1 = re.sub(r"\^.([ACTGNacgtn.,])", r"\1", pos_bases_original).upper()    
                pos_bases_sub_2 = re.sub("[.,]",ref_base, pos_bases_sub_1).upper()
                self.dict_pileup[line_number]= pos_bases_sub_2
            else:
                pos_bases_sub = re.sub("[.,]",ref_base, pos_bases_original).upper()
                self.dict_pileup[line_number]= pos_bases_sub


    # This function calculates the nucleotide diversity of the Pileup object
	def nucleotide_diversity(self):
		# First we create a list to hold the nucleotide diversity of each base position
        nuc_diversity_list=[]

        # Iterate over each value in the dict_pileup (representing the polymorphisms of each position)
        for pos in self.dict_pileup.values():
                
            # Count the number of each bases in the string
            A = pos.upper().count("A")
            T = pos.upper().count("T")
            C = pos.upper().count("C")
            G = pos.upper().count("G")
            reads = A + C + T + G
            Di = (float((A*C)+(A*G)+(A*T)+(C*G)+(C*T)+(G*T)))/(float((reads**2 - reads)/2))
                
            # Store the value nucleotide diversity in this position to the list
            nuc_diversity_list.append(Di)

        # Calculate average nucleotide diversity buy adding all the values in the list
        # and dividing by the number of base positions (number of itens in the dict_pileup)
        avg_nuc_diversity = sum(nuc_diversity_list)/len(self.dict_pileup.values())
        return(round(avg_nuc_diversity,4))

    # This function counts the number of polymorphic sites 
    def polymorphic_sites(self, value_type, cutoff):
    	
    	number_polymorphic_sites = 0

    	for key, value in dict_pileup.items():
    		A = value.upper().count("A")
        	T = value.upper().count("T")
        	C = value.upper().count("C")
        	G = value.upper().count("G")

        	list_n_bases = [A,T,C,G]
        	list_n_bases.sort

        	if list_n_bases[2] > cutoff:    
                    number_polymorphic_sites += 1

        return(number_polymorphic_sites)
'''

class Pileup(object):     '''This class is in charge of dealing with pileup
files, parsing it and calculating # diversity indices like nucleotide diversity
'''     def __init__(self, pileup_file, min_reads_per_position=0):
        
        # Create the dictionary to store each base data
        self.dict_pileup = {}

        #Load the file and read lines
        self.pileup_file = open(pileup_file, "r").readlines()
        
        # Iterate over the lines of the pileup file storing the information 
        # of each base position  
        for n,line in enumerate(self.pileup_file):
            line_elements =  line.split("\t")
            ref= line_elements[0]
            line_number = line_elements[1]
            ref_base = line_elements[2]
            pos_depth= line_elements[3]
            base_qualities= line_elements[5]
            pos_bases_original = line_elements[4]
            if n < min_reads_per_position:
                pos_bases_sub_1 = re.sub(r"\^.([ACTGNacgtn.,])", r"\1", pos_bases_original).upper()    
                pos_bases_sub_2 = re.sub("[.,]",ref_base, pos_bases_sub_1).upper()
                self.dict_pileup[line_number]= pos_bases_sub_2
            else:
                pos_bases_sub = re.sub("[.,]",ref_base, pos_bases_original).upper()
                self.dict_pileup[line_number]= pos_bases_sub
    
    def nucleotide_diversity(self):
            # dict_pileup = parsing_pileup(input_mpileup)
            nuc_diversity_list=[]

            for pos in self.dict_pileup.values():
                A = pos.upper().count("A")
                T = pos.upper().count("T")
                C = pos.upper().count("C")
                G = pos.upper().count("G")
                reads = A + C + T + G
                Di = (float((A*C)+(A*G)+(A*T)+(C*G)+(C*T)+(G*T)))/(float((reads**2 - reads)/2))
                #return(Di)
                #nuc_diversity_pos = nuc_diversity_site(pos)
                nuc_diversity_list.append(Di)

            avg_nuc_diversity = sum(nuc_diversity_list)/len(self.dict_pileup.values())
            #print(round(avg_nuc_diversity,4))
            return(round(avg_nuc_diversity,4))

    def polymorphic_sites(self, value_type, cutoff):

            number_polymorphic_sites = 0

            if value_type == "percentage":
            	for key, value in self.dict_pileup.items():
	                A = value.upper().count("A")
	                T = value.upper().count("T")
	                C = value.upper().count("C")
	                G = value.upper().count("G")

	                list_n_bases = [A,T,C,G]
	                list_n_bases.sort

	                #if list_n_bases[2] < sum(list_n_bases)*(cutoff/100):
	                if list_n_bases[2] > cutoff:    
	                    number_polymorphic_sites += 1

	        if value_type == "number_reads":
            	for key, value in self.dict_pileup.items():
	                A = value.upper().count("A")
	                T = value.upper().count("T")
	                C = value.upper().count("C")
	                G = value.upper().count("G")

	                list_n_bases = [A,T,C,G]
	                list_n_bases.sort

	                #if list_n_bases[2] < sum(list_n_bases)*(cutoff/100):
	                if list_n_bases[2] > cutoff:    
	                    number_polymorphic_sites += 1




            return(number_polymorphic_sites)







# These are python functions that may be used on the main workflow
def filtering_clusters_by_percentage(input1, minimun_perc, output1):

	file_in = open (input1, "r")
	fasta_in = file_in.readlines()
	file_out = open (output1, "w")
	number_reads = 0
	list_sizes = []

	for number, line in enumerate(fasta_in):
		if line[0]==">":
			n_search = re.search(";size=(\d+)", line)
			n = int(n_search.group(1))
			number_reads += n


	for number, line in enumerate(fasta_in):
		if line[0]==">":
			n_search = re.search(";size=(\d+)", line)
			n = int(n_search.group(1))
			if n >= (number_reads * float(minimun_perc)):
				file_out.write(line.rstrip() + str(round(float((n/number_reads)*100),2)) + "%" + "\n")
				file_out.write(fasta_in[number+1])


def SAM_to_fasta(SAM_file, fasta_file):

	SAM = open (SAM_file, "r")
	SAM_lines = SAM.readlines()
	fasta_file = open (fasta_file, "w")

	for line in SAM_lines:
		if line[0] != "@":
			line_elements = line.split("\t")
			fasta_file.write(">" + str(line_elements[0] + "\n"))
			fasta_file.write(line_elements[9]+ "\n")



# This function converts fastq files into fasta files
def fastq_to_fasta(fastq_file, fasta_file):
	startTime = datetime.now()

	# Open input fastq file and create output fasta file 
	fastq1 = open (fastq_file, "r")
	fastq = fastq1.readlines()
	fasta = open (fasta_file, "w")
	
	print("Number of lines in the file " + str(len(fastq)))

	# Create lists with the number of the lines carrying the name of reads and bases 
	list_name_lines = [] # 0
	list_bases_lines = [] # 1
	 
	for a in range(0,(len(fastq)+1),4):
		list_name_lines.append(a)
	
	for b in range(1,(len(fastq)+1),4):
		list_bases_lines.append(b)
	
	# Iterate through each line in the file and copy only the ones that represent the
	# name and the bases of each read to the new fasta file. Also replacing the "@" for 
	# the ">" before each read name 	
	for number, line in enumerate(fastq):
		if number in list_name_lines:
			mod1 = re.sub('^@','>', line)
			fasta.write(mod1)
			
		elif number in list_bases_lines:
			fasta.write(line)
	# Print the time it took to convert fastq to fasta
	print ("Time to convert fastq to fasta: ",datetime.now() - startTime)


# This function converts multiline fasta into single line
def MultiLine_fasta_to_SingleLine(input1, output1):

	fasta = open (input1, "r")
	fasta_ML = fasta.readlines()
	fasta_SL = open (output1, "w")

	for n,line in enumerate(fasta_ML):
		if line.startswith(">"):
			fasta_ML[n] =  "\n" + line.rstrip() 
		else:
			fasta_ML[n]=line.rstrip()

	data = ('|'.join(fasta_ML))

	data = data[1:]
	
	x = re.sub(";\|","\n",data)
	y = re.sub("\|","",x)
	fasta_SL.write(y+"\n")

# This function transform clusters back to single reads
def clusters_to_single_reads(input1, output1):

	fasta = open (input1, "r")
	fasta_ML = fasta.readlines()
	fasta = open (output1, "w")

	for n,line in enumerate(fasta_ML):
		if line.startswith(">"):
			count_match = re.search("size=(\d+)", line)
			count_number = int(count_match.group(1))
			read_name = line[:line.find(";size")]
			
			for i in range(0,count_number):
				fasta.write(read_name + "_" + str(i)+"\n")
				fasta.write(fasta_ML[n+1])

