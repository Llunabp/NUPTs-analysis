import sys
import argparse

from subprocess import run


import pyfastx
import statistics


def parse_arguments():
    #Create the arguments and the help menssage to introduce data to the program
    desc = "Found the not aligned parts of reads and if exist multiple reads"
    parser = argparse.ArgumentParser(description=desc)

    help_aligned_reads_file = "Aligned reads to analyze"
    parser.add_argument("--aligned" ,
                        "-a", type=str,
                        help=help_aligned_reads_file,
                        required=True)
    help_acceptation_value = "Minimun read query coverage"
    parser.add_argument("--acceptation" ,
                        "-v", type=float,
                        help=help_acceptation_value,
                        required=False,
                        default = 1)
    help_reference_guided = "Use a nuclear reference genome to detect NU(M/P)Ts"
    parser.add_argument("--reference", 
                        "-r",
                        type=str, help=help_reference_guided,
                        required=False,
                        default="")    
    help_reads_sequences_file = "Sequences of reads to analyze"
    parser.add_argument("--sequences" ,
                        "-s", type=str,
                        help=help_reads_sequences_file,
                        required=True)
    help_seq_aligned_file = "File with sequences of reads without propper query"
    parser.add_argument("--file" ,
                        "-f", type=str,
                        help=help_seq_aligned_file,
                        required=True)
    help_blast_results= "File with blast results"
    parser.add_argument("--blast",
                        "-b", type=str,
                        help=help_blast_results,
                        required=True)
    help_e_value_filter= "E value that should pass the blast results"
    parser.add_argument("--e_value",
                        "-e", type=float,
                        help=help_e_value_filter,
                        default=1e-50)
    help_number_threads= "Number of threads"
    parser.add_argument("--threads",
                        "-t", type=int,
                        help=help_number_threads,
                        default = 3)
    help_length_filter= "Length aligned that should pass the blast results"
    parser.add_argument("--length",
                        "-l", type=int,
                        help=help_length_filter,
                        default=1000)
    help_median_filter= "Limit to agroup in order to do median"
    parser.add_argument("--median",
                        "-m", type=int,
                        help=help_median_filter,
                        default=15)
    help_cloro_limit= "Limit to group reads of same cloroplast"
    parser.add_argument("--cloro",
                        "-c", type=int,
                        help=help_cloro_limit,
                        default=100)
    help_min_reads= "Minimun reads to the median"
    parser.add_argument("--min",
                        "-i", type=int,
                        help=help_min_reads,
                        default=10)
   
   
    return parser


def get_options():
    #Use arguments in order to introduce the necessary data to the program
    parser = parse_arguments()
    options = parser.parse_args()
    alignments_fhand = Path(options.aligned)
    acceptation_value = options.acceptation
    sequences_file = options.sequences
    seq_aligned_file = Path(options.file)
    blast_results = Path(options.blast)
    e_value = options.e_value
    number_threads = options.threads
    length_filter = options.length
    median_filter = options.median
    cloro_limit = options.cloro
    min_reads = options.min
    if options.reference:
        reference_genome = Path(options.reference)
    else:
        reference_genome = False

    return {"alignments_fhand" : alignments_fhand,
            "acceptation_value": acceptation_value,
            "reference_genome": reference_genome,
            "sequences_file" : sequences_file,
            "seq_aligned_file" : seq_aligned_file,
            "blast_results" : blast_results,
            "e_value" : e_value,
            "number_threads" : number_threads,
            "length_filter" : length_filter,
            "median_filter" : median_filter,
            "cloro_limit" : cloro_limit,
            "min_reads" : min_reads
            }


def get_reads_alignments_info(reads_fhand):
#Get info of diferent reads from an archive in paf format
    reads_alignments_info = {}
    for line in reads_fhand:
        if line:
            line = line.rstrip()
            line = line.split()
            read_name = line[0]
            read_length = int(line[1])
            relative_strand = line[4]
            target_seq_name = line[5]
            target_positions = int(line[7]), int(line[8])
            aligned_read_positions = (int(line[2]),int(line[3]), relative_strand)
            total_alignment =  int(line[3]) - int(line[2])
            if read_name in reads_alignments_info:
                if read_length != reads_alignments_info[read_name][0]['length']:
                    msg = "Read have differents lengths {}".format(read_name)
                    raise RuntimeError(msg)
                reads_alignments_info[read_name].append({'length' : read_length, 
                                                         'positions' : aligned_read_positions, 
                                                         'target_seq_name' : target_seq_name, 
                                                         'target_positions' : target_positions,
                                                         'total_aligned_length': total_alignment})
            
            else:
                reads_alignments_info[read_name] = [{'length' : read_length, 
                                                    'positions' : aligned_read_positions, 
                                                    'target_seq_name' : target_seq_name, 
                                                    'target_positions' : target_positions,
                                                    'total_aligned_length': total_alignment}]
    return reads_alignments_info

def get_reads_with_multiple_matches(reads_alignments_info):
#Agroup reads that are fragmented
    return [read for read, read_data in reads_alignments_info.items() if len(read_data) > 1]


def calculate_reads_query_coverage(alignments_info):
#Calculate reads query coverage (aligned part respect total)
    coverages = {}
    for read, alignments in alignments_info.items():
        read_length = alignments[0]["length"]
        total_nucleotides = 0
        for alignment in alignments:
            total_nucleotides += alignment["total_aligned_length"]
        coverage = float(total_nucleotides / read_length)
        coverages[read] = coverage
    return coverages

def get_reads_under_coverage(reads_coverage, acceptation_value = 0.9):
#Obtain reads that their query coverage is under a determinate acceptation value
    return [read_name for read_name, coverage in reads_coverage.items() if coverage <= acceptation_value]


def filter_reads_by_name(alignments_info, readnames):
#Obtain the rest of the information of the reads that not presente propper query coverage
    return {readname: alignments_info[readname] for readname in readnames}


def get_the_sequence_of_filtered_reads(reads_sequences,filtered_reads,not_aligned_file):
    with open(not_aligned_file, 'w') as out_fhand:     
        for read in filtered_reads: 
            if read in reads_sequences:
                sequence = "{} \n".format(reads_sequences[read].seq)
                header = ">{} \n".format(read)
                out_fhand.write(header)
                out_fhand.write(sequence)
                out_fhand.flush()
    return not_aligned_file


def get_alignments_with_blast (seq_aligned_file,reference_genome,num_threads, results_blast):
    cmd = ['blastn','-query', str(seq_aligned_file), '-db', str(reference_genome),'-num_threads', str(num_threads), '-outfmt', "\"6 std slen qlen\" > {}".format(results_blast)]
    results_blast = run(" ".join(cmd), shell=True,
                capture_output=True)
    return results_blast


def pass_evalue_filter_to_blast(results_blast,e_value_filter,acceptable_blast_fpath):
    with open(results_blast, 'r') as blast_fhand:
        with open(acceptable_blast_fpath, 'w') as acceptable_blast_fhand:
            for line in blast_fhand:
                line = line.rstrip()
                line = line.split()
                evalue = float(line[10])
                if evalue  <= e_value_filter or evalue == 0.0:
                    acceptable_blast_fhand.write("{} \n".format(" ".join(line)))


def pass_length_filter_to_blast(length_filter,acceptable_blast_fpath,final_blast_fpath):
    with open(acceptable_blast_fpath, 'r') as blast_fhand:
        with open(final_blast_fpath, 'w') as final_blast_fhand:
            for line in blast_fhand:
                line = line.rstrip()
                line = line.split()
                total_length = int(line[13])
                length = int(line[3])
                region = line[1]
                if ((length + length_filter) >= total_length) or region == 'NC_000932.1':
                    final_blast_fhand.write("{} \n".format(" ".join(line)))


def get_blast_info(blast_fpath, include=[], ignore=[], organelle_hits=False):
    blast_info = {}
    with open(blast_fpath, "r") as blast_fhand:
        for line in blast_fhand:
           if line:
                line = line.rstrip()
                line = line.split()
                read_name = line[0]
                subject_name = line[1]
                subject_start = int(line[8])
                subject_end = int(line[9])
                query_start = int(line[6])
                query_end = int(line[7])
                query_length = int(line[-1])
                if subject_end < subject_start:
                    subject_start, subject_end = subject_end, subject_start
                    strand = "-"
                else:
                    strand = "+"
                if subject_name not in include and include:
                    continue
                if subject_name in ignore and ignore:
                    continue
                check_name = "{}_0".format(read_name)
                if check_name not in blast_info:
                    blast_info[check_name] = {"subject_start": subject_start, 
                                             "subject_end": subject_end,
                                             "query_start": query_start,
                                             "query_end": query_end,
                                             "strand": strand,
                                             "subject_name": subject_name,
                                             "length": query_length,
                                             "read_name": read_name}
                elif organelle_hits:
                    if (blast_info[check_name]["query_end"] - query_start) < 100 and (blast_info[check_name]["query_end"] - query_start) > 0:
                        blast_info[check_name]["query_end"] = query_end
                    elif (query_end - blast_info[check_name]["query_start"]) < 100 and (query_end - blast_info[check_name]["query_start"]) > 0:
                        blast_info[check_name]["query_start"] = query_start
                    elif query_end < blast_info[check_name]["query_start"] or query_start > blast_info[check_name]["query_end"]:
                        blast_info["{}_1".format(read_name)] = {"subject_start": subject_start, 
                                                                "subject_end": subject_end,
                                                                "query_start": query_start,
                                                                "query_end": query_end,
                                                                "strand": strand,
                                                                "subject_name": subject_name,
                                                                "length": query_length,
                                                                "read_name": read_name}

                             
    return blast_info

def get_possible_NUPTs_reads(blast_info):
    possible_NUPTs = {}
    for read, features in blast_info.items():
        for i in range(0,len(features)-1):
            region = features[i]['region']
            positions = features[i]['positions']
            genome_positions = features[i]['genome_positions']
        if region == 'NC_000932.1' and (read not in possible_NUPTs):
            start_genome = genome_positions[0]
            end_genome = genome_positions[1]
            direction = 'forward'
            if start_genome > end_genome:
                direction = 'reverse'
                genome_positions = (end_genome,start_genome)
            possible_NUPTs[read] = {'positions' : positions, 'genome_positions' : genome_positions, 'cloro_direct' : direction}
    return possible_NUPTs

def get_positions_of_NUPTs_in_genome(possible_NUPTs, blast_info):
    NUPTs_reads_positions={}
    cont = 0
    for read, features in blast_info.items():
        chromosome_positions = {}
        if read in possible_NUPTs:
            cont += 1
            for i in range(0,len(features)-1):
                region = features[i]['region']
                start_position_NUPT = int(possible_NUPTs[read]['positions'][0])
                end_position_NUPT = int(possible_NUPTs[read]['positions'][1])

                start_position_genome = int(features[i]['genome_positions'][0])
                end_position_genome = int(features[i]['genome_positions'][1])
                direction = 'forward'
                
                if start_position_genome > end_position_genome:
                    start_position_genome, end_position_genome = end_position_genome, start_position_genome
                    start_position_NUPT, end_position_NUPT = end_position_NUPT, start_position_NUPT
                    direction = 'reverse'
                    
                NUPT_start = start_position_genome + start_position_NUPT
                NUPT_end = int(end_position_genome) - int(end_position_NUPT)
                NUPT_position = (NUPT_start,NUPT_end)
                if NUPT_start > NUPT_end:
                    NUPT_position = (NUPT_end,NUPT_start)
                if region not in chromosome_positions and region != 'NC_000932.1' and region != 'NC_037304.1': 
                    chromosome_positions[region] = [NUPT_position]
                elif region != 'NC_000932.1' and region != 'NC_037304.1':
                    chromosome_positions[region].append (NUPT_position)
            NUPTs_reads_positions[read]={'cloroplast' : possible_NUPTs[read]['genome_positions'], 'cloro_direct' : possible_NUPTs[read]['cloro_direct'], 'region': region, 'chromosome' : chromosome_positions, 'direction' : direction}
    return NUPTs_reads_positions

def get_file_with_results(NUPTs_reads_positions,NUPTs_file):
    with open(NUPTs_file, 'w') as NUPTs_fhand: 
        first_line = ' Read Cloro_Start Cloro_End Chromosome Genom_Start Genom_End Direction\n'
        NUPTs_fhand.write(first_line)
        for read,positions in NUPTs_reads_positions.items():
            Read = read
            Cloro_Start = positions['cloroplast'][0]
            Cloro_End = positions['cloroplast'][1]
            cloro_direct = positions['cloro_direct']
            Direction = positions['direction']
            chromosomes = positions['chromosome']
            for Chromosome, genom_positions in chromosomes.items():
                for pair in genom_positions:
                    Genom_Start = pair[0]
                    Genom_End= pair[1]
                    line = '{} {} {} {} {} {} {} {} \n'.format(Read,Cloro_Start,Cloro_End,cloro_direct,Chromosome,Genom_Start,Genom_End,Direction)
                    NUPTs_fhand.write(line)
                    NUPTs_fhand.flush()
    return NUPTs_file

def get_reads_that_aligne_same_chloroplast(cloro_limit,NUPTs_reads_positions):
    reads_grouped_by_cloro = {}
    read_list = []
    for Read in NUPTs_reads_positions:
            Cloro_Start = int(NUPTs_reads_positions[Read]['cloroplast'][0])
            Cloro_End = int(NUPTs_reads_positions[Read]['cloroplast'][1])
            range_values = (Cloro_Start,Cloro_End)
            if Read not in read_list:
                for read, positions in NUPTs_reads_positions.items():
                    if read not in read_list:
                        cloro_start = int(positions['cloroplast'][0])
                        cloro_end = int(positions['cloroplast'][1])
                        if cloro_start >= (Cloro_Start - cloro_limit) and cloro_end <= (Cloro_End + cloro_limit):
                            read_list.append(read)
                            genome_positions = NUPTs_reads_positions[read]['chromosome']
                            list_genom_positions = []
                            count = 0
                            for chromosome,position in genome_positions.items():
                                count +=1
                                for pair in position:
                                    list_genom_positions.append(pair)
                            for i in range(0,len(list_genom_positions)-1):
                                genom_positions = list_genom_positions[i]
                                if range_values not in reads_grouped_by_cloro:
                                    reads_grouped_by_cloro[range_values] = [{'read_name':read, 'cloro_positions':positions['cloroplast'], 'genom_positions':genom_positions}]
                                else:
                                    reads_grouped_by_cloro[range_values].append({'read_name':read, 'cloro_positions':positions['cloroplast'], 'genom_positions':genom_positions})
    return reads_grouped_by_cloro

def get_reads_groups_by_genome_proximity(median_filter,reads_grouped_by_cloro, min_reads):
    NUPTs_positions = {}
    cont = 1
    read_list =[]
    for range_values, reads in reads_grouped_by_cloro.items():
        i = len(reads) -1
        for pat in range(0,i):
            Start_genom = reads[pat]['genom_positions'][0]
            End_genom = reads[pat]['genom_positions'][1]
            group_genome_start = []
            group_genome_end = []
            group_reads = []
            read_list =[]
            for a in range(0,len(reads)-1):
               
                read = reads[a]['read_name']
                if read not in read_list:
                    start_genom = reads[a]['genom_positions'][0]
                    end_genom = reads[a]['genom_positions'][1]
                    if start_genom >= (Start_genom - median_filter) and end_genom <= (End_genom + median_filter):
                        read_list.append(read)
                        group_genome_start.append(start_genom)
                        group_genome_end.append(end_genom)
                        group_reads.append(read)
            if len(group_genome_start) >= min_reads: 
                
                name = 'NUPT_{}'.format(cont)
                NUPTs_positions[name]= {'reads': group_reads, 'genom_start': group_genome_start,'genom_end': group_genome_end, 'cloroplast': range_values}
                cont +=1
    return NUPTs_positions 


def write_final_NUPTs_in_file(final_file,NUPTs_final):
    with open(final_file, 'w') as final_fhand:
        first_line = 'NUPT  cloro_start  cloro_end  region  genome_start  genome_end \n'
        final_fhand.write(first_line)
        for NUPT,features in NUPTs_final.items():
            cloro_start = features['cloro_start']
            cloro_end = features['cloro_end']
            region = features['region']
            genome_start = features['genome_start']
            genome_end = features['genome_end']
            line_NUPT = '{}  {}  {}  {}  {}  {}  \n'.format(NUPT,cloro_start,cloro_end,region,genome_start,genome_end)
            final_fhand.write(line_NUPT)
            final_fhand.flush()
    return final_file

def get_nupts_positions(nuclear_hits, sorted_organelle_hits):
        #FEATURES -> COSAS DEL CLOROPLASTOS 
        #NUCLEAR HITS -> COSAS DEL NUCLE
        NUPT_reads = {}
        for readname, features in sorted_organelle_hits.items():
            check_readname = "{}_0".format(readname[:-2])
            if check_readname in nuclear_hits:
                nuclear_hit = nuclear_hits[check_readname]
                nupt_start = 0
                nupt_end = 0
                cloro_start = features["subject_start"]
                cloro_end = features["subject_end"]
                if features["strand"] == "-" and nuclear_hit["strand"] == "-":
                    nupt_start = nuclear_hit["subject_start"] + features["length"] - features["query_end"]
                    nupt_end = nupt_start + (features["query_end"] - features["query_start"])
                    
                    #print(readname, features["strand"], nuclear_hit["strand"], nupt_start, nupt_end)
                if features["strand"] == "+" and nuclear_hit["strand"] == "+":
                    nupt_start = nuclear_hit["subject_start"] + (features["query_start"])
                    nupt_end = nupt_start + ((features["query_end"] - features["query_start"]))
                    
                    #print(readname, features["strand"], nuclear_hit["strand"], nupt_start, nupt_end)
                if features["strand"] == "+" and nuclear_hit["strand"] == "-":
                     nupt_start = nuclear_hit["subject_start"] + (features["length"] - features["query_end"])
                     nupt_end = nupt_start + (features["query_end"] - features["query_start"]) 
                     
                if features["strand"] == "-" and nuclear_hit["strand"] == "+":
                     nupt_end = nuclear_hit["subject_start"] + (features["query_end"])
                     nupt_start = nupt_end - (features["query_end"] - features["query_start"])
                     
                name_to_insert = readname
                if name_to_insert in NUPT_reads:
                    name_to_insert = "{}_1".format(readname[:-2])
                NUPT_reads[name_to_insert]= {'read' : readname, 
                             'nupt_start' : nupt_start, 
                             'nupt_end' : nupt_end, 
                             'cloro_start' : cloro_start, 
                             'cloro_end': cloro_end,
                             "chrom": nuclear_hit["subject_name"]}
        return NUPT_reads


def group_reads_of_same_NUPT(NUPT_reads,median_filter,cloro_limit, min_reads):
    p_cloro_start = 0
    p_cloro_end = 0
    p_nupt_start = 0
    p_nupt_end = 0
    Group_positions ={}
    NUPT_start_positions = []
    NUPT_end_positions = []
    readnames = []
    count = 0
    chrom = ""
    for readname, values in NUPT_reads.items():
        if not p_cloro_start:
            NUPT_start_positions.append(values["nupt_start"])
            NUPT_end_positions.append(values["nupt_end"])
            p_cloro_start = values["cloro_start"]
            p_cloro_end = values["cloro_end"]
            chrom = values["chrom"]
            readnames = [readname]
        else:
            if abs(NUPT_start_positions[-1] - values["nupt_start"]) <= median_filter and abs(NUPT_end_positions[-1] - values["nupt_end"]) <= median_filter:
                NUPT_start_positions.append(values["nupt_start"])
                NUPT_end_positions.append(values["nupt_end"])
                p_cloro_start = values["cloro_start"]
                p_cloro_end = values["cloro_end"]
                chrom = values["chrom"]
                readnames.append(readname)
            else:
                name = "Group_{}".format(count)
                Group_positions[name] = {"nupt_starts": NUPT_start_positions.copy(),
                                            "nupt_ends": NUPT_end_positions.copy(),
                                            "organelle_start": p_cloro_start,
                                            "organelle_end": p_cloro_end,
                                            "nuclear": chrom,
                                            "reads": readnames.copy()}
                count += 1
                
                chrom = values["chrom"]
                NUPT_start_positions = [values["nupt_start"]]
                NUPT_end_positions = [values["nupt_end"]]
                p_cloro_start = values["cloro_start"]
                p_cloro_end = values["cloro_end"]
                readnames = [readname]
               
    return Group_positions

def order_group_by_start_position(Group_positions):
    start_positions = {}
    for group, values in Group_positions.items():
        start_position = values["nupt_starts"][0]
        if start_position not in start_positions:
            start_positions[start_position] = group
        else:
            start_position += 1
            start_positions[start_position] = group

        ordered_starts = start_positions.keys()
        ordered_starts = sorted(ordered_starts)
    ordered_groups= []
    for starts in ordered_starts:
        ordered_groups.append(start_positions[starts])
    Group_positions_ordered = {}
    for Group in ordered_groups:
        info = Group_positions[Group]
        Group_positions_ordered[Group] = info
    return Group_positions_ordered
        
        
def join_true_NUPTs (Group_positions_ordered, median_filter,min_reads):
    p_start_organ = 0
    count = 0
    NUPTs_positions = {}
    joined_groups_starts =[]
    joined_groups_end =[]
    list_reads = []
    for group,value in Group_positions_ordered.items():
        starts_positions = value['nupt_starts']
        ends_positions = value['nupt_ends']
        start_organ = value['organelle_start']
        end_organ = value['organelle_end']
        chrom = value['nuclear']
        reads = value['reads']
        
        if not p_start_organ:
            for start in starts_positions:
                joined_groups_starts.append(start)
            for end in ends_positions:
                joined_groups_end.append(end)
            for read in reads:
                list_reads.append(read)

            p_start_organ = value['organelle_start']
            p_end_organ = value['organelle_end']
            p_chrom = value['nuclear']
            
        elif abs(joined_groups_starts[-1]-starts_positions[0]) <=  median_filter and abs(joined_groups_end[-1]-ends_positions[0]) <=  median_filter:     
            if p_start_organ == start_organ and p_end_organ == end_organ and chrom == p_chrom:
                for start in starts_positions:
                    joined_groups_starts.append(start)
                for end in ends_positions:
                    joined_groups_end.append(end)
                for read in reads:
                    list_reads.append(read)
            p_start_organ = value['organelle_start']
            p_end_organ = value['organelle_end']
            p_chrom = value['nuclear'] 
                   
        elif len(joined_groups_starts) >= min_reads:
            name = 'NUPT_{}'.format(count)
            NUPTs_positions[name] = {'nupt_starts' : joined_groups_starts, 
                                     'nupt_ends' : joined_groups_end,
                                     'organelle_start' : p_start_organ, 
                                     'organelle_end' : p_end_organ,
                                     'nuclear' : chrom,
                                     'reads' : list_reads}
            count +=1
            joined_groups_starts =[]
            joined_groups_end =[]
            for start in starts_positions:
                joined_groups_starts.append(start)
            for end in ends_positions:
                joined_groups_end.append(end)
            p_start_organ = value['organelle_start']
            p_end_organ = value['organelle_end']
            p_chrom = value['nuclear'] 
            list_reads = []
            for read in reads:
                list_reads.append(read)

        else:  
            joined_groups_starts =[]
            joined_groups_end =[]
            for start in starts_positions:
                joined_groups_starts.append(start)
            for end in ends_positions:
                joined_groups_end.append(end)
            p_start_organ = value['organelle_start']
            p_end_organ = value['organelle_end']
            p_chrom = value['nuclear']
            list_reads = []
            for read in reads:
                list_reads.append(read)
          

        
    return NUPTs_positions

        

def obtain_final_info_NUPTs (NUPT_positions):
    NUPT_info = {}
    for NUPT,values in NUPT_positions.items():
        starts_position =values['nupt_starts']
        ends_position =values['nupt_ends']
        reads = values['reads']
        start_position = int(statistics.median(starts_position))
        end_position = int(statistics.median(ends_position))
        positions_nuclear = (start_position,end_position)
        positions_organelle = (values['organelle_start'], values['organelle_end'])
        NUPT_info[NUPT] = {'nuclear' : values['nuclear'], 'positions_nuclear' : positions_nuclear, 'positions_organelle' : positions_organelle, 'reads': reads}
    return NUPT_info

def write_results_table(NUPTs_info, results_fhand):
    line = ["#NUPT_ID", "CHROM_ID","NUCLEAR_START", "NUCLEAR_END", "ORGANELLE_START", "ORGANELLE_END", "READS_ID", "\n"]
    results_fhand.write("\t".join(line))
    results_fhand.flush()
    for nupt, features in NUPTs_info.items():
        nupt_id = nupt
        chrom_id = features["nuclear"]
        nuclear_start = str(features["positions_nuclear"][0])
        nuclear_end = str(features["positions_nuclear"][1])
        organelle_start = str(features["positions_organelle"][0])
        organelle_end = str(features["positions_organelle"][1])
        reads = ",".join(features["reads"])
        line = [nupt_id, chrom_id, nuclear_start, nuclear_end, organelle_start, organelle_end, reads, "\n"]
        results_fhand.write("\t".join(line))
        results_fhand.flush()

def main():
    options = get_options()
    alignments_fpath = options['alignments_fhand']
    acceptation_value = options["acceptation_value"]
    reads_sequences = pyfastx.Fastq(options['sequences_file'])
    seq_aligned_file = options["seq_aligned_file"]
    results_blast = options['blast_results']
    e_value_filter = options['e_value']
    length_filter = options['length_filter']
    acceptable_blast_fpath = "filtered_blast_evalue_{}.tab".format(e_value_filter)
    final_blast_fpath = "filtered_blast_evalue_{}_{}.tab".format(e_value_filter, length_filter)
    reference_genome = options['reference_genome']
    num_threads = options['number_threads']
    min_reads = options['min_reads']
    median_filter = options['median_filter']
    cloro_limit = options['cloro_limit']
    with open(alignments_fpath) as alignments_fhand:
        alignments_info = get_reads_alignments_info(alignments_fhand)

    #get_reads_with_multiple_matches(alignments_info)
    coverages = calculate_reads_query_coverage(alignments_info) #ok
    reads_under_coverage = get_reads_under_coverage(coverages, acceptation_value=acceptation_value) #ok
    filtered_reads = filter_reads_by_name(alignments_info, reads_under_coverage) #ok
    get_the_sequence_of_filtered_reads(reads_sequences,filtered_reads,seq_aligned_file) #ok
    get_alignments_with_blast (seq_aligned_file,reference_genome,num_threads, results_blast) #ok
    pass_evalue_filter_to_blast(results_blast,e_value_filter,acceptable_blast_fpath) #ok
    pass_length_filter_to_blast(length_filter,acceptable_blast_fpath,final_blast_fpath) #ok
    organelle_hits = get_blast_info(final_blast_fpath, include="NC_000932.1", organelle_hits=True)
    sorted_organelle_hits = dict(sorted(organelle_hits.items(), key=lambda k_v: k_v[1]["subject_start"]))
    nuclear_hits = get_blast_info(final_blast_fpath, ignore=["NC_037304.1", "NC_000932.1"])
    nupts_reads = get_nupts_positions(nuclear_hits, sorted_organelle_hits)
    Group_positions = group_reads_of_same_NUPT(nupts_reads,median_filter,cloro_limit, min_reads)
    Group_positions_ordered = order_group_by_start_position(Group_positions)
    NUPTs_positions = join_true_NUPTs (Group_positions_ordered, median_filter,min_reads)
    NUPTs_info = obtain_final_info_NUPTs(NUPTs_positions)
    print(NUPTs_info)
    with open("NUPTs_results.table.test.tsv", "w") as results_fhand:
        write_results_table(NUPTs_info, results_fhand)


if __name__ == '__main__':
    main()