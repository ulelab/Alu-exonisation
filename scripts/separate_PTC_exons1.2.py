'''
Created on Jul 7, 2015

@author: Nejc Haberman

Script will receive a fasta file with a header of regulated exon ID, regulated exon position relative to CDS start, all exonic positions relative to CDS start position 
and a FASTA sequence of CDS.

As a result it will find a premature termination codon (PTC) in regulated exons

(4) Will exons introduce a PTC (premature termination codon)?
Predict frame of the exon based on the dominant isoform of the gene and YES if UAG,UAA,UGA within this frame (it contains a PTC)
If it doesnt have PTC, then check if exon length not dividable by 3. If not, then write YES (it will cause a frame-shift and thereby likely to introduce PTC). If it is dividable by 3, then write NO.

'''


import sys

def get_exon_coding_position(reg_exon_pos, exon_positions):
    reg_exon_start = int(reg_exon_pos.rsplit(':')[0])
    reg_exon_end = int(reg_exon_pos.rsplit(':')[1])
    distance = 0    #intialized start position of regulated exon on CDS 
    for i in range(0,exon_positions.__len__()): 
        exon_start = int(exon_positions[i].rsplit(':')[0])  #we count with 0
        exon_end = int(exon_positions[i].rsplit(':')[1])
        exon_length = exon_end - exon_start
        
        if i == 0 and exon_start < 0:   #in case the first exon is overlapping CDS start
            exon_length = exon_end - 1  # we count with 1 now that's why -1
        if reg_exon_start >= exon_end:
            distance += exon_length
            distance += 1
        else:
            if exon_start <= reg_exon_start:    #in case we have an exon that overlaps with the regulated exon
                distance += reg_exon_start - exon_start
            return distance
    return distance
            
# this script will return exonic length in CDS; input is array of exonic positions, regulated_exon and CDS length
def get_exonic_lenght(exon_positions, reg_exon_pos, cds_length):
    reg_exon_start = int(reg_exon_pos.rsplit(':')[0])
    reg_exon_end = int(reg_exon_pos.rsplit(':')[1])
    total_exonic_length = 0
    overlap = False
    for i in range(0,exon_positions.__len__()): 
        exon_start = int(exon_positions[i].rsplit(':')[0])  #we count with 0
        exon_end = int(exon_positions[i].rsplit(':')[1])
        exon_length = exon_end - exon_start + 1
        
        if i == 0 and exon_start < 0:   #in case the first exon is overlapping CDS start
            exon_length = exon_end - 1  # we count with 1 now that's why -1
            
        if reg_exon_start <= exon_start and reg_exon_end >= exon_end:   # in case regulated exon fully overlap the exon
            exon_start = reg_exon_start
            exon_end = reg_exon_end
            overlap = True
            
        if reg_exon_end >= exon_start and reg_exon_end <= exon_end: #in case regulated exon overlaps upstream we extend it 
            exon_start = reg_exon_start
            if reg_exon_end >= exon_end:    #we also look at the end
                exon_end = reg_exon_end
            overlap = True
        
        if reg_exon_start >= exon_start and reg_exon_start <= exon_end: #in case regulated exon overlaps downstream we extend it 
            exon_end = reg_exon_end
            overlap = True 
        
        if exon_end <= cds_length:
            total_exonic_length += exon_length
        elif exon_start < cds_length:
            exon_length = cds_length - exon_start + 1
            total_exonic_length += exon_length
            
    if overlap == False:
        total_exonic_length += reg_exon_end - reg_exon_start + 1
    return total_exonic_length



def get_total_number_of_exons_in_cds(exon_positions, cds_length):
    exon_num = 0
    for i in range(0,exon_positions.__len__()): 
        exon_start = int(exon_positions[i].rsplit(':')[0])  #we count with 0
        exon_end = int(exon_positions[i].rsplit(':')[1])
        if exon_start < cds_length:
            exon_num += 1
    return exon_num

def get_number_of_exons_upstream_reg_exon(reg_exon_pos, exon_positions):
    reg_exon_start = int(reg_exon_pos.rsplit(':')[0])
    reg_exon_end = int(reg_exon_pos.rsplit(':')[1])
    exon_num = 0
    for i in range(0,exon_positions.__len__()): 
        exon_start = int(exon_positions[i].rsplit(':')[0])  #we count with 0
        exon_end = int(exon_positions[i].rsplit(':')[1])
        if reg_exon_start > exon_end:
            exon_num += 1
    return exon_num

# script will return a distance from PTC to the downstream junction (end of the second exon)
def get_downstream_junction_position(ptc, reg_exon_pos, exon_positions):
    reg_exon_start = int(reg_exon_pos.rsplit(':')[0])
    reg_exon_end = int(reg_exon_pos.rsplit(':')[1])
    reg_exon_length = reg_exon_end - reg_exon_start
    PTC_exon = False    # we search for the downstream exon from regulated exon
    downstream_junction_distance = -1
    for i in range(0,exon_positions.__len__()): 
        exon_start = int(exon_positions[i].rsplit(':')[0])  #we count with 0
        exon_end = int(exon_positions[i].rsplit(':')[1])
        if PTC_exon:
            exon_length = exon_end - exon_start
            downstream_junction_distance = reg_exon_length - ptc + exon_length
            return downstream_junction_distance
        if reg_exon_start > exon_end:
            if reg_exon_end < exon_end: #then they overlap and regulated_exon end is shorter then exon end so it will be extended
                reg_exon_length = exon_end - reg_exon_start
            PTC_exon = True
    return downstream_junction_distance
            
# script will receive a sequence and motid and return a list of all positions with containing that motif
def get_motif_positions(seq, motif, positions):
    motif_pos = seq.find(motif)
    if motif_pos == -1:
        return positions
    else:
        positions.append(motif_pos)
    pos = motif_pos
    while motif_pos != -1:
        motif_pos = seq[pos+1:].find(motif) #we search for the same motif downstream fro mthe previous one
        if motif_pos != -1:
            pos = motif_pos + pos + 1
            positions.append(pos)
    return positions


def filter(fname_in, fname_out):
    fin = open(fname_in, "rt")
    fout = open(fname_out, "w")
    line = fin.readline()
    fout.write("reg_exon_id\tstatus\ttotal_exonic_distance\tPTC_pos\tdistance_between_PTC_EJC\tPTC_seq\treg_exon_seq\tremaining_exonic_CDS_lenght\ttotal_exonic_CDS_length\tnum_usptream_exons\tnum_downstream_exons\tdistance_from_PTC_to_second_EJC\n")    #header
    while line:
        PTC = False
        info = line.rstrip('\n').replace('>','')
        info_col = info.rsplit('|')
        reg_exon_id = info_col[0]
        reg_exon_pos = info_col[1]  #regulated exon position
        reg_exon_start = int(reg_exon_pos.rsplit(':')[0]) - 1  #we count with 0
        reg_exon_end = int(reg_exon_pos.rsplit(':')[1])
        exon_positions = info_col[2].rsplit(',')    #all exon positions on CDS
        
        seq = fin.readline().rstrip('\n')
        seq = str(seq).upper()
        cds_length = seq.__len__()
        distance = get_exon_coding_position(reg_exon_pos, exon_positions)   #get the length of all exons from the start of CDS to regulated exon position
        total_exonic_length = get_exonic_lenght(exon_positions,reg_exon_pos, cds_length)
            
        num_total_exon_number = get_total_number_of_exons_in_cds(exon_positions, cds_length)
        num_upstream_exons_from_regulated = get_number_of_exons_upstream_reg_exon(reg_exon_pos, exon_positions)
        num_downstream_exons_from_regulated = num_total_exon_number - num_upstream_exons_from_regulated
            
        if distance % 3 == 0:   # if position is dividable by 3 look for TAG, TGA, TAA (either of them) and if codon position relative to regulated exons start is dividable by 3 as well then YES
            exon_seq = seq[reg_exon_start:reg_exon_end]
            ptc_pos = []
            ptc_pos = get_motif_positions(exon_seq, "TAG", ptc_pos)
            ptc_pos = get_motif_positions(exon_seq, "TGA", ptc_pos)
            ptc_pos = get_motif_positions(exon_seq, "TAA", ptc_pos)
            
            if ptc_pos.__len__() > 0:
                ptc_pos.sort()
                for i in range(0,ptc_pos.__len__()):
                    if (ptc_pos[i]+3) % 3 == 0:
                        remaining_CDS_lenght = total_exonic_length - distance - ptc_pos[i]+1
                        downstream_junction_distance = get_downstream_junction_position(ptc_pos[i]+1, reg_exon_pos, exon_positions)
                        fout.write(reg_exon_id + '\t'+ "PTC" + '\t' + str(distance) + '\t' + str(reg_exon_start+ptc_pos[i]+1) + '\t' + str(exon_seq.__len__() - ptc_pos[i]) + '\t' + exon_seq[ptc_pos[i]:(ptc_pos[i]+3)] + '\t' + exon_seq + '\t' + str(remaining_CDS_lenght) + '\t' + str(total_exonic_length) + '\t' + str(num_upstream_exons_from_regulated) + '\t' + str(num_downstream_exons_from_regulated) + '\t' + str(downstream_junction_distance) + '\n')
                        #print(reg_exon_id + '\t'+ "PTC" + '\t' + str(distance) + '\t' + str(reg_exon_start+ptc_pos[i]+1) + '\t' + exon_seq[ptc_pos[i]:(ptc_pos[i]+3)] + '\t' + exon_seq + '\n')
                        PTC = True
                        break
            
        elif (distance+1) % 3 == 0:   #  if not, if +1 one make it dividable by 3 then codon position to relative start also needs to be + 1 devidable by 3
            exon_seq = seq[reg_exon_start:reg_exon_end]
            ptc_pos = []
            ptc_pos = get_motif_positions(exon_seq, "TAG", ptc_pos)
            ptc_pos = get_motif_positions(exon_seq, "TGA", ptc_pos)
            ptc_pos = get_motif_positions(exon_seq, "TAA", ptc_pos)
            
            if ptc_pos.__len__() > 0:
                ptc_pos.sort()
                for i in range(0,ptc_pos.__len__()):
                    if (ptc_pos[i]+2) % 3 == 0: # +2 !!!! we shifted the reading frame
                        remaining_CDS_lenght = total_exonic_length - distance - ptc_pos[i]+1
                        downstream_junction_distance = get_downstream_junction_position(ptc_pos[i]+1, reg_exon_pos, exon_positions)
                        fout.write(reg_exon_id + '\t'+ "PTC+1" + '\t' + str(distance) + '\t' + str(reg_exon_start+ptc_pos[i]+1) + '\t' + str(exon_seq.__len__() - ptc_pos[i]) + '\t' + exon_seq[ptc_pos[i]:(ptc_pos[i]+3)] + '\t' + exon_seq + '\t' + str(remaining_CDS_lenght) + '\t' + str(total_exonic_length) + '\t' + str(num_upstream_exons_from_regulated) + '\t' + str(num_downstream_exons_from_regulated) + '\t' + str(downstream_junction_distance) + '\n')
                        #print(reg_exon_id + '\t'+ "PTC+1" + '\t' + str(distance) + '\t' + str(reg_exon_start+ptc_pos[i]+1) + '\t' + exon_seq[ptc_pos[i]:(ptc_pos[i]+3)] + '\t' + exon_seq + '\n')
                        PTC = True
                        break
                
        elif (distance+2) % 3 == 0:   #if it's +2 .....
            exon_seq = seq[reg_exon_start:reg_exon_end]
            ptc_pos = []
            ptc_pos = get_motif_positions(exon_seq, "TAG", ptc_pos)
            ptc_pos = get_motif_positions(exon_seq, "TGA", ptc_pos)
            ptc_pos = get_motif_positions(exon_seq, "TAA", ptc_pos)
            
            if ptc_pos.__len__() > 0:
                ptc_pos.sort()
                for i in range(0,ptc_pos.__len__()):
                    if (ptc_pos[i]+1) % 3 == 0: # +1 !!!! we shifted the reading frame
                        remaining_CDS_lenght = total_exonic_length - distance - ptc_pos[i]+1
                        downstream_junction_distance = get_downstream_junction_position(ptc_pos[i]+1, reg_exon_pos, exon_positions)
                        fout.write(reg_exon_id + '\t'+ "PTC+2" + '\t' + str(distance) + '\t' + str(reg_exon_start+ptc_pos[i]+1) + '\t' + str(exon_seq.__len__() - ptc_pos[i]) + '\t' + exon_seq[ptc_pos[i]:(ptc_pos[i]+3)] + '\t' + exon_seq + '\t' + str(remaining_CDS_lenght) + '\t' + str(total_exonic_length) + '\t' + str(num_upstream_exons_from_regulated) + '\t' + str(num_downstream_exons_from_regulated) + '\t' + str(downstream_junction_distance)  + '\n')
                        #print(reg_exon_id + '\t'+ "PTC+2" + '\t' + str(distance) + '\t' + str(reg_exon_start+ptc_pos[i]+1) + '\t' + exon_seq[ptc_pos[i]:(ptc_pos[i]+3)] + '\t' + exon_seq + '\n')
                        PTC = True
                        break
                    
        if PTC == False:
            if (reg_exon_end - reg_exon_start) % 3 == 0:   # and then we check if the regulated exon length is devidable by 3 (it's a good exon) if not it's a bady
                fout.write(reg_exon_id + '\t'+ "GOOD" + '\t' + str(distance) + '\t' + "NA" + '\t' + "NA" + '\t' + "NA" + '\t' + exon_seq + '\t' + "NA"  + '\t' + str(total_exonic_length) + '\t' + str(num_upstream_exons_from_regulated) + '\t' + str(num_downstream_exons_from_regulated) + '\t' + "NA" + '\n')
            else:        
                fout.write(reg_exon_id + '\t'+ "BAD" + '\t' + str(distance) + '\t' + "NA" + '\t' + "NA" + '\t' + "NA" + '\t' + exon_seq + '\t' + "NA"  + '\t' + str(total_exonic_length) + '\t' + str(num_upstream_exons_from_regulated) + '\t' + str(num_downstream_exons_from_regulated) + '\t' + "NA" + '\n')
        
        line = fin.readline()
    fout.close()
    fin.close()


'''

filter("/media/nebo/SAMSUNG/UCL-backup/2015.04.27@ALUs-Jan/CDSs-PTC-regulated-exons/hg19-cds-gene_symbols-longest-CDS2-regulated_exons-positions.fasta","/media/nebo/SAMSUNG/UCL-backup/2015.04.27@ALUs-Jan/CDSs-PTC-regulated-exons/test-results.tab")

'''
if sys.argv.__len__() == 3:
    fname_in = sys.argv[1]
    fname_out = sys.argv[2]
    filter(fname_in, fname_out)
else:
    print("You need two arguments to run the script")

    
    
    