#!/usr/bin/env python
'''
Created on 29 March 2016

@author: Igor Ruiz de los Mozos

Check the best 3SS score in 2nt upstream and 2nt downstream

NOTE copy the folder of splicemodels from MAXSCAN3SS to the running folder

It will print the best 3´ss score on the defined 3´ss two nucleotides upstream and two nucleotides downstream

Usage:
        Usage: \t python get_best_3SS.py bedIN bedOUT tempDir


'''
#Genome to get fasta
hg19='/media/igor/DATA/UCL/Genomes/Human/ucsc.hg19.fasta'


import sys
import operator
import subprocess




def runUNIXCommands(commands):
    #print commands
    try:
        retcode = subprocess.check_call(commands, shell=True)
        if retcode < 0:
            print >>sys.stderr, "Child was terminated by signal", -retcode
        else:
            #print >>sys.stderr, "Child returned", retcode
            pass
    except OSError, e:
        print >>sys.stderr, "Execution failed:", e

runUNIXCommands("cp -r " + "/home/igor/Programs/MaxEntScan/splicemodels" + " .")

def flank_positions(line_in, left_shift, right_shift):
    line = line_in
    left_shift = int(left_shift)
    right_shift = int(right_shift)
    while line:
        col = line.rstrip('\n').rsplit('\t')
        if line.__len__() > 6:
            chr = col[0]
            start = col[1]
            end = col[2]
            strand = col[5]
            shift_start = shift_end = 0
            if strand == '+':
                shift_start = int(start) - int(left_shift)
                shift_end = int(end) + int(right_shift)
            else:
                shift_end = int(end) + int(left_shift)
                shift_start = int(start) - int(right_shift)
            if shift_start >= 0:
                line_out = (chr + '\t' + str(shift_start) + '\t' + str(shift_end) + '\t' + col[3] + '\t' + col[0] + ':' + col[1] + ':' + col[2] + '\t' + col[5] + '\n')
                #line_out = (chr + '\t' + str(shift_start) + '\t' + str(shift_end) + '\t' + col[5] + '\t' + '' + '\t' + col[5] + '\n')
        return line_out
        #line = fin.readline()  ### Just read the first line




def get_best_3SS(fin_fname, fout_name, genome_dir):
    fin = open(fin_fname, "rt")
    fout = open(fout_fname, "w")
    line = fin.readline()
    vc1 = vc2 = vc3 = vc4 = vc5 = v_error = 0           # Initialize variables to count each position


    while line:
        col = line.rstrip('\n').rsplit('\t')
        chr = col[0]
        start = col[1]
        end = col[2]
        strand = col[5]
        runUNIXCommands("mkdir -p " + genome_dir)
        v1 = v2 = v3 = v4 = v5 = 0                      # Initialize variables to get the best 3SS score
        
        ######## original 3SS ########

        temp_out_v1 = (genome_dir + chr + "_" + start + "_" + end +"_v1.bed")
        out = open(temp_out_v1, "w")
        out.write(flank_positions(line, "19", "3"))   # 19nt upstream 3 downstream
        out.close()
        temp_out2_v1 = temp_out_v1.replace('_v1.bed', '')
        runUNIXCommands('bedtools ' + ' getfasta -s -fi '+ hg19 + ' -bed ' + temp_out2_v1 + '_v1.bed -fo ' + temp_out2_v1 + '_v1.fasta' )
        runUNIXCommands("perl /home/igor/Programs/MaxEntScan/score3.pl " + temp_out2_v1 + "_v1.fasta > " + temp_out2_v1 +"_score_v1.fasta" )
        fin_v1 = open(temp_out2_v1 +"_score_v1.fasta", "rt")
        line_v1 = fin_v1.readline()
        fin_v1.close()
        toq_v1 = line_v1.rstrip('\n').rsplit('\t')
        if toq_v1.__len__() > 1:
            v1 = toq_v1[1]
            #print "SS score v1    " + str(v1)


        ######## 3SS 1nt upstream ########

        temp_out_v2 = (genome_dir + chr + "_" + start + "_" + end +"_v2.bed")
        out = open(temp_out_v2, "w")
        out.write(flank_positions(line, "20", "2"))
        out.close()
        temp_out2_v2 = temp_out_v2.replace('_v2.bed', '')
        runUNIXCommands('bedtools ' + ' getfasta -s -fi '+ hg19 + ' -bed ' + temp_out2_v2 + '_v2.bed -fo ' + temp_out2_v2 + '_v2.fasta' )
        runUNIXCommands("perl /home/igor/Programs/MaxEntScan/score3.pl " + temp_out2_v2 + "_v2.fasta > " + temp_out2_v2 +"_score_v2.fasta" )
        fin_v2 = open(temp_out2_v2 +"_score_v2.fasta", "rt")
        line_v2 = fin_v2.readline()
        fin_v2.close()
        toq_v2 = line_v2.rstrip('\n').rsplit('\t')
        if toq_v2.__len__() > 1:
            v2 = float(toq_v2[1])
            #print "SS score v2  " + str(v2)


        ######## 3SS 2nt upstream ########

        temp_out_v3 = (genome_dir + chr + "_" + start + "_" + end +"_v3.bed")
        out = open(temp_out_v3, "w")
        out.write(flank_positions(line, "21", "1"))
        out.close()
        temp_out2_v3 = temp_out_v3.replace('_v3.bed', '')
        runUNIXCommands('bedtools ' + ' getfasta -s -fi '+ hg19 + ' -bed ' + temp_out2_v3 + '_v3.bed -fo ' + temp_out2_v3 + '_v3.fasta' )
        runUNIXCommands("perl /home/igor/Programs/MaxEntScan/score3.pl " + temp_out2_v3 + "_v3.fasta > " + temp_out2_v3 +"_score_v3.fasta" )
        fin_v3 = open(temp_out2_v3 +"_score_v3.fasta", "rt")
        line_v3 = fin_v3.readline()
        fin_v3.close()
        toq_v3 = line_v3.rstrip('\n').rsplit('\t')
        if toq_v3.__len__() > 1:
            v3 = float(toq_v3[1])
            #print "SS score v3  " + str(v3)



        ######## 3SS 1nt downstream ########

        temp_out_v4 = (genome_dir + chr + "_" + start + "_" + end +"_v4.bed")
        out = open(temp_out_v4, "w")
        out.write(flank_positions(line, "18", "4"))
        out.close()
        temp_out2_v4 = temp_out_v4.replace('_v4.bed', '')
        runUNIXCommands('bedtools ' + ' getfasta -s -fi '+ hg19 + ' -bed ' + temp_out2_v4 + '_v4.bed -fo ' + temp_out2_v4 + '_v4.fasta' )
        runUNIXCommands("perl /home/igor/Programs/MaxEntScan/score3.pl " + temp_out2_v4 + "_v4.fasta > " + temp_out2_v4 +"_score_v4.fasta" )
        fin_v4 = open(temp_out2_v4 +"_score_v4.fasta", "rt")
        line_v4 = fin_v4.readline()
        fin_v4.close()
        toq_v4 = line_v4.rstrip('\n').rsplit('\t')
        if toq_v4.__len__() > 1:
            v4 = float(toq_v4[1])
            #print "SS score v4  " + str(v4)


        
        ######## 3SS 2nt downstream ########

        temp_out_v5 = (genome_dir + chr + "_" + start + "_" + end +"_v5.bed")
        out = open(temp_out_v5, "w")
        out.write(flank_positions(line, "17", "5"))
        out.close()
        temp_out2_v5 = temp_out_v5.replace('_v5.bed', '')
        runUNIXCommands('bedtools ' + ' getfasta -s -fi '+ hg19 + ' -bed ' + temp_out2_v5 + '_v5.bed -fo ' + temp_out2_v5 + '_v5.fasta' )
        runUNIXCommands("perl /home/igor/Programs/MaxEntScan/score3.pl " + temp_out2_v5 + "_v5.fasta > " + temp_out2_v5 +"_score_v5.fasta" )
        fin_v5 = open(temp_out2_v5 +"_score_v5.fasta", "rt")
        line_v5 = fin_v5.readline()
        fin_v5.close()
        toq_v5 = line_v5.rstrip('\n').rsplit('\t')
        if toq_v5.__len__() > 1:
            v5 = toq_v5[1]
            #print "SS score v5  " + str(v5)


        runUNIXCommands("rm -r " + genome_dir)   # Remove the temporal directory


        # BedFile testing and writing
        # write only if any vx differ from 0
        #  if vx = 0               1) bedtools get fasta is out of scope OR
        #                          2) MAXSCANscore has fail because there are NNNNs in the sequence

        if v1 != 0 and v2 != 0 and v3 != 0 and v4 != 0 and v5 != 0:

            v1 = float(v1)
            v2 = float(v2)
            v3 = float(v3)
            v4 = float(v4)
            v5 = float(v5)

            # float values to dictionary and afterwards get the key of the MAX. value
        
            SSscores = {'v1': (v1), 'v2': (v2), 'v3': (v3), 'v4': (v4), 'v5': (v5)}

            best3SSscore = max(SSscores.iteritems(), key=operator.itemgetter(1))[0]
            #print  "Best 3SS score _improved   " + str(best3SSscore)  + "   " + str(SSscores[best3SSscore]) + "\n"
            #print best3SSscore


            if best3SSscore == "v5":        # 2 downstream
                fout.write(flank_positions(line, "-2", "2"))
                #print "V5   " + flank_positions(line, "-2", "2")
                #print "2 downstream"
                vc5 += 1

            elif best3SSscore == "v4":      # 1 downstream
                fout.write(flank_positions(line, "-1", "1"))
                #print "V4   " + flank_positions(line, "-1", "1")
                #print "1 downstream"
                vc4 += 1

            elif best3SSscore == "v1":      # Same 3SS
                fout.write(line)
                #print "V1   " + line
                #print "same 3SS"
                vc1 += 1

            elif best3SSscore == "v2":      # 1 upstream
                fout.write(flank_positions(line, "1", "-1"))
                #print "V2   " + flank_positions(line, "1", "-1")
                #print "1 upstream"
                vc2 +=1

            elif best3SSscore == "v3":      # 2 upstream
                fout.write(flank_positions(line, "2", "-2"))
                #print "V3   " + flank_positions(line, "2", "-2")
                #print "2 upstream"
                vc3 += 1
            else:
                print "Error:"


        else:
            print "########### bed out of range or with NNNNNNNNNNN"
            v_error += 1

        
        ### "NEXT LINE "
        line = fin.readline()


    total = vc1 + vc2 + vc3 + vc4 + vc5
    print "2 downstream " + str(vc5)
    print "1 downstream " + str(vc4)
    print "same 3SS     " + str(vc1)
    print "1 upstream   " + str(vc2)
    print "2 upstream   " + str(vc3)
    print "Total        " + str(total)
    print "BED REJECTED " + str(v_error)




if sys.argv.__len__() == 4:
    fin_fname = sys.argv[1]
    fout_fname = sys.argv[2]
    genome_dir = sys.argv[3]
    get_best_3SS(fin_fname, fout_fname, genome_dir)
else:
    #print str(sys.argv.__len__())
    print "error:\t4 \n Usage: \t python get_best_3SS.py bedIN bedOUT tempDir"



'''

fin_fname='/Users/Igor/Desktop/LiftOver_bedPositions/All_Aluexons_3SS_sort.bed'
#fin_fname='/Users/Igor/Desktop/LiftOver_bedPositions/test3SS_2.bed'
fout_fname='/Users/Igor/Desktop/LiftOver_bedPositions/All_Aluexons_3SS_corrected.bed'
genome_dir='/Users/Igor/Desktop/LiftOver_bedPositions/TEMP/'

'''