from urllib.request import urlopen  # Python 3 version of urllib2
import pprint  # importing pretty printing
from Bio import Entrez  # importing Entrez
from Bio import SeqIO  # importing SeqIO
from Bio.SeqUtils import GC

# QUESTION 1
print('*********QUESTION 1********************************************************************************************')
print('------------------HBA1 SEQUENCE ID, NAME, DESCRIPTION----------------------------------------------------------')
Entrez.email = 'im1247@nyu.edu'    # inserting email
HBA1_document = Entrez.efetch(db="nucleotide", rettype="gb", id="NM_000558")    # getting HBA1
open_doc1 = open("NM_000558.fasta", 'w')
gbseq1 = SeqIO.read(HBA1_document, "genbank")
SeqIO.write(gbseq1, open_doc1, "fasta")
open_doc1.close()
# Print the sequence ID, name, and description of these sequence records.
print(gbseq1)   # this is printing out the description of HBA1
print('-------------------SEQUENCE OF HBA1----------------------------------------------------------------------------')
pprint.pprint(gbseq1.seq)  # pretty printing out the sequence
print('-------------------HBA2 SEQUENCE ID, NAME, AND DESCRIPTION-----------------------------------------------------')
HBA2_document = Entrez.efetch(db="nucleotide", rettype="gb", id="NM_000517")  # getting the nucleotide sequence for HBA2
open_doc2 = open("NM_000517.fasta", 'w')
gbseq = SeqIO.read(HBA2_document, "genbank")
SeqIO.write(gbseq, open_doc2, "fasta")
HBA1_document.close()
HBA2_document.close()
open_doc2.close()
print(gbseq)
print('---------------------SEQUENCE OF HBA2--------------------------------------------------------------------------')
pprint.pprint(gbseq.seq)    # pretty printing of sequence

# QUESTION 2
print('*********QUESTION 2********************************************************************************************')
response = urlopen('https://raw.githubusercontent.com/izzykayu/HW_Week3_Bioinformatics/master/seq_id_list.fasta')
out = open('q2_seq_id_list.fasta', 'w')  # out will be the file that the html fasta will be written to
fasta = response.read().decode("utf-8", "ignore")  # allows ls_orchid fasta information be retrieved from url link
# print(fasta)  # printing to see how it looks
out.write(fasta)  # this writes the html fasta file onto fasta file in my directory
out.close()  # closes the file
# from  a text file with sequence IDs --> making a python list
record_id_list = []  # creating a list to store the sequence name and lengths that will be written to a file
with open('q2_seq_id_list.fasta', "r") as inf:
    for line in inf:
        values = line.split()
        record_id_list.append(values)

flattened = [val for sublist in record_id_list for val in sublist]  # removing extra brackets from record_id_list
print(flattened)  # printing seq id list

record_list = []
for i in flattened:
    q2_seq = Entrez.efetch(db="nucleotide", rettype="gb", id=i)  # getting the nucleotide sequence for HBA2
    q2gbseq = SeqIO.read(q2_seq, "genbank")
    len_str = str(len(q2gbseq.seq))
    # print(q2gbseq.seq)
    # print(len_str)
    record = 'sequence name:      ' + q2gbseq.id + '\nlength of sequence: ' + len_str
    print(record)
    record_list.append(record)    # this prints sequence name and length of sequence

the_file = open('list_genbank_sequences.txt', 'w')  # creating a file with the list of sequence names and lengths
for item in record_list:
    the_file.write("%s\n" % item)  # this indicates the format


the_file.close()  # to ensure file saves

# QUESTION 3
print('*********QUESTION 3********************************************************************************************')
# url of ecogene.fasta
# Import a large set of sequences from a Fasta file: ecogene.fasta
ecogene = urlopen('http://fenyolab.org/presentations/Introduction_Biostatistics_Bioinformatics_2015/pdf/lecture5/ecogene.fasta')
out_f = open('ecogene.fasta', 'w')  # out will be the file that the html fasta will be written to
ecogene_fasta = ecogene.read().decode("utf-8", "ignore")  # allows ls_orchid fasta information be retrieved from url link
# print(fasta)  # printing to see how it looks
out_f.write(ecogene_fasta)  # this writes the html fasta file onto ecogene fasta file in my directory
out_f.close()  # closes the file
name_list = []  # list of ids
length_list_all = []  # list of all lengths as numbers
length_list_records = []  # list to record sequences with lengths over 300
GC_list_records = []    # list to record sequences with GC > 60
for eco_seq_record in SeqIO.parse('ecogene.fasta', 'fasta'):
    length_ = len(eco_seq_record.seq)
    length_list_all.append(length_)
    name_ = eco_seq_record.id
    name_list.append(name_)
    if length_ >= 300:  # getting list of sequences for with length equal to or over 300
        length_list_records.append(eco_seq_record.seq)
    if GC(eco_seq_record.seq) > 60:
        GC_list_records.append(eco_seq_record.seq)

# A. How many sequences are in this file?
print('NUMBER OF SEQ IN FILE:  ', len(name_list))
# 4294 sequences in this file

print("RECORD OF FIRST SEQUENCE IN ECOGENE FASTA FILE:")
# B. printing sequence ID, name, and description for the first Sequence Record.
first_record = next(SeqIO.parse("ecogene.fasta", "fasta"))
print(first_record)
# this record is much less detailed than the metabase genbank records in question 1:
# ID: eschColi_K12_refSeq_b0001
# Name: eschColi_K12_refSeq_b0001
# Description: eschColi_K12_refSeq_b0001 range=chr:190-255 5'pad=0 3'pad=0 strand=+ repeatMasking=none
# Number of features: 0
# Seq('ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAAC...TGA', SingleLetterAlphabet())

# C. What is the total length of all of the sequences in the file (just the DNA, not
# headers)
print("TOTAL LENGTH OF ALL SEQ IN FILE:  ", sum(length_list_all))
#  LENGTH OF ALL SEQ IN FILE:   4130063

# D. Make a new FASTA file with just the sequences >= 300 bp in length
length_file = open('length_over_300_sequences.txt', 'w')  # creating a file with the list of seq w/ lengths >= 300
for s in length_list_records:   # s stands for sequence
    length_file.write("%s\n" % s)  # this indicates the format


length_file.close()  # to ensure file saves
# E. Make a new FASTA file with just the sequences with %GC > 60
GC_file = open('GC_over_60_sequences.txt', 'w')  # creating a file with the list of GC > 60 sequences
for s in GC_list_records:   # s stands for sequence
    GC_file.write("%s\n" % s)  # this indicates the format


GC_file.close()     # to ensure file saves
