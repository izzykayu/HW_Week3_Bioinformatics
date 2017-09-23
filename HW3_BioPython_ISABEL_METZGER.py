from urllib.request import urlopen  # Python 3 version of urllib2
import pprint  # importing pretty printing
from Bio import Entrez  # importing Entrez
from Bio import SeqIO  # importing SeqIO

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

# handle = Entrez.esearch(db="nucleotide",term="E coli[Orgn] AND RefSeq", idtype="refseq", rettype="gb", retmode="text", retmax='1000')
# record = Entrez.read(handle)
# print(record["Count"])  # count of all records
# eco_id_list = record["IdList"]  #IDs of E coli
# print('------------\n', eco_id_list)

# for i in eco_id_list:
#     eco_seq = Entrez.efetch(db="nucleotide", rettype="gb", id=i)  # getting the nucleotide sequence for HBA2
#     gbseq = SeqIO.read(eco_seq, "genbank")
#     print(gbseq)