
import random
from Bio import SeqIO
import os
import glob
import pandas as pd
from random import randint


class Base(object):
    A = 'A'
    C = 'C'
    G = 'G'
    T = 'T'
    N = 'N'

def pick_random_sequence(fasta_path):

    records_dict = SeqIO.index(fasta_path, "fasta")
    print "pick random number between 0 and {}\n\n".format(len(records_dict))
    random_index = random.randint(0, len(records_dict))
    i = 0
    random_title = ''
    for title in records_dict.keys():
        if i == random_index:
            random_title = title
            break
        i += 1
    print "random_title {}".format(random_title)
    records_dict.close()
    return random_title


def get_sequence_by_title(fasta_path, title):
    records = SeqIO.index(fasta_path, "fasta")
    try:
        res =  records[title].seq
        res = res.__str__()
    except KeyError:
        res = None
    records.close()
    return res



def pick_random_sequnce_from_fasta(fasta_dir):
    os.chdir(fasta_dir)
    files = glob.glob("*.fasta")
    title=None
    sequences = []
    regions = []
    for file in files:
        region = -1
        for i in range(1, 5 + 1):
            if str(i) in file:
                region = i
        if title is None:
            title = pick_random_sequence(file)
        sequence = get_sequence_by_title(file, title)
        if sequence:
            sequences.append(sequence)
            regions.append(region)
    print "Found {} regions for title {}".format(len(sequences), title)
    return regions, sequences, title

def get_reads_from_sequence(sequence, read_len=126):
    read = "".join(sequence[:read_len])
    print read

    base_dict = {Base.A: Base.T,
                 Base.G: Base.C,
                 Base.C: Base.G,
                 Base.T: Base.A,
                 Base.N: Base.N}

    end_of_seq = sequence[-126:]
    # print end_of_seq
    reversed_end_of_seq = end_of_seq[::-1]
    # print reversed_end_of_seq
    comp_rev_end_of_seq = "".join([base_dict[b] for b in reversed_end_of_seq])
    print '\n\nrev:'
    print comp_rev_end_of_seq
    print 'end of seq'
    print sequence[-126:]
    print 'read:'
    print read
    print 'begin sequence:'
    print sequence[:126]
    return read, comp_rev_end_of_seq, read+end_of_seq


def change_sequence(sequence, bases_to_change, read_len):
    changed_indexes = []
    for i in range(0, bases_to_change):
        index = randint(0, 2*read_len -1 )
        while index in changed_indexes:
            index = randint(0, 2*read_len -1 )
        changed_indexes.append(index)
        sequence[index].replace(sequence[index], change_base(sequence[index]))
    print "change indexes = {}".format(changed_indexes)
    return sequence


def change_base(old_base):
    bases_dict = {0:'A', 1: 'C', 2: 'G', 3:'T'}
    new_base = bases_dict[randint(0, 3)]
    while new_base == old_base:
        new_base  = bases_dict[randint(0, 3)]
    print "change base {} -> {}".format(old_base, new_base)
    return new_base


def get_random_reads(fasta_dir, read_len, amount, bases_to_change_for_region=0, region_to_change=0):
    regions, seqs, title = pick_random_sequnce_from_fasta(fasta_dir)
    reads = []
    reversed_reads = []
    sequences = []
    i = 0
    for seq in seqs:
        if i < region_to_change:
            seq = change_sequence(seq, bases_to_change_for_region, read_len)
            i += 1
            print("change the reference title = {}".format(title))
        read, rev_read, sequence= get_reads_from_sequence(seq, read_len)
        reads += amount*[read]
        reversed_reads += amount*[rev_read]
        sequences.append(sequence)
    return regions, reads, reversed_reads, sequences, title

def create_mock_fastq_from_fastq_file(fasta_dir, read_len, fastq, new_fastq_path, new_rev_fastq_path,
                                      list_of_percents=[100], ref_amount = 1000, ref_to_change=0, regions_to_change=0,
                                      bases_to_change=0, write_to_fasta=False):
    reads = []
    rev_reads = []
    seqs_for_logs = []
    priors = []
    titles = []
    regions = []
    i = 0
    for percent in list_of_percents:
        amount = int(round((ref_amount/5)*percent/100))
        if i<ref_to_change:
            region_part, reads_part, rev_reads_part, sequences, title = get_random_reads(fasta_dir, read_len, amount, bases_to_change, regions_to_change)
            i += 1
        else:
            region_part, reads_part, rev_reads_part, sequences, title = get_random_reads(fasta_dir, read_len, amount)
        reads += reads_part
        rev_reads += rev_reads_part
        seqs_for_logs += sequences
        priors += [float(percent)/100]*len(sequences)
        titles += [title]*len(sequences)
        regions += region_part

    print "create mock data len = {}".format(len(seqs_for_logs))
    data =pd.DataFrame.from_dict({'sequence': seqs_for_logs, 'prior': priors, 'id': titles, 'region': regions})
    data.to_csv("/home/vered/EMIRGE/EMIRGE-data/mock_reads/expected_res.csv", index=False)

    if write_to_fasta:
        write_reads_to_fasta(reads, rev_reads, new_fastq_path, new_rev_fastq_path)
    else:
        write_reads_to_fastq(reads, rev_reads, fastq, new_fastq_path, new_rev_fastq_path)


def write_reads_to_fasta(reads, rev_reads, new_fastq_path, new_rev_fastq_path):
    print "reads = {} = {}".format(len(reads), len(rev_reads))
    with open(new_fastq_path, 'w') as new_fa:
        with open(new_rev_fastq_path, 'w') as new_rev_fa:
            for i in range(0, len(reads)):
                new_fa.write(">" + str(i) + "\n")
                new_fa.write(reads[i] + "\n")
                new_fa.write("\n")

                new_rev_fa.write(">" + str(i) + "\n")
                new_rev_fa.write(rev_reads[i] + "\n")
                new_rev_fa.write("\n")
    print "Done write_reads_to_fasta"


def write_reads_to_fastq(reads, rev_reads, fastq, new_fastq_path, new_rev_fastq_path):
    i = 0
    r_ix = 0
    print "reads = {} = {}".format(len(reads), len(rev_reads))
    with open(fastq, 'r') as ref_fastq:
        with open(new_fastq_path, 'w') as new_fastq:
            with open(new_rev_fastq_path, 'w') as new_rev_fastq:
                for line in ref_fastq:
                    if r_ix < len(reads) or (r_ix == len(reads) and (i%4 != 0)):
                        if i % 4 == 1:
                            new_fastq.write(reads[r_ix] + '\n')
                            new_rev_fastq.write(rev_reads[r_ix] + '\n')
                            r_ix += 1
                            # print r_ix
                        else:
                            new_fastq.write(line)
                            new_rev_fastq.write(line)
                        i += 1
                    else:
                        break
    print "\nwrote {} reads to fastq file".format(r_ix)

def create_mock_fasta(old_fasta_path, new_fasta_path, amount):
    os.chdir(old_fasta_path)
    old_fasta_files = glob.glob("*.fasta")
    for file in old_fasta_files:
        new_file_path = os.path.join(new_fasta_path, file)
        file_path = os.path.join(old_fasta_path, file)
        with open(file_path, 'r') as f:
            lines = f.readlines()
            lines = [lines[i] for i in range(amount*7)]
            with open(new_file_path, "w") as f1:
                f1.writelines(lines)


def get_percents_list(amount_of_refs=15, max_percent=20):
    list_of_percents = []
    for i in range(0, amount_of_refs-1):
        remainder = (100 - sum(list_of_percents))/2
        remainder = remainder if remainder > 0 else 0
        percent = randint(1, min(remainder, max_percent ))
        list_of_percents.append(percent)
    print("get_percents_list -> {}".format(list_of_percents))
    curr_sum = sum(list_of_percents)
    list_of_percents.append(100 - curr_sum)
    for p in list_of_percents:
        if p < 0 :
            raise Exception("Percent must be positive number {}".format(p))
    return list_of_percents

def create_mock_reads():
    """
    using the full fasta DB as reference
    :return:
    # """
    # new_fastq = "/home/vered/EMIRGE/EMIRGE-data/mock_reads/RDB53_CATTGACG_L007_R1_001.fastq"
    # new_rev_fastq = "/home/vered/EMIRGE/EMIRGE-data/mock_reads/RDB53_CATTGACG_L007_R2_001.fastq"
    new_fastq = "/home/vered/EMIRGE/EMIRGE-data/mock_reads/RDB53_CATTGACG_L007_R1_001.fa"
    new_rev_fastq = "/home/vered/EMIRGE/EMIRGE-data/mock_reads/RDB53_CATTGACG_L007_R2_001.fa"
    fastq = "/home/vered/EMIRGE/EMIRGE-data/RDB53_CATTGACG_L007_R1_001.fastq"
    fasta_dir = "/home/vered/EMIRGE/EMIRGE-data"
    read_len = 126

    percents = get_percents_list(15)
    sum = 0
    for l in percents:
        sum += l
    if abs(sum - 100) > 0.001:
        raise Exception("Not 100% -->{}%".format(sum))

    create_mock_fastq_from_fastq_file(fasta_dir, read_len, fastq, new_fastq, new_rev_fastq, list_of_percents=percents,
                                      ref_amount=130000, ref_to_change=0, regions_to_change=0, bases_to_change=0, write_to_fasta=True)

def create_mock_reads_from_mock_refs():
    mock_refs = '/home/vered/EMIRGE/EMIRGE-data/mock_fasta/'
    new_fastq = "/home/vered/EMIRGE/EMIRGE-data/mock_reads/RDB53_CATTGACG_L007_R1_001.fastq"
    new_rev_fastq = "/home/vered/EMIRGE/EMIRGE-data/mock_reads/RDB53_CATTGACG_L007_R2_001.fastq"
    fastq = "/home/vered/EMIRGE/EMIRGE-data/RDB53_CATTGACG_L007_R1_001.fastq"
    read_len = 126

    percents = get_percents_list()
    sum = 0
    for l in percents:
        sum += l
    if abs(sum - 100) > 0.001:
        raise Exception("Not 100% -->{}%".format(sum))

    create_mock_fastq_from_fastq_file(mock_refs, read_len, fastq, new_fastq, new_rev_fastq, list_of_percents=percents,
                                      ref_amount=5000, ref_to_change=1, regions_to_change=5, bases_to_change=1)

def create_mock_refs():
    old_fasta_path = '/home/vered/EMIRGE/EMIRGE-data/'
    new_fasta_path = '/home/vered/EMIRGE/EMIRGE-data/mock_fasta/'
    amount = 100
    create_mock_fasta(old_fasta_path, new_fasta_path, amount)

if __name__ == '__main__':
   create_mock_reads()
   # create_mock_refs()
   # create_mock_reads_from_mock_refs()
   # create_mock_reads_from_mock_refs()


# #mock reads - 01/08/2017, changed reference
# change base G -> A
# change indexes = [69]
# change the reference title = 4029628.03