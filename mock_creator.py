import sys
from optparse import OptionParser, OptionGroup, SUPPRESS_HELP

import random
from Bio import SeqIO
import os
import glob
import pandas as pd
from random import randint
import shutil


class Base(object):
    A = 'A'
    C = 'C'
    G = 'G'
    T = 'T'
    N = 'N'

    base_dict = {A: T,
                 G: C,
                 C: G,
                 T: A,
                 N: N}


class FastaFile(object):
    def __init__(self, directory, file_name, total_num_of_regions=None, region=None):
        self.path = os.path.join(directory, file_name)
        self.file_name = file_name
        self.region = region
        self.total_num_of_regions = total_num_of_regions

    def initialize(self):
        if self.region == None:
            if self.total_num_of_regions is None:
                raise Exception("Missing 'total number of regions'")
            self.region = self.get_region_from_file_name()

    def get_region_from_file_name(self):
        region = -1
        for i in range(1, self.total_num_of_regions + 1):
            if str(i) in self.file_name:
                region = i
        return region

class MockBacterium(object):
    def __init__(self,
                 reference_id,
                 frequency,
                 read_length,
                 total_amount_of_reads_per_region,
                 bases_to_change_per_region=0,
                 regions_to_change=0
                 ):
        self.bases_to_change_per_region = bases_to_change_per_region
        self.regions_to_change = regions_to_change
        self.id = reference_id
        self.new_id = reference_id
        self.frequency = frequency
        self.amount_of_reads_per_region = self.frequency*total_amount_of_reads_per_region/100
        self.read_length = read_length
        self.amplified_regions = []
        self.full_sequences = []
        self.sequences = []
        self.reads = []
        self.pair_reads = []
        self.regions=[]

    def initialize(self, fasta_files):
        if self.regions_to_change > 0:
            self.new_id += '#'
        changed_regions = 0
        for file in fasta_files:
            sequence = self.get_sequence_by_title(file.path)
            if sequence is None: #the bacteria wasn't amplified in this region
                continue
            if sequence:
                if changed_regions < self.regions_to_change:
                    sequence = self.change_sequence(sequence[:self.read_length] + sequence[-1 * self.read_length:])
                    changed_regions += 1
                self.sequences.append(sequence[:self.read_length] + sequence[-1 * self.read_length:])
                self.regions.append(file.region)
                self.full_sequences.append(sequence)
                self.amplified_regions.append(file.region)
                self.reads.append("".join(sequence[:self.read_length]))
                self.pair_reads.append(self.extract_paired_read_from_sequence(sequence))

    def extract_paired_read_from_sequence(self, sequence):
        end_of_seq = sequence[-1 * self.read_length:]
        reversed_end_of_seq = end_of_seq[::-1]
        pair_read = "".join([Base.base_dict[b] for b in reversed_end_of_seq])

        return pair_read

    def get_sequence_by_title(self, fasta_path):

        records = SeqIO.index(fasta_path, "fasta")
        try:
            res = records[self.id].seq
            res = res.__str__()
        except KeyError:
            print "Could not find the id = {} in fasta = {}".format(self.id, fasta_path)
            res = None
        records.close()

        return res

    def change_sequence(self, sequence):
        changed_indexes = []
        for i in range(0, self.bases_to_change_per_region):
            index = randint(0, 2 * self.read_length - 1)
            while index in changed_indexes:
                index = randint(0, 2 * self.read_length - 1)
            changed_indexes.append(index)
            if index < len(sequence)-1:
                sequence = sequence[:index] + self.change_base(sequence[index]) + sequence[index+1:]
            else:
                sequence = sequence[:index] + self.change_base(sequence[index])
        print "change indexes = {}".format(changed_indexes)
        return sequence

    def change_base(self, old_base):
        bases_dict = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
        new_base = bases_dict[randint(0, 3)]
        while new_base == old_base:
            new_base = bases_dict[randint(0, 3)]
        print "change base {} -> {}".format(old_base, new_base)
        return new_base

    def get_all_reads(self):
        reads = []
        for read in self.reads:
            reads += int(self.amount_of_reads_per_region)*[read]
        return reads

    def get_all_paired_reads(self):
        pair_reads = []
        for read in self.pair_reads:
            pair_reads += int(self.amount_of_reads_per_region) * [read]
        return pair_reads


def pick_random_bacterias(fasta_path, total_num_of_bacteria):
    random_ixs = []
    random_ids = []
    records_dict = SeqIO.index(fasta_path, "fasta")
    # print "pick_random_bacterias: pick random number between 0 and {}\n\n".format(len(records_dict))
    for i in range(total_num_of_bacteria):
        random_index = random.randint(0, len(records_dict))
        while random_index in random_ixs:
            random_index = random.randint(0, len(records_dict))
        random_ixs.append(random_index)

    i = 0
    for id in records_dict.keys():
        if i in random_ixs:
            random_ids.append(id)
            random_ixs.remove(i)
            if len(random_ixs) == 0:
                break
        i += 1

    records_dict.close()
    return random_ids


def get_fasta_files(fasta_dir, total_num_of_regions):
    fasta_files = []
    os.chdir(fasta_dir)
    files = glob.glob("*.fasta")
    for file in files:
        fasta_file = FastaFile(fasta_dir, file, total_num_of_regions)
        fasta_file.initialize()
        fasta_files.append(fasta_file)
    return fasta_files


def create_mock_bacterium(bacterium_frequencies,
                          fasta_dir,
                          total_num_of_regions,
                          read_length,
                          total_amount_of_reads_per_region,
                          bacterium_to_change,
                          regions_to_change,
                          bases_to_change,
                          keep_original=False):
    fasta_files = get_fasta_files(fasta_dir, total_num_of_regions)
    bacterias_ids = pick_random_bacterias(fasta_files[0].path, len(bacterium_frequencies))
    bacterium =[]
    for i in range(len(bacterium_frequencies)):
        if bacterium_to_change > 0:
            bacteria = MockBacterium(bacterias_ids[i], bacterium_frequencies[i], read_length, total_amount_of_reads_per_region, bases_to_change, regions_to_change)
            bacterium_to_change -= 1
            if keep_original:
                print "Keeping original bacteria {}. changed frequency = {}, original frequency = {}".format(bacterias_ids[i],
                                                                                                             bacterium_frequencies[i],
                                                                                                             bacterium_frequencies[i+1])
                original_bacteria = MockBacterium(bacterias_ids[i], bacterium_frequencies[i+1], read_length, total_amount_of_reads_per_region)
                original_bacteria.initialize(fasta_files)
                bacterium.append(original_bacteria)
                i += 1
        else:
            bacteria = MockBacterium(bacterias_ids[i], bacterium_frequencies[i], read_length, total_amount_of_reads_per_region)
        bacteria.initialize(fasta_files)
        bacterium.append(bacteria)

    print "Done create_mock_bacterium"
    return bacterium


def write_bacterium_to_fastq(mock_bacterium,
                             base_fastq_path,
                             mock_fastq_path1,
                             mock_fastq_path2):
    reads = []
    paired_reads = []
    for bacteria in mock_bacterium:
        reads.extend(bacteria.get_all_reads())
        paired_reads.extend(bacteria.get_all_paired_reads())
    with open(base_fastq_path, 'r') as base_fastq:
        with open(mock_fastq_path1, 'w') as mock_fastq1:
            with open(mock_fastq_path2, 'w') as mock_fastq2:
                read_ix = 0
                fastq_ix = 0
                for line in base_fastq:
                    if read_ix >= len(reads) and fastq_ix % 4 == 0:
                        break
                    if fastq_ix % 4 == 1:
                        mock_fastq1.write(reads[read_ix] + '\n')
                        mock_fastq2.write(paired_reads[read_ix] + '\n')
                        read_ix += 1
                    else:
                        mock_fastq1.write(line)
                        mock_fastq2.write(line)
                    fastq_ix += 1

    print "Done write_bacterium_to_fastq"


def write_bacterium_to_fasta(mock_bacterium,
                             mock_fasta_path1,
                             mock_fasta_path2,
                             min_frequency):
    with open(mock_fasta_path1, 'w') as mock_fasta1:
        with open(mock_fasta_path2, 'w') as mock_fasta2:
            for bacteria in mock_bacterium:
                sequences = bacteria.full_sequences*int(100*bacteria.frequency/min_frequency)
                for seq in sequences:
                    mock_fasta1.write(">{}\n".format(bacteria.id))
                    mock_fasta1.write(seq + '\n')
                    mock_fasta1.write('\n')

    print "Done write_bacterium_to_fasta"


def write_mock_bacterium_to_csv(mock_bacterium,
                                csv_path="/home/vered/EMIRGE/EMIRGE-data/mock_reads/expected_res.csv"):
    sequences = []
    frequencies = []
    ids = []
    regions = []
    sum_frequency = sum([b.frequency for b in mock_bacterium])
    for bacteria in mock_bacterium:
        for seq in bacteria.sequences:
            sequences.append(seq)
        regions += bacteria.regions
        frequencies += len(bacteria.regions)*[float(bacteria.frequency)/(100*sum_frequency)]
        ids += len(bacteria.regions)*[bacteria.new_id]

    data = pd.DataFrame.from_dict({'sequence': sequences, 'prior': frequencies, 'id': ids, 'region': regions})
    data.to_csv(csv_path, index=False)
    print "Done write_mock_bacterium_to_csv"


def create_mock_fasta(source_fasta_path, new_fasta_path, amount):
    os.chdir(source_fasta_path)
    source_fasta_files = glob.glob("*.fasta")
    reference_ids = pick_random_bacterias(os.path.join(source_fasta_path, source_fasta_files[0]), amount)
    for file in source_fasta_files:
        new_file_path = os.path.join(new_fasta_path, file)
        with open(new_file_path, "w") as f1:
            file_path = os.path.join(source_fasta_path, file)
            records = SeqIO.index(file_path, "fasta")
            for id in reference_ids:
                try:
                    res = records[id].seq
                    sequence = res.__str__()
                    f1.writelines(">{}\n{}\n\n".format(id, sequence))
                except Exception:
                    pass


def get_frequency_list(amount_of_refs=15, max_percent=20, min_frequency=0.5, isSorted=False):
    list_of_frequency = []
    for i in range(0, amount_of_refs-1):
        remainder = int((100 - sum(list_of_frequency))*min_frequency)
        remainder = remainder if remainder > 0 else 0
        percent = randint(1, min(remainder, max_percent )/min_frequency)*min_frequency  # accuracy of 0.5
        list_of_frequency.append(percent)
    curr_sum = sum(list_of_frequency)
    list_of_frequency.append(100 - curr_sum)
    for p in list_of_frequency:
        if p < 0 :
            raise Exception("Percent must be positive number {}".format(p))
    total = 0
    for l in list_of_frequency:
        total += l
    if abs(total - 100) > 0.001:
        raise Exception("frequency sum is not 100% -->{}%".format(total))
    if isSorted:
        list_of_frequency.sort()
    print("get_frequency_list -> {}".format(list_of_frequency))
    return list_of_frequency


def create_mock_reads_directory(mock_reads_dir):
    if os.path.exists(mock_reads_dir):
        shutil.rmtree(mock_reads_dir)
    os.mkdir(mock_reads_dir)


def create_mock_reads(fasta_directory="/home/vered/EMIRGE/EMIRGE-data/",
                      mock_reads_dir='/home/vered/EMIRGE/EMIRGE-data/mock_reads/keep_original_1_change_each_region/',
                      base_fastq="/home/vered/EMIRGE/EMIRGE-data/RDB53_CATTGACG_L007_R1_001.fastq",
                      total_number_of_bacterias=15,
                      min_frequency=0.5,
                      keep_original=True,
                      sort_frequency=False,
                      bacterium_to_change=1,
                      regions_to_change=1,
                      bases_to_change=1,
                      total_amout_of_reads = 130000
                      ):
    """
    using the full fasta DB as reference
    :return:
    # """

    create_mock_reads_directory(mock_reads_dir)

    base_read_file_name = 'mockReads_R'
    mock_fastq1 = os.path.join(mock_reads_dir, base_read_file_name + "1.fastq")
    mock_fastq2 = os.path.join(mock_reads_dir, base_read_file_name + "2.fastq")
    mock_fa1 = os.path.join(mock_reads_dir, base_read_file_name + "1.fa")
    mock_fa2 = os.path.join(mock_reads_dir, base_read_file_name + "2.fa")
    expected_data_path = os.path.join(mock_reads_dir, "expected_res.csv")

    # mock_refs = '/home/vered/EMIRGE/EMIRGE-data/mock_fasta/'
    # fasta_dir = mock_refs

    total_num_of_regions = 5
    total_amount_of_reads_per_region = total_amout_of_reads/total_num_of_regions
    read_len = 126

    frequencies = get_frequency_list(total_number_of_bacterias, min_frequency=min_frequency, isSorted=sort_frequency)
    mock_bacterium = create_mock_bacterium(frequencies,
                                           fasta_directory,
                                           total_num_of_regions,
                                           read_len,
                                           total_amount_of_reads_per_region,
                                           bacterium_to_change=bacterium_to_change,
                                           regions_to_change=regions_to_change,
                                           bases_to_change=bases_to_change,
                                           keep_original=keep_original)

    write_bacterium_to_fasta(mock_bacterium, mock_fa1, mock_fa2, min_frequency)
    write_bacterium_to_fastq(mock_bacterium, base_fastq, mock_fastq1, mock_fastq2)
    write_mock_bacterium_to_csv(mock_bacterium, expected_data_path)


def create_mock_refs():
    old_fasta_path = '/home/vered/EMIRGE/EMIRGE-data/'
    new_fasta_path = '/home/vered/EMIRGE/EMIRGE-data/mock_fasta/'
    amount = 5000
    create_mock_fasta(old_fasta_path, new_fasta_path, amount)


def get_pair(seq):
    reversed_end_of_seq = seq[::-1]
    pair_read = "".join([Base.base_dict[b] for b in reversed_end_of_seq])
    return pair_read


def fix_expected_res(expected_path):
    df = pd.DataFrame.from_csv(expected_path, index_col=None)
    priors_sum = sum(df.drop_duplicates('id')['prior'])
    df.prior = df.prior.apply(lambda r : float(r)/priors_sum)
    reads = df.sequence.apply(lambda r: r[:126])


    pair_reads = df.sequence.apply(lambda r: get_pair(r[-126:]))

    df.sequence = reads + pair_reads

    df.to_csv(expected_path)

def merge_fasta_paired_reads_with_reads(fa1_path, fa2_path, output_path):
    with open(fa1_path, 'r') as fa1:
        with open(fa2_path, 'r') as fa2:
            with open(output_path, 'w') as output:
                row_idx = 0
                for (fa1_line, fa2_line) in zip(fa1, fa2):
                    if row_idx%3==1:
                        output.write(fa1_line[:-1] + fa2_line[:-1] + '\n')
                    else:
                        output.write(fa1_line)
                    row_idx += 1


def split_fastq_with_errors(fastq, fastq1_path, fastq2_path, read_len=126):
    with open(fastq, 'r') as fastq:
        with open(fastq1_path, 'w') as fastq1:
            with open(fastq2_path, 'w') as fastq2:
                row_idx = 0
                for line in fastq:
                    if row_idx%4==1 or row_idx%4==3 :
                        fastq1.write(line[:read_len] + '\n')
                        fastq2.write(line[-1*read_len:])
                    else:
                        fastq1.write(line)
                        fastq2.write(line)
                    row_idx += 1


def fix_randomreads_output_undo_reverse_complete(fastq_path, fixed_fastq_path):
    with open(fastq_path, 'r') as fastq:
        with open(fixed_fastq_path, 'w') as fixed_fastq:
            to_fix = False
            line_ix = 0
            for line in fastq:
                new_line = line
                if line_ix%4 == 0:
                    if 'chr1_1' in line:
                        to_fix = True
                    else:
                        to_fix = False
                if to_fix:
                    if line_ix%4 == 1:
                        pair = get_pair(line[:-1]) # without the '\n
                        new_line = pair + '\n'
                    if line_ix%4 == 3:
                        new_line = line[:-1][::-1] + '\n'
                fixed_fastq.write(new_line)
                line_ix += 1

def test_random_reads_without_pair(fixed_fastq_path):

    with open(fixed_fastq_path, 'r') as fixed_fastq:
        line_ix = 0
        for line in fixed_fastq:
            prev_line = ''
            if line_ix%4 == 1:
                if prev_line != '':
                    if not prev_line in line:
                        print 'Not the same'
                prev_line = line
            line_ix += 1

def randomreads_to_emrige_input(fastq, fixed_fastq, fastq1_path, fastq2_path, read_len=126):
    fix_randomreads_output_undo_reverse_complete(fastq, fixed_fastq)
    split_fastq_with_errors(fixed_fastq, fastq1_path, fastq2_path)


def parse_arguments(argv = sys.argv[1:]):
    """
    command line interface to emirge

    """
    parser = OptionParser("create mock data")

    # REQUIRED
    group_reqd = OptionGroup(parser, "Required flags",
                             "These flags are all required to run EMIRGE, and may be supplied in any order.")

    group_reqd.add_option("-f", dest="fasta_directory",
                      type="string",
                      help="path to the directory contain fasta file for each region")
    group_reqd.add_option("-m", dest="mock_reads_directory",
                      type="string",
                      help="path to the directory that will contain the results fastq\fa files")
    group_reqd.add_option("-r", dest="base_fastq_path",
                          type="string",
                          help="path to the fastq file that will serve as reference fastq for the mock fastq")
    group_reqd.add_option("-f", dest="total_number_of_bacterias",
                      type="int",
                      help="total number of bacterias in the mixure")
    parser.add_option_group(group_reqd)
    (options, args) = parser.parse_args(argv)
    return options

if __name__ == '__main__':
    # fix_expected_res('~/EMIRGE/EMIRGE-data/mock_reads/mock_reads_one_change_3_in_each_region_low_frequency/expected_res.csv')
    # options = parse_arguments()
    # create_mock_reads(options.fasta_directory, options.mock_reads_directory, options.base_fastq_path, options.total_number_of_bacterias)
    create_mock_reads()
    # create_mock_refs()
   # create_mock_reads_from_mock_refs()
   #  fix_randomreads_output_undo_reverse_complete('/home/vered/EMIRGE/EMIRGE-data/mock_read_errors_test/random_reads.fastq', '/home/vered/EMIRGE/EMIRGE-data/mock_read_errors_test/random_reads_fixed.fastq')
   #  test_random_reads_without_pair('/home/vered/EMIRGE/EMIRGE-data/mock_read_errors_test/random_reads_fixed.fastq')
   # merge_fasta_paired_reads_with_reads('/home/vered/EMIRGE/EMIRGE-data/mock_read_errors_test/mockReads_R1.fa',
   #                                     '/home/vered/EMIRGE/EMIRGE-data/mock_read_errors_test/mockReads_R2.fa' ,
   #                                     '/home/vered/EMIRGE/EMIRGE-data/mock_read_errors_test/randomreads_low.fa')

    # randomreads_to_emrige_input('/home/vered/EMIRGE/EMIRGE-data/mock_read_errors_test/randomreads_low.fastq',
    #                             '/home/vered/EMIRGE/EMIRGE-data/mock_read_errors_test/randomreads_low_fixed.fastq',
    #                         '/home/vered/EMIRGE/EMIRGE-data/mock_read_errors_test/randomreads_low_R1.fastq',
    #                         '/home/vered/EMIRGE/EMIRGE-data/mock_read_errors_test/randomreads_low_R2.fastq')