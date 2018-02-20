import sys
from optparse import OptionParser, OptionGroup, SUPPRESS_HELP

import random
from Bio import SeqIO
import os
import glob
import pandas as pd
from random import randint
import shutil

from emirge_utills import *


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


def get_fasta_files(fasta_dir, total_num_of_regions):
    fasta_files = []
    os.chdir(fasta_dir)
    files = glob.glob("*.fasta")
    for file in files:
        fasta_file = FastaFile(fasta_dir, file, total_num_of_regions)
        fasta_file.initialize()
        fasta_files.append(fasta_file)
    return sorted(fasta_files)


def create_original_bacterium(bacterias_ids,
                              fasta_dir,
                              total_num_of_regions,
                              read_length):
    fasta_files = get_fasta_files(fasta_dir, total_num_of_regions)
    bacterium =[]
    for i in range(len(bacterias_ids)):
        bacteria = MockBacterium(bacterias_ids[i], 0, read_length, 1, 0, 0)
        bacteria.initialize(fasta_files)
        bacterium.append(bacteria)

    print "Done create_mock_bacterium"
    return bacterium



def get_mock_df(mock_csv_path):
    df = pd.DataFrame.from_csv(mock_csv_path, index_col=None)
    df['id'] = df['id'].apply(lambda id: id.split("#")[0])
    return df


def compare_mock_to_original(fasta_directory,
                             mock_csv_path):
    """
    using the full fasta DB as reference,
    :return:
    # """

    mock_df = get_mock_df(mock_csv_path)

    total_num_of_regions = 5
    read_len = 126
    mock_bacterium = create_original_bacterium(mock_df['id'].unique().tolist(),
                                               fasta_directory,
                                               total_num_of_regions,
                                               read_len)
    for mock_bacteria in mock_bacterium:
        curr_bacteria_df = mock_df[mock_df.id == mock_bacteria.id]
        if len(curr_bacteria_df) > 5:
            print "test {}".format(mock_bacteria.id)
            curr_bacteria_df = curr_bacteria_df.sort_values(by=['prior'], ascending=False)
            curr_bacteria_df = curr_bacteria_df.drop_duplicates('region')
        if(len(curr_bacteria_df) == 0):
            print("ERROR in id = {}".format(mock_bacteria.id))
            continue
        mock_seqs = curr_bacteria_df['sequence'].tolist()
        for mock_seq, org_seq in zip(mock_seqs, mock_bacteria.sequences):
            changes = [ "id = {} -  {}!={}".format(ix, base_mock, base_org)
                        for base_mock, base_org, ix in zip(mock_seq, org_seq, range(len(org_seq)))
                        if base_mock!= base_org]
            if len(changes) != 0:
                print "Bacteria id = {}".format(mock_bacteria.id)
                print changes

def get_command_line_arguments_parser():
    """
    command line interface to emirge

    """
    USAGE = \
        """usage: %prog FASTA_DIR MOCK_CSV_PATH [options]

        Compare mock bacterium to the original becterium in the Fasta DB
        """

    parser = OptionParser(USAGE)
    return parser


def main(argv = sys.argv[1:]):
    define_logger(logging.DEBUG)
    parser = get_command_line_arguments_parser()

    # ACTUALLY PARSE ARGS
    (options, args) = parser.parse_args(argv)

    # minimal sanity checking of input
    if len(args) != 2:
        parser.error(
            "FASTA_DIR and MOCK_CSV_PATH are required, and all other options should have a flag associated with them (options without flags: %s)" % args)
        return

    fasta_dir = os.path.abspath(args[0])
    mock_csv_path = os.path.abspath(args[1])

    compare_mock_to_original(fasta_dir,
                             mock_csv_path)


if __name__ == "__main__":
    main()
