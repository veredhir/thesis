#!/usr/bin/env python2
from optparse import OptionParser, OptionGroup
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

DEFAULT_SIMILARITY_RATE = 0.8


def get_cmd_arguments_parser():
    USAGE = \
        """usage: %prog --fa <reference_fasta_path> [options]

        Art Illumina fastq output contains the sequences and the reverse compliment sequences.
        This script fix the output using the reference fasta file.
        The output is fastq file without reverse compliment sequences.


        usage examples:

        1. fixing fasta files - each sequence will have unique id:
           %prog --fa <reference_fasta_path> --fix

        2. fix fastq files:
           %prog --fa <reference_fasta_path> --fq1 <fastq1_path> --fq2 <fastq1_path>

        """

    parser = OptionParser(USAGE)

    group_required = OptionGroup(parser, "Required parameters")

    group_required.add_option("--fa", dest="reference_fasta_path",
                          type="string", default=None,
                          help="reference fasta path")
    group_required.add_option("--fq1", dest="fastq1_path",
                             type="string", default=None,
                          help="first fastq file to fix")
    group_required.add_option("--fq2", dest="fastq2_path",
                              type="string", default=None,
                              help="second fastq file to fix")
    group_required.add_option("--fix", dest="fix_fasta_id",
                              action="store_true", default=False,
                              help="change the fasta file, so each sequence will have unique id. "
                                  "one should have to run art_illumina on the fixed fasta file.")

    parser.add_option_group(group_required)
    return parser


def fix_records(fq1_record, fq2_record, fa_record, min_similarity_rate=DEFAULT_SIMILARITY_RATE):
    counter_fq1 = 0
    for base_fa, base_fq1 in zip(str(fa_record.seq), str(fq1_record.seq)):
        if base_fa == base_fq1:
            counter_fq1 += 1

    similarity_rate1 = float(counter_fq1)/len(fq1_record.seq)
    if similarity_rate1 > min_similarity_rate:
        return fq1_record, fq2_record
    else:
        return fq2_record, fq1_record


def validate_records(fq1_record, fa_record, min_similarity_rate=DEFAULT_SIMILARITY_RATE):
    counter_fq1 = 0
    for base_fa, base_fq1 in zip(str(fa_record.seq), str(fq1_record.seq)):
        if base_fa == base_fq1:
            counter_fq1 += 1

    similarity_rate1 = float(counter_fq1)/len(fq1_record.seq)
    if similarity_rate1 > min_similarity_rate:
        return True
    else:
        return False


def fix_fastq_using_fasta(ref_fasta_path, fastq1_path, fastq2_path,):
    fixed_fq1_records = []
    fixed_fq2_records = []
    fa_generator = SeqIO.parse(ref_fasta_path, "fasta")
    fa_record = fa_generator.next()

    counter=0
    for fq1_record, fq2_record in zip(SeqIO.parse(fastq1_path, "fastq"), SeqIO.parse(fastq2_path, "fastq")):
        while not fa_record.id == fq1_record.id.partition('-')[0]:
            # print "fa = {}, fq1 = {}".format(fa_record.id , fq1_record.id.partition('-')[0])
            fa_record = fa_generator.next()
        fixed_records = fix_records(fq1_record, fq2_record, fa_record)
        fixed_fq1_records.append(fixed_records[0])
        fixed_fq2_records.append(fixed_records[1])
        counter+=1

    print "----------------\nReads length = {}\n".format(counter)


    SeqIO.write(fixed_fq1_records, fastq1_path, "fastq")
    SeqIO.write(fixed_fq2_records, fastq2_path, "fastq")


def validte_output(ref_fasta_path, fastq1_path, max_mismatch_allowed=0.05):
    fa_generator = SeqIO.parse(ref_fasta_path, "fasta")
    fa_record = fa_generator.next()

    fastq_len = len(list(SeqIO.parse(fastq1_path, "fastq")))
    fasta_len = len(list(SeqIO.parse(ref_fasta_path, "fasta")))
    print("fasta len = {} fastq len = {}".format(fasta_len, fastq_len))

    mismatch_ctr = 0
    for fq1_record, fq2_record in zip(SeqIO.parse(fastq1_path, "fastq"), SeqIO.parse(fastq1_path, "fastq")):
        while not fa_record.id == fq1_record.id.partition('-')[0]:
            fa_record = fa_generator.next()
        if not validate_records(fq1_record, fa_record):
            mismatch_ctr += 1
    mismatch = float(mismatch_ctr)/fastq_len
    if mismatch > max_mismatch_allowed:
        print("FAILED FIXING ART ILLUMINA OUTPUT [{}/{}, Mismach allowed = {}]".format(mismatch_ctr, fastq_len, max_mismatch_allowed))
    else:
        print("FIX art Illumina [{} <= {}]".format(mismatch, max_mismatch_allowed))


def update_fasta_ids(fasta_path):
    prev_record = None
    fasta_records=[]
    ctr=0
    for fa_record in SeqIO.parse(fasta_path, "fasta"):
        if prev_record:
            if prev_record.id == fa_record.id:
                new_id = str(fa_record.id) + str(ctr)
                fa_record.id = new_id
                fa_record.description = new_id
                ctr += 1
            else:
                prev_record=fa_record
                ctr = 0
        else:
            prev_record = fa_record
        fasta_records.append(fa_record)
    SeqIO.write(fasta_records, fasta_path, "fasta")
    print("done update_fasta_ids!")


def main(argv = sys.argv[1:]):
    parser = get_cmd_arguments_parser()
    (options, args) = parser.parse_args(argv)
    if not options.reference_fasta_path:
        parser.error("missing required arguments")

    if options.fix_fasta_id:
        update_fasta_ids(options.reference_fasta_path)
    else:
        validte_output(options.reference_fasta_path, options.fastq1_path)
        fix_fastq_using_fasta(options.reference_fasta_path, options.fastq1_path, options.fastq2_path)
        validte_output(options.reference_fasta_path, options.fastq1_path)

if __name__ == "__main__":
    main()