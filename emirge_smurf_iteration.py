import glob
import os
import os.path
from __builtin__ import staticmethod
from ast import literal_eval
from gmpy2 import popcount
from itertools import izip
import math

import numpy as np
import pandas as pd
from Bio import SeqIO


from emirge_headers import *
from emirge_const import quals
from emirge_utills import *
from seqBIn import *


class Primer(object):

    def __init__(self):
        self.all = {}
        self.len = 0
        self.histogram = {}

    def add_primer(self, region, primers_list):
        self.all.update({region: primers_list})
        for primer in primers_list:
            self.histogram.update({primer: 0})
        self.len += 1

    def print_histogram(self):
        logging.debug("\n\n=============================")
        logging.debug("PRIMERS HISTOGRAM:")
        for primer, counter in self.histogram.iteritems():
            logging.debug("{} --> {}".format(primer, counter))
        logging.debug("\n=============================")
        logging.debug("\n")

    def update_primers_using_histogram(self, limit=20):
        self.print_histogram()
        for i, primers in self.all.iteritems():
            new_primers = []
            for primer in primers:
                if self.histogram[primer] > limit:
                    new_primers.append(primer)
            logging.info("Region {}: new primers = {}, old primers = {}".format(i, len(new_primers), len(primers)))
            self.all[i] = new_primers


class EmirgeThresholds():
    """hold the thresholds for emirge computation"""
    def __init__(self,
                 max_changed_bases_rate_for_split,
                 min_coverage_for_split,
                 min_minor_prob_for_split,
                 min_similar_bases_rate_for_merge,
                 max_priors_diff_for_stability_test):
        self.max_changed_bases_rate_for_split = max_changed_bases_rate_for_split
        self.min_coverage_for_split = min_coverage_for_split
        self.min_minor_prob_for_split = min_minor_prob_for_split
        self.min_similar_bases_rate_for_merge = min_similar_bases_rate_for_merge
        self.max_priors_diff_for_stability_test = max_priors_diff_for_stability_test


class EmirgePaths():
    def __init__(self, working_dir):
        self.reference = os.path.join(working_dir, "reference_db.csv")
        self.full_reference_name = "full_reference_db.csv"
        self.full_reference = ""
        self.current_state = os.path.join(working_dir, "curr_state.csv")
        self.final_results = os.path.join(working_dir, "final_results.csv")
        self.mapping = os.path.join(working_dir, "mapping.csv")
        self.unique_ref_to_ref = os.path.join(working_dir, "unique_ref_id_to_ref_id.csv")
        self.posteriors = os.path.join(working_dir, "posteriors.csv")
        self.read_quals = os.path.join(working_dir, "reads_db.csv")


class EmirgeIteration(object):

    def __init__(self, working_dir,
                 prev_emirge_iteration=None,
                 fastq_path=None,
                 reversed_fastq_path=None,
                 primers_path=None,
                 fasta_path=None,
                 read_len=None,
                 number_of_regions=5,
                 update_weight_using_the_reads=False,
                 # thresholds:
                 max_changed_bases_rate_for_split=0.0025,  # (2.5%)
                 min_coverage_for_split=0.005,  # (0.5%)
                 min_minor_prob_for_split=0.08,
                 min_similar_bases_rate_for_merge=0.999,  # (98%)
                 max_priors_diff_for_stability_test=0.005,
                 allow_split=False,
                 debug_mode=True):  # (0.5%)
        """
        :param prev_emirge_iteration:
        :param working_dir:
        :param fastq_path: path to the fastq file
        :param reversed_fastq_path: path to the reversed fastq file
        :param primers_path: path to csv file containing table with the regions as header,
                             each column contains the primers for the specific region.
        :param fasta_path: path to directory containing the fasta files for each region
        :param read_len: read length
        :param reference_path: contain the reference path in the new binary format. if not exists - None
        :param max_changed_bases_rate_for_split: (# of changed bases in the split reference)/( read length)
        :param min_coverage_for_split:  (% of reads mapped to the split reference - estimated)
        :param min_minor_prob_for_split:  (the second best prob N in index i, if higher then threshold - will be candidate for split)
        :param min_similar_bases_rate_for_merge: (# of similar bases in 2 candidates) / (read length)
        :param max_priors_diff_for_stability_test: % priors diff between previous and current iteration
        :param number_of_regions:
        :param update_weight_using_the_reads: update the reference db with the number
               of the region that were amplified according to the reads.
        """

        # try to use weight calculate after getting the relevant references using the actual primers:
        self.number_of_regions = number_of_regions
        self.primers = Primer()
        self.paths = EmirgePaths(working_dir)
        self.th = EmirgeThresholds(max_changed_bases_rate_for_split,
                                   min_coverage_for_split,
                                   min_minor_prob_for_split,
                                   min_similar_bases_rate_for_merge,
                                   max_priors_diff_for_stability_test)
        self.ref_format = ReferenceFormat()
        self.reads_full_data_format = ReadsFullDataFormat()
        self.prev_priors_for_stability_test = None

        # First EMIRGE iteration:
        if prev_emirge_iteration is None:
            self.init_pre_process(primers_path,
                                  fastq_path,
                                  reversed_fastq_path,
                                  read_len,
                                  fasta_path,
                                  update_weight_using_the_reads,
                                  allow_split,
                                  debug_mode)
        else:
            self.init(prev_emirge_iteration)

    @time_it
    def init_pre_process(self,
                         primers_path,
                         fastq_path,
                         reversed_fastq_path,
                         read_len,
                         fasta_path,
                         update_weight_using_the_reads,
                         allow_split,
                         debug_mode):
        self.paths.full_reference = os.path.join(fasta_path, self.paths.full_reference_name)
        self.debug_mode = debug_mode
        self.iteration_index = 0
        logging.info("\n\n ITERATION {}".format(self.iteration_index))
        logging.info("number of regions = {}".format(self.number_of_regions))
        self.read_len = read_len
        self.allow_split = allow_split
        self.initialize_primers(primers_path)
        reads_df = self.prepare_reads(fastq_path, reversed_fastq_path)
        self.primers.update_primers_using_histogram()
        if not self.is_processed_reference_file_exists():
            self.prepare_references(fasta_path)
            self.get_unique_amplified_references()
        else:
            full_ref_df = pd.DataFrame.from_csv(self.paths.full_reference, index_col=None)
            full_ref_df.to_csv(self.paths.reference, index=False)

        self._find_mapping(reads_df)
        self._find_initial_weight(update_weight_using_the_reads)
        self.calc_initial_priors()


    @time_it
    def init(self, prev_emirge_iteration):
        """
        Initialize class using the previous iteration.
        :param prev_emirge_iteration: type EmirgeIteration
        """
        self.iteration_index = prev_emirge_iteration.iteration_index + 1
        logging.info("\n\n ITERATION {}".format(self.iteration_index))

        self.debug_mode = prev_emirge_iteration.debug_mode
        self.number_of_regions = prev_emirge_iteration.number_of_regions
        self.allow_split = prev_emirge_iteration.allow_split
        self.read_len = prev_emirge_iteration.read_len
        self.paths.read_quals = prev_emirge_iteration.paths.read_quals
        self.th = prev_emirge_iteration.th
        curr_state_df = pd.DataFrame.from_csv(prev_emirge_iteration.paths.current_state, index_col=None)
        curr_state_df = curr_state_df[[CurrentStateFormat.Reference_id,
                                       CurrentStateFormat.Region,
                                       CurrentStateFormat.Weight,
                                       CurrentStateFormat.Priors] +
                                       CurrentStateFormat.Bases.all]

        curr_state_df.to_csv(self.paths.current_state, index=False)
        self.prev_priors_for_stability_test = curr_state_df[[CurrentStateFormat.Reference_id,
                                                             CurrentStateFormat.Priors]].drop_duplicates()

        reference_df = curr_state_df[
            curr_state_df[CurrentStateFormat.Priors] > 0.1 * self.th.max_priors_diff_for_stability_test]
        logging.info("ref columns = {}".format(reference_df.columns))
        reference_df = reference_df[[CurrentStateFormat.Reference_id,
                                      CurrentStateFormat.Region,
                                      CurrentStateFormat.Weight] +
                                     Base.all].drop_duplicates()
        reference_df.to_csv(self.paths.reference, index=False)
        posteriors_df = pd.DataFrame.from_csv(prev_emirge_iteration.paths.posteriors, index_col=None)
        posteriors_df.to_csv(self.paths.posteriors, index=False)
        reads_df = pd.DataFrame.from_csv(prev_emirge_iteration.paths.mapping,
                                         index_col=None)[[  MappingForamt.Region,
                                                            MappingForamt.Group_id,
                                                            MappingForamt.Count] +
                                                            Base.all].drop_duplicates()
        self._find_mapping(reads_df)


    def is_processed_reference_file_exists(self):
        if os.path.exists(self.paths.full_reference):
            return True
        else:
            return False


    @time_it
    def _get_ref_for_region(self, fasta_path, read_length, region_ix):
        logging.info("Start processing fasta, path = %s", fasta_path)

        data_dicts = []
        for record in SeqIO.parse(fasta_path, "fasta"):
            # logging.debug("Start processing fasta, path = %s", fasta_path)
            sequence = record.seq.__str__()
            title = record.id.__str__()
            record_dict = self.ref_format.get_ref_dict(title, region_ix, sequence_to_bin(
                sequence[0:read_length] + sequence[-1 * read_length:]))
            data_dicts.append(record_dict)

        logging.debug("Created dict for region = {}".format(region_ix))
        ref_df = pd.DataFrame(data_dicts, columns=ReferenceFormat.temp_all)

        return ref_df

    @staticmethod
    def _get_region_from_fasta_name(file_name, size):
        for i in range(1, size+1):
            if str(i) in file_name:
                return i
        raise Exception("Invalid file name: {}, doesn't contain ragion id".format(file_name))

    def _get_region_using_primer(self, read, minimum_match_score=0.8):
        match_size = 15
        read_primer = read[:match_size]
        best_score = minimum_match_score
        region = None

        for i, primers in self.primers.all.iteritems():
            score = max(similar(read_primer, primer[:match_size]) for primer in primers)
            if score >= best_score:
                best_score = score
                region = i

        # update the primers histogram from the reads:
        if region is not None:
            for primer in self.primers.all[region]:
                if similar(read_primer, primer[:match_size]) == best_score:
                    self.primers.histogram[primer] += 1

        return region

    @time_it
    def _get_mapped_references(self):
        mapping_df = pd.DataFrame.from_csv(self.paths.mapping, index_col=None)
        full_ref_df = pd.DataFrame.from_csv(self.paths.reference, index_col=None)

        # Merging - taking only mapped references in mapped regions
        # e.g after conversation with Noam and Gary at 14-03-2017, do not assume all regions were sampled correctly.
        mapped_reference = mapping_df[[MappingForamt.Ref_id, MappingForamt.Region]].drop_duplicates()
        curr_ref_df = mapped_reference.merge(full_ref_df,
                                             on=[HeadersFormat.Unique_Ref_id, HeadersFormat.Region],
                                             how='left')
        logging.info("Mapped references: {}/{}".format(
            len(curr_ref_df.drop_duplicates(CurrentStateFormat.Reference_id)), len(full_ref_df.drop_duplicates(CurrentStateFormat.Reference_id))))

        return curr_ref_df

    @time_it
    def calc_initial_priors(self):
        """
        calc prior for each reference, e.g: P(Si) for each reference Si.
        P(Si) = (sum_j(P(rj|Si)/Wi))/ sum_k(sum_j(P(rj|Sk)/Wk)
            where P(rj|Si) is equal to 1 if the read j mapped to Si and zero otherwise.
            if, for example, rj is mapped with the same probability to Si and Sj then, P(ri|Si) = P(ri|Sk) = 0.5.
        """
        mapping_df = pd.DataFrame.from_csv(self.paths.mapping, index_col=None)

        # get all the mapped reference and all their data:
        curr_state_without_priors = self._get_mapped_references()

        # x reads in the group, the group mapped to y refs -> the mapp weight is 1/y, the mapped reads counter is x/y
        mapping_df['mapped_reads_counter'] = mapping_df[MappingForamt.Count]*mapping_df[MappingForamt.Mapp_weight]
        # sum how many weighted reads mapped to each reference.
        ref_with_counter_df = mapping_df.groupby(MappingForamt.Ref_id)['mapped_reads_counter'].sum().reset_index()

        # add the reference weight to each reference
        curr_state_with_priors = ref_with_counter_df.merge(curr_state_without_priors, on=HeadersFormat.Unique_Ref_id, how="left")

        curr_state_with_priors[CurrentStateFormat.Priors] = curr_state_with_priors.apply(
            lambda row: float(row['mapped_reads_counter'])/ float(row[CurrentStateFormat.Weight]), axis=1)

        # sum over the priors (on the reference without the regions duplication.)
        priors_weight = curr_state_with_priors[[CurrentStateFormat.Reference_id, CurrentStateFormat.Priors]].drop_duplicates()[[CurrentStateFormat.Priors]].sum(0)

        # normalize by the priors overall weight.
        curr_state_with_priors[CurrentStateFormat.Priors] = curr_state_with_priors.apply(
            lambda row: float(row[CurrentStateFormat.Priors]) / priors_weight, axis=1)
        curr_state_columns = [CurrentStateFormat.Reference_id, CurrentStateFormat.Region, CurrentStateFormat.Weight, CurrentStateFormat.Priors] + CurrentStateFormat.Bases.all
        curr_state_with_priors = curr_state_with_priors[curr_state_columns]
        curr_state_with_priors.to_csv(self.paths.current_state)



    @time_it
    def _find_initial_weight(self, use_weight_using_reads):
        """
        In the first iteration, get only reference which mapped to all existing regions of the reference.
        :param reads_df:
        :return:
        """
        full_ref_df = pd.DataFrame.from_csv(self.paths.reference, index_col=None)
        mapping_df = pd.DataFrame.from_csv(self.paths.mapping, index_col=None)

        all_mapped_refs = mapping_df[[MappingForamt.Ref_id, MappingForamt.Region]].drop_duplicates()
        if use_weight_using_reads:
            ref_with_weight = all_mapped_refs.groupby([MappingForamt.Ref_id]).count().reset_index()
        else:
            ref_with_weight = full_ref_df.groupby([MappingForamt.Ref_id]).count().reset_index()
        ref_with_weight.rename(columns={MappingForamt.Region: CurrentStateFormat.Weight}, inplace=True)
        ref_with_weight = ref_with_weight.drop_duplicates()[[CurrentStateFormat.Weight, CurrentStateFormat.Reference_id]]
        refs_without_weight = full_ref_df[ReferenceFormat.Bases.all + [ReferenceFormat.Ref_Id,
                                                                       ReferenceFormat.Region]]
        # get the actual weight for each mapped reference
        ref_with_weight = pd.merge(refs_without_weight,
                                   ref_with_weight,
                                   on=ReferenceFormat.Ref_Id,
                                   how='right')
        ref_with_weight.to_csv(self.paths.reference, index=False)


    @time_it
    def _find_mapping(self, unique_reads_df, full_ref_df=None):
        if full_ref_df is None:
            full_ref_df = pd.DataFrame.from_csv(self.paths.reference, index_col=None)
        ref_df_copy = full_ref_df.copy()

        for base in Base.all:
            full_ref_df[base].update(full_ref_df[base].apply(lambda r: int(r)))
            unique_reads_df[base] = unique_reads_df[base].apply(lambda r: int(r))

        unique_reads_for_region_groups = unique_reads_df.groupby(HeadersFormat.Region)
        references_grouped_by_regions = full_ref_df.groupby(HeadersFormat.Region)

        mapping_dfs = []
        for region, reads_df in unique_reads_for_region_groups:
            reads_df.reset_index(inplace=True)
            try:
                ref_df = references_grouped_by_regions.get_group(region)
            except:
                logging.info("Skip region {}, no reference match to the region.".format(region))
                continue
            ref_df.reset_index(inplace=True)
            # Find the number of matches between each read - ref couple:
            dfs = []
            chunk_size = 1000
            logging.info("Mapping for region = {}, reads size = {}, amount of chunks = {}".format(region, len(reads_df), len(ref_df)/chunk_size))
            chunks = int(len(ref_df)/chunk_size) + 1
            for chunk in np.array_split(ref_df, chunks):
                reads_and_refs_chunk_df = pd.DataFrame.merge(chunk, reads_df, on=HeadersFormat.Region)
                a = np.bitwise_and(reads_and_refs_chunk_df['A_x'], reads_and_refs_chunk_df['A_y'])
                c = np.bitwise_and(reads_and_refs_chunk_df['C_x'], reads_and_refs_chunk_df['C_y'])
                g = np.bitwise_and(reads_and_refs_chunk_df['G_x'], reads_and_refs_chunk_df['G_y'])
                t = np.bitwise_and(reads_and_refs_chunk_df['T_x'], reads_and_refs_chunk_df['T_y'])
                reads_and_refs_chunk_df['Score'] = np.bitwise_or(np.bitwise_or(a, c), np.bitwise_or(g, t))
                reads_and_refs_chunk_df['Score'] = reads_and_refs_chunk_df['Score'].apply(lambda r: popcount(r))

                reads_and_refs_chunk_df = reads_and_refs_chunk_df[reads_and_refs_chunk_df['Score'] >=
                                                                  (reads_and_refs_chunk_df.groupby(HeadersFormat.Group_id)['Score'].transform(max)- 2)]
                dfs.append(reads_and_refs_chunk_df)

            logging.info("Done scoring chunks")
            reads_and_refs_df = pd.concat(dfs, ignore_index=True)
            max_scores_series = reads_and_refs_df.groupby(HeadersFormat.Group_id)['Score'].transform(max)
            reads_and_refs_df = reads_and_refs_df[(reads_and_refs_df['Score'] == max_scores_series )]

            rename_base_dict = {}
            for base in Base.all:
                rename_base_dict.update({base + '_y': base}) # keep only the reads bases

            reads_and_refs_df.rename(columns=rename_base_dict, inplace=True)
            reads_and_refs_df[MappingForamt.Mapp_weight] = 0
            reads_and_refs_df = reads_and_refs_df[MappingForamt.full_header]
            mapping_dfs.append(reads_and_refs_df)

            logging.info("Found mapping for region = {}".format(region))

        mapping_df = pd.concat(mapping_dfs)
        mapping_df = self.merge_similar_mapped_reference(mapping_df.copy(), ref_df_copy)

        # Find to how many refs each read group was mapped(approximate)
        mapping_df[MappingForamt.Mapp_weight] = 1 / mapping_df.groupby(HeadersFormat.Group_id)[
            MappingForamt.Count].transform('count').apply(float)

        if abs((sum(unique_reads_df.Count) - int(sum(mapping_df.Count*mapping_df.Mapp_weight)))) > 1:
            raise Exception("MAPPING ERROR: the count is wrong! unique_reads = {}, counter = {}".format(sum(unique_reads_df.Count), sum(mapping_df.Count*mapping_df.Mapp_weight)))

        mapping_df.to_csv(self.paths.mapping, index=False)
        return mapping_df

    @time_it
    def merge_similar_mapped_reference(self, mapping_df, references):
        mapped_references = mapping_df[[MappingForamt.Ref_id, MappingForamt.Region]].merge(references, on=[MappingForamt.Ref_id, MappingForamt.Region], how='left')
        ref_count_before = len(mapped_references[MappingForamt.Ref_id].drop_duplicates())
        mapped_references[ReferenceFormat.Original_Id] = mapped_references[CurrentStateFormat.Reference_id]

        UNIQUE_SEQUENCE = 'unique_sequence'
        UNIQUE_SEQUENCES_GROUP = 'unique_sequence_group'
        UNIQUE_SEQUENCES_GROUP_ID = 'unique_sequences_group_id'
        mapped_references[UNIQUE_SEQUENCE] = mapped_references.groupby([ReferenceFormat.Region] + Base.all).grouper.group_info[0]

        # find for each reference the groups ids of his regions
        # example -> ref1reg1 -> groupA, ref1reg2 -> groupB --> ref1 --> groupA, groupB
        ref_to_unique_groups = pd.DataFrame({UNIQUE_SEQUENCES_GROUP: mapped_references.groupby(ReferenceFormat.Original_Id)[
            UNIQUE_SEQUENCE].apply(list)}).reset_index()
        ref_to_unique_groups[UNIQUE_SEQUENCES_GROUP] = ref_to_unique_groups[UNIQUE_SEQUENCES_GROUP].astype(str)

        # find unique references (should be the same in all regions
        ref_to_unique_groups[UNIQUE_SEQUENCES_GROUP_ID] = \
            ref_to_unique_groups.groupby(UNIQUE_SEQUENCES_GROUP).grouper.group_info[0]

        ref_to_unique_groups = ref_to_unique_groups[[UNIQUE_SEQUENCES_GROUP_ID, ReferenceFormat.Original_Id]]
        if self.iteration_index == 0:
            cols = CurrentStateFormat.Bases.all + [ReferenceFormat.Region, ReferenceFormat.Original_Id]
        else:
            cols = CurrentStateFormat.Bases.all + [ReferenceFormat.Region, ReferenceFormat.Original_Id, CurrentStateFormat.Weight]


        new_refs = mapped_references.merge(ref_to_unique_groups,
                                    on=ReferenceFormat.Original_Id,
                                    how='right')
        # map all the original ids to unique one, original a -> unique 1, original b -> unique 1
        original_ref_id_to_unique_id = new_refs[[ReferenceFormat.Original_Id, UNIQUE_SEQUENCES_GROUP_ID]]
        new_refs = new_refs.drop_duplicates([UNIQUE_SEQUENCES_GROUP_ID, ReferenceFormat.Region])
        ref_id_to_unique_id = new_refs[[ReferenceFormat.Original_Id, UNIQUE_SEQUENCES_GROUP_ID]]
        new_refs = new_refs[cols]
        new_refs.rename(columns={ReferenceFormat.Original_Id: ReferenceFormat.Ref_Id}, inplace=True)
        new_refs_ids = new_refs[[ReferenceFormat.Ref_Id]].drop_duplicates()
        new_refs = new_refs_ids.merge(references, on=ReferenceFormat.Ref_Id, how='left')
        ref_count_after = len(new_refs[MappingForamt.Ref_id].drop_duplicates())
        logging.info("UNIQUE REFERENCES: before = {}, after = {}".format(ref_count_before, ref_count_after))
        new_refs.to_csv(self.paths.reference, index=False)

        # only the unique ids will appear in this df: original a -> unique 1
        ref_id_to_unique_id.rename(columns={ReferenceFormat.Original_Id:ReferenceFormat.Ref_Id}, inplace=True)
        ref_id_to_original_id = ref_id_to_unique_id.merge(original_ref_id_to_unique_id, on=UNIQUE_SEQUENCES_GROUP_ID)

        # all the original will be in this list: original 1 -> original 1, original 2 -> original 1
        ref_id_to_original_id = ref_id_to_original_id[[ReferenceFormat.Original_Id, ReferenceFormat.Ref_Id]]
        ref_id_to_original_id.drop_duplicates(inplace=True)

        mapping_df.rename(columns={ReferenceFormat.Ref_Id: ReferenceFormat.Original_Id}, inplace=True)
        mapping_df = mapping_df.merge(ref_id_to_original_id, on=ReferenceFormat.Original_Id, how='left')
        mapping_df.drop_duplicates([MappingForamt.Region, MappingForamt.Ref_id, MappingForamt.Group_id], inplace=True)
        return mapping_df[MappingForamt.full_header]

    @time_it
    def prepare_references(self, fasta_path):
        os.chdir(fasta_path)

        df = pd.DataFrame(columns=ReferenceFormat.temp_all)
        df.to_csv(self.paths.full_reference, index=False, header=True)

        files = glob.glob("*.fasta")
        for file in files:
            region_ix = self._get_region_from_fasta_name(file, len(files))
            ref_df = self._get_ref_for_region(file, self.read_len, region_ix)
            ref_df.to_csv(self.paths.full_reference, index=False, header=False, mode='a')


    @time_it
    def get_unique_amplified_references(self):
        """
        Find the amplified references using the primers found from the reads.
        Work only on the first iteration
        """
        reference_df = pd.DataFrame.from_csv(self.paths.full_reference, index_col=None)

        for base in Base.all:
            # get only the primers bits
            logging.info("Base: {}".format(base))
            reference_df[base + '_primer'] = reference_df[base].apply(lambda r: int(r) >> (2*self.read_len - 15))

        reference_for_region_groups = reference_df.groupby(HeadersFormat.Region)

        amplified_references_for_regions = {}
        for region, ref_df_in_region in reference_for_region_groups:
            ref_df = ref_df_in_region.reset_index().copy()
            primers = self.primers.all[region]
            high_score_refs_for_prior = []
            for primer in primers:
                bin_primer = sequence_to_bin(primer[:15])
                a = np.bitwise_and(ref_df[Base.A + '_primer'], int(bin_primer[0]))
                c = np.bitwise_and(ref_df[Base.C + '_primer'], int(bin_primer[1]))
                g = np.bitwise_and(ref_df[Base.G + '_primer'], int(bin_primer[2]))
                t = np.bitwise_and(ref_df[Base.T + '_primer'], int(bin_primer[3]))
                # ref_df.insert(len(ref_df.columns), 'Score', np.bitwise_or(np.bitwise_or(a, c), np.bitwise_or(g, t)))
                score_as_bin = np.bitwise_or(np.bitwise_or(a, c), np.bitwise_or(g, t))
                ref_df.loc[:, 'Score'] = score_as_bin.apply(lambda r: popcount(r))
                high_score_refs_for_prior += ref_df[ref_df['Score'] > 10][ReferenceFormat.Original_Id].tolist()
                logging.info('Found {} refs in region {}'.format(len(ref_df[ref_df['Score'] > 10]), region))
            amplified_references_for_regions.update({region: high_score_refs_for_prior})

        relevant_refs = len(reference_df)*[False]
        for region, ref in amplified_references_for_regions.iteritems():
            relevant_refs = relevant_refs | \
                            ((reference_df[ReferenceFormat.Original_Id].isin(ref)) &
                             (reference_df[ReferenceFormat.Region] == region))
        reference_df = reference_df[relevant_refs]
        reference_df = reference_df[[ReferenceFormat.Original_Id, ReferenceFormat.Region] + ReferenceFormat.Bases.all]
        # find duplication of regions:
        # example: ref1reg1 -> groupA, ref1reg2 -> groupB, ref2reg1 -> groupC, ref2reg2 -> groupB
        reference_df[ReferenceFormat.Group_id] = reference_df.groupby([ReferenceFormat.Region] + ReferenceFormat.Bases.all).grouper.group_info[0]
        # find for each reference the groups ids of his regions
        # example -> ref1reg1 -> groupA, ref1reg2 -> groupB --> ref1 --> groupA, groupB
        ref_to_unique_groups = pd.DataFrame({'uniques_groups': reference_df.groupby(ReferenceFormat.Original_Id)[
            ReferenceFormat.Group_id].apply(list)}).reset_index()
        ref_to_unique_groups['uniques_groups'] = ref_to_unique_groups['uniques_groups'].astype(str)

        # find unique references (should be the same in all regions
        ref_to_unique_groups[UnqiueRefToRefFormat.Unique_id] = ref_to_unique_groups.groupby('uniques_groups').grouper.group_info[0]

        map_ref_to_unique_ref = ref_to_unique_groups[[UnqiueRefToRefFormat.Unique_id,
                                                      ReferenceFormat.Original_Id]].drop_duplicates()
        map_ref_to_unique_ref.to_csv(self.paths.unique_ref_to_ref, index=False)

        unique_ref_df = reference_df.merge(ref_to_unique_groups,
                                          on=ReferenceFormat.Original_Id,
                                          how='right')
        unique_ref_df = unique_ref_df[ReferenceFormat.header]
        unique_ref_df.drop_duplicates(inplace=True)

        # initialize current run reference file
        unique_ref_df.to_csv(self.paths.reference, index=False)

        # save the reference DB in cvs format, for future use:
        unique_ref_df.to_csv(self.paths.full_reference, index=False)

        logging.info("Found {}/{} unique reference".format(
            len(unique_ref_df.drop_duplicates(ReferenceFormat.Ref_Id)), len(reference_df.drop_duplicates(ReferenceFormat.Original_Id))))
        return unique_ref_df

    @time_it
    def prepare_reads(self, fastq_path, reversed_fastq_path, unique_group_size_fraction=10000):
        logging.info("Start processing: path = %s, reversed path = %s", fastq_path, reversed_fastq_path)
        base_comp_dict = {Base.A: Base.T,
                          Base.G: Base.C,
                          Base.C: Base.G,
                          Base.T: Base.A,
                          Base.N: Base.N}
        data_dicts=[]
        read_index = 0
        all_reads_counter = 0
        for (record, reversed_record) in izip(SeqIO.parse(fastq_path, "fastq"),
                                              SeqIO.parse(reversed_fastq_path, "fastq")):
            read = record.seq.__str__()
            read_quals = record.letter_annotations['phred_quality']
            reversed_read = reversed_record.seq.__str__()[::-1]
            reversed_read_quals = reversed_record.letter_annotations['phred_quality'][::-1]
            quality = read_quals + reversed_read_quals
            second_read = ''.join(map(lambda base: base_comp_dict[base], reversed_read))
            bin = sequence_to_bin(read + second_read)
            region = self._get_region_using_primer(read)
            if region is not None:
                record_dict = self.reads_full_data_format.get_dict(read_index, region, bin, quality)
                data_dicts.append(record_dict)
                read_index += 1
            all_reads_counter += 1
        logging.info("Reads with matched region found = {} out of {} reads".format(read_index, all_reads_counter))
        reads_df = pd.DataFrame(data_dicts, columns=self.reads_full_data_format.all)
        reads_df[ReadsFullDataFormat.Quals] = reads_df[ReadsFullDataFormat.Quals].astype(str)
        if self.debug_mode:
            reads_df.to_csv(self.paths.read_quals, columns=self.reads_full_data_format.all, index=False)

            # reads_df = pd.DataFrame.from_csv(self.paths.read_quals, index_col=None)
        reads_df[ReadsFormat.Group_id] = reads_df.groupby(
            [self.reads_full_data_format.Region] + self.reads_full_data_format.Bases.all).grouper.group_info[0]

        unique_reads_df = pd.DataFrame({MappingForamt.Count: reads_df.groupby(
            [self.reads_full_data_format.Region,
             self.reads_full_data_format.Group_id] +
             self.reads_full_data_format.Bases.all).size()}).reset_index()

        reads_full_size = len(reads_df)
        unique_group_minimum_size = reads_full_size / unique_group_size_fraction
        unique_reads_df = unique_reads_df.loc[unique_reads_df[MappingForamt.Count] > unique_group_minimum_size]
        groups_list = unique_reads_df[MappingForamt.Group_id].tolist()
        reads_df = reads_df[reads_df[ReadsFormat.Group_id].isin(groups_list)]
        reads_df[[ReadsFormat.Id, ReadsFormat.Group_id, ReadsFormat.Quals]].to_csv(self.paths.read_quals, index=False)
        number_of_reads = len(reads_df)
        logging.info("Dropped {}/{} reads (too low frequency)".
                     format(reads_full_size - number_of_reads, reads_full_size))
        return unique_reads_df

    @time_it
    def calc_prob_N(self):
        """
            Pr(N=n) -> Pr(S[i] = N) for each base in each sequence, for each N {A, C, G, T}

            The result is 4 list for each sequence. ProbN_A, ProbN_C, ProbN_G, ProbN_T
            Saved into the CurrentStateTable.

            If read or sequence is new this round (not seen at t-1), then
            there is no Pr(S|R) from previous round,
            so we 1 instead. (same probability for each Pr(S|R)

           get_unique_amplified_references'prob N) If initial iteration, all reads and seqs are new, so all calcs
            for Pr(N=n) use the prior as weighting factor instead of
            previous round's posterior.
        """
        # Get data:
        bases = [Base.A, Base.C, Base.G, Base.T]
        current_state_df, mapping_df, reads_df, posteriors_df = self.get_dfs_for_calc_probN(bases)
        logging.info("length of dataframes: reads = {}, mapping = {}, current_state_df = {}"
                     .format(len(reads_df), len(mapping_df), len(current_state_df)))
        reads_data_columns, ref_data_columns, prob_success_data_columns, prob_failure_data_columns, rename_dict = self.get_columns_for_calc_probN(bases)

        prob_n_full_dict = {Base.A:[], Base.C:[], Base.G:[], Base.T:[]}
        mapping_grouped_by_ref = mapping_df.groupby(CurrentStateFormat.Reference_id)

        i=0
        for ref_group_id, mapping_ref_df in mapping_grouped_by_ref:
            mapping_ref_df = mapping_ref_df.reset_index()
            ref_full_data = mapping_ref_df.merge(current_state_df, on=[CurrentStateFormat.Reference_id, CurrentStateFormat.Region], how='left', suffixes=('_read', '_ref'))
            # logging.info("len of mapping_ref_df.merge(current_state_df) = {}".format(len(ref_full_data)))
            ref_full_data = ref_full_data.merge(reads_df, on=MappingForamt.Group_id, how='left')
            # logging.info("len of merge of (mapping, reads, curr_state) = {}".format(len(ref_full_data)))

            if posteriors_df is not None:
                ref_full_data = ref_full_data.merge(posteriors_df, on = [ReadsFullDataFormat.Id, CurrentStateFormat.Reference_id], how='left')
                # logging.info("len of merge of all df = {}, len(posteriors) = {}".format(len(ref_full_data), len(posteriors_df)))

            prob_n_for_base = {}
            prob_success = ref_full_data[prob_success_data_columns].copy()
            prob_success.rename(columns=rename_dict, inplace=True)
            prob_fail = ref_full_data[prob_failure_data_columns].copy()
            prob_fail.rename(columns=rename_dict, inplace=True)
            mapp_weight = ref_full_data[MappingForamt.Mapp_weight]

            if self.iteration_index == 0:
                posteriors = mapp_weight
            else:
                posteriors = ref_full_data[HeadersFormat.Posterior]

            for base in bases:
                reads = ref_full_data[reads_data_columns[base]].rename(columns=rename_dict)
                # for each base we calculate P(nk = base):
                prob_n_for_base.update({base: (reads.multiply(prob_success) + (1 - reads).multiply(prob_fail))
                                       .multiply(posteriors, axis='index')})

                # refs = ref_full_data[ref_data_columns[base]].rename(columns=rename_dict)

                # prob_n_for_base.update({base: (reads.multiply(refs.multiply(prob_success)) +
                #                                (1 - reads.multiply(refs)).multiply(prob_fail))
                #                                .multiply(posteriors, axis='index')})

                # Calculate P(Nk)for each region of the current reference:
                prob_n_for_base[base][HeadersFormat.Region] = ref_full_data[HeadersFormat.Region]
                prob_n_for_base[base][HeadersFormat.Posterior] = posteriors
                prob_n_for_base[base] = prob_n_for_base[base].groupby(HeadersFormat.Region).apply(sum)

                # Divide each sequence (reference + region) in the sum of the posteriors in the specific sequence.
                prob_n_for_base[base] = prob_n_for_base[base].divide(prob_n_for_base[base][HeadersFormat.Posterior], axis=0)

                # "groupby" turn the region into the index:
                prob_n_for_base[base][HeadersFormat.Region] = prob_n_for_base[base].index
                prob_n_for_base[base][MappingForamt.Ref_id] = ref_group_id
                prob_n_full_dict[base] += [prob_n_for_base[base]]


        for base in bases:
            prob_n_full_dict[base] = pd.concat(prob_n_full_dict[base], ignore_index=True)
            # prob_n_full_dict[base].to_csv("/tmp/prob_" + base + ".csv", index=False)

        return prob_n_full_dict

    def get_dfs_for_calc_probN(self, bases):
        """
        get the relevant dataframe for the calculation of prob N.
        :param bases:
        :return: current_state_df, mapping_df, reads_df, posteriors_df
        """
        reads_df = pd.DataFrame.from_csv(self.paths.read_quals, index_col=None)
        reads_df[ReadsFormat.Quals] = reads_df[ReadsFormat.Quals].apply(literal_eval)

        mapping_df = pd.DataFrame.from_csv(self.paths.mapping, index_col=None)

        current_state_df = pd.DataFrame.from_csv(self.paths.current_state, index_col=None)

        posteriors_df = None
        if (self.iteration_index > 0):
            posteriors_df = pd.DataFrame.from_csv(self.paths.posteriors, index_col=None)

        for base in bases:
            mapping_df[base] = mapping_df[base].apply(int).apply(lambda r: bin(r)[2:].zfill(2 * self.read_len))
            current_state_df[base] = current_state_df[base].apply(int).apply(
                lambda r: bin(r)[2:].zfill(2 * self.read_len))
            current_state_df[CurrentStateFormat.ProbN + base] = 0

        for i in range(2 * self.read_len):
            reads_df[str(i) + "_prob_success"] = reads_df[ReadsFormat.Quals].apply(
                lambda r: quals.ONE_MINUS_P[int(r[i])])
            reads_df[str(i) + "_prob_fail"] = reads_df[ReadsFormat.Quals].apply(lambda r: quals.P_DIV_3[int(r[i])])
            for base in Base.all:
                mapping_df[str(i) + base] = mapping_df[base].apply(lambda r: int(r[i]))
                current_state_df[str(i) + base] = current_state_df[base].apply(lambda r: int(r[i]))

        return current_state_df, mapping_df, reads_df, posteriors_df

    def get_dfs_for_calc_likelihoods(self, bases):
        reads_df = pd.DataFrame.from_csv(self.paths.read_quals, index_col=None)
        reads_df[ReadsFormat.Quals] = reads_df[ReadsFormat.Quals].apply(literal_eval)

        mapping_df = pd.DataFrame.from_csv(self.paths.mapping, index_col=None)

        # current_state_df = pd.DataFrame.from_csv(self.paths.current_state, index_col=None)

        posteriors_df = None
        if (self.iteration_index > 0):
            posteriors_df = pd.DataFrame.from_csv(self.paths.posteriors, index_col=None)

        for base in bases:
            mapping_df[base] = mapping_df[base].apply(int).apply(lambda r: bin(r)[2:].zfill(2 * self.read_len))
            # current_state_df[base] = current_state_df[base].apply(int).apply(
            #     lambda r: bin(r)[2:].zfill(2 * self.read_len))

        for i in range(2 * self.read_len):
            reads_df[str(i) + "_prob_success"] = reads_df[ReadsFormat.Quals].apply(
                lambda r: quals.ONE_MINUS_P[int(r[i])])
            reads_df[str(i) + "_prob_fail"] = reads_df[ReadsFormat.Quals].apply(lambda r: quals.P_DIV_3[int(r[i])])
            for base in Base.all:
                mapping_df[str(i) + base] = mapping_df[base].apply(lambda r: int(r[i]))
                # current_state_df[str(i) + base] = current_state_df[base].apply(lambda r: int(r[i]))

        return mapping_df, reads_df, posteriors_df

    def get_columns_for_calc_probN(self, bases):

        reads_data_columns = {Base.A: [], Base.C: [], Base.G: [], Base.T: []}
        ref_data_columns = {Base.A: [], Base.C: [], Base.G: [], Base.T: []}
        prob_success_data_columns = []
        prob_failure_data_columns = []
        rename_dict = {}

        for i in range(2 * self.read_len):
            for base in bases:
                reads_data_columns[base] += [str(i) + base + "_read"]
                ref_data_columns[base] += [str(i) + base + "_ref"]
                rename_dict.update({str(i) + base + "_read": str(i), str(i) + base + "_ref": str(i)})
            prob_success_data_columns += [str(i) + "_prob_success"]
            prob_failure_data_columns += [str(i) + "_prob_fail"]
            rename_dict.update({str(i) + "_prob_success": str(i), str(i) + "_prob_fail": str(i)})

        return reads_data_columns, ref_data_columns, prob_success_data_columns, prob_failure_data_columns, rename_dict

    def get_columns_for_calc_likelihoods(self, bases):

        reads_data_columns = {Base.A: [], Base.C: [], Base.G: [], Base.T: []}
        # ref_data_columns = {Base.A: [], Base.C: [], Base.G: [], Base.T: []}
        prob_success_data_columns = []
        prob_failure_data_columns = []
        rename_dict = {}
        prob_n_columns = []
        for i in range(2 * self.read_len):
            for base in bases:
                reads_data_columns[base].append(str(i) + base)
                # ref_data_columns[base] += [str(i) + base + "_ref"]
                rename_dict.update({str(i) + base : str(i)})
            prob_success_data_columns.append(str(i) + "_prob_success")
            prob_failure_data_columns.append(str(i) + "_prob_fail")
            prob_n_columns.append(str(i))
            rename_dict.update({str(i) + "_prob_success": str(i), str(i) + "_prob_fail": str(i)})

        return prob_n_columns, reads_data_columns, prob_success_data_columns, prob_failure_data_columns, rename_dict

    @time_it
    def calc_likelihoods(self, prob_n_dict):
        """
        Calc the P(r|S)
        this is the first part of the calculation of P(S|r)
        P(r|S) = Mul_k(sum_n(P(b_k|n)P(n)))
        e.g:
            exp ( sum over all bases in each read of
                    log(
                        sum over {A, C, G, T} of P(S[x] == N)*P(r[x] == N) were S is the reference which the read mapped to.
                        )
                )

        e^(log_a + log_b) = e^(log_a)*e^(log_b) = a*b (faster way)
        :return:
        """

        bases = [Base.A, Base.C, Base.G, Base.T]
        mapping_df, reads_df, posteriors_df = self.get_dfs_for_calc_likelihoods(bases)
        prob_n_columns, reads_data_columns, prob_success_data_columns, prob_failure_data_columns, rename_dict = self.get_columns_for_calc_likelihoods(bases)

        likelihoods = []
        chunk_size = 1000
        chunks = int(len(reads_df) / chunk_size) + 1
        logging.info("\nLIKELIHOOD: chunk size = {}, number of chunks = {}, mapping size = {}, reads size = {}".
                     format(chunk_size, chunks, len(mapping_df), len(reads_df)))
        for chunk in np.array_split(reads_df, chunks):
            chunk_reads_full_data = pd.DataFrame.merge(chunk, mapping_df, on=HeadersFormat.Group_id, how='left')
            res_dict = {Base.A: [], Base.C: [], Base.G: [], Base.T: []}
            for base in bases:
                chunk_full_data_df = chunk_reads_full_data.merge(prob_n_dict[base],
                                                                 on=[MappingForamt.Ref_id, HeadersFormat.Region],
                                                                 how='left')
                reads = chunk_full_data_df[reads_data_columns[base]]
                reads.rename(columns=rename_dict, inplace=True)
                prob_n = chunk_full_data_df[prob_n_columns]
                prob_success = chunk_full_data_df[prob_success_data_columns].rename(columns=rename_dict)
                prob_fail = chunk_full_data_df[prob_failure_data_columns].rename(columns=rename_dict)

                reads_prob = reads.multiply(prob_success) + (1-reads).multiply(prob_fail)
                res_dict[base] = prob_n.multiply(reads_prob)
                size_for_debug = len(chunk_full_data_df)

            res = res_dict[Base.A] + res_dict[Base.C] + res_dict[Base.G] + res_dict[Base.T]
            res = res[prob_n_columns].applymap(np.log)
            likelihood = pd.DataFrame({HeadersFormat.Likelihood: res.apply(sum, axis=1)}).apply(np.exp)
            likelihood[MappingForamt.Ref_id] = chunk_reads_full_data[MappingForamt.Ref_id]
            likelihood[ReadsFullDataFormat.Id] = chunk_reads_full_data[ReadsFullDataFormat.Id]
            likelihoods.append(likelihood)

        # logging.info("likelihood, before concat")
        likelihood_df = pd.concat(likelihoods, ignore_index=True)
        return likelihood_df


    @time_it
    def calc_posteriors(self, prob_n_dict):
        """
        P(S|r) = P(r|S)P(s) / sum_i(P(r|Si)P(Si))
        :param prob_n_dict:
        :return:
        """
        likelihood_df = self.calc_likelihoods(prob_n_dict)
        priors_df = pd.DataFrame.from_csv(self.paths.current_state, index_col=None)[
            [CurrentStateFormat.Reference_id, CurrentStateFormat.Priors]].drop_duplicates()

        posteriors = likelihood_df.merge(priors_df, on=HeadersFormat.Unique_Ref_id, how='left')

        # P(r|S)P(s)
        posteriors[PosteriorsFormat.Posterior] = posteriors.\
            apply(lambda r: float(r[PosteriorsFormat.Likelihood]) * float(r[CurrentStateFormat.Priors]), axis=1)

        # denominator_df --> sum_i(P(r|Si)P(Si))
        denominator_df = pd.DataFrame({'denominator': posteriors.groupby(PosteriorsFormat.Read_id)[
            PosteriorsFormat.Posterior].sum()}).reset_index()

        posteriors = posteriors.merge(denominator_df, how='left')

        # P(S|r) = P(r|S)P(s) / sum_i(P(r|Si)P(Si))
        posteriors[PosteriorsFormat.Posterior] = posteriors.apply(
            lambda r: r[PosteriorsFormat.Posterior] / r['denominator'] if r['denominator'] != 0 else 0, axis=1)
        posteriors = posteriors.fillna(0)

        posteriors = posteriors[[PosteriorsFormat.Posterior,
                                 PosteriorsFormat.Read_id,
                                 PosteriorsFormat.Ref_id]]

        posteriors.to_csv(self.paths.posteriors, index=False)

    @time_it
    def calc_priors(self):
        """
            calc the priors P(s_i) = sum_j(P(s_i|rj)/ sum_i, j(P(s_i|r_j)
            where P(s|r) is the posteriors.
        """
        posteriors_df = pd.DataFrame.from_csv(self.paths.posteriors, index_col=None)
        new_priors = posteriors_df.groupby(PosteriorsFormat.Ref_id)[PosteriorsFormat.Posterior].sum().reset_index()

        # normalize in the reference # regions:
        curr_state_df = pd.DataFrame.from_csv(self.paths.current_state, index_col=None)
        ref_weight = curr_state_df[[CurrentStateFormat.Reference_id, CurrentStateFormat.Weight]].drop_duplicates()
        new_priors = ref_weight.merge(new_priors, on=CurrentStateFormat.Reference_id, how='right')
        new_priors[CurrentStateFormat.Priors] = new_priors.\
            apply(lambda r: r[PosteriorsFormat.Posterior]/ r[CurrentStateFormat.Weight], axis=1)

        posteriors_sum = new_priors.Priors.sum()
        new_priors[CurrentStateFormat.Priors] = new_priors.apply(lambda r: r[CurrentStateFormat.Priors] / posteriors_sum, axis=1)

        new_priors = new_priors[[CurrentStateFormat.Priors, CurrentStateFormat.Reference_id]].drop_duplicates()

        # get privius priors values:
        old_priors = curr_state_df[[CurrentStateFormat.Reference_id, CurrentStateFormat.Priors]].drop_duplicates()

        # Take only only the prior values which where not calculated this round.
        old_priors = old_priors[(~old_priors[CurrentStateFormat.Reference_id].isin(new_priors[PosteriorsFormat.Ref_id]))]
        priors_df = pd.concat([old_priors, new_priors], ignore_index=True)

        # update curr state with the new priors values.
        cols = [CurrentStateFormat.Reference_id,
                CurrentStateFormat.Region,
                CurrentStateFormat.Weight] + CurrentStateFormat.Bases.all
        curr_state_df = curr_state_df[cols]
        curr_state_df = curr_state_df.merge(priors_df)

        curr_state_df.to_csv(self.paths.current_state, index=False)

    def initialize_primers(self, primers_path):
        primers_df = pd.DataFrame.from_csv(primers_path, index_col=None)
        for region in primers_df.columns:
            self.primers.add_primer(int(region), primers_df[region].fillna('').tolist())
        logging.info("Primers were update with {} regions".format(len(primers_df.columns)))

    def _get_new_posteriors_and_priors(self, posteriors_df, minors, refs):
        if len(minors) > 0:
            minor_fraction_avg = np.mean(minors)
        else:
            minor_fraction_avg = 0

        ref_ids = list(set([str(r.get(CurrentStateFormat.Reference_id)) for r in refs]))
        ref_ids = sorted(ref_ids, key = len, reverse=False)
        if len(ref_ids) > 1:
            for ref in refs:
                if str(ref.get(CurrentStateFormat.Reference_id)) == ref_ids[0]:
                    ref[CurrentStateFormat.Priors] = ref[CurrentStateFormat.Priors]*(1 - minor_fraction_avg)
                else:
                    ref[CurrentStateFormat.Priors] = ref[CurrentStateFormat.Priors] * minor_fraction_avg

        new_posteriors = []
        for _, posterior in posteriors_df.iterrows():
            posterior_dict = {PosteriorsFormat.Ref_id: ref_ids[0],
                              PosteriorsFormat.Read_id: posterior[PosteriorsFormat.Read_id],
                              PosteriorsFormat.Posterior: posterior[PosteriorsFormat.Posterior] * (
                              1 - minor_fraction_avg)}
            new_posteriors += [posterior_dict]
            if minor_fraction_avg and len(ref_ids) > 1:
                posterior_dict = {PosteriorsFormat.Ref_id: ref_ids[1],
                                  PosteriorsFormat.Read_id: posterior[PosteriorsFormat.Read_id],
                                  PosteriorsFormat.Posterior: posterior[PosteriorsFormat.Posterior] * (minor_fraction_avg)}
                new_posteriors += [posterior_dict]

        return new_posteriors, refs


    def _get_new_references(self, test, reference_suffix):
        # FIND BEST MATCH REFERENCES
        full_probs_df = test.copy()
        best_match = full_probs_df.copy()
        best_binary_sequence = full_probs_df['best_sequence'].apply(sequence_to_bin)
        best_match[Base.A] = best_binary_sequence.apply(lambda r: r[0])
        best_match[Base.C] = best_binary_sequence.apply(lambda r: r[1])
        best_match[Base.G] = best_binary_sequence.apply(lambda r: r[2])
        best_match[Base.T] = best_binary_sequence.apply(lambda r: r[3])
        # if split -> change the priors to prior*(1-avg_minor) else prior stay the same.
        best_match[CurrentStateFormat.Priors].update(best_match[CurrentStateFormat.Priors]*(1-best_match['is_split']*best_match['avg_minor']))

        # FIND SECOND BEST REFERENCES
        second_best_match = full_probs_df[full_probs_df['is_split'] == True]
        if len(second_best_match) > 0:
            logging.info('************************')
            logging.info("splitting {} references".format(len(second_best_match.drop_duplicates(CurrentStateFormat.Reference_id))))
            logging.info("splitting {}".format(second_best_match[CurrentStateFormat.Reference_id].drop_duplicates().tolist()))
            logging.info('************************\n')
            second_best_match[CurrentStateFormat.Reference_id].update(
                second_best_match[CurrentStateFormat.Reference_id].apply(lambda r: r + reference_suffix))
            second_best_binary_sequence = second_best_match['2-best_sequence'].apply(sequence_to_bin)
            second_best_match[Base.A] = second_best_binary_sequence.apply(lambda r: r[0])
            second_best_match[Base.C] = second_best_binary_sequence.apply(lambda r: r[1])
            second_best_match[Base.G] = second_best_binary_sequence.apply(lambda r: r[2])
            second_best_match[Base.T] = second_best_binary_sequence.apply(lambda r: r[3])
            second_best_match[CurrentStateFormat.Priors].update(
                second_best_match[CurrentStateFormat.Priors]*second_best_match['avg_minor'])

        new_refs = pd.concat([best_match, second_best_match], ignore_index=True)[CurrentStateFormat.Bases.all +
                                                                                [CurrentStateFormat.Reference_id,
                                                                                 CurrentStateFormat.Priors,
                                                                                 CurrentStateFormat.Region,
                                                                                 CurrentStateFormat.Weight]]
        logging.info("References: {}".format(new_refs[CurrentStateFormat.Reference_id].drop_duplicates().tolist()))
        return new_refs

    def _get_new_posterior(self, max_prob_n_full_data, full_posteriors, reference_suffix):
        reference_minors = max_prob_n_full_data[['is_split', 'avg_minor', CurrentStateFormat.Reference_id]].drop_duplicates()
        # we need only reference from the original posteriors dataframe
        posteriors = reference_minors.merge(full_posteriors, on=CurrentStateFormat.Reference_id, how='right')
        posteriors[PosteriorsFormat.Posterior] = posteriors[PosteriorsFormat.Posterior]*(1-posteriors['is_split']*posteriors['avg_minor'])
        # Create the new posteriors (of the splitted references)
        new_posteriors = posteriors[posteriors['is_split']==True]
        new_posteriors[CurrentStateFormat.Reference_id].update(new_posteriors[CurrentStateFormat.Reference_id].apply(lambda r: r + reference_suffix))
        new_posteriors[PosteriorsFormat.Posterior].update(new_posteriors[PosteriorsFormat.Posterior] * new_posteriors['avg_minor'])

        posteriors_df = pd.concat([posteriors, new_posteriors], ignore_index=True)[[PosteriorsFormat.Posterior,
                                                                                   PosteriorsFormat.Read_id,
                                                                                   PosteriorsFormat.Ref_id]]

        logging.info("Posteriors - References: {}".format(posteriors_df[PosteriorsFormat.Ref_id].drop_duplicates().tolist()))
        return posteriors_df.copy()

    @time_it
    def update_references(self, probN_dict):
        curr_state_df = pd.DataFrame.from_csv(self.paths.current_state, index_col=None)
        full_posteriors = pd.DataFrame.from_csv(self.paths.posteriors, index_col=None)

        # create table of the probabilities: the columns contains the the probablity of each base in each index - 0A, 0C, 0G, 0T, 1A, etc
        #                                    the rows are the references and the regions
        prob_n_ac = pd.merge(probN_dict[Base.A], probN_dict[Base.C], on=[HeadersFormat.Unique_Ref_id, HeadersFormat.Region], suffixes=(Base.A, Base.C))
        prob_n_gt = pd.merge(probN_dict[Base.G], probN_dict[Base.T], on=[HeadersFormat.Unique_Ref_id, HeadersFormat.Region], suffixes=(Base.G, Base.T))
        prob_n_full = pd.merge(prob_n_ac, prob_n_gt, on=[HeadersFormat.Unique_Ref_id, HeadersFormat.Region])

        base_ix_columns = []
        bases_columns = []
        second_bases_columns = []
        second_probs_columns = []

        for i in range(0, 2*self.read_len):
            base_ix_columns.append(str(i))  # the probability of the most probable base in index i
            bases_columns.append('base_' + str(i))  # the most probable base to be in in index i
            second_bases_columns.append('2-base_' + str(i))  # the second most probable base to be in in index i
            second_probs_columns.append('2-' + str(i)) # the probability of the second most probable base in index i

        for base_ix in base_ix_columns:
            ix_cols = [base_ix + Base.A, base_ix + Base.C, base_ix + Base.G, base_ix + Base.T]
            # prob_n_for_ix: row for each reference-region pair
            prob_n_for_ix = prob_n_full[ix_cols]

            #  take the second best base, if the probability is larger than the the split threshold
            prob_n_full[['2-' + base_ix, base_ix]] = pd.DataFrame(np.sort(prob_n_for_ix.values)[:, -2:])
            prob_n_full['2-' + base_ix] = prob_n_full['2-' + base_ix].apply(lambda r: r if r > self.th.min_minor_prob_for_split else 0)

            arank = prob_n_for_ix.apply(np.argsort, axis=1)
            argsort_dict = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
            prob_n_full[['2-base_' + base_ix, 'base_' + base_ix]] = arank[[2, 3]].applymap(lambda x: argsort_dict[x])
            # select best probability if the second is too low.
            prob_n_full['2-base_' + base_ix] = prob_n_full.apply(lambda r: r['base_' + base_ix] if r['2-' + base_ix]==0 else r['2-base_' + base_ix], axis=1)


        prob_n_full['best_sequence'] = prob_n_full[bases_columns].apply(lambda r: "".join(r.tolist()), axis=1)
        prob_n_full['2-best_sequence'] = prob_n_full[second_bases_columns].apply(lambda r: "".join(r.tolist()), axis=1)
        prob_n_full['second_best_counter'] = prob_n_full[second_probs_columns].apply(np.count_nonzero, axis=1)

        prob_n_full[second_probs_columns] = prob_n_full[second_probs_columns].replace(0, np.NaN)
        prob_n_full['minor'] = prob_n_full[second_probs_columns].mean(axis=1)
        prob_n_full['minor'].fillna(0, inplace=True)
        prob_n_full[second_probs_columns].fillna(0, inplace=True)

        prob_n_relevant = prob_n_full[[CurrentStateFormat.Reference_id, CurrentStateFormat.Region, 'best_sequence','2-best_sequence' ,'second_best_counter', 'minor']]

        max_prob_n_full_data = prob_n_relevant.merge(curr_state_df, on=[HeadersFormat.Unique_Ref_id, HeadersFormat.Region])

        logging.info("NAN: minor = {}, Prior = {}, Weight = {}".format(max_prob_n_full_data['minor'].isnull().values.sum(),
                                                                       max_prob_n_full_data[CurrentStateFormat.Priors].isnull().values.sum(),
                                                                       max_prob_n_full_data[CurrentStateFormat.Weight].isnull().values.sum()))
        max_prob_n_full_data['expected_coverage_minor'] = max_prob_n_full_data['minor'] * max_prob_n_full_data[CurrentStateFormat.Priors] * self.number_of_regions / max_prob_n_full_data[CurrentStateFormat.Weight]

        if self.allow_split and self.iteration_index > 0:
            # 1. The rate of changed bases in the split reference: (how many bases will change)/(length of read)
            # 2. The min rate to create a new split reference: (prior*minor) * (# of regions/ reference weight)

            max_prob_n_full_data['is_split'] = (max_prob_n_full_data['second_best_counter'] / float(2*self.read_len) < self.th.max_changed_bases_rate_for_split) \
                                               & (max_prob_n_full_data['expected_coverage_minor'] >= self.th.min_coverage_for_split)
            max_prob_n_full_data['is_split'].update(
                max_prob_n_full_data.groupby(CurrentStateFormat.Reference_id)['is_split'].transform(any).reset_index()['is_split'])
            max_prob_n_full_data['avg_minor'] = max_prob_n_full_data.groupby(CurrentStateFormat.Reference_id)['minor'].transform(np.mean).reset_index()['minor']
        else:
            max_prob_n_full_data['is_split'] = False
            max_prob_n_full_data['avg_minor'] = 0

        split_reference_suffix = self.iteration_index*math.pow(10, -1*self.iteration_index - 2)
        new_refs = self._get_new_references(max_prob_n_full_data, split_reference_suffix)
        new_posteriors = self._get_new_posterior(max_prob_n_full_data, full_posteriors, split_reference_suffix)

        new_refs.to_csv(self.paths.current_state, index=False)
        new_posteriors.to_csv(self.paths.posteriors, index=False)
        logging.info("{}/{} references found".format(len(new_refs.drop_duplicates(ReferenceFormat.Ref_Id)), len(curr_state_df.drop_duplicates(CurrentStateFormat.Reference_id))))

        return new_refs

    @time_it
    def merge_references(self, curr_refs, curr_posteriors):
        """
        merge reference that have similar sequences in all their amplified regions
        :param curr_refs:
        :param curr_posteriors:
        :return: the merged references and the merged posteriors
        """
        curr_refs_len = len(curr_refs.drop_duplicates(CurrentStateFormat.Reference_id))

        # find duplication of regions:
        # example #1:
        # Reference #1: region 1 --> 'AA' (id =1)
        #               region 2 --> 'AC' (id =2)
        # Reference #2: region 1 --> 'CC' (id =3)
        #               region 2 --> 'AC' (id =4)
        curr_refs[ReferenceFormat.Original_Id] = curr_refs[CurrentStateFormat.Reference_id]

        # example #1:
        # Reference #1: region 1 --> 'AA' (unique id =1)
        #               region 2 --> 'AC' (unique id =2)
        # Reference #2: region 1 --> 'CC' (unique id =3)
        #               region 2 --> 'AC' ****(unique id =2)****
        curr_refs = self.get_unique_group_for_similar_sequences(curr_refs)

        # find for each reference the groups ids of his regions
        # example #2:
        # --> the unique groups of reference #1: [unique id 1, unique id 2]
        # --> the unique groups of reference #2: [unique id 3, unique id 2]
        # --> the unique groups of reference #3: [unique id 3, unique id 2]
        ref_to_unique_groups = pd.DataFrame({'uniques_groups': curr_refs.groupby(ReferenceFormat.Original_Id)[
            'Group'].apply(list)}).reset_index()
        ref_to_unique_groups['uniques_groups'] = ref_to_unique_groups['uniques_groups'].astype(str)

        # find unique references (should be the same in all regions
        # example #3:
        # --> reference #1: original id = 1, unique id = 1
        # --> reference #2: original id = 2, unique id = 2
        # --> reference #3: original id = 3, unique id = 2 (the same unique id as reference #2)
        ref_to_unique_groups[UnqiueRefToRefFormat.Unique_id] = \
                ref_to_unique_groups.groupby('uniques_groups').grouper.group_info[0]

        ref_to_unique_groups = ref_to_unique_groups[[UnqiueRefToRefFormat.Unique_id, ReferenceFormat.Original_Id]]
        curr_refs = curr_refs[CurrentStateFormat.Bases.all + [CurrentStateFormat.Region,
                                                              CurrentStateFormat.Weight,
                                                              CurrentStateFormat.Priors,
                                                              ReferenceFormat.Original_Id]]

        # add the unique id to the curr_refs
        unique_ref_df = curr_refs.merge(ref_to_unique_groups,
                                          on=ReferenceFormat.Original_Id,
                                          how='right')

        unique_ref_df_with_priors = unique_ref_df[[UnqiueRefToRefFormat.Unique_id,
                                                   CurrentStateFormat.Priors]].drop_duplicates()

        # using the new unique id to calc the merged priors.
        unique_ref_df_with_priors = unique_ref_df_with_priors.groupby(UnqiueRefToRefFormat.Unique_id)[CurrentStateFormat.Priors]\
                                                             .sum().reset_index()

        # merge the ref data with the priors, using the merged id
        unique_ref_df = unique_ref_df[CurrentStateFormat.Bases.all + [CurrentStateFormat.Region,
                                                                      CurrentStateFormat.Weight,
                                                                      UnqiueRefToRefFormat.Unique_id,
                                                                      ReferenceFormat.Original_Id]]
        unique_ref_df = unique_ref_df.merge(unique_ref_df_with_priors, on=UnqiueRefToRefFormat.Unique_id, how='left')

        # keep only one instance of each merged group:
        unique_ref_df.drop_duplicates([UnqiueRefToRefFormat.Unique_id, CurrentStateFormat.Region], inplace=True)

        # keep only the original id
        merged_refs = unique_ref_df[CurrentStateFormat.Bases.all + [CurrentStateFormat.Region,
                                                                    CurrentStateFormat.Weight,
                                                                    CurrentStateFormat.Priors,
                                                                    ReferenceFormat.Original_Id]]
        merged_refs.rename(columns={ReferenceFormat.Original_Id: CurrentStateFormat.Reference_id}, inplace=True)

        logging.info("Found {}/{} unique reference".format(
            len(unique_ref_df.drop_duplicates(CurrentStateFormat.Reference_id)), curr_refs_len))

        merged_posteriors = curr_posteriors.merge(ref_to_unique_groups,
                                                  left_on=PosteriorsFormat.Ref_id,
                                                  right_on=CurrentStateFormat.Reference_id,
                                                  how='left')
        merged_posteriors = merged_posteriors[[PosteriorsFormat.Posterior, PosteriorsFormat.Ref_id, PosteriorsFormat.Read_id]]
        merged_posteriors[PosteriorsFormat.Posterior] = merged_posteriors.groupby([PosteriorsFormat.Ref_id,
                                                                                   PosteriorsFormat.Read_id])[PosteriorsFormat.Posterior].transform(sum)
        merged_posteriors.drop_duplicates(inplace=True)

        return merged_refs, merged_posteriors

    def get_unique_group_for_similar_sequences(self, old_reference, min_similarity_rate_for_merge=None):
        if min_similarity_rate_for_merge is None:
            min_similarity_rate_for_merge = self.th.min_similar_bases_rate_for_merge

        references_scored = self.get_similarity_score_for_each_reference_pair(old_reference)
        similar_references_all_data = references_scored[references_scored['Score'] >= min_similarity_rate_for_merge]

        similar_references_map = get_matched_references(similar_references_all_data)
        similar_references_and_region = get_similar_references_and_region(similar_references_map)

        non_similar_references_map = self.get_non_similar_references(old_reference, similar_references_and_region)

        all_sequences = pd.concat([non_similar_references_map, similar_references_map])

        # find group id for all similar sequences in the same region
        all_sequences['Group'] = all_sequences.groupby([CurrentStateFormat.Reference_id + "_x",
                                                        CurrentStateFormat.Region]).grouper.group_info[0]

        reference_with_group = pd.concat([all_sequences[[CurrentStateFormat.Reference_id + "_x",
                                                         CurrentStateFormat.Region,
                                                         'Group']].
                                         drop_duplicates().
                                         rename(columns={CurrentStateFormat.Reference_id + "_x": CurrentStateFormat.Reference_id}),
                                          all_sequences[[CurrentStateFormat.Reference_id + "_y",
                                                         CurrentStateFormat.Region,
                                                         'Group']].
                                         drop_duplicates().
                                         rename(columns={CurrentStateFormat.Reference_id + "_y": CurrentStateFormat.Reference_id})])

        reference_with_group = reference_with_group.drop_duplicates()

        # add the data from the old_reference df
        all_data = reference_with_group.merge(old_reference, on=[CurrentStateFormat.Reference_id, CurrentStateFormat.Region])

        return all_data


    def get_non_similar_references(self, all_reference_df, similar_reference):
        """
        get only the references from all_reference_df that not in similar_reference (not as ref_a and not as ref_b)
        :param all_reference_df:
        :param similar_reference:
        :return:
        """
        all_reference_df_merged = all_reference_df.merge(similar_reference,
                                                   on=[CurrentStateFormat.Reference_id, CurrentStateFormat.Region],
                                                   how='left', indicator=True)
        not_similar_references = all_reference_df_merged[all_reference_df_merged['_merge'] == 'left_only'][
                        [CurrentStateFormat.Reference_id,
                         CurrentStateFormat.Region]]
        not_similar_references[CurrentStateFormat.Reference_id + "_x"] = not_similar_references[CurrentStateFormat.Reference_id]
        not_similar_references[CurrentStateFormat.Reference_id + "_y"] = not_similar_references[CurrentStateFormat.Reference_id]
        not_similar_references = not_similar_references[[CurrentStateFormat.Reference_id + "_y",
                                                         CurrentStateFormat.Reference_id + "_x",
                                                         CurrentStateFormat.Region]]
        return not_similar_references

    def get_similarity_score_for_each_reference_pair(self, reference_df):
        for base in Base.all:
            reference_df[base].update(reference_df[base].apply(lambda r: int(r)))

        all_sequences = reference_df.merge(reference_df.copy(), on=CurrentStateFormat.Region, how='left')
        a = np.bitwise_and(all_sequences['A_x'], all_sequences['A_y'])
        c = np.bitwise_and(all_sequences['C_x'], all_sequences['C_y'])
        g = np.bitwise_and(all_sequences['G_x'], all_sequences['G_y'])
        t = np.bitwise_and(all_sequences['T_x'], all_sequences['T_y'])
        all_sequences['Score'] = np.bitwise_or(np.bitwise_or(a, c), np.bitwise_or(g, t))
        all_sequences['Score'].update(all_sequences['Score'].apply(lambda r: popcount(r) / float(self.read_len * 2)))
        # take only the similar sequences:
        similar_sequences = all_sequences[all_sequences['Score'] != 1]
        # take only one instance of the merged couple:
        similar_sequences = similar_sequences[similar_sequences[CurrentStateFormat.Priors + "_x"] >= similar_sequences[CurrentStateFormat.Priors + "_y"]]
        similar_sequences = similar_sequences.rename(columns={CurrentStateFormat.Priors + "_x": CurrentStateFormat.Priors})
        logging.info("merge_similar_sequences: len = {}".format(len(similar_sequences)))

        return similar_sequences[[CurrentStateFormat.Reference_id + '_x',
                                  CurrentStateFormat.Reference_id + '_y',
                                  CurrentStateFormat.Region,
                                  CurrentStateFormat.Priors,
                                  'Score']]

    @time_it
    def is_stable_state(self, new_reference_df):
        if self.iteration_index <= 0:
            return False

        new_df = new_reference_df[[CurrentStateFormat.Reference_id, CurrentStateFormat.Priors]].drop_duplicates()
        merged = self.prev_priors_for_stability_test.merge(new_df, on=CurrentStateFormat.Reference_id, how='outer')
        merged.fillna(0, inplace=True)
        merged['Prior_diff'] = merged[CurrentStateFormat.Priors + "_x"] - merged[CurrentStateFormat.Priors + "_y"]
        merged['is_stable'] = merged['Prior_diff'].apply(lambda r:
                                                         True if abs(r) < self.th.max_priors_diff_for_stability_test else False)
        if False in merged['is_stable'].tolist():
            logging.info("Priors are not stable")
            return False
        else:
            logging.info("Stable state!")
            new_reference_df['Sequence'] = new_reference_df.apply(lambda r: bin_to_sequence(r[CurrentStateFormat.Bases.A],
                                                                            r[CurrentStateFormat.Bases.C],
                                                                            r[CurrentStateFormat.Bases.G],
                                                                            r[CurrentStateFormat.Bases.T],
                                                                            2 * self.read_len), axis=1)
            new_reference_df[[CurrentStateFormat.Reference_id,
                      CurrentStateFormat.Region,
                      CurrentStateFormat.Priors,
                      'Sequence']].to_csv(self.paths.final_results, index=False)
            return True


    @time_it
    def do_iteration(self):
        prob_n_dict = self.calc_prob_N()
        self.calc_posteriors(prob_n_dict)
        self.calc_priors()
        new_reference_df =  self.update_references(prob_n_dict)
        is_stable = self.is_stable_state(new_reference_df)
        return is_stable


def get_matched_references(test):
    """
    allowing each reference a to map to multiple references
    allowing each reference b to be mapped by only one reference
    :param similar_sequences: df - ref_x ref_y prior score region etc
    :return: df: ref_x ref_y region
    """
    logging.debug("Found {} sequences".format(len(test)))
    similar_sequences = test.copy()
    # allowing each reference b to be mapped by only one reference
    similar_sequences.sort_values(CurrentStateFormat.Priors, inplace=True, ascending=False)
    similar_sequences.drop_duplicates([CurrentStateFormat.Region, CurrentStateFormat.Reference_id + "_y"], inplace=True)

    # allowing each reference a to map to multiple references: ref_a should not contains any reference from ref_b
    def get_unique_similar_groups(df):
        to_list = df[CurrentStateFormat.Reference_id + "_y"]
        df = df[~df[CurrentStateFormat.Reference_id + "_x"].isin(to_list)]
        return df

    similar_sequences = similar_sequences.groupby(CurrentStateFormat.Region).apply(get_unique_similar_groups)
    if similar_sequences.empty:
        return pd.DataFrame(columns=[CurrentStateFormat.Reference_id + "_x",
                                     CurrentStateFormat.Reference_id + "_y",
                                     CurrentStateFormat.Region])
    logging.debug("Found {} similar sequences in single region".format(len(similar_sequences)))
    return similar_sequences[[CurrentStateFormat.Reference_id + "_x",
                              CurrentStateFormat.Reference_id + "_y",
                              CurrentStateFormat.Region]]


def get_similar_references_and_region(similar_sequences):
    ref_list = pd.concat([similar_sequences[[CurrentStateFormat.Reference_id + "_x",
                                             CurrentStateFormat.Region]].
                         drop_duplicates().
                         rename(columns={CurrentStateFormat.Reference_id + "_x": CurrentStateFormat.Reference_id}),
                          similar_sequences[[CurrentStateFormat.Reference_id + "_y",
                                             CurrentStateFormat.Region]].
                         drop_duplicates().
                         rename(columns={CurrentStateFormat.Reference_id + "_y": CurrentStateFormat.Reference_id})])
    return ref_list


def get_emirge_iteration():
    data_path = "/home/vered/EMIRGE/EMIRGE-data/"
    reads_fastq_path = data_path + "RDB53_CATTGACG_L007_R1_001.fastq"
    reversed_reads_fastq_path = data_path + "RDB53_CATTGACG_L007_R2_001.fastq"
    working_dir = data_path + "testing"
    reference_path = working_dir + "/full_reference_db.csv"
    fasta_path = data_path
    read_len = 126
    emirge_iteration = EmirgeIteration(working_dir, reads_fastq_path, reversed_reads_fastq_path,
                                   fasta_path, read_len, reference_path=reference_path)
    return emirge_iteration

def get_emirge_iteration_mock():
    data_path = "/home/vered/EMIRGE/EMIRGE-data/mock_data/"
    reads_fastq_path = data_path + "reads_mock1.fastq"
    reversed_reads_fastq_path = data_path + "reads_mock2.fastq"
    working_dir = "/home/vered/EMIRGE/EMIRGE-data/mock_tests"
    # reference_path = working_dir + "/full_reference_db.csv"
    fasta_path = data_path
    read_len = 126
    emirge_iteration = EmirgeIteration(working_dir, reads_fastq_path, reversed_reads_fastq_path,
                                   fasta_path, read_len)
    return emirge_iteration


def main():
    define_logger(logging.INFO)
    iteration = get_emirge_iteration_mock()
    iteration.do_iteration()

if __name__ == '__main__':
    main()