

import numpy as np
from emirge_headers import *
from emirge_utills import *
import os
import matplotlib.pyplot as plt
import pandas as pd
from pandas.tools.plotting import table
from optparse import OptionParser, OptionGroup
import sys

VALIDATION_THRESHOLD_THRESHOLD = 0.00001
MAX_ALLOWED_MISMATCH_THRESHOLD = 10


class Header():
    ref_id = 'ref'
    prior = 'prior'
    sequence = 'sequence'
    region = 'region'
    weight = 'weight'
    all = [ref_id, prior, sequence, region, weight]



def validate_priors(df, threshold=VALIDATION_THRESHOLD_THRESHOLD):
    sum_prior = sum(df.drop_duplicates(Header.ref_id)[Header.prior])
    df.prior = df.prior.apply(lambda r: r / sum_prior)

    sum_prior = sum(df.drop_duplicates(Header.ref_id)[Header.prior])

    logging.debug("sum of priors is {}".format(sum_prior))
    if abs(sum_prior - 1) > threshold:
        raise Exception("sum of prior is not 1")


def calc_mismatch_score(seq_t, seq_e):
    score = 0
    for i in range(len(seq_t)):
        if (seq_t[i] != seq_e[i]) & (seq_e[i] != 'N'):
            score += 1
    return score


def calc_dfs_mismatch_score(test_df, expected_df, mismatch_threshold):
    """
    :param test_df:
    :param expected_df:
    :param mismatch_threshold:
    :return: merged df with score col contains the number of mismatch between the sequences
    """
    merged = pd.merge(test_df, expected_df, on=Header.region, suffixes=('_t', '_e'))
    merged['mismatch_score'] = merged.apply(lambda r: calc_mismatch_score(r[Header.sequence + '_t'], r[Header.sequence + '_e']), axis=1)
    merged = merged[merged['mismatch_score'] < mismatch_threshold][['mismatch_score',
                                                           Header.ref_id + '_t',
                                                           Header.ref_id + '_e',
                                                           Header.region,
                                                           Header.prior + '_t',
                                                           Header.prior + '_e',
                                                           Header.weight + '_e',
                                                           Header.weight + '_t']]
    merged['mismatch_score'] = merged.groupby([Header.ref_id + '_e', Header.ref_id + '_t'])['mismatch_score'].transform('sum')
    merged = merged.drop_duplicates([Header.ref_id + '_t', Header.ref_id + '_e'])
    merged['weight_diff'] = merged[Header.weight + '_e'] - merged[Header.weight + '_t']
    merged['prior_diff_rel_e'] = (merged[Header.prior + '_e'] - merged[Header.prior + '_t']) / merged[Header.prior + '_e']
    merged['prior_diff_rel_e'] = merged['prior_diff_rel_e'].apply(lambda r: abs(r))
    merged['prior_diff'] = merged[Header.prior + '_e'] - merged[Header.prior + '_t']
    merged['prior_diff'] = merged['prior_diff'].apply(lambda r: abs(r))
    merged['diff_score'] = merged['weight_diff'] + 2*merged['mismatch_score']

    return merged


def compare_specific_reference(actual_res_path, expected_res_path, result_dir, test_name,  max_allowed_mismatch, ref_expected_id=None):
    test_df = get_test_df(actual_res_path)
    expected_df = get_expected_df(expected_res_path)

    match_df = get_expected_test_map_df(test_df, expected_df, max_allowed_mismatch)
    match_df = match_df[[Header.ref_id + "_e", Header.ref_id + "_t"]].drop_duplicates()
    test_df = test_df.rename(columns={Header.ref_id: Header.ref_id + '_t', Header.sequence: Header.sequence + '_t', Header.prior: Header.prior + '_t'})
    expected_df = expected_df.rename(
        columns={Header.ref_id: Header.ref_id + '_e', Header.sequence: Header.sequence + '_e', Header.prior: Header.prior + '_e'})
    scored_merge_df = pd.merge(test_df, match_df, on=[Header.ref_id + '_t'])
    scored_merge_df = pd.merge(expected_df, scored_merge_df, on=[Header.region, Header.ref_id + '_e'])
    scored_merge_df['mismatch_score'] = scored_merge_df.apply(
        lambda r: calc_mismatch_score(r[Header.sequence + '_t'], r[Header.sequence + '_e']), axis=1)

    if not ref_expected_id:
        ref_expected_ids = scored_merge_df[Header.ref_id + '_e'].drop_duplicates().tolist()
    else:
        ref_expected_ids = list(ref_expected_id)

    for ref_expected_id in ref_expected_ids:
        mapped_data = scored_merge_df[(scored_merge_df[Header.ref_id + '_e'] == ref_expected_id)]
        presented_data = mapped_data[[Header.ref_id + '_e', Header.ref_id + '_t', Header.region, 'prior_e', 'prior_t', 'mismatch_score' ]].copy()
        presented_data = presented_data.sort([Header.ref_id + '_e', Header.ref_id + '_t', 'prior_t', Header.region],ascending=False )
        presented_data.index = range(1, len(presented_data) + 1)
        presented_data['prior_e'] = presented_data['prior_e'].apply(lambda val: "{0:.2f}%".format(val * 100))
        presented_data['prior_t'] = presented_data['prior_t'].apply(lambda val: "{0:.2f}%".format(val * 100))

        ax = plt.subplot(111, frame_on=False)  # no visible frame
        ax.xaxis.set_visible(False)  # hide the x axis
        ax.yaxis.set_visible(False)  # hide the y axis

        table(ax, presented_data, loc='center')  # where df is your data frame
        # table(ax, presented2, loc='center')

        res_name = 'emirge_smurf_'+ test_name + '_reference_id_' + str(ref_expected_id) + '.png'
        im_path = os.path.join(result_dir, res_name)

        plt.tight_layout()
        plt.savefig(im_path, bbox_inches='tight')
        plt.clf()
        logging.info("saving results to: {}".format(im_path))


def test_how_many_mismatch_in_each_reference(test_df, expected_df, match_df, max_allowed_mismatch, original_fasta_path=None):
    match_df.sort_values([Header.prior + '_t'], axis=0, ascending=False, inplace=True)
    match_df = match_df[[Header.ref_id +"_e", Header.ref_id +"_t"]].drop_duplicates(Header.ref_id +"_e")
    test_df = test_df.rename(columns={Header.ref_id: Header.ref_id + '_t', Header.sequence: Header.sequence + '_t'})
    expected_df = expected_df.rename(columns={Header.ref_id: Header.ref_id + '_e', Header.sequence: Header.sequence + '_e'})
    scored_merge_df = pd.merge(test_df, match_df, on=[Header.ref_id + '_t'])
    scored_merge_df = pd.merge(expected_df, scored_merge_df, on=[Header.ref_id + '_e', Header.region])
    scored_merge_df['mismatch_score'] = scored_merge_df.apply(
        lambda r: calc_mismatch_score(r[Header.sequence + '_t'], r[Header.sequence + '_e']), axis=1)

    scored_groupby_ref_expected = scored_merge_df.groupby(Header.ref_id + '_e')

    perfect_match=0
    non_perfect_match=0
    for ref_e, df in scored_groupby_ref_expected:
        curr_perfect_match = True
        mismatch_data = df.apply(
                lambda r: ["ref_t = {}\{}[{}], {}!={} |".format(r.ref_t, r.region, i,  r.sequence_e[i], r.sequence_t[i])
                           for i in range(len(r.sequence_e)) if (r.sequence_e[i] != r.sequence_t[i]) & (r.sequence_e[i]!='N') ],
                axis=1)
        print "\nRef expected: {}:".format(ref_e)
        for mismatch_target in mismatch_data:
            if len(mismatch_target) != 0:
                curr_perfect_match=False
                print "".join(mismatch_target)
        if curr_perfect_match:
            perfect_match += 1
        else:
            non_perfect_match += 1
    print "\n perfect match/ matches {}/{}".format(perfect_match, perfect_match+non_perfect_match)


def test_mock_reads_for_similarity(expected_res_path):
    expected_df = get_expected_df(expected_res_path)
    validate_priors(expected_df)

    scored_merge_df = pd.merge(expected_df, expected_df, on=Header.region, suffixes=('_t', '_e'))
    scored_merge_df['mismatch_score'] = scored_merge_df.apply(
        lambda r: calc_mismatch_score(r[Header.sequence + '_t'], r[Header.sequence + '_e']), axis=1)

    mapped_data = scored_merge_df[
        (scored_merge_df[Header.ref_id + '_e'] != scored_merge_df[Header.ref_id + '_t']) & (scored_merge_df['mismatch_score'] < 5)]

    similar_ref = mapped_data.groupby(Header.ref_id + '_e')
    for g_id, g in similar_ref:
        print "\nref_id = {}:".format(g_id)
        matched_references = g[Header.ref_id + '_t'].drop_duplicates()
        for ref_t in matched_references:
            ref_t_data  = g[g[Header.ref_id + '_t']==ref_t]
            print ref_t_data[['mismatch_score', Header.ref_id + '_t', Header.region]]


def get_expected_test_map_df(test_df, expected_df, max_allowed_mismatch):
    validate_priors(test_df)
    validate_priors(expected_df)

    scored_merge_df = calc_dfs_mismatch_score(test_df, expected_df, max_allowed_mismatch)
    # scored_merge_df = scored_merge_df1[scored_merge_df1['mismatch_score']==0]
    # best_match_weight_df = scored_merge_df.groupby(Header.ref_id + '_e')['weight_diff'].min().reset_index()
    # best_match_df = best_match_weight_df.groupby(Header.ref_id + '_e')['diff_score'].min().reset_index()
    best_match_df = scored_merge_df.groupby([Header.ref_id + '_t'])['diff_score'].min().reset_index()

    match_full_data_df = best_match_df.merge(scored_merge_df, on=['diff_score', Header.ref_id + '_t'], how='inner')
    match_full_data_df = match_full_data_df.drop_duplicates(Header.ref_id + '_t')
    match_full_data_df['actual_prior'] = match_full_data_df.groupby('ref_e')['prior_t'].transform('sum')
    match_full_data_df['actual_prior_diff'] = match_full_data_df['actual_prior'] - match_full_data_df[
        Header.prior + '_e']
    match_full_data_df['actual_prior_diff_rel'] = match_full_data_df['actual_prior_diff'] / match_full_data_df[
        Header.prior + '_e']
    match_full_data_df['count'] = match_full_data_df.groupby('ref_e')['ref_t'].transform('count')

    if abs(sum(match_full_data_df.drop_duplicates('ref_e')['actual_prior']) - 1) > 0.001:
        raise Exception("ERROR! could not find map unique map between the test and the expected results")
    return match_full_data_df


def present_compared_data(match_full_data_df, name, results_dir):
    match_full_data_df = match_full_data_df.drop_duplicates('ref_e')
    presented_data = match_full_data_df[[Header.ref_id + '_e', 'count', 'actual_prior', 'prior_e', 'actual_prior_diff']].copy()
    presented_data.rename(columns={'count': '# of refs',
                                   'actual_prior': 'actual frequency',
                                   'prior_e': 'expected frequency',
                                   'actual_prior_diff':'frequency difference'}, inplace=True)
    presented_data.index = range(1, len(presented_data) + 1)
    presented_data.loc[len(presented_data)+1] = [None, sum(presented_data['# of refs']), sum(presented_data['actual frequency']), sum(presented_data['expected frequency']), None]
    presented_data['expected frequency'] = presented_data['expected frequency'].apply(lambda val: "{0:.2f}%".format(val * 100))
    presented_data['actual frequency'] = presented_data['actual frequency'].apply(lambda val: "{0:.2f}%".format(val * 100))
    presented_data['frequency difference'] = presented_data['frequency difference'].apply(lambda val: "{0:.2f}%".format(val * 100))
    presented_data.index = range(1, len(presented_data)) + ['sum']

    ax = plt.subplot(111, frame_on=False )  # no visible frame
    ax.xaxis.set_visible(False)  # hide the x axis
    ax.yaxis.set_visible(False)  # hide the y axis

    table(ax, presented_data, loc='center')  # where df is your data frame
    # table(ax, presented2, loc='center')

    res_name = name + '_emirge_smurf_compare.png'
    plt.tight_layout()

    im_path = os.path.join(results_dir, res_name)

    plt.savefig(im_path, bbox_inches='tight')
    plt.clf()
    logging.info("saved results to {}".format(im_path))


def compare(test_df, expected_df, results_dir, name, max_allowed_mismatch):
    match_full_data_df = get_expected_test_map_df(test_df, expected_df, max_allowed_mismatch)

    present_compared_data(match_full_data_df, name, results_dir)
    test_how_many_mismatch_in_each_reference(test_df, expected_df, match_full_data_df, max_allowed_mismatch)


def get_presented_data(match_full_data_df):
    presented_data = match_full_data_df[['ref_e', 'count', 'prior_e', 'actual_prior_diff']]
    presented_data.rename(columns={'count': '# of refs',
                                   'actual_prior': 'actual frequency',
                                   'prior_e': 'expected frequency',
                                   'actual_prior_diff':'freq diff'}, inplace=True)
    presented_data.index = range(1, len(presented_data) + 1)
    presented_data.loc[len(presented_data)+1] = ["-", sum(presented_data['# of refs']), sum(presented_data['expected frequency']), 0]
    presented_data['expected frequency'] = presented_data['expected frequency'].apply(lambda val: "{0:.2f}%".format(val * 100))
    presented_data['freq diff'] = presented_data['freq diff'].apply(lambda val: "{0:.2f}%".format(val * 100))
    presented_data.index = range(1, len(presented_data)) + ['sum']
    return presented_data


def get_test_df(path):
    """
    :param path: path to 'final_results.csv produced by emirge_smurf.py
    :return: df hold the final results
    """
    df = pd.DataFrame.from_csv(path, index_col=None)
    df = df.rename(columns = {'Sequence': Header.sequence,
                              HeadersFormat.Region: Header.region,
                              HeadersFormat.Priors: Header.prior,
                              'Unique_Reference_id': Header.ref_id})
    df[Header.weight] = df.groupby(Header.ref_id)[Header.region].transform('count')
    return df[Header.all]


def get_expected_df(path):
    """
    :param path: path to
    :return: path to expected_res.csv produced by mock_creator.py
    """
    df = pd.DataFrame.from_csv(path, index_col=None)
    df = df.rename(columns={'sequence': Header.sequence,
                            'region': Header.region,
                            'prior': Header.prior,
                            'id': Header.ref_id})
    df[Header.weight] = df.groupby(Header.ref_id)[Header.region].transform('count')
    return df[Header.all]


def main_compare_length():
    res_dir = "/home/vered/EMIRGE/EMIRGE-data/mock_for_noam_test/results"
    expected_dir = "/home/vered/EMIRGE/EMIRGE-data/mock_for_noam_test/"
    STATIC = "_static.csv"
    STATIC_WEIGHT = "_static_weight.csv"
    WEIGHT = "_weight.csv"
    BASIC = ".csv"
    FINAL_RES = "final_results_"
    indexes = ['5', '10', '15']
    indexes = ['15']

    for i in indexes:
        static = get_test_df(os.path.join(res_dir, FINAL_RES + i + STATIC))
        static_weight = get_test_df(os.path.join(res_dir, FINAL_RES + i + STATIC_WEIGHT))
        weight = get_test_df(os.path.join(res_dir, FINAL_RES + i + WEIGHT))
        basic = get_test_df(os.path.join(res_dir, FINAL_RES + i + BASIC))
        expected = get_expected_df(os.path.join(expected_dir, "mock_" + i + "seq/reads/expected_res.csv"))

        static_compare = get_presented_data(compare(static, expected, i + "_static"))
        static_weight_compare = get_presented_data(compare(static_weight, expected, i + "_static_weight"))
        weight_compare = get_presented_data(compare(weight, expected, i + "_weight"))
        basic_compare = get_presented_data(compare(basic, expected, i + "_basic"))

        full_comparison = pd.merge(pd.merge(static_compare, static_weight_compare, on=['ref_e', 'expected frequency'],
                                            suffixes=(" s", " s+w")),
                                   pd.merge(weight_compare, basic_compare, on=['ref_e', 'expected frequency'],
                                            suffixes=(" w", "-")))
        suffixes = [" s+w", " w", " s", "-"]
        full_comparison = full_comparison[
            ["# of refs" + s for s in suffixes] + ['expected frequency'] + ['freq diff' + s for s in
                                                                            suffixes]]  # of match references

        fig = plt.figure(figsize=(10, 4), dpi=300)
        ax = fig.add_subplot(111, frame_on=False)  # no visible frame
        ax.xaxis.set_visible(False)  # hide the x axis
        ax.yaxis.set_visible(False)  # hide the y axis

        the_table = table(ax, full_comparison, loc='center')  # where df is your data frame
        the_table.set_fontsize(18)
        # the_table.scale(3, 3)

        res_name = i + '_emirge_smurf_full_compare.png'
        # plt.show()
        plt.tight_layout()
        plt.savefig(os.path.join('/home/vered/EMIRGE/EMIRGE-data/', res_name), bbox_inches='tight')

        plt.clf()
        fig.clear()


def main_compare_basic(expected_res_path, actual_res_path, working_dir, test_name, max_allowed_mismatch):
    test_df = get_test_df(actual_res_path)
    expected_df = get_expected_df(expected_res_path)

    compare(test_df, expected_df, working_dir, test_name, max_allowed_mismatch)



def get_command_line_arguments_parser():
    USAGE = \
        """usage: %prog EXPECTED_RES_PATH ACTUAL_RES_PATH WORKING_DIR [required_options] [options]

        Compare the expected results to the actual results of emirge_smurf
        """

    parser = OptionParser(USAGE)

    group_compare = OptionGroup(parser, "Optional flags",
                             "These flags are all optional.")

    group_compare.add_option("-r", dest="reference_id",
                          type="string", default=None,
                          help="reference id to compare")
    group_compare.add_option("--ra", dest="reference_id_all",
                             action="store_true", default=False,
                             help="compare all the references")
    group_compare.add_option("-s", dest="similarity",
                          action="store_true", default=False,
                          help="Test the similarity between the reads in the mixture")
    group_compare.add_option("-n", dest="test_name",
                             type="string", default="test",
                             help="The test name")
    group_compare.add_option("-m", dest="max_allowed_mismatch",
                             type="int", default=10,
                             help="The test name")

    parser.add_option_group(group_compare)
    return parser


def main(argv = sys.argv[1:]):
    define_logger(logging.DEBUG)
    parser = get_command_line_arguments_parser()

    # ACTUALLY PARSE ARGS
    (options, args) = parser.parse_args(argv)

    # minimal sanity checking of input
    if len(args) != 3:
        parser.error(
            "EXPECTED_RES_PATH ACTUAL_RES_PATH WORKING_DIR are required, and all options except DIR should have a flag associated with them (options without flags: %s)" % args)
        return

    expected_res_path = os.path.abspath(args[0])
    actual_res_path = os.path.abspath(args[1])
    working_dir = os.path.abspath(args[2])


    if options.similarity:
        test_mock_reads_for_similarity(expected_res_path)
    elif options.reference_id:
        compare_specific_reference(actual_res_path, expected_res_path, working_dir, options.test_name, options.max_allowed_mismatch, options.reference_id)
    elif options.reference_id_all:
        compare_specific_reference(actual_res_path, expected_res_path, working_dir, options.test_name, options.max_allowed_mismatch)
    else:
        main_compare_basic(expected_res_path, actual_res_path, working_dir, options.test_name, options.max_allowed_mismatch)



if __name__ == "__main__":
    main()



