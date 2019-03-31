

import numpy as np
from emirge_headers import *
from emirge_utills import *
import os
import matplotlib.pyplot as plt
import pandas as pd
# from pandas.tools.plotting import table
from optparse import OptionParser, OptionGroup
import glob
import sys


VALIDATION_THRESHOLD_THRESHOLD = 0.00001
EXPECTED_RES_FILE_NAME = "expected_res.csv"
ACTUAL_RES_FILE_NAME = "emirge_smurf_WFalseSTrue.csv"

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

    df = df[Header.all]
    validate_priors(df)

    return df



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
    df = df[Header.all]
    validate_priors(df)

    return df


def is_expected_contain_actual(merged, actual_df, expected_id, actual_id):
    """
    :param merged:
    :param actual_df:
    :param expected_id:
    :param actual_id:
    :return:
    """
    common_regions = len(merged[(merged[Header.ref_id + '_e'] == expected_id) & (merged[Header.ref_id + '_a'] == actual_id)])
    actual_length = len(actual_df[actual_df[Header.ref_id] == actual_id])
    return common_regions == actual_length


def is_actual_equals_expected(merged, expected_df, expected_id, actual_id):
    """
    :param merged:
    :param actual_df:
    :param expected_id:
    :param actual_id:
    :return:
    """
    common_regions = len(merged[(merged[Header.ref_id + '_e'] == expected_id) & (merged[Header.ref_id + '_a'] == actual_id)])
    expected_length = len(expected_df[expected_df[Header.ref_id] == expected_id])
    return common_regions == expected_length



def calc_table(expected_path, actual_path):
    expected_df = get_expected_df(expected_path)
    actual_df = get_test_df(actual_path)
    unique_actual_ref = actual_df[Header.ref_id].unique().tolist()
    unique_expected_ref = expected_df[Header.ref_id].unique().tolist()
    score_table_contain = pd.DataFrame(columns=unique_actual_ref, index=unique_expected_ref)
    score_table_eqaul = pd.DataFrame(columns=unique_actual_ref, index=unique_expected_ref)
    merged = expected_df.merge(actual_df, on=[Header.region, Header.sequence],how="inner", suffixes=('_e', '_a'))
    common_regions_df = merged.groupby([Header.ref_id+"_e", Header.ref_id + '_a'])[Header.region].count().reset_index()
    common_regsion_df = common_regions_df.merge(
        merged[[Header.ref_id+"_e", Header.ref_id + '_a', Header.weight + '_e', Header.weight + '_a']],
        on=[Header.ref_id+"_e", Header.ref_id + '_a'],
        how="inner")
    common_regsion_df['is_contains'] = (common_regsion_df[Header.weight + '_a'] == common_regsion_df[Header.region])
    common_regsion_df['is_equals'] = (common_regsion_df[Header.weight + '_e'] == common_regsion_df[Header.region])

    score_table_eqaul = pd.pivot_table(common_regsion_df,
                                       values='is_equals',
                                       index=[Header.ref_id + '_e'],
                                       columns=[Header.ref_id + '_a'],
                                       aggfunc=np.sum,
                                       fill_value=0)
    score_table_contain = pd.pivot_table(common_regsion_df,
                                         values='is_contains',
                                         index=[Header.ref_id + '_e'],
                                         columns=[Header.ref_id + '_a'],
                                         aggfunc=np.sum,
                                         fill_value=0)
    # logging.debug("expected results length={}, actual={}".format(len(expected_df), len(actual_df)))
    # for expected_ref in unique_expected_ref:
    #     for actual_ref in unique_actual_ref:
    #         score_table_contain.set_value(expected_ref,
    #                                       actual_ref,
    #                                       is_expected_contain_actual(merged, actual_df, expected_ref, actual_ref))
    #         score_table_eqaul.set_value(expected_ref,
    #                                     actual_ref,
    #                                     is_actual_equals_expected(merged, expected_df, expected_ref, actual_ref))

    actual_priors = actual_df.drop_duplicates(Header.ref_id)[Header.prior]
    actual_priors.index = actual_df.drop_duplicates(Header.ref_id)[Header.ref_id]
    expected_priors = expected_df.drop_duplicates(Header.ref_id)[Header.prior]
    expected_priors.index = expected_df.drop_duplicates(Header.ref_id)[Header.ref_id]
    return score_table_contain, score_table_eqaul, actual_priors, expected_priors


def calc_recall(score_table, expected_priors):
    """
    Recall  = sum(I_Ei * f_Ei)
    :param score_table:
    :param expected_priors:
    :return: recall
    """
    expected_indicator = (score_table > 0).any(axis=1)
    recall = sum(expected_indicator.mul(expected_priors))
    logging.info( "recall = {}".format(recall))
    # logging.info(score_table[expected_indicator == False].index)
    return recall


def calc_precision(score_table, actual_priors):
    """
    Precision = sum(I_Ai * f_Ai)
    :param score_table:
    :param actual_priors:
    :return: precision
    """
    actual_indicator = (score_table > 0).any()
    precision = sum(actual_indicator.mul(actual_priors))
    logging.info( "precision = {}".format(precision))
    return precision

    # logging.info( score_table.transpose()[actual_indicator == False].index)


def get_cmd_arguments(argv = sys.argv[1:]):
    USAGE = \
        """usage: %prog WORKING_DIR [required_options] [options]

        Compare the expected results to the actual results of emirge_smurf
        """

    parser = OptionParser(USAGE)

    # group_compare = OptionGroup(parser, "Optional flags",
    #                          "These flags are all optional.")
    #
    # group_compare.add_option("-r", dest="reference_id",
    #                       type="string", default=None,
    #                       help="reference id to compare")
    #
    # parser.add_option_group(group_compare)

    (options, args) = parser.parse_args(argv)

    if len(args) != 1:
        parser.error(
            "WORKING_DIR is required, and all options except DIR should have a flag associated with them (options without flags: %s)" % args)
        return

    return options, args


class Statistics(object):
    def __init__(self, data):
        self.mean = np.mean(data)
        self.std = np.std(data)
        self.min = min(data)
        self.max = max(data)


class SingleExperimentData(object):
    def __init__(self, changed_bacterias, changed_bases):
        self.changed_bacterias=changed_bacterias
        self.changed_bases=changed_bases
        self.precisions_similarity=[]
        self.recall_similarity = []
        self.precisions_contains = []
        self.recall_contains=[]

    def add_data(self, precisions_similarity, recall_similarity, precisions_contains, recall_contains):
        self.precisions_similarity.append(precisions_similarity)
        self.precisions_contains.append(precisions_contains)
        self.recall_contains.append(recall_contains)
        self.recall_similarity.append(recall_similarity)

    @staticmethod
    def get_exp_key(changed_bacterias, changed_bases):
        return "bctr{}bases{}".format(changed_bacterias, changed_bases)


    def get_my_exp_key(self):
        return SingleExperimentData.get_exp_key(self.changed_bacterias, self.changed_bases)


    def is_changed_bacterias(self, changed_bacterias):
        if self.changed_bacterias == changed_bacterias:
            return True
        return False

    def is_changed_bases(self, changed_bases):
        if self.changed_bases == changed_bases:
            return True
        return False


class ExperimentsData(object):
    def __init__(self):
        self.data={}

    def add_experiment(self, experiment):
        self.data.update({experiment.get_my_exp_key(): experiment})

    def add_data(self, precisions_similarity, recall_similarity, precisions_contains, recall_contains, changed_bacterias, changed_bases):
        if precisions_similarity is not None and precisions_similarity is not None:
            self.data[SingleExperimentData.get_exp_key(changed_bacterias, changed_bases)].\
                add_data(precisions_similarity, recall_similarity, precisions_contains, recall_contains)

    def get_experiments(self, changed_bacterias):
        res = []
        for exp in self.data.values():
            if exp.is_changed_bacterias(changed_bacterias):
                res.append(exp)
        return res



def main():
    define_logger(logging.DEBUG)

    # ACTUALLY PARSE ARGS
    (options, args) = get_cmd_arguments()

    working_dir = os.path.abspath(args[0])

    os.chdir(working_dir)
    files = glob.glob("test*")
    # for file in files:
    os.chdir(working_dir)
    # os.chdir(file)

    experiments=ExperimentsData()
    experiments.add_experiment(SingleExperimentData(1, 1))
    experiments.add_experiment(SingleExperimentData(1, 10))
    experiments.add_experiment(SingleExperimentData(1, 20))
    experiments.add_experiment(SingleExperimentData(3, 1))
    experiments.add_experiment(SingleExperimentData(3, 10))
    experiments.add_experiment(SingleExperimentData(3, 20))

    for test_dir in glob.glob("exp*"):
        try:
            print "test_dir = {}".format(test_dir)
            test_id = int(test_dir.split("exp")[1])
            changed_bacterias=0
            changed_bases=0
            if test_id <= 300:
                changed_bacterias=1
                if test_id<=100:
                    changed_bases=1
                elif test_id<=200:
                    changed_bases=10
                else:
                    changed_bases=20
            else:
                changed_bacterias=3
                if test_id <= 400:
                    changed_bases = 1
                elif test_id <= 500:
                    changed_bases = 10
                else:
                    changed_bases = 20

            actual_path = os.path.join(test_dir, "test_{}".format(test_id), ACTUAL_RES_FILE_NAME)
            expected_path = os.path.join(test_dir, "test_{}".format(test_id), EXPECTED_RES_FILE_NAME)

            contain_table, similarity_table, actual_priors, expected_priors = calc_table(expected_path, actual_path)

            experiments.add_data(calc_precision(similarity_table, actual_priors),
                                 calc_recall(similarity_table, expected_priors),
                                 calc_precision(contain_table, actual_priors),
                                 calc_recall(contain_table, expected_priors),
                                 changed_bacterias, changed_bases)
        except Exception as ex:
            print ex


    bacterias1_exp=experiments.get_experiments(1)
    bacterias3_exp=experiments.get_experiments(3)
    bacterias1_exp.sort(key=lambda x: x.changed_bases)
    bacterias3_exp.sort(key=lambda x: x.changed_bases)

    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(7, 9))
    x_ix = range(len(bacterias1_exp))
    x_ticklabels = [exp.changed_bases for exp in bacterias1_exp]

    print [Statistics(exp.precisions_similarity).mean for exp in bacterias3_exp]

    axes[0].errorbar(x_ix, [Statistics(exp.precisions_similarity).mean for exp in bacterias3_exp],
                     yerr=[Statistics(exp.precisions_similarity).std for exp in bacterias3_exp],
                     fmt = 'o', color='g', label='3 changed bacterias')
    axes[0].errorbar(x_ix, [Statistics(exp.precisions_similarity).mean for exp in bacterias1_exp],
                     yerr=[Statistics(exp.precisions_similarity).std for exp in bacterias1_exp],
                     fmt='o', color='b', label='1 changed bacterias')
    axes[0].set_xticks(x_ix)
    axes[0].set_xticklabels(x_ticklabels)
    axes[0].set_title('Evaluation - test results bacteria is SIMILAR to the actual one')
    axes[0].set_xlabel("number of bases changed in each bacteria")
    axes[0].set_ylabel("Precision")
    legend = axes[0].legend()

    axes[1].errorbar(x_ix, [Statistics(exp.recall_similarity).mean for exp in bacterias3_exp],
                     yerr=[Statistics(exp.recall_similarity).std for exp in bacterias3_exp],
                     fmt='o', color='g', label='3 changed bacterias')


    axes[1].errorbar(x_ix, [Statistics(exp.recall_similarity).mean for exp in bacterias1_exp],
                     fmt='o', color='b', yerr=[Statistics(exp.recall_similarity).std for exp in bacterias1_exp],
                     label='1 changed bacterias')
    axes[1].set_xticks(x_ix)
    axes[1].set_xticklabels(x_ticklabels)
    axes[1].set_title('Evaluation - test results bacteria is SIMILAR to the actual one')
    axes[1].set_xlabel("number of bases changed in each bacteria")
    axes[1].set_ylabel("Recall")
    legend = axes[0].legend()

    plt.tight_layout()

    plt.show()
    plt.savefig(os.path.join(working_dir, "evaluation_similarity_test.png"), bbox_inches='tight')



    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(7, 9))
    x_ix = range(len(bacterias1_exp))
    x_ticklabels = [exp.changed_bases for exp in bacterias1_exp]

    axes[0].errorbar(x_ix, [Statistics(exp.precisions_contains).mean for exp in bacterias3_exp],
                     yerr=[Statistics(exp.precisions_contains).std for exp in bacterias3_exp],
                     fmt='o', color='g', label='3 changed bacterias')
    axes[0].errorbar(x_ix, [Statistics(exp.precisions_contains).mean for exp in bacterias1_exp],
                     yerr=[Statistics(exp.precisions_contains).std for exp in bacterias1_exp],
                     fmt='o', color='b', label='1 changed bacterias')
    axes[0].set_xticks(x_ix)
    axes[0].set_xticklabels(x_ticklabels)
    axes[0].set_title('Evaluation - test results bacteria is CONAINS in the actual one')
    axes[0].set_xlabel("number of bases changed in each bacteria")
    axes[0].set_ylabel("Precision")
    legend = axes[0].legend()

    axes[1].errorbar(x_ix, [Statistics(exp.recall_contains).mean for exp in bacterias3_exp],
                     yerr=[Statistics(exp.recall_contains).std for exp in bacterias3_exp],
                     fmt='o', color='g', label='3 changed bacterias')

    axes[1].errorbar(x_ix, [Statistics(exp.recall_contains).mean for exp in bacterias1_exp],
                     yerr=[Statistics(exp.recall_contains).std for exp in bacterias1_exp],
                     fmt='o', color='b', label='1 changed bacterias')
    axes[1].set_xticks(x_ix)
    axes[1].set_xticklabels(x_ticklabels)
    axes[1].set_title('Evaluation - test results bacteria is CONAINS in the actual one')
    axes[1].set_xlabel("number of bases changed in each bacteria")
    axes[1].set_ylabel("Recall")
    legend = axes[0].legend()

    plt.tight_layout()

    plt.show()
    plt.savefig(os.path.join(working_dir, "evaluation_conains_test.png"), bbox_inches='tight')


if __name__ == "__main__":
    main()



