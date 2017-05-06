# import sys
# from time import strftime
from io import StringIO
import pandas
from pybedtools import BedTool
from sys_tool import find_directory
from abstract_feature_util import AbstractFeatureUtil


# def __read_metadata(src):
#     """
#         designed to read `csnps.ucsc_metadata.txt` and `rsnps.ucsc_metadata.txt`, which have headers inside
#     """
#
#     try:
#         # col.names = c("name",	"chrom", "chromStart", "chromEnd", rep("NULL", 16))
#         # colClasses = c("character", "factor", "numeric", "numeric", rep("NULL", 16))
#         metadata = pandas.read_csv(src, sep='\t')
#     # names=["name", "chrom", "chromStart", "chromEnd"]
#     # dtype={"name" : str, "chrom" : str, "chromStart" : np.int, "chromEnd" : np.int}
#     except IOError:
#         sys.exit("Could not open file " + src)
#
#     return metadata[["name", "chrom", "chromStart", "chromEnd"]]
#
#
# def __read_dhs_uniform_peak(src=None):
#     """
#         designed to read `UniformDnaseIHS`, which has NO header inside
#     """
#
#     if src is None:
#         # use the default one
#         directory = find_directory("dhs")
#         src = directory + "/UniformDnaseIHS"
#
#     try:
#         # col.names = c("chrom", "chromStart", "chromEnd", "NULL", "score", "NULL")
#         dhs_uni_pk = pandas.read_csv(src, sep='\t', header=None, usecols=[0, 1, 2, 4],
#                                      names=["chrom", "chromStart", "chromEnd", "score"])
#     # dtype={"chrom" : str, "chromStart" : np.int, "chromEnd" : np.int, "score" : np.int}
#     except IOError:
#         sys.exit("Could not open file " + src)
#
#     return dhs_uni_pk


# def __amend_chrom_start(metadata):
#     """
#         Some cSNPs' lengths are not 1-bp. E.g a deletion can be 0-bp while a insertion can be 12-bp.
#         Add a new column `amendedChromStart` which equals floor((chromStart + chromEnd)/2), which would not affect
#             the 1-bp-long SNPs
#     """
#
#     metadata.loc[:, "amendedChromStart"] = numpy.floor(
#         metadata.loc[:, "chromStart"].add(metadata.loc[:, "chromEnd"]).divide(2))


# def __extract_dhs_score(metadata, dhs_uni_pk):
#     """
#         Match every SNP to DHSs by `chrom`, `chromStart`, `chromEnd`;
#
#         E.g. A matched group for SNP `rs1` is like below
#
#             name	score
#             rs1		[1,2,3]
#
#         metadata - a dataframe of SNPs' `name`s, `chrom`s, `chromStart`s and `chromEnd`s
#         dhs_uni_pk - a dataframe of DHSs' `chrom`s, `chromStart`s, `chromEnd`s and `score`s
#     """
#     dhs_score = group_match(metadata, "chrom", "name", "amendedChromStart", dhs_uni_pk, "chrom", "chromStart",
#                             "chromEnd", "score")
#     dhs_score.loc[:, "uniformDhsScore"] = dhs_score.loc[:, "score"].apply(sum)
#     dhs_score.loc[:, "uniformDhsCount"] = dhs_score.loc[:, "score"].apply(len)
#
#     return dhs_score


# def __write_dhs_score(dhs_score, dest):
#     """
#         Write the matched DHS score to a file; columns are separated by tabs ("\t")
#
#         dhs_score - a dataframe that contains the matching result
#         filename - the name of the output file
#     """
#
#     dhs_score.to_csv(dest, sep='\t', header=True, columns=["name", "uniformDhsScore", "uniformDhsCount"])


# def extract_dhs_score(src, dest):
#
#     print(strftime("[%Y-%m-%d %H:%M:%S]") + "extract_dhs_score: reading UCSC metadata...")
#     metadata = __read_metadata(src)
#     __amend_chrom_start(metadata)
#
#     print(strftime("[%Y-%m-%d %H:%M:%S]") + "extract_dhs_score: reading DHS unique peaks...")
#     dhs_uni_peaks = __read_dhs_uniform_peak()
#
#     print(strftime("[%Y-%m-%d %H:%M:%S]") + "extract_dhs_score: matching...")
#     dhs_score = __extract_dhs_score(metadata, dhs_uni_peaks)
#     __write_dhs_score(dhs_score, dest)
#
#     print(strftime("[%Y-%m-%d %H:%M:%S]") + "extract_dhs_score: complete!")

class DhsUtil(AbstractFeatureUtil):
    def __init__(self):
        self._dhs_bed_obj = None

    @property
    def dhs_bed_obj(self):
        # Lazy Initialization
        if self._dhs_bed_obj is None:
            dhs_bed_fn = "{dir}/{fn}".format(dir=self.src_data_dir, fn=self.src_data_fn)
            dhs_bed_obj = BedTool(dhs_bed_fn)

            self._dhs_bed_obj = dhs_bed_obj

        return self._dhs_bed_obj

    def get_feat(self, _input):
        """
        :param _input: the SNP data frame in BED format
        :return:
        """
        snp_dfm = _input

        # rsnp_bed_obj = BedTool(rsnp_bed_fn)
        snp_bed_obj = BedTool(snp_dfm.to_string(index=False, header=False, index_names=False), from_string=True)

        # loj == left outer join.
        # intersect(`loj=True`) will return -1 dhsScore if no match
        # Do not use loj here; will fillna(0) later
        intx = snp_bed_obj.intersect(self.dhs_bed_obj, wb=True)
        intx_dfm = pandas.read_table(StringIO(str(intx)), header=None,
                                     names=['snpChrom', 'snpChromStart', 'snpChromEnd', 'snpName',
                                            'dhsChrom', 'dhsChromStart', 'dhsChromEnd', 'dhsName', 'dhsScore',
                                            'dhsStrand'])
        dhs_score_dfm = intx_dfm.groupby(by='snpName')['dhsScore'].agg({"uniformDhsScore": sum,
                                                                        "uniformDhsCount": len}).reset_index()

        snp_dhs = snp_dfm.merge(dhs_score_dfm, how='left', left_on='name', right_on='snpName')
        snp_dhs.fillna(0, inplace=True)

        return snp_dhs.loc[:, ["name", "uniformDhsScore", "uniformDhsCount"]]

    def save_temp(self, _result):
        _result.to_csv(self.temp_dest, sep='\t', header=True, index=False,
                       columns=["name", "uniformDhsScore", "uniformDhsCount"])


if __name__ == '__main__':
    dhs_util = DhsUtil()
    dhs_util.src_data_dir = find_directory("dhs")
    dhs_util.src_data_fn = "UniformDnaseIHS"

    rsnp_bed_fn = "{}/RSNP_50kb.bed".format(find_directory('dhs'))
    rsnp_dfm = pandas.read_table(rsnp_bed_fn, header=None,
                                 names=['chrom', 'chromStart', 'chromEnd', 'name'])

    dhs_util.temp_dest = 'foo_rsnp.txt'
    dhs_util.extract(_input=rsnp_dfm)

    csnp_bed_fn = "{}/CSNP_50kb.bed".format(find_directory('dhs'))
    csnp_dfm = pandas.read_table(csnp_bed_fn, header=None,
                                 names=['chrom', 'chromStart', 'chromEnd', 'name'])

    dhs_util.temp_dest = 'foo_csnp.txt'
    dhs_util.extract(_input=csnp_dfm)


    # dhs_bed_fn = find_directory("dhs") + "/UniformDnaseIHS"
    # rsnp_bed_fn = "{}/RSNP_50kb.bed".format(find_directory('dhs'))
    # rsnp_dfm = pandas.read_table(rsnp_bed_fn, header=None,
    #                              names=['chrom', 'chromStart', 'chromEnd', 'name'])
    #
    # dhs_bed_obj = BedTool(dhs_bed_fn)
    # # rsnp_bed_obj = BedTool(rsnp_bed_fn)
    # rsnp_bed_obj = BedTool(rsnp_dfm.to_string(index=False, header=False, index_names=False), from_string=True)
    #
    # # loj == left outer join.
    # # intersect(`loj=True`) will return -1 dhsScore if no match
    # # Do not use loj here; will fillna(0) later
    # intx = rsnp_bed_obj.intersect(dhs_bed_obj, wb=True)
    # intx_dfm = pandas.read_table(StringIO(str(intx)), header=None,
    #                              names=['snpChrom', 'snpChromStart', 'snpChromEnd', 'snpName',
    #                                     'dhsChrom', 'dhsChromStart', 'dhsChromEnd', 'dhsName', 'dhsScore', 'dhsStrand'])
    # dhs_score_dfm = intx_dfm.groupby(by='snpName')['dhsScore'].agg({"uniformDhsScore": sum,
    #                                                                 "uniformDhsCount": len}).reset_index()
    #
    # rsnp_dhs = rsnp_dfm.merge(dhs_score_dfm, how='left', left_on='name', right_on='snpName')
    # rsnp_dhs.fillna(0, inplace=True)
    #
    # print(rsnp_dhs)

    # rsnp_dfm = pandas.read_table("{}/RSNP_50kb.bed".format(find_directory('dhs')), header=None,
    #                              names=['chrom', 'chromStart', 'chromEnd', 'name'])
    #
    # csnp_dfm = pandas.read_table("{}/CSNP_50kb.bed".format(find_directory('dhs')), header=None,
    #                              names=['chrom', 'chromStart', 'chromEnd', 'name'])

    # __amend_chrom_start(rsnp_dfm)
    # __amend_chrom_start(csnp_dfm)

    # print(rsnp_dfm.loc[rsnp_dfm['amendedChromStart'] != rsnp_dfm['chromStart']].index.tolist())
    # print(csnp_dfm.loc[csnp_dfm['amendedChromStart'] != csnp_dfm['chromStart']].index.tolist())

    # print(all(rsnp_dfm['amendedChromStart'] == rsnp_dfm['chromStart']))
    # print(all(csnp_dfm['amendedChromStart'] == csnp_dfm['chromStart']))


