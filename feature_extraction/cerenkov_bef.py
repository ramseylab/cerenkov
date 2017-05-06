import os
import pandas
import metadata_util
import tss_util
import dhs_util
import genome_seg_util
import gerp_util
import sys_tool
import argparse

# read the CERENKOV_DIR from the environment, but default to the directory where the script is located
CERENKOV_DIR = os.getenv('CERENKOV_DIR', '.')


class BefMaker:
    def __init__(self, rsid_path, artifact_names):
        self.rsid = rsid_path
        self.artifact = artifact_names

    def metadata(self):
        metadata_util.extract_metadata(src=self.rsid,
                                       dest_csv=self.artifact["metadata"],
                                       dest_bed=self.artifact["bed"])

    def tss_dist(self):
        tss_util.extract_min_tss_dist(src=self.rsid, dest=self.artifact["tssDist"])

    def dhs_score(self):
        dhs_util.extract_dhs_score(src=self.artifact["metadata"], dest=self.artifact["dhsScore"])

    def genome_seg_annot(self):
        genome_seg_util.extract_genome_seg_annot(src=self.rsid,
                                                 dest=self.artifact["genomeSeg"])

    def gerp(self):
        gerp_util.extract_gerp(src_bed=self.artifact["bed"], dest=self.artifact["gerp"])

    def dna_shape_and_gc_content(self):
        r_file = "getFlankingSequenceFeatures.R"

        _input = self.artifact["metadata"]
        _output = self.artifact["seq"]

        sys_tool.run_r_script(r_file, [_input, _output])

    def eqtl_pvalue(self):
        r_file = "ExtractEqtlPvalue_v2.R"

        eqtl_folder = sys_tool.find_directory("eqtl")

        _input = self.artifact["metadata"]
        _output = self.artifact["eqtlPvalue"]

        sys_tool.run_r_script(r_file, [eqtl_folder, _input, _output])

    def merge_artifacts(self):
        print("[cerenkov_bef] merging...")

        _metadata = pandas.read_csv(self.artifact["metadata"], sep='\t')
        _tss_dist = pandas.read_csv(self.artifact["tssDist"], sep='\t')
        _uniform_dhs = pandas.read_csv(self.artifact["dhsScore"], sep='\t')
        _eqtl_pvalue = pandas.read_csv(self.artifact["eqtlPvalue"], sep='\t')
        _sequence = pandas.read_csv(self.artifact["seq"], sep='\t')
        _chmmsw = pandas.read_csv(self.artifact["genomeSeg"], sep='\t')
        _gerp = pandas.read_csv(self.artifact["gerp"], sep='\t')

        _all = _metadata.merge(_tss_dist, on='name', how='left', left_index=True). \
            merge(_uniform_dhs.loc[:, ["name", "uniformDhsScore", "uniformDhsCount"]],
                  on='name', how='left', left_index=True). \
            merge(_eqtl_pvalue, on='name', how='left', left_index=True). \
            merge(_sequence, on='name', how='left', left_index=True). \
            merge(_chmmsw.loc[:, ["name", "ch1Name", "ch2Name", "ch3Name", "ch4Name", "ch5Name", "ch6Name",
                                  "sw1Name", "sw2Name", "sw3Name", "sw4Name", "sw5Name", "sw6Name"]],
                  on='name', how='left', left_index=True). \
            merge(_gerp, on='name', how='left', left_index=True)

        _all.to_csv(self.artifact["merged"], sep='\t', header=True, index=False)

    def make_artifacts(self):
        print("[cerenkov_bef] extracting metadata...")
        self.metadata()
        print("[cerenkov_bef] extracting minimum TSS distances...")
        self.tss_dist()
        print("[cerenkov_bef] extracting gerp scores...")
        self.gerp()
        print("[cerenkov_bef] extracting segment annotation...")
        self.genome_seg_annot()
        print("[cerenkov_bef] extracting eQTL p-values...")
        self.eqtl_pvalue()
        print("[cerenkov_bef] extracting DNAShape and GC content...")
        self.dna_shape_and_gc_content()
        print("[cerenkov_bef] extracting DNase I HS uniform peak scores...")
        self.dhs_score()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract BEF features.", allow_abbrev=False)

    parser.add_argument('-w', '--workspace', dest='workspace', type=str, required=True, help="the workspace path")
    parser.add_argument('-p', '--artifact-prefix', dest='artifact_prefix', type=str, required=True,
                        help="prefix of artifact files")
    parser.add_argument('-r', '--rsid', dest='rsid', type=str, required=True,
                        help="path of the rsid file")

    args = parser.parse_args()

    _ws_dir = CERENKOV_DIR + "/" + args.workspace  # Workspace directory
    _af_prefix = args.artifact_prefix
    _rsid_path = args.rsid

    print("[cerenkov_bef] workspace: {}; use rsid: {}; artifact prefix: {};".
          format(_ws_dir, _rsid_path, _af_prefix))

    _artifact_ext = dict(
        bed=".bed",
        metadata=".ucsc_metadata.txt",
        tssDist=".tss_dist.txt",
        dhsScore=".dhs_uni_pk.txt",
        genomeSeg=".chmmsw.txt",
        gerp=".gerp.txt",
        eqtlPvalue=".eqtl_pvalue.txt",
        seq=".sequence_features.txt",
        merged=".ALL.txt"
    )

    # prepend workspace path and prefix to these artifact names
    prelude = _ws_dir + '/' + _af_prefix
    _artifact_names = {k: prelude + v for k, v in _artifact_ext.items()}

    bef_maker = BefMaker(rsid_path=_rsid_path, artifact_names=_artifact_names)

    bef_maker.make_artifacts()
    bef_maker.merge_artifacts()
