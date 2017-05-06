"""
Created on Aug 3, 2016

@author: ramseylab
"""

import pandas
from sqlalchemy import create_engine
from sqlalchemy.engine.url import URL
from sqlalchemy.pool import NullPool
import chrom_tool as ct
import numbers
import ast
import re
from allele_tool import parse_gb_alleles_col, parse_gb_allele_freqs_col
from genome_browser_tool import bin_from_range


class GenomeBrowserClient:
    # Make sure you have created such a user in MySQL
    # For MySQL 5.6
    #   GRANT SELECT PRIVILEGES ON hg19.* To 'bud'@'localhost' IDENTIFIED BY 'earth';
    #   GRANT SELECT PRIVILEGES ON hgmd_pro.* To 'bud'@'localhost' IDENTIFIED BY 'earth';
    # For MySQL 5.7
    #   CREATE USER 'bud'@'localhost' IDENTIFIED BY 'earth';
    #   GRANT SELECT ON hg19.* TO 'bud'@'localhost';
    #   GRANT SELECT ON hgmd_pro.* TO 'bud'@'localhost';
    #   FLUSH PRIVILEGES;
    __db_url = dict(
        local_hg19=dict(
            drivername='mysql+pymysql',
            host='localhost',
            port='3306',
            username='bud',
            password='earth',
            database='hg19',
            query={'charset': 'utf8'}
        ),

        remote_hg19=dict(
            drivername='mysql+pymysql',
            host='genome-mysql.cse.ucsc.edu',
            port='3306',
            username='genome',
            password='',
            database='hg19',
            query={'charset': 'utf8'}
        ),
    )

    def __init__(self, config_key):
        # db = create_engine('mysql://bud:earth@localhost:3306/hg19') # require module `MySQLdb`
        #   default dialect is 'mysql+mysql-python'
        #   `MySQLdb` is a fork of MySQL-python with added support for Python 3
        #   See http://docs.sqlalchemy.org/en/latest/core/engines.html#mysql

        # db = create_engine('mysql+pymysql://bud:earth@localhost:3306/hg19') # require module `PyMySQL`

        # For `poolclass`, see http://stackoverflow.com/a/8705750

        self.db = create_engine(URL(**GenomeBrowserClient.__db_url[config_key]), poolclass=NullPool)
        self.conn = self.db.connect()

        # Subtraction between integer values, where one is of type UNSIGNED, produces an unsigned result by default.
        #   If the difference is negative, an error results because it must be unsigned.
        # Coordinates are unsigned int. We'll use subtraction between coordinates to get TSS distances,
        #   so we must enable this mode.
        self.conn.execute("SET sql_mode = 'NO_UNSIGNED_SUBTRACTION'")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.conn.close()
        self.db.dispose()

    def select_tss_dist(self, chrom, chrom_start, chrom_end):
        # CASE g.strand
        #     WHEN '+' THEN s.chromStart - g.txStart
        #     WHEN '-' THEN g.txEnd - s.chromStart
        #         # Yao edited so that 'upstream' is always a negative distance, 2016.04.20
        #     # WHEN '+' THEN g.txStart - s.chromStart
        #     # WHEN '-' THEN s.chromStart - g.txEnd
        #         # Steve edited so that 'upstream' is always a positive distance, 2015.10.08
        #         # commented out by Yao, 2016.04.20
        #     # WHEN '-' THEN g.txEnd - s.chromStart  # Satpreet's original code 2015.08.31
        query_template = '''
                         SELECT
                             CASE strand
                                 WHEN '+' THEN {chromStart} - txStart
                                 WHEN '-' THEN txEnd - {chromStart}
                             END as tssDistance
                         FROM
                             ensGene
                         WHERE
                             chrom = "{chrom}" {binCondtion}
                         ORDER BY
                             abs(tssDistance)
                         '''
        _bin = bin_from_range(chrom_start, chrom_end)
        bin_condition = "AND bin between {} and {}".format(_bin - 1, _bin + 1)

        query = query_template.format(chromStart=chrom_start, chrom=chrom, binCondtion=bin_condition)
        result = self.conn.execute(query).fetchall()

        # If no TSS matched in nearby bins, search across all bins
        # I.e. calculate distances to all TSS on the SNP's chromosome
        if len(result) == 0:
            query = query_template.format(chromStart=chrom_start, chrom=chrom, binCondtion="")
            result = self.conn.execute(query).fetchall()

        # `result` is a list of one-element tuples. Change it into a list of integers
        return [r[0] for r in result]

    def compute_tss_dist(self, rsid, adjacent_bins=1):

        if not isinstance(adjacent_bins, numbers.Integral):
            raise ValueError("Required: 'adjacent_bins' must be an int. Actually: adjacent_bins == {}".
                             format(adjacent_bins))
        else:
            if adjacent_bins > 0:
                # search within a bin window
                bin_condition = "g.bin between s.bin - {} and s.bin + {}".format(adjacent_bins, adjacent_bins)
            elif adjacent_bins == 0:
                # search within the same bin
                bin_condition = "g.bin = s.bin"
            elif adjacent_bins == -1:
                # search across all bins
                bin_condition = ""
            else:
                raise ValueError("Invalid 'adjacent_bins' value: " + adjacent_bins)

        # if binCondition is not empty, concat conditions with AND
        condition_op = " AND" if bin_condition else ""

        snps = (", ".join("'" + x + "'" for x in rsid))
        chroms = (", ".join("'" + x + "'" for x in ct.REGULAR_CHR))

        # query = '''
        #         SELECT
        #             s.name, s.chrom, s.chromStart, s.chromEnd,
        #             g.name as tssGene, g.txStart, g.txEnd, g.strand,
        #         CASE g.strand
        #             WHEN '+' THEN s.chromStart - g.txStart
        #             WHEN '-' THEN g.txEnd - s.chromStart
        #                 # Yao edited so that 'upstream' is always a negative distance, 2016.04.20
        #             # WHEN '+' THEN g.txStart - s.chromStart
        #             # WHEN '-' THEN s.chromStart - g.txEnd   # Steve edited so that 'upstream' is always a positive distance, 2015.10.08
        #                 # commented out by Yao Yao, 2016.04.20
        #             # WHEN '-' THEN g.txEnd - s.chromStart  # Satpreet's original code 2015.08.31
        #         END as tssDistance
        #         FROM
        #             snp146 s
        #         LEFT OUTER JOIN
        #             ensGene g
        #         ON
        #             g.bin = s.bin # Speeds up JOINs # Satpreet's original code 2015.08.31
        #             # g.bin between s.bin - 2 and s.bin + 2 # Yao edited so that no NA values output, 2016.08.03
        #             AND g.chrom = s.chrom
        #         WHERE
        #             s.name IN  ( 'rs9264942', 'rs9267551', 'rs9277535', 'rs9282699' )
        #         ORDER BY
        #             name, abs(tssDistance)
        #         '''

        query = '''
                SELECT
                    s.name, s.chrom, s.chromStart, s.chromEnd,
                    g.name as tssGene, g.txStart, g.txEnd, g.strand,
                CASE g.strand
                    WHEN '+' THEN s.chromStart - g.txStart
                    WHEN '-' THEN g.txEnd - s.chromStart
                END as tssDistance
                FROM
                    snp146 s
                    INNER JOIN
                        ensGene g
                    ON
                ''' + bin_condition + condition_op + '''
                        g.chrom = s.chrom
                WHERE
                    s.name IN ( ''' + snps + ''')
                    AND s.chrom IN ( ''' + chroms + ''')
                ORDER BY
                    name, abs(tssDistance)
                '''

        rows = self.conn.execute(query)

        df = pandas.DataFrame(rows.fetchall())
        df.columns = rows.keys()

        return df

    def identify_genome_seg(self, rsid):
        snps = (", ".join("'" + x + "'" for x in rsid))
        chroms = (", ".join("'" + x + "'" for x in ct.REGULAR_CHR))

        query = '''
                SELECT s.name, s.chrom, s.chromStart, s.chromEnd,
                    SUBSTRING_INDEX(GROUP_CONCAT(DISTINCT ch1.name), ',', 1) as ch1Name,
                    SUBSTRING_INDEX(GROUP_CONCAT(DISTINCT ch2.name), ',', 1) as ch2Name,
                    SUBSTRING_INDEX(GROUP_CONCAT(DISTINCT ch3.name), ',', 1) as ch3Name,
                    SUBSTRING_INDEX(GROUP_CONCAT(DISTINCT ch4.name), ',', 1) as ch4Name,
                    SUBSTRING_INDEX(GROUP_CONCAT(DISTINCT ch5.name), ',', 1) as ch5Name,
                    SUBSTRING_INDEX(GROUP_CONCAT(DISTINCT ch6.name), ',', 1) as ch6Name,
                    SUBSTRING_INDEX(GROUP_CONCAT(DISTINCT sw1.name), ',', 1) as sw1Name,
                    SUBSTRING_INDEX(GROUP_CONCAT(DISTINCT sw2.name), ',', 1) as sw2Name,
                    SUBSTRING_INDEX(GROUP_CONCAT(DISTINCT sw3.name), ',', 1) as sw3Name,
                    SUBSTRING_INDEX(GROUP_CONCAT(DISTINCT sw4.name), ',', 1) as sw4Name,
                    SUBSTRING_INDEX(GROUP_CONCAT(DISTINCT sw5.name), ',', 1) as sw5Name,
                    SUBSTRING_INDEX(GROUP_CONCAT(DISTINCT sw6.name), ',', 1) as sw6Name
                FROM
                    snp146 s
                    LEFT OUTER JOIN
                        wgEncodeAwgSegmentationChromhmmGm12878 ch1
                    ON
                        s.bin = ch1.bin
                        AND s.chromStart BETWEEN ch1.chromStart AND (ch1.chromEnd - 1)
                        AND s.chrom = ch1.chrom
                    LEFT OUTER JOIN
                        wgEncodeAwgSegmentationChromhmmH1hesc ch2
                    ON
                        s.bin = ch2.bin
                        AND s.chromStart BETWEEN ch2.chromStart AND (ch2.chromEnd - 1)
                        AND s.chrom = ch2.chrom
                    LEFT OUTER JOIN
                        wgEncodeAwgSegmentationChromhmmHelas3 ch3
                    ON
                        s.bin = ch3.bin
                        AND s.chromStart BETWEEN ch3.chromStart AND (ch3.chromEnd - 1)
                        AND s.chrom = ch3.chrom
                    LEFT OUTER JOIN
                        wgEncodeAwgSegmentationChromhmmHepg2 ch4
                    ON
                        s.bin = ch4.bin
                        AND s.chromStart BETWEEN ch4.chromStart AND (ch4.chromEnd - 1)
                        AND s.chrom = ch4.chrom
                    LEFT OUTER JOIN
                        wgEncodeAwgSegmentationChromhmmHuvec ch5
                    ON
                        s.bin = ch5.bin
                        AND s.chromStart BETWEEN ch5.chromStart AND (ch5.chromEnd - 1)
                        AND s.chrom = ch5.chrom
                    LEFT OUTER JOIN
                        wgEncodeAwgSegmentationChromhmmK562 ch6
                    ON
                        s.bin = ch6.bin
                        AND s.chromStart BETWEEN ch6.chromStart AND (ch6.chromEnd - 1)
                        AND s.chrom = ch6.chrom
                    LEFT OUTER JOIN
                        wgEncodeAwgSegmentationSegwayGm12878 sw1
                    ON
                        s.bin = sw1.bin
                        AND s.chromStart BETWEEN sw1.chromStart AND (sw1.chromEnd - 1)
                        AND s.chrom = sw1.chrom
                    LEFT OUTER JOIN
                        wgEncodeAwgSegmentationSegwayH1hesc sw2
                    ON
                        s.bin = sw2.bin
                        AND s.chromStart BETWEEN sw2.chromStart AND (sw2.chromEnd - 1)
                        AND s.chrom = sw2.chrom
                    LEFT OUTER JOIN
                        wgEncodeAwgSegmentationSegwayHelas3 sw3
                    ON
                        s.bin = sw3.bin
                        AND s.chromStart BETWEEN sw3.chromStart AND (sw3.chromEnd - 1)
                        AND s.chrom = sw3.chrom
                    LEFT OUTER JOIN
                        wgEncodeAwgSegmentationSegwayHepg2 sw4
                    ON
                        s.bin = sw4.bin
                        AND s.chromStart BETWEEN sw4.chromStart AND (sw4.chromEnd - 1)
                        AND s.chrom = sw4.chrom
                    LEFT OUTER JOIN
                        wgEncodeAwgSegmentationSegwayHuvec sw5
                    ON
                        s.bin = sw5.bin
                        AND s.chromStart BETWEEN sw5.chromStart AND (sw5.chromEnd - 1)
                        AND s.chrom = sw5.chrom
                    LEFT OUTER JOIN
                        wgEncodeAwgSegmentationSegwayK562 sw6
                    ON
                        s.bin = sw6.bin
                        AND s.chromStart BETWEEN sw6.chromStart AND (sw6.chromEnd - 1)
                        AND s.chrom = sw6.chrom
                    WHERE
                        s.name IN ( ''' + snps + ''')
                        AND s.chrom IN ( ''' + chroms + ''')
                    GROUP BY s.name, s.chrom, s.chromStart, s.chromEnd
                ORDER BY s.name
                '''

        rows = self.conn.execute(query)

        df = pandas.DataFrame(rows.fetchall())
        df.columns = rows.keys()

        return df

    def fetch_coord(self, rsid):
        snps = ", ".join("'{}'".format(x) for x in rsid)
        chroms = ", ".join("'{}'".format(x) for x in ct.REGULAR_CHR)
        clazz = "'single'"

        query = '''
            select name, chrom, chromStart, chromEnd
            from snp146
            where name IN ( ''' + snps + ''') AND
                  chrom IN (''' + chroms + ''') AND
                  class = ''' + clazz

        rows = self.conn.execute(query)

        df = pandas.DataFrame(rows.fetchall())
        df.columns = rows.keys()

        return df

    def fetch_alleles(self, rsid):
        snps = ", ".join("'{}'".format(x) for x in rsid)
        chroms = ", ".join("'{}'".format(x) for x in ct.REGULAR_CHR)
        clazz = "'single'"

        query = '''
            select name, chrom, strand, alleles, alleleFreqs
            from snp146
            where name IN ( ''' + snps + ''') AND
                  chrom IN (''' + chroms + ''') AND
                  class = ''' + clazz

        rows = self.conn.execute(query)

        df = pandas.DataFrame(rows.fetchall())
        df.columns = rows.keys()

        # These 2 columns are longblob in mysql
        # Sqlalchemy will read them as bytes not strings
        df.loc[:, "alleles"] = df["alleles"].apply(parse_gb_alleles_col)
        df.loc[:, "alleleFreqs"] = df["alleleFreqs"].apply(parse_gb_allele_freqs_col)

        return df

    def fetch_metadata(self, rsid):
        snps = ", ".join("'" + x + "'" for x in rsid)

        # Use '%%' to escape '%' in python
        query = '''
                SELECT
                    s.name, s.chrom, s.chromStart, s.chromEnd, s.alleles, s.strand,
                    s.class, s.alleleFreqs,
                    dh.score as masterDhsScore, dh.sourceCount as masterDhsCount,
                    tf.tfName, tf.tfCount,
                    pc.sumData/pc.validCount as phastCons
                FROM
                    snp146 s
                LEFT OUTER JOIN
                    wgEncodeAwgDnaseMasterSites dh
                ON
                    dh.bin = s.bin # Speeds up JOINs
                    AND s.chromStart BETWEEN dh.chromStart AND dh.chromEnd - 1
                    AND dh.chrom = s.chrom
                LEFT OUTER JOIN
                    phyloP46wayPlacental pc
                ON
                    pc.bin = s.bin # Speeds up JOINs
                    AND s.chromStart BETWEEN pc.chromStart AND pc.chromEnd - 1
                    AND pc.chrom = s.chrom
                LEFT OUTER JOIN
                (
                    SELECT
                        s.name as snp, GROUP_CONCAT(tf.name) as tfName, COUNT(tf.name) as tfCount
                    FROM
                        snp146 s
                    LEFT OUTER JOIN
                        wgEncodeRegTfbsClusteredV3 tf
                    ON
                        tf.bin = s.bin # Speeds up JOINs
                        AND s.chromStart BETWEEN tf.chromStart AND tf.chromEnd - 1
                        AND tf.chrom = s.chrom
                    WHERE
                        s.name IN ( ''' + snps + ''')
                        AND tf.name != 'POLR2A'
                    GROUP BY s.name
                ) tf
                ON
                    tf.snp = s.name
                WHERE
                    s.name IN  ( ''' + snps + ''')
                '''

        rows = self.conn.execute(query)

        df = pandas.DataFrame(rows.fetchall())
        df.columns = rows.keys()

        # These 4 columns are longblob in mysql
        # Sqlalchemy will read them as bytes not strings
        # df.loc[:, "refNCBI"] = df["refNCBI"].apply(lambda x: x.decode("utf-8") if x != b'' else numpy.nan)
        # df.loc[:, "refUCSC"] = df["refUCSC"].apply(lambda x: x.decode("utf-8") if x != b'' else numpy.nan)
        df.loc[:, "alleles"] = df["alleles"].apply(parse_gb_alleles_col)
        df.loc[:, "alleleFreqs"] = df["alleleFreqs"].apply(parse_gb_allele_freqs_col)

        return df

    def in_protein_coding_exon(self, chrom, chrom_start):
        """
        Steve:

        > Any SNP that is located in a protein-coding exon should be excluded,
        > regardless of the strand on which the gene is located.
        """
        query = '''
                SELECT cdsStart, cdsEnd, exonStarts, exonEnds
                FROM ensGene
                WHERE chrom="{chrom}"
                    AND cdsStart != cdsEnd
                    AND cdsStart <= {chromStart}
                    AND cdsEnd > {chromStart}
                '''.format(chrom=chrom, chromStart=chrom_start)

        rows = self.conn.execute(query)

        for row in rows.fetchall():
            # parse the comma-separated string into a tuple
            # E.g. '117292896,117304741,' => (117292896,117304741)
            exon_starts = ast.literal_eval(row[2].decode("utf-8"))
            exon_ends = ast.literal_eval(row[3].decode("utf-8"))

            for start, end in zip(exon_starts, exon_ends):
                if start <= chrom_start < end:
                    return True

        return False

    def in_nki_lad(self, chrom, chrom_start):
        """
        NKI: Netherlands Cancer Institute
        LAD: Lamina Associated Domains
        """
        query = '''
                SELECT count(1)
                FROM laminB1Lads
                WHERE chrom="{chrom}"
                    AND chromStart <= {chromStart}
                    AND chromEnd > {chromStart}
                '''.format(chrom=chrom, chromStart=chrom_start)

        rows = self.conn.execute(query)

        # fetchone() returns a one-element tuple
        return rows.fetchone()[0]

    def select_vista_enhancer(self, chrom, chrom_start):
        """
        NKI: Netherlands Cancer Institute
        LAD: Lamina Associated Domains
        """
        query = '''
                SELECT name, score
                FROM vistaEnhancers
                WHERE chrom="{chrom}"
                    AND chromStart <= {chromStart}
                    AND chromEnd > {chromStart}
                '''.format(chrom=chrom, chromStart=chrom_start)

        rows = self.conn.execute(query)

        # fetchall() returns list of tuples as [(vh_name, vh_score)]
        return rows.fetchall()

    def select_chia_pet_cluster(self, table_name):
        query = '''
                SELECT name
                FROM {tableName}
                '''.format(tableName=table_name)

        rows = self.conn.execute(query)

        # Each `row` is a one-element tuple, and the element is a string like
        #   "chr1:839717..840790-chr1:872838..874070,6".
        # Split each string by ":", "::", "-" and ","
        # The last part after splitting is unused, so `[:-1]` is applied here
        data = list(map(lambda row: re.split(":|\.\.|-|,", row[0])[:-1], rows.fetchall()))
        df = pandas.DataFrame(data)
        df.columns = ['b1_chrom', 'b1_chromStart', 'b1_chromEnd', 'b2_chrom', 'b2_chromStart', 'b2_chromEnd']

        for col in ['b1_chromStart', 'b1_chromEnd', 'b2_chromStart', 'b2_chromEnd']:
            df.loc[:, col] = pandas.to_numeric(df[col])

        return df

# if __name__ == '__main__':
#     with GenomeBrowserClient('local_hg19') as gb_client:
#         print(gb_client.select_vista_enhancer('chr1', 3190581))
