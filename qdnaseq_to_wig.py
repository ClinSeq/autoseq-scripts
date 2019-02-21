#!/usr/bin/env python

from collections import OrderedDict
import pandas as pd
import begin
import logging

@begin.start(auto_convert=True)
def generate_wigs(qdnaseq_filename, copynumber_wig_filename, readcount_wig_filename):
    logging.basicConfig(level=logging.DEBUG)

    col_types = OrderedDict({
        "chromosome": object,
        "start": int,
        "end": int,
        "bases": float,
        "gc": float,
        "mappability": float,
        "blacklist": float,
        "residual": object,
        "use": object,
        "readcount": int,
        "copynumber": object,
        "segmented": object,
    })

    df = pd.read_table(qdnaseq_filename, dtype = col_types)
    df_nonnull = df.loc[~(df["segmented"].isnull()),:]

    curr_chrom = None
    curr_pos = -1
    idx = 0
    with open(copynumber_wig_filename, 'w') as copynumber_wigfile, open(readcount_wig_filename, 'w') as readcount_wigfile:
        for index, row in df_nonnull.iterrows():
            if row["chromosome"] != curr_chrom or row["start"] != curr_pos + 15000:
                idx += 1
                curr_chrom = row["chromosome"]
                print >> copynumber_wigfile, "fixedStep chrom={} start={} step=15000 span=15000".format(str(curr_chrom), str(row["start"]))
                print >> readcount_wigfile, "fixedStep chrom={} start={} step=15000 span=15000".format(str(curr_chrom), str(row["start"]))
            print >> copynumber_wigfile, "{}".format(row["copynumber"])
            print >> readcount_wigfile, "{}".format(row["readcount"])
            curr_pos = row["start"]
