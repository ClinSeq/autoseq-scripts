#!/usr/bin/env python
import sys
import click
import pandas as pd
from functools import partial


# Adjust segment medians from the adjusted bin values:
def calc_new_log2(segments_row, bins=None):
    bins_idx = (bins.chromosome == segments_row.chromosome) & \
        (bins.start >= segments_row.start) & \
        (bins.end <= segments_row.end)
    return bins.log2[bins_idx].median(skipna=True)


@click.command()
@click.option('--input-cnr', type=click.Path(exists=True), help='Input cnr file', required=True)
@click.option('--input-cns', type=click.Path(exists=True), help='Input cns file', required=True)
@click.option('--input-reference', type=click.Path(exists=True), help='Input reference file, for generating fixed cns and cnr files', required=True)
@click.option('--output-cnr', default='CNVkit_fixed.cnr', help='Output cnr file')
@click.option('--output-cns', default='CNVkit_fixed.cns', help='Output cns file')
def main(input_cnr, input_cns, input_reference, output_cnr, output_cns):
    """
Fix an input CNV-kit cnr and cns file using the specified table of reference data.
"""

    # Load the sample:
    bins = pd.read_table(input_cnr)
    segments = pd.read_table(input_cns)

    # Load the reference data for fixing the samples:
    reference_data = pd.read_table(input_reference)

    background_idx = (bins.gene == "Background")
    autosomal_idx = ~(bins.chromosome.isin(["X", "Y"]))

    # Median center target and antitarget bins:
    median_log2_bg = bins.log2[background_idx & autosomal_idx].median(skipna=True)
    median_log2_not_bg = bins.log2[~background_idx & autosomal_idx].median(skipna=True)
    bins.loc[background_idx, "log2"] = bins.log2[background_idx] - median_log2_bg
    bins.loc[~background_idx, "log2"] = bins.log2[~background_idx] - median_log2_not_bg

    # Adjust bin values based on healthy donor reference data:
    bins = bins.merge(reference_data)
    bins = bins.sort_values(["chromosome", "start"])
    bins.log2 = bins.log2 - bins.referenceMedian

    segments.log2 = segments.apply(partial(calc_new_log2, bins=bins), axis=1)

    # Write output files:
    bins.iloc[:, :6].to_csv(output_cnr, sep="\t", index=None)
    segments.to_csv(output_cns, sep="\t", index=None)


if __name__ == "__main__":
    sys.exit(main())