#!/usr/bin/env python

import vcf
import begin
import logging

class InvalidInputVcfException(Exception):
    pass


@begin.start(auto_convert=True)
def extract_somatic_afs(vcf_filename, output="output.bedGraph"):
    logging.basicConfig(level=logging.DEBUG)

    vcf_reader = vcf.Reader(filename=vcf_filename)

    # Extract the names of all tumor (i.e. non-normal) samples:
    non_normal_samples = [sample for sample in vcf_reader.samples if sample.split("-")[3] != "N"]

    # Currently just cope with input files containing exactly one tumor sample:
    if len(non_normal_samples) != 1:
        raise InvalidInputVcfException("Input VCF must include exactly one tumor (i.e. non-normal) sample.")

    non_normal_sample = non_normal_samples[0]

    with open(output, 'w') as output_file:
        for record in vcf_reader:
            tumor_genotype = record.genotype(non_normal_sample)
            if tumor_genotype.data.DP > 0:
                allelic_fraction = tumor_genotype.data.AO/tumor_genotype.data.DP
                print >> output_file, "{}\t{}\t{}\t{}".format(record.CHROM, record.POS, record.POS, allelic_fraction)