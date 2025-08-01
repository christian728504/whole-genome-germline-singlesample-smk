import subprocess
from pathlib import Path
import shlex
import csv
import argparse
import sys

CMD = (
    """
    VerifyBamID \
    --Verbose \
    --NumPC 4 \
    --NumThread {threads} \
    --Output {output_self_sm} \
    --BamFile {input_bam} \
    --Reference {ref_fasta} \
    --UDPath {contamination_sites_ud} \
    --MeanPath {contamination_sites_mu} \
    --BedPath {contamination_sites_bed} \
    """
)

def parse_and_run(args):
    return_code = 0
    
    threads = args.threads
    output_self_sm = args.output_self_sm
    input_bam = args.input_bam
    ref_fasta = args.ref_fasta
    contamination_sites_ud = args.contamination_sites_ud
    contamination_sites_bed = args.contamination_sites_bed
    contamination_sites_mu = args.contamination_sites_mu
    contamination_underestimation_factor = args.contamination_underestimation_factor

    cmd = CMD.format(
        threads=threads,
        output_self_sm=output_self_sm.rstrip(".selfSM"),
        input_bam=input_bam,
        ref_fasta=ref_fasta,
        contamination_sites_ud=contamination_sites_ud,
        contamination_sites_bed=contamination_sites_bed,
        contamination_sites_mu=contamination_sites_mu,
    )
    split_cmd = shlex.split(cmd)

    print(f"Running command:\n\n{cmd}")

    subprocess.check_call(split_cmd)

    print("Verifying output")
    
    with Path(output_self_sm).open("r") as selfSM:
        reader = csv.DictReader(selfSM, delimiter='\t')
        i = 0
        for row in reader:
            if float(row["FREELK0"])==0 and float(row["FREELK1"])==0:
              # a zero value for the likelihoods implies no data. This usually indicates a problem rather than a real event.
              # if the bam isn't really empty, this is probably due to the use of a incompatible reference build between
              # vcf and bam.
              print("Found zero likelihoods. Bam is either very-very shallow, or aligned to the wrong reference (relative to the vcf).", file=sys.stderr)
              return_code = 1
        print(float(row["FREEMIX"])/contamination_underestimation_factor, file=sys.stderr)
        i = i + 1
        # there should be exactly one row, and if this isn't the case the format of the output is unexpectedly different
        # and the results are not reliable.
        if i != 1:
            print("Found %d rows in .selfSM file. Was expecting exactly 1. This is an error"%(i), file=sys.stderr)
            return_code = 2
    return return_code

def main(args):
    result = parse_and_run(args)
    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("threads", type=int, help="Number of threads")
    parser.add_argument("output_self_sm", help="Output dot selfSM file")
    parser.add_argument("input_bam", help="Input bam file")
    parser.add_argument("ref_fasta", help="Reference fasta file")
    parser.add_argument("contamination_sites_ud", help="UD path")
    parser.add_argument("contamination_sites_bed", help="Bed path")
    parser.add_argument("contamination_sites_mu", help="Mu path")
    parser.add_argument("contamination_underestimation_factor", type=float, help="Underestimation factor")
    args = parser.parse_args()
    result = main(args)
    if result != 0:
        print(f"An error occurded. Exiting with code {result}")
    sys.exit(result)
