#!/usr/bin/env python


import subprocess
import argparse
import pandas as pd
import os


def find_primers(assembly_file, primer_5, primer_3, output_dir, max_mismatches=5, amplicon_file=None):
    """
    Find primers in an assembly, extract the porA sequence, and optionally compare it to a set of amplicons.
    """
    # Extract assembly ID from filename
    assembly_id = os.path.splitext(os.path.basename(assembly_file))[0]

    # Define file paths
    locate_5_output = os.path.join(output_dir, f"seqkit_locate_5start_{assembly_id}.tsv")
    locate_3_output = os.path.join(output_dir, f"seqkit_locate_3end_{assembly_id}.tsv")
    bed_file = os.path.join(output_dir, f"{assembly_id}.bed")
    subseq_output = os.path.join(output_dir, f"{assembly_id}_porA.fa")
    wgs_out = os.path.join(output_dir, f"{assembly_id}_porA_wgs_match.txt")

    for mismatch in range(0, max_mismatches + 1):
        # Run locate for 5' and 3' primers up to max mismatches 
        run_locate_with_mismatch(primer_5, locate_5_output, mismatch, assembly_file)
        run_locate_with_mismatch(primer_3, locate_3_output, mismatch, assembly_file)

        # Check primer hits
        locate_5_df, locate_5_warning = check_primer_hits(locate_5_output)
        locate_3_df, locate_3_warning = check_primer_hits(locate_3_output)

        # Handle warnings and retry with increased mismatch
        if locate_5_warning:
            print(f"5' primer issue in {assembly_file}: {locate_5_warning}. Retrying with mismatch={mismatch + 1}")
            continue
        if locate_3_warning:
            print(f"3' primer issue in {assembly_file}: {locate_3_warning}. Retrying with mismatch={mismatch + 1}")
            continue

        # Validate primer pair and ensure they are on the same contig
        if locate_5_df is not None and locate_3_df is not None:
            seqID_5, end_5 = locate_5_df.iloc[0][["seqID", "end"]]
            seqID_3, start_3 = locate_3_df.iloc[0][["seqID", "start"]]

            if seqID_5 == seqID_3:
                # Create BED file with primer positions
                with open(bed_file, "w") as bed:
                    bed.write(f"{seqID_5}\t{end_5}\t{start_3 - 1}\n")

                # Extract subsequence based on BED file
                subprocess.run(["seqkit", "subseq", "--quiet", "--bed", bed_file, "-o", subseq_output, assembly_file])
                print(f"Subsequence extracted to {subseq_output}")

                # Compare extracted sequence to amplicons (if provided)
                if amplicon_file:
                    run_blast(amplicon_file, subseq_output, wgs_out)
                else:
                    print("No amplicon file provided.")
                return  # Exit early after successful processing

            else:
                print(f"Primers not on the same contig in {assembly_file}. Retrying with mismatch={mismatch + 1}")

    print(f"Max mismatches reached. No valid primer pair found in {assembly_file}.")


def run_locate_with_mismatch(primer, output_file, mismatch, assembly_file):
    """
    Run seqkit locate with a specified number of mismatches.
    """
    try:
        subprocess.run(
            ["seqkit", "locate", "-p", primer, "-i", "-o", output_file, "-m", str(mismatch), assembly_file],
            check=True, capture_output=True, text=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Error running seqkit locate for primer '{primer}' on file '{assembly_file}': {e}")


def check_primer_hits(output_file):
    """
    Check the results of seqkit locate for primer hits.
    """
    try:
        df = pd.read_csv(output_file, sep="\t")
        if len(df) == 0:
            return None, "No primers found."
        if len(df) > 1:
            return None, "Multiple primers found."
        return df, None
    except pd.errors.EmptyDataError:
        return None, "No primers found."


def run_blast(amplicon_file, subseq_output, wgs_output):
    """
    Run BLAST to compare an extracted porA sequence to a set of amplicons.
    Write the assembly ID and BLAST results if criteria are met; otherwise, write the assembly ID and 'NA'.
    """
    # Extract assembly ID from subseq_output file name
    assembly_id = os.path.basename(subseq_output).split("_porA")[0]

    blastn_cmd = [
        "blastn",
        "-subject", amplicon_file,
        "-query", subseq_output,
        "-outfmt", "6 qaccver saccver slen pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-perc_identity", "100",
        "-culling_limit", "1",
        "-max_hsps", "1"
    ]


    try:
        # Run BLAST and capture its output
        with subprocess.Popen(blastn_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True) as proc:
            for line in proc.stdout:
                cols = line.strip().split("\t")
                if len(cols) < 12:
                    continue

                qaccver, saccver, slen, pident, length = cols[0], cols[1], int(cols[2]), float(cols[3]), int(cols[4])

                # Check if the hit meets the criteria
                if pident == 100.0 and slen == length:
                    # Write the result to the valid hits file
                    with open(wgs_output, "w") as vf:
                        vf.write(f"{assembly_id}\t{qaccver}\t{saccver}\n")
                    return  # Exit early as culling_limit ensures one result
            else:
                # This executes if no line in stdout meets the criteria
                with open(wgs_output, "w") as vf:
                    vf.write(f"{assembly_id}\tNA\tNA\n")

    except subprocess.CalledProcessError as e:
        # Handle BLAST failure
        with open(wgs_output, "w") as vf:
            vf.write(f"{assembly_id}\tNA\tNA\n")


def main():
    """
    Main function to parse command-line arguments and process assemblies.
    """
    parser = argparse.ArgumentParser(description="Find porA primers in assemblies and extract the porA sequence.")
    parser.add_argument("--ass", required=True, help="Directory containing assemblies")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--amplicons", help="File of amplicons to compare to extracted porA sequences")
    args = parser.parse_args()

    # Collect all .fasta or .fa files in the provided directory
    assembly_files = [
        os.path.join(args.ass, f) for f in os.listdir(args.ass)
        if f.endswith(".fasta") or f.endswith(".fa")
    ]

    # Primers
    primer_5 = "CCACAATTATGGTTAGCTTA"
    primer_3 = "CTCTCCAAAACTTAACTTCTCA"

    # Process each assembly file
    for assembly_file in assembly_files:
        find_primers(assembly_file, primer_5, primer_3, args.out, amplicon_file=args.amplicons)


if __name__ == "__main__":
    main()
