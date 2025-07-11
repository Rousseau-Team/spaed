# Functions to fetch domains sequences from whole sequence file

import os
import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('pred_path', help='Path to spaed prediction file. Make sure entries have same name as in fasta.')
    parser.add_argument('--fasta_path', help='Path to fasta file or folder containing fasta files. Make sure headers and/or file names correspond to entries in prediction file.')
    parser.add_argument('--output_path', help="Path to output file.", default="predicted_domain_seqs.faa")
    parser.add_argument('--mode', help='Whether to get sequences for predicted domains, linkers, or disordered regions. options: domains, linkers, disordered.', default="domains")
    args = parser.parse_args()


    pred_path = args.pred_path
    fasta_path = args.fasta_path
    output_path = args.output_path
    mode = args.mode

    return pred_path, fasta_path, output_path, mode


# Load prediction file
def load_predictions(filename):
    df = pd.read_csv(filename, index_col= 0)
    df.index = df.index.astype(str)
    df.linkers = df.linkers.astype(str)
    df.disordered = df.disordered.astype(str)
    df.domains = df.domains.astype(str)
    return df

# Search for fasta with protID
def search_for_fasta(directory, name):
    for root, dirs, files in os.walk(directory, topdown=True):
        for e in files+dirs:
            if os.path.splitext(e)[0] == name:
                if os.path.splitext(e)[1] in [".faa", ".fasta", ".fa"]:
                    return os.path.join(root, e)

    print(f"No fasta file found for {name} with extension .fa, .faa or .fasta")
    return None


# Load sequence(s) from fasta
def load_fasta(fasta_path, prot):
    if prot != "":
        fasta_path = search_for_fasta(fasta_path, prot)
        if fasta_path is None: return None

    seqs = dict()
    with open(fasta_path) as f:
        lines = f.readlines()
        for i in range(0, len(lines)):
            s=lines[i].strip()
            if s == "": pass
            elif s[0] == '>':
                key=s[1:].split()[0]
            else:
                seqs[key] = s
    return seqs

# Get sequence from dict of {protID: seq}
def get_seq(fasta_seqs, prot):
    try: seq = fasta_seqs[prot]
    except: print(f"{prot} not found in fasta file.")
    return seq


# Fetch domains
def subset_dom_seq(preds, prot, seq, mode="domains"):
    if preds.loc[prot, mode] == "": return []

    dom_seqs = []
    for dom in preds.loc[prot, mode].split(";"):
        start_ind = int(dom.split("-")[0])-1
        end_ind = int(dom.split("-")[1])-1

        dom_seq = seq[start_ind:end_ind+1]
        dom_seqs.append(dom_seq)

    return dom_seqs


# Save domains for all proteins in fasta
def save_to_fasta(dom_seqs, prot, output_path):
    with open(output_path, "a") as f:
        for i, dom in enumerate(dom_seqs):
            header=f">{prot}_{i}"
            f.write(header); f.write("\n")
            f.write(dom); f.write("\n")


# Launch all
def fetch_domains(pred_path, fasta_path, output_path="predicted_domain_seqs.faa", mode="domains"):
    preds = load_predictions(pred_path)

    if os.path.isfile(fasta_path):
        try: fasta_seqs = load_fasta(fasta_path, prot="")
        except: print(f"Fasta file not formatted correctly: {fasta_path}")

    if os.path.isfile(output_path):
        print("Output file already exists and will be overwritten.")
        os.remove(output_path)

    for prot in preds.index:
        try:
            #fetch sequence from fasta
            if os.path.isfile(fasta_path):
                seq = get_seq(fasta_seqs, prot)
            else:
                fasta_seq = load_fasta(fasta_path, prot)
                if fasta_seq is None: continue
                seq = get_seq(fasta_seq, prot)

            #subset domains from sequence
            dom_seqs = subset_dom_seq(preds, prot, seq, mode)

            #save to output fasta file
            save_to_fasta(dom_seqs, prot, output_path)

        except:
            print(f"Sequence {prot} skipped. Please investigate.")

if __name__ == '__main__':
    pred_path, fasta_path, output_path, mode = parse_args()
    fetch_domains(pred_path, fasta_path, output_path, mode)

