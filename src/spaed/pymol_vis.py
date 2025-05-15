import os
import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('pred_path', help='Path to spaed prediction file. Make sure entries have same name as in pdb filename(s).')
    parser.add_argument('pdb_path', help='Path to pdb file or folder containing pdb files. Make sure name of file(s) correspond to entries in prediction file.')
    parser.add_argument('--output_folder', help="Path to output file. (default ./pymol_vis)", default="pymol_vis")
    parser.add_argument('--output_type', help='Whether to save pymol session (.pse) or png of structure. (options [pse png both]).', default="both")
    args = parser.parse_args()


    pred_path = args.pred_path
    pdb_path = args.pdb_path
    output_folder = args.output_folder
    output_type = args.output_type

    return pred_path, pdb_path, output_folder, output_type



def vis_structure(pdb, sample_name, pred, output_folder, output_type='png'):
    #Load file
    cmd.feedback('disable', 'all', 'actions')
    cmd.feedback('disable', 'all', 'results')
    cmd.reinitialize()
    cmd.load(pdb)
    cmd.orient()
    cmd.center()
    cmd.zoom('center', 35)

    #load domain pred
    try:
        doms = str(pred.loc[sample_name, "domains"])
        disordered = str(pred.loc[sample_name, "disordered"])
        linkers = str(pred.loc[sample_name, "linkers"])
    except:
        print(f"{sample_name} not found. Make sure name.pdb corresponds to name.json (pae file).")
        return

    #show cartoon
    cmd.hide('all')
    cmd.show('cartoon')

    #color
    if linkers != "":
        for link in linkers.split(";"):
            cmd.color_deep('gray80', f"resi " + link, 0)

    if disordered != "":
        for dis in disordered.split(";"):
            cmd.color_deep('red', f"resi " + dis, 0)

    cols = ["cyan", "magenta", "forest", "slate", "yellow", "orange", "deepteal", "hotpink", "green", "purple", "salmon", "lightorange"]
    for i, dom in enumerate(doms.split(";")):
        cmd.color_deep(cols[i], f"resi {dom}", 0)

    if (output_type == "png") | (output_type == "both"):
        cmd.set('ray_trace_mode', 1)
        cmd.bg_color("white")
        cmd.ray(1280,720)

        #save file
        cmd.png(os.path.join(output_folder, f'{sample_name}.png'), dpi=300)

    if (output_type == "pse") | (output_type == "both"):
        cmd.save(os.path.join(output_folder, f'{sample_name}.pse'))

    if (output_type != "pse") & (output_type != "png") & (output_type != "both"):
        print("Output type incorrectly specified. Can either be 'png', 'pse' or 'both'.")



def pymol_vis(pdb_path, pred_path, output_folder, output_type='png'):
    from pymol import cmd

    #load prediction file
    assert  os.path.isfile(pred_path), "Specified domain prediction file was not found or is not a file."
    preds = pd.read_csv(pred_path, index_col=0)

    #make sure output folder exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    ##fetch sequence from fasta
    if os.path.isfile(pdb_path): #only one structure to process
        sample_name = os.path.basename(pdb_path).replace(".pdb", "")
        vis_structure(pdb_path, sample_name, preds, output_folder, output_type)

    else: #multiple structures to process
        for name in os.listdir(pdb_path):
            if name.endswith(".pdb"):
                sample_name = name.replace(".pdb", "")
                vis_structure(os.path.join(pdb_path, name), sample_name, preds, output_folder, output_type)

def main():
    pred_path, pdb_path, output_folder, output_type = parse_args()
    pymol_vis(pdb_path, pred_path, output_folder, output_type)

if __name__ == '__main__':
    #pred_path, pdb_path, output_folder, output_type = parse_args()
    #pymol_vis(pdb_path, pred_path, output_folder, output_type)
    main()
