import pandas as pd
import seaborn as sb
import os
import itertools

# Various coloring helper functions 
def color_sig(val, cutoff=0.1):
    color = "#87CEFA" if val < cutoff else "#FF6347"
    return "background-color: {}".format(color)


def color_reg_type(val):
    if val > 0:
        color = "#e0fde0"
    elif val < 0:
        color = "#FCE5E5"
    else:
        color = "#FFFFFF"
    return "background-color: {}".format(color)


def compile_clip_degs(de_data, other_peak_data, input_data, mapping_path, outdir):
    """Compiles gene expression data and clipdb/other CLIP peak data into a single color coded document/report
    args
        de_data
        other_peak_data
        input_data
        mapping_path
    output:
        outdir/compiled_clip.tsv
        outdir/compiled_clip.xls
    """
    base_table = pd.read_csv(input_data, sep="\t")
    # Add the supplementary data to the front first
    for path in (os.path.join(other_peak_data, i) for i in os.listdir(other_peak_data)):
        peak_table = pd.read_csv(path, sep="\t")
        # Use the filename to distinguish between each
        dataset_name = os.path.splitext(os.path.basename(path))[0]
        peak_table = peak_table[["ens_gene", "peak_counts"]]
        peak_table = peak_table.rename(columns={"peak_counts": "{}_peak_counts".format(dataset_name)})
        base_table = base_table.merge(peak_table, how="outer")

    # Add in all the DE gene data
    for path in (os.path.join(de_data, i) for i in os.listdir(de_data)):
        de_table = pd.read_csv(path, sep="\t")
        de_table = de_table[["target_id", "qval", "b"]]
        # Use the filename to distinguish between each
        dataset_name = os.path.splitext(os.path.basename(path))[0]
        de_table = de_table.rename(columns={"target_id": "ens_gene",
                                            "qval": "{}_qval".format(dataset_name),
                                            "b": "{}_LEFC".format(dataset_name)})
        # print(de_table.head())
        base_table = base_table.merge(de_table, how="outer", on="ens_gene")

    # Merge in the ens_id symbol mappings for ease of reading
    mappings = pd.read_csv(mapping_path, sep="\t")
    mappings.columns = ["ens_gene", "ens_trans", "symbol"]
    mappings = mappings[["ens_gene", "symbol"]].drop_duplicates()
    base_table = mappings.merge(base_table, how="right")

    # Different columns for effectove coloring
    qval_cols = [i for i in base_table.columns if i.endswith("qval")]
    lefc_cols = [i for i in base_table.columns if i.endswith("LEFC")]
    peak_counts = [i for i in base_table.columns if i.endswith("peak_counts")]
    # fill targets
    targets = {i: 0 for i in itertools.chain(peak_counts, ["agreement", "agreement_prop"])}
    print(targets)

    # Fill NA values with int-negatives (0 for foldchange, 1 for qvals)
    base_table = base_table.fillna(targets)\
                           .fillna({i: 0 for i in lefc_cols})\
                           .fillna({i: 1 for i in qval_cols})

    view = base_table.style.applymap(color_sig, subset=qval_cols)
    view = view.applymap(color_reg_type, subset=lefc_cols)
    cmap = sb.light_palette("green", as_cmap=True)
    view = view.background_gradient(cmap, subset=["agreement", "agreement_prop"] + peak_counts)
    
    # write tsv
    base_table.to_csv(os.path.join(outdir, "compiled_clip.tsv"), sep="\t", index=False)
    # write excel
    writer = pd.ExcelWriter(os.path.join(outdir, "compiled_clip.xlsx"))
    view.to_excel(writer, "compiled_clip", engine="openpyxl")
    writer.save()
    # write... html?
    # view.to_html("compiled_clip.html")
    return view
