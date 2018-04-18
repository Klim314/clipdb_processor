"""
Assumes that data is a clipDB bulk download in an input dir
CLIPDB dir must be in folder input/clipdb and must contain:
    (GENE_CIMS, GENE_CITS, GENE_PARalyzer, GENE_Piranha) directories
    and an annotation db containing filename-> sample_name mappings

TODO: Automate annotation collection?
setup directory
"""
import matplotlib
matplotlib.use("Agg")
import os
import shutil
import logging
import pandas as pd
import seaborn as sb
from clip_analyze import ClipFinder, ClipAnalyzer
from clip_compile import compile_clip_degs
from collections import Counter

logging.basicConfig(level=logging.DEBUG)


if "runtype" not in config:
    # print("Please select a runtype: init")
    config["runtype"] = "default"

base_dirs = ["temp/compiled", "input", "output/bedfiles/piranha", "output/bedfiles/specialized",
             "output/peak_summaries/piranha", "output/peak_summaries/specialized"]

if config["runtype"] == "initialize" or config["runtype"] == "init":
    output = base_dirs
elif config["runtype"] == "default":
    # Get the gene designator
    dat = [i for i in os.listdir("input") if i.endswith("Piranha")]
    gene_sym = dat[0].split("_")[0]
    specialized = [i for i in os.listdir("input") if i.startswith(gene_sym) and not i.endswith("Piranha")]
    piranha = gene_sym + "_Piranha"
    bedfile_list = os.listdir(os.path.join("input", piranha))
    print("Specialized: ", specialized)
    print("Piranha: ", piranha)
    print("Bedfiles: ", bedfile_list)

    # Perform integrity check
    df = pd.read_csv(os.path.join("input", "annotations.tsv"), sep="\t")
    expected_bedfiles = df["file_name"].tolist()
    expected_files = [os.path.splitext(i)[0] for i in expected_bedfiles]
    # Check that all files are present
    # Annotations are currently used to control what will be processed. Maybe change this
    assert(all(i in bedfile_list for i in expected_bedfiles))

    # Check that all specialized files are present
    expected_specialized = []
    for i in specialized:
        expected_specialized.extend(os.listdir(os.path.join("input", i)))

    assert(all(i in bedfile_list for i in expected_specialized))

    # output = [expand("/".join(["temp", "compiled", "specialized", "{bedfile}"]), bedfile=expected_bedfiles),
    #           expand("/".join(["temp", "compiled", "piranha", "{bedfile}"]), bedfile=expected_bedfiles)]
    # output = [expand("/".join(["output", "peak_summaries", "specialized", "{file}.tsv"]), file=expected_files),
    #           expand("/".join(["output", "peak_summaries", "piranha", "{file}.tsv"]), file=expected_files)]

    output = ["output/specialized/report.txt", "output/piranha/report.txt"]
    print(output)

rule all:
    input: output

rule prepare_dirs:
    output: base_dirs
    run:
        for i in output:
            if not os.path.exists(i):
                os.makedirs(i)

"""
Compilaton of raw data to a temporary directory for later, easier handling
"""
rule compile_data:
    input: ["/".join(["input", i]) for i in specialized]
    output: expand("/".join(["temp", "compiled", "specialized", "{bedfile}"]), bedfile=expected_bedfiles) 
    run:
        print("Compiling bedfiles fron: {}".format(input))
        holder = []
        for root_path in input:
            holder.extend([os.path.join(root_path, i) for i in os.listdir(root_path)])
        print(holder)
        for bedfile_path in holder:
            outpath = "temp/compiled/specialized/{}".format(os.path.basename(bedfile_path))
            print("Copying from {} to {}".format(bedfile_path, outpath))
            shutil.copy(bedfile_path, outpath)


input_data = ["input/{}_Piranha"]
input_data = [i.format(gene_sym) for i in input_data]

rule compile_piranha_data:
    input: input_data
    output: expand("/".join(["temp", "compiled", "piranha", "{bedfile}"]), bedfile=expected_bedfiles)
    run:
        holder = []
        for root_path in input:
            holder.extend([os.path.join(root_path, i) for i in os.listdir(root_path)])
        for bedfile_path in holder:
            outpath = "temp/compiled/piranha/{}".format(os.path.basename(bedfile_path))
            shutil.copy(bedfile_path, outpath)


print("SUmmarize peaks\n", expand("output/peak_summaries/{{dataset}}/{expected_files}.tsv", expected_files=expected_files))
"""
Creates a summary of gene peaks, peak counts and exon peaks for an individual clip-bed
"""
rule summarize_peaks:
    input: expand("temp/compiled/{{dataset}}/{expected_bedfiles}", expected_bedfiles=expected_bedfiles)
    output:
        summaries = "output/peak_summaries/{dataset}",
        files = expand("output/peak_summaries/{{dataset}}/individual/{expected_files}.tsv", expected_files=expected_files),
        summary = "output/peak_summaries/{dataset}/summary.tsv"
    params:
        test = 1
    run:
        # annotation = os.path.join("input", "annotations.tsv")
        cf = ClipFinder({"mouse": "/mnt/d/work/sc/data/ensembl_mmus_2018_01_03_ucsc-chroms.bed",
                         "human": "/mnt/d/work/sc/data/ensembl_hsapiens_2018_01_03_ucsc-chroms.bed"},
                        {"mouse": "/mnt/d/work/sc/data/ensembl_mmus_exon_2018_01_03_ucsc-chroms.bed",
                         "human": "/mnt/d/work/sc/data/ensembl_hsapiens_exon_2018_01_03_ucsc-chroms.bed"},
                        annotations="input/annotations.tsv",
                        crossmap="/mnt/d/work/sc/data/ensembl_human_mappings_2018_01_05.tsv")
        cf.load_mappings({"mouse": "/mnt/d/work/sc/data/ensembl_mmusculus_mappings_2017_11_16.tsv",
                          "human": "/mnt/d/work/sc/data/ensembl_human_2018_01_05.tsv"})
        print(input)
        # Generate individual peak-count tsvs for each input ped
        for in_bed, out_tsv in zip(input, output.files):
            in_bed, out_tsv = str(in_bed), str(out_tsv)
            print(in_bed, out_tsv)
            cf.process(in_bed, out_tsv, "mouse")

        # Collated the peak-count tsvs and compile into easily understandable data
        # Note. Final vector depends on
        base = pd.DataFrame({"ens_gene": [],
                             "peak_counts": [],
                             "category": []})
        for summ_data_path in (str(i) for i in output.files):
            logging.debug(summ_data_path)
            df = pd.read_csv(summ_data_path, sep="\t")
            df["category"] = summ_data_path
            base = pd.concat([base, df])

        res = base[base["peak_counts"] > 0].groupby("ens_gene").category.count()
        res = pd.DataFrame({"ens_gene": res.index,
                            "agreement": res.values,
                            "agreement_prop": res.values / len(output.files)})
        res.to_csv(output.summary, sep="\t", index=False)

rule integrate_expression_data:
    input:
        peak_summary = "output/peak_summaries/{dataset}/summary.tsv",
        de_data = "input/exp_data",
        other_peak_data = "input/supplementary_peaks",
        mapping_path = "/mnt/d/work/sc/data/ensembl_mmusculus_mappings_2017_11_16.tsv"
    output:
        "output/{dataset}/clip_summary/compiled_clip.tsv",
        "output/{dataset}/clip_summary/compiled_clip.xlsx"
    params:
        outdir = "output/{dataset}/clip_summary"
    run:
        compiled, view = compile_clip_degs(input.de_data, input.other_peak_data,
                                           input.peak_summary, input.mapping_path,
                                           params.outdir)
        # Agreement distribution plots
        agreement = [int(i) for i in compiled["agreement"]]
        plt = sb.distplot(agreement)
        plt.figure.savefig("output/{}/clip_summary/agreement-plot.png".format(wildcards.dataset), kde=False)
        matplotlib.pyplot.clf()

        # Remove the zeros
        agreement = [int(i) for i in compiled["agreement"] if i != 0]
        plt = sb.distplot(agreement)
        plt.figure.savefig("output/{}/clip_summary/agreement-plot_nozeros.png".format(wildcards.dataset), kde=False)




"""
Analysis of clip-peak data and the creation of the following plots
    PCA
    Clustering (maybe don't use the PCs for this)

"""
rule analyze_data:
    input: expand("output/peak_summaries/{{dataset}}/individual/{expected_files}.tsv", expected_files=expected_files)
    output:
        "output/{dataset}/plots/pca",
        # "output/{dataset}/plots/cluster"
    params:
        indir = "output/peak_summaries/{dataset}/individual",
        outdir = "output/{dataset}/plots"
    run:
        ca = ClipAnalyzer("input/annotations.tsv")
        ca.load_data(params.indir)
        ca.generate_all_plots(params.outdir)

rule compile_report:
    input:
        "output/{dataset}/plots/pca",
        # "output/{dataset}/plots/cluster",
        "output/peak_summaries/{dataset}/summary.tsv",
        "output/{dataset}/clip_summary/compiled_clip.tsv",
        "output/{dataset}/clip_summary/compiled_clip.xlsx"
    output:
        "output/{dataset}/report.txt"
    run:
        with open(str(output), "w") as f:
            f.write("DONE")
