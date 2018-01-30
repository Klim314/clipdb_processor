#! /usr/bin/env python
"""
TODO: Handle issues with multi-species mapping Currently maps human directly to mouse. THATS BAD
"""
import matplotlib

matplotlib.use("Agg")

import os
import pandas as pd
import numpy as np
import seaborn as sns

from sklearn.decomposition import PCA
from collections import defaultdict
from pybedtools import BedTool
from pprint import pprint
from matplotlib import pyplot as plt


class ClipFinder:
    """Processes a CLIP-peak bedfile, calculating peak counts and exon peak counts for each gene found within the bedfile
    """

    def __init__(self, ref_genes, ref_exons, annotations, crossmap):
        """
        gene_bed: bedfile containing gene ranges. name field must be gene identifier
        exon_bed: bedfile containing exon ranges.
        annotations: annotations for the given dataset
        crossmap: cross species mappings

        """
        self.ref_gene = {species: BedTool(path) for species, path in ref_genes.items()}
        self.ref_exon = {species: BedTool(path) for species, path in ref_exons.items()}
        self.annotations = pd.read_csv(annotations, sep="\t")

        self.mappings = dict()
        self.crossmap = self._load_crossmap(crossmap)

    def _get_species(self, string):
        # Gets the species from a bedfile/tsv/samplename string
        sample_name = os.path.basename(string)
        sample_name = os.path.splitext(sample_name)[0]
        print("sample_name: ", sample_name)
        species = self.annotations[self.annotations["sample_name"] == sample_name]["species"]
        try:
            species = species.tolist()[0]
        except:
            print(sample_name)
            print(self.annotations["sample_name"])
            raise
        return species

    def _get_species_ens(self, ens_id):
        if ens_id.startswith("ENSG"):
            species = "human"
        elif ens_id.startswith("ENSMUSG"):
            species = "mouse"
        else:
            species = None
        return species

    def _load_crossmap(self, crossmap):
        """From an ensembl biomart dump containing various species ENS_GENE_IDS
           creates a table mapping the various species ENS_GENE_IDS to each other
        """
        df = pd.read_csv(crossmap, sep="\t")
        rename_targets = dict()
        for col in df.columns:
            if col == "Gene name" or col == "symbol":
                rename_targets[col] = "symbol"
                continue
            data = df[col]
            # Find the first non-blank entry in the column
            for i in data:
                if i and not pd.isnull(i):
                    ens_id = i
                    break
            # Map ENS_ID_PREFIX to species and add to rename_dict
            print("Identifying for {}".format(ens_id))
            species = self._get_species_ens(ens_id)
            if not species:
                print("Unable to handle data with code {}".format(ens_id))
            rename_targets[col] = species
        print(rename_targets)
        return df.rename(columns=rename_targets).fillna('')


    def load_mapping(self, mapping, species):
        """
        Loads the gene mapping from an ensembl biomart tsv
        """
        df = pd.read_csv(mapping, sep="\t")
        mapper = {"Gene stable ID": "ens_gene",
                  "Gene name": "symbol",
                  "Gene start (bp)": "start",
                  "Gene end (bp)": "end",
                  "Strand": "strand",
                  "Chromosome/scaffold name": "chr",
                  "Exon region start (bp)": "start",
                  "Exon region end (bp)": "end",
                  "Transcript stable ID": "ens_trans"}
        df.columns = [mapper[i] for i in df.columns]
        df = df[["ens_gene", "symbol"]].drop_duplicates
        self.mappings[species] = df

    def load_mappings(self, mappings):
        """
        loads multiple mapping tsvs for specific species
        mappings: (mapping_path, species)
        """
        for species, path in mappings.items():
            self.load_mapping(path, species)

    def count_peaks(self, sample_clip_bed):
        """Analyzes a sample bed, identifying genes contained within and the respective number of counts
        args: 
            sample_clip_bed: string. Path to a sample clip bedfile to be analyzed
        returns
        peak_counts dict keyed by by gene identifier with contents (peak_counts, exon_peak_counts)
        """
        sample_bed = BedTool(sample_clip_bed)
        # Identify the species associated with the bedfile
        species = self._get_species(sample_clip_bed)

        # Perform bedtools intersects to identify regions of interest for the specific species
        gene_bed = self.ref_gene[species].intersect(sample_bed)
        exon_bed = self.ref_exon[species].intersect(sample_bed)
        print("{} gene peaks, {} exon_peaks".format(len(gene_bed), len(exon_bed)))

        peak_counts = defaultdict(lambda: list((0, 0)))
        # Count genes
        for interval in gene_bed:
            peak_counts[interval.name][0] += 1
        for interval in exon_bed:
            peak_counts[interval.name][1] += 1
        return peak_counts

    def export(self, peak_data, outpath, project=None):
        """Exports calculated peak data to outpath

        if projected, drops those without mappings
        """
        species = self._get_species_ens(next(iter(peak_data)))
        try:
            self.mappings[species]
        except:
            print("No mapping found for species: {}, call the load_mapping method to load".format(species))
            raise
        keys = sorted(peak_data)
        ids, peak_counts, exon_counts = [], [], []
        for ens_id in keys:
            peak_count, exon_count = peak_data[ens_id]
            ids.append(ens_id)
            peak_counts.append(peak_count)
            exon_counts.append(exon_count)
        # Export as pandas
        df = pd.DataFrame({"ens_gene": ids,
                           "peak_counts": peak_counts,
                           "exon_counts": exon_counts})

        # if projection is desired
        if project and project != species:
            print("projecting from {} to {}".format(species, project))
            cross = self.crossmap[[species, project]].rename(columns={species: "ens_gene",
                                                                      project: "ens_gene_projected"})
            df = df.merge(cross, how="left", on="ens_gene")
            df = df.rename(columns={"ens_gene": "junk", "ens_gene_projected": "ens_gene"})
        df = df[["ens_gene", "peak_counts", "exon_counts"]]
        total = len(df)
        # print(df["ens_gene"].head())
        # print(df["ens_gene"].head() != "")
        df = df[df["ens_gene"].notnull()]
        df = df[df["ens_gene"] != ""]
        dropped = total - len(df)
        print("Exporting to: {}, {} entries dropped due to lack of mappers".format(outpath, dropped))
        df.to_csv(outpath, sep="\t", index=False)
        return 0

    def process(self, sample_clip_bed, outpath, project=None):
        print("Processing: ", sample_clip_bed, outpath)
        peaks = self.count_peaks(sample_clip_bed)
        self.export(peaks, outpath, project=project)


class Sample:
    def __init__(self, annotation_row, dataset_dir):
        self.name = annotation_row["sample_name"]
        self.tech = annotation_row["technology"]
        self.species = annotation_row["species"]
        self.tissue = annotation_row["tissue"]
        self.disease = annotation_row["disease"]
        self.material = annotation_row["material"]
        self.batch = annotation_row["batch"]
        self.peak_data = pd.read_csv(os.path.join(dataset_dir,
                                                  os.path.splitext(annotation_row["file_name"])[0] + ".tsv"),
                                     sep="\t")

    def generate_arrays(self, index_map):
        self.gene_peak_arr = [0 for i in index_map]
        self.exon_peak_arr = [0 for i in index_map]
        for index, ens_gene, gene_peaks, exon_peaks in self.peak_data.itertuples():
            self.gene_peak_arr[index_map[ens_gene]] = gene_peaks
            self.exon_peak_arr[index_map[ens_gene]] = exon_peaks



class ClipAnalyzer:
    """Performs the quick analysis of the peak data
    """
    def __init__(self, annotations):
        self.annotations = pd.read_csv(annotations, sep="\t")

    def load_data(self, dataset_dir):
        # check that the annotations and targeted dataset are identical
        # All targets correct
        annotation_tsvs = [os.path.splitext(os.path.basename(file_name))[0] + ".tsv" for file_name in self.annotations["file_name"]]
        try:
            assert(sorted(annotation_tsvs) == sorted((i for i in os.listdir(dataset_dir) if i.endswith(".tsv"))))
        except:
            pprint(sorted(annotation_tsvs))
            pprint(sorted(os.listdir(dataset_dir)))
            raise
        self.sample_data = dict()
        for index, row in self.annotations.iterrows():
            sample = Sample(row, dataset_dir)
            self.sample_data[sample.name] = sample

        # get all genes involved in our dataset
        genes = set()
        for sample_name, sample in self.sample_data.items():
            print(sample.peak_data.columns)
            for ens_gene in sample.peak_data["ens_gene"]:
                genes.add(ens_gene)

        # Generate the arrays for the various samples
        index_map = {item: index for index, item in enumerate(sorted(genes))}
        for sample_name, sample in self.sample_data.items():
            sample.generate_arrays(index_map)

    def PCA_data(self, data_type="gene_peaks"):
        sample_names = sorted([i for i in self.sample_data])
        samples = [self.sample_data[name] for name in sample_names]
        if data_type == "gene_peaks":
            data = np.array([sample.gene_peak_arr for sample in samples])
        elif data_type == "exon_peaks" or data_type == "exon":
            data = np.array([sample.exon_peak_arr for sample in samples])
        pca = PCA(n_components=2)
        pcs = pca.fit(data).transform(data)
        pc1 = [i[0] for i in pcs]
        pc2 = [i[1] for i in pcs]
        pca_data = pd.DataFrame({"sample_name": [sample.name for sample in samples],
                                 "pc1": pc1,
                                 "pc2": pc2})
        print(pca_data)
        pca_data = pca_data.merge(self.annotations, how="left")
        return pca_data

    def PCA_plot(self, pca_data, outdir):
        if not os.path.exists(outdir):
            print("creating directory {} for pca".format(outdir))
            os.makedirs(outdir)
        for hue in ["batch", "species", "tissue"]:
            plt = sns.lmplot("pc1", "pc2", pca_data, hue=hue, fit_reg=False)
            plt.savefig(os.path.join(outdir, "{}.png".format(hue)))
        return 0

    def generate_PCA_plots(self, outdir):
        pca_data = self.PCA_data()
        self.PCA_plot(pca_data, outdir)
        return 0

    def generate_cluster_plots(self, outdir):
        if not os.path.exists(outdir):
            print("creating directory {} for cluster".format(outdir))
            os.makedirs(outdir)
        # Load the PCA data just for categories
        pca_data = self.PCA_data()

        sample_names = sorted([i for i in self.sample_data])
        samples = [self.sample_data[name] for name in sample_names]
        cluster_data = np.array([sample.gene_peak_arr for sample in samples])
        for category in ["batch", "species", "tissue"]:
            color_map = dict(zip(pca_data[category], sns.color_palette("Paired", len(pca_data[category]))))
            colors = list(pca_data[category].map(color_map))
            plt = sns.clustermap(cluster_data, xticklabels=[], yticklabels=list(pca_data[category]), row_colors=colors)
            plt.savefig(os.path.join(outdir, "{}.png".format(category)))
        return 0

    def generate_all_plots(self, plot_outdir):
        print("GENERATING PLOTS IN: {}".format(plot_outdir))
        self.generate_PCA_plots(os.path.join(plot_outdir, "pca"))
        self.generate_cluster_plots(os.path.join(plot_outdir, "cluster"))


if __name__ == "__main__":
    pass
    # cf.load_mapping("/mnt/d/work/sc/data/ensembl_mmusculus_mappings_2017_11_16.tsv")
    # peaks = cf.analyze("/mnt/d/work/sc/workflows/runs/clipdb_fus/input/FUS_Piranha/FUS_E-MTAB-1223-ERR208898.bed")
    # cf.export(peaks, "/mnt/d/work/temp/tempout.tsv")
    # pca_data = ca.PCA_data()

    # cluster_data = list(zip(pca_data["pc1"], pca_data["pc2"]))
    # print(cluster_data)
    # category = "batch"
    # color_map = dict(zip(pca_data[category], sns.color_palette("hus1", len(pca_data[category]))))
    # colors = list(pca_data[category].map(color_map))
    # print(colors)
    # plt = sns.clustermap(cluster_data, xticklabels=["pc1", "pc2"], yticklabels=list(pca_data[category]), row_colors=colors)
    # plt.savefig(os.path.join("temp", "cluster.png"))



