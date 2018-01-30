#! /usr/bin/env python

from clip_analyze import ClipFinder

cf = ClipFinder({"mouse": "/mnt/d/work/sc/data/ensembl_mmus_2018_01_03_ucsc-chroms.bed",
                 "human": "/mnt/d/work/sc/data/ensembl_hsapiens_2018_01_03_ucsc-chroms.bed"},
                {"mouse": "/mnt/d/work/sc/data/ensembl_mmus_exon_2018_01_03_ucsc-chroms.bed",
                 "human": "/mnt/d/work/sc/data/ensembl_hsapiens_exon_2018_01_03_ucsc-chroms.bed"},
                annotations="/mnt/d/work/sc/workflows/runs/clipdb_fus/input/annotations.tsv",
                crossmap="/mnt/d/work/sc/data/ensembl_human_mappings_2018_01_05.tsv")
cf.load_mappings({"mouse": "/mnt/d/work/sc/data/ensembl_mmusculus_mappings_2017_11_16.tsv",
                  "human": "/mnt/d/work/sc/data/ensembl_human_2018_01_05.tsv"})
cf.load_mappings({"mouse": "/mnt/d/work/sc/data/ensembl_mmusculus_mappings_2017_11_16.tsv",
                  "human": "/mnt/d/work/sc/data/ensembl_human_2018_01_05.tsv"})
# peaks = cf.count_peaks("/mnt/d/work/sc/workflows/runs/clipdb_fus/input/FUS_Piranha/FUS_E-MTAB-1223-ERR208898.bed")
# cf.export(peaks, "/mnt/d/work/temp/tempout.tsv")
# cf.export(peaks, "/mnt/d/work/temp/tempout_human.tsv", project="human")
peaks = cf.count_peaks("/mnt/d/work/sc/workflows/runs/clipdb_fus/temp/compiled/specialized/FUS_SRX029328.bed")
cf.export(peaks, "temp/exported_peaks.tsv")

cf.process("/mnt/d/work/sc/workflows/runs/clipdb_fus/temp/compiled/specialized/FUS_SRX029328.bed", "temp/processed_test.tsv", "mouse")