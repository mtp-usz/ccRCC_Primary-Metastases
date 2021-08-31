#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
#setwd("/Users/abdullah/science/usz/research/projects/kidneyPriMet/")
library(ComplexHeatmap)

read.file <- function(fileName, row.names = 1) {
  data = read.table(fileName, header = T, sep="\t", na.strings = "-")
  classes = sapply(data, class)
  data = read.table(fileName, header = T, sep="\t", colClasses = classes, row.names = row.names, na.strings = "-")
  return(data)
}

#file.all      = "kidney_priMet_oncoprint"
file.3perc    = "kidney_priMet_oncoprint_svnSplit_top3perc"
file.3percVus = "kidney_priMet_oncoprint_svnSplit_vus_top3perc"
file.snv      = "kidney_priMet_oncoprint_svnSplit"
file.phenoPri = "samplesOrderedByImmunoOfPrimary.txt"
file.phenoMet = "samplesOrderedByImmunoOfFirstMetastase.txt"

#data.all      = read.file(paste0(file.all, ".tsv"))
data.3perc    = read.file(paste0(file.3perc, ".tsv"))
data.3percVUS = read.file(paste0(file.3percVus, ".tsv"))
data.snv      = read.file(paste0(file.snv, ".tsv"))
data.phenoPri = read.file(file.phenoPri)
data.phenoMet = read.file(file.phenoMet)

col = c("P_SNV" = "green",
        "P_AMP" = "red",
        "P_LOSS" = "blue",
        "P_FUSION" = "purple",
        "P_Missense" = "green",
        "P_Nonsense" = "springgreen3",
        "P_Frameshift" = "orange",
        "P_Indel" = "darkorange3",
        "P_Splicesite" = "yellowgreen",
        "P_Promoter" = "darkolivegreen",
        "M1_1_SNV" = "green",
        "M1_1_AMP" = "red",
        "M1_1_LOSS" = "blue",
        "M1_1_FUSION" = "purple",
        "M1_1_Missense" = "green",
        "M1_1_Nonsense" = "springgreen3",
        "M1_1_Frameshift" = "orange",
        "M1_1_Indel" = "darkorange3",
        "M1_1_Splicesite" = "yellowgreen",
        "M1_1_Promoter" = "darkolivegreen",
        "M2_1_SNV" = "green",
        "M2_1_AMP" = "red",
        "M2_1_LOSS" = "blue",
        "M2_1_FUSION" = "purple",
        "M2_1_Missense" = "green",
        "M2_1_Nonsense" = "springgreen3",
        "M2_1_Frameshift" = "orange",
        "M2_1_Indel" = "darkorange3",
        "M2_1_Splicesite" = "yellowgreen",
        "M2_1_Promoter" = "darkolivegreen",
        "M2_2_SNV" = "green",
        "M2_2_AMP" = "red",
        "M2_2_LOSS" = "blue",
        "M2_2_FUSION" = "purple",
        "M2_2_Missense" = "green",
        "M2_2_Nonsense" = "springgreen3",
        "M2_2_Frameshift" = "orange",
        "M2_2_Indel" = "darkorange3",
        "M2_2_Splicesite" = "yellowgreen",
        "M2_2_Promoter" = "darkolivegreen",
        "desert" = "black",
        "excluded" = "blue",
        "excluded/inflamed" = "blue",
        "inflamed" = "red",
        "inflamed/excluded" = "red",
        "unknown" = "grey")

alter_fun = list(
  background = function(x, y, w, h) {
    grid.polygon(
      unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
      unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
      gp = gpar(fill = "grey", col = "white"))
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
      unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
      gp = gpar(fill = "grey", col = "white"))
  },
  P_SNV = function(x, y, w, h) {
    grid.polygon(
      unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
      unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
      gp = gpar(fill = col["P_SNV"], col = "white"))
  },
  M1_1_SNV = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
      unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
      gp = gpar(fill = col["M1_1_SNV"], col = "white"))
  },
  M2_1_SNV = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x), 
      unit.c(y + 0.5*h, y - 0.5*h, y),
      gp = gpar(fill = col["M2_1_SNV"], col = "white"))
  },
  M2_2_SNV = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x, x - 0.5*w), 
      unit.c(y + 0.5*h, y, y + 0.5*h),
      gp = gpar(fill = col["M2_2_SNV"], col = "white"))
  },
  P_AMP = function(x, y, w, h) {
    grid.polygon(
      unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
      unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
      gp = gpar(fill = col["P_AMP"], col = "white"))
  },
  M1_1_AMP = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
      unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
      gp = gpar(fill = col["M1_1_AMP"], col = "white"))
  },
  M2_1_AMP = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x), 
      unit.c(y + 0.5*h, y - 0.5*h, y),
      gp = gpar(fill = col["M2_1_AMP"], col = "white"))
  },
  M2_2_AMP = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x, x - 0.5*w), 
      unit.c(y + 0.5*h, y, y + 0.5*h),
      gp = gpar(fill = col["M2_2_AMP"], col = "white"))
  },
  P_LOSS = function(x, y, w, h) {
    grid.polygon(
      unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
      unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
      gp = gpar(fill = col["P_LOSS"], col = "white"))
  },
  M1_1_LOSS = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
      unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
      gp = gpar(fill = col["M1_1_LOSS"], col = "white"))
  },
  M2_1_LOSS = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x), 
      unit.c(y + 0.5*h, y - 0.5*h, y),
      gp = gpar(fill = col["M2_1_LOSS"], col = "white"))
  },
  M2_2_LOSS = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x, x - 0.5*w), 
      unit.c(y + 0.5*h, y, y + 0.5*h),
      gp = gpar(fill = col["M2_2_LOSS"], col = "white"))
  },
  M1_1_FUSION = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
      unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
      gp = gpar(fill = col["M1_1_FUSION"], col = "white"))
  },
  M2_1_FUSION = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x), 
      unit.c(y + 0.5*h, y - 0.5*h, y),
      gp = gpar(fill = col["M2_1_FUSION"], col = "white"))
  },
  M2_2_FUSION = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x, x - 0.5*w), 
      unit.c(y + 0.5*h, y, y + 0.5*h),
      gp = gpar(fill = col["M2_2_FUSION"], col = "white"))
  },
  P_Missense = function(x, y, w, h) {
    grid.polygon(
      unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
      unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
      gp = gpar(fill = col["P_Missense"], col = "white"))
  },
  M1_1_Missense = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
      unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
      gp = gpar(fill = col["M1_1_Missense"], col = "white"))
  },
  M2_1_Missense = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x), 
      unit.c(y + 0.5*h, y - 0.5*h, y),
      gp = gpar(fill = col["M2_1_Missense"], col = "white"))
  },
  M2_2_Missense = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x, x - 0.5*w), 
      unit.c(y + 0.5*h, y, y + 0.5*h),
      gp = gpar(fill = col["M2_2_Missense"], col = "white"))
  },
  P_Nonsense = function(x, y, w, h) {
    grid.polygon(
      unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
      unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
      gp = gpar(fill = col["P_Nonsense"], col = "white"))
  },
  M1_1_Nonsense = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
      unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
      gp = gpar(fill = col["M1_1_Nonsense"], col = "white"))
  },
  M2_1_Nonsense = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x), 
      unit.c(y + 0.5*h, y - 0.5*h, y),
      gp = gpar(fill = col["M2_1_Nonsense"], col = "white"))
  },
  M2_2_Nonsense = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x, x - 0.5*w), 
      unit.c(y + 0.5*h, y, y + 0.5*h),
      gp = gpar(fill = col["M2_2_Nonsense"], col = "white"))
  },
  P_Frameshift = function(x, y, w, h) {
    grid.polygon(
      unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
      unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
      gp = gpar(fill = col["P_Frameshift"], col = "white"))
  },
  M1_1_Frameshift = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
      unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
      gp = gpar(fill = col["M1_1_Frameshift"], col = "white"))
  },
  M2_1_Frameshift = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x), 
      unit.c(y + 0.5*h, y - 0.5*h, y),
      gp = gpar(fill = col["M2_1_Frameshift"], col = "white"))
  },
  M2_2_Frameshift = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x, x - 0.5*w), 
      unit.c(y + 0.5*h, y, y + 0.5*h),
      gp = gpar(fill = col["M2_2_Frameshift"], col = "white"))
  },
  P_Indel = function(x, y, w, h) {
    grid.polygon(
      unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
      unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
      gp = gpar(fill = col["P_Indel"], col = "white"))
  },
  M1_1_Indel = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
      unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
      gp = gpar(fill = col["M1_1_Indel"], col = "white"))
  },
  M2_1_Indel = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x), 
      unit.c(y + 0.5*h, y - 0.5*h, y),
      gp = gpar(fill = col["M2_1_Indel"], col = "white"))
  },
  M2_2_Indel = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x, x - 0.5*w), 
      unit.c(y + 0.5*h, y, y + 0.5*h),
      gp = gpar(fill = col["M2_2_Indel"], col = "white"))
  },
  P_Splicesite = function(x, y, w, h) {
    grid.polygon(
      unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
      unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
      gp = gpar(fill = col["P_Splicesite"], col = "white"))
  },
  M1_1_Splicesite = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
      unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
      gp = gpar(fill = col["M1_1_Splicesite"], col = "white"))
  },
  M2_1_Splicesite = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x), 
      unit.c(y + 0.5*h, y - 0.5*h, y),
      gp = gpar(fill = col["M2_1_Splicesite"], col = "white"))
  },
  M2_2_Splicesite = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x, x - 0.5*w), 
      unit.c(y + 0.5*h, y, y + 0.5*h),
      gp = gpar(fill = col["M2_2_Splicesite"], col = "white"))
  },
  P_Promoter = function(x, y, w, h) {
    grid.polygon(
      unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
      unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
      gp = gpar(fill = col["P_Promoter"], col = "white"))
  },
  M1_1_Promoter = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
      unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
      gp = gpar(fill = col["M1_1_Promoter"], col = "white"))
  },
  M2_1_Promoter = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x), 
      unit.c(y + 0.5*h, y - 0.5*h, y),
      gp = gpar(fill = col["M2_1_Promoter"], col = "white"))
  },
  M2_2_Promoter = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x, x - 0.5*w), 
      unit.c(y + 0.5*h, y, y + 0.5*h),
      gp = gpar(fill = col["M2_2_Promoter"], col = "white"))
  }
)

# Legend + title for Oncoprint
column_title = "OncoPrint for Kidney Primary-Metastses Cohort"
heatmap_legend_param = list(title = "Alternations", at = c("P_SNV", "P_AMP", "P_LOSS", "P_FUSION", "M1_1_SNV", "M1_1_AMP", "M1_1_LOSS", "M1_1_FUSION"), 
                            labels = c("Primary Short Variance", "Primary Amplification", "Primary Homozygous Loss", "Primary Rearrangement", "Metastasis Short Variance", "Metastasis Amplification", "Metastasis Homozygous Loss", "Metastasis Rearrangement"))
heatmap_legend_param.snv = list(title = "Alternations", at = c("P_Missense", "P_Nonsense", "P_Splicesite", "P_Promoter", "P_Indel", "P_Frameshift", "P_AMP", "P_LOSS", "P_FUSION", 
                                                               "M1_1_Missense", "M1_1_Nonsense", "M1_1_Splicesite", "M1_1_Promoter", "M1_1_Indel", "M1_1_Frameshift", "M1_1_AMP", "M1_1_LOSS", "M1_1_FUSION"), 
                            labels = c("Primary Missense", "Primary Nonsense", "Primary Splice Site", "Primary Promoter", "Primary InDel", "Primary Frameshift", "Primary Amplification", "Primary Homozygous Loss", "Primary Rearrangement",
                                       "Metastasis Missense", "Metastasis Nonsense", "Metastasis Splice Site", "Metastasis Promoter", "Metastasis InDel", "Metastasis Frameshift", "Metastasis Amplification", "Metastasis Homozygous Loss", "Metastasis Rearrangement"))

# plot Oncoprints
#pdf(paste0(file.all, ".pdf"), width=20, height=10)
#oncoPrint(data.all,
#          alter_fun = alter_fun,
#          col = col, 
#          column_title = column_title,
#          show_column_names = TRUE,
#          column_order = row.names(data.phenoPri),
#          heatmap_legend_param = heatmap_legend_param,
#          top_annotation = HeatmapAnnotation(MetPheno = data.phenoMet[colnames(data.all),], 
#                                             PriPheno = data.phenoPri[colnames(data.all),],
#                                             column_barplot = anno_oncoprint_barplot(), 
#                                             col=list(MetPheno = col, PriPheno = col)))
#dev.off()

pdf(paste0(file.snv, "_sortByImmunoOfPrimary.pdf"), width=20, height=10)
oncoPrint(data.snv,
          alter_fun = alter_fun,
          col = col,
          column_order = row.names(data.phenoPri),
          column_title = column_title,
          heatmap_legend_param = heatmap_legend_param.snv,
          show_column_names = TRUE,
          top_annotation = HeatmapAnnotation(MetPheno = data.phenoMet[colnames(data.snv),], 
                                             PriPheno = data.phenoPri[colnames(data.snv),], 
                                             cbar = anno_oncoprint_barplot(), 
                                             col=list(MetPheno = col, PriPheno = col)))
dev.off()

pdf(paste0(file.snv, "_sortByImmunoOfFirstMetastase.pdf"), width=20, height=10)
oncoPrint(data.snv,
          alter_fun = alter_fun,
          col = col,
          column_order = row.names(data.phenoMet),
          column_title = column_title,
          show_column_names = TRUE,
          heatmap_legend_param = heatmap_legend_param.snv,
          top_annotation = HeatmapAnnotation(MetPheno = data.phenoMet[colnames(data.snv),], 
                                             PriPheno = data.phenoPri[colnames(data.snv),], 
                                             cbar = anno_oncoprint_barplot(), 
                                             col=list(MetPheno = col, PriPheno = col))) 
dev.off()

# add Fusion 
col = c(col, "P_FUSION" = "purple")
alter_fun = c(alter_fun, 
               P_FUSION = function(x, y, w, h) {
                    grid.polygon(
                    unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
                    unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
                    gp = gpar(fill = col["P_FUSION"], col = "white"))
                }
              )
heatmap_legend_param.snv$at = c(heatmap_legend_param.snv$at,"P_FUSION")
heatmap_legend_param.snv$labels = c(heatmap_legend_param.snv$labels, "Primary Rearrangement")

pdf(paste0(file.3percVus, "_sortByImmunoOfPrimary.pdf"), width=20, height=10)
oncoPrint(data.3percVUS,
          alter_fun = alter_fun,
          col = col,
          column_order = row.names(data.phenoPri),
          column_title = column_title,
          heatmap_legend_param = heatmap_legend_param.snv,
          show_column_names = TRUE,
          top_annotation = HeatmapAnnotation(MetPheno = data.phenoMet[colnames(data.3percVUS),], 
                                             PriPheno = data.phenoPri[colnames(data.3percVUS),], 
                                             cbar = anno_oncoprint_barplot(), 
                                             col=list(MetPheno = col, PriPheno = col)))
dev.off()
