# Files in `Code` folder

00_experiments.R
: Input is `fit.eb`.
: Output: Venn diagrams (png) comparing all conditions

00_covariances.R
: Input is `fit.eb`.
: Output: graphs in Rdata format

01_cancer_rate.R
: IN: `samplelist.tsv`.
: Out: Plot of breast cancer rate given diet, given age:diet and combinig ICR

02_read_and_normalize.R
: in: `intermitent_sdrf_new.tsv`+ array data in cel format.
: This is the file we developped together in July 2021
: Out: `full.Data` with raw expressions, `norm.Data` with normalized expr, `mouse.Data` with normalized data for mouse only. All `RData` files have the same `features` table.

03_diet_age_tissue_model.R
: Input is ~~`norm.Data`~~, `mouse.Data`.
: Merge labels for ICRR+ICRRF -> ICR. Fits `expr ~ diet:age:tissue + 0` and does contrasts vs BASELINE
: Output is `fit.eb`

04_pca.R
: Input is `full.Data`, `mouse.Data`.
: Principal Component Analysis of both (raw and normalized data) cases.

05_test_correlations.R
: Input is `mouse.Data`.
: Correlation between conditions/tissues for the subset of animals with samples in all tissues. Only Blood and MFP are highly correlated
: out: plots

05.5-heatmap.R
: Input is `mouse.Data`.
: All conditions, all miRNA, no contrast, no DE.

06_count_DE.R
: Input is `fit.eb`.
: List of topTable for each contrast. Count nrow on each.
: Out: `DE.table` (list of data frames, one for each contrast), `DE.number` (nrow for each data frame), `gnames` (list of vectors, with names of each DE miRNA in each condition), `gnames.unique` (vector of miRNA names, combining all `genomes`)

07_venn_table.R
: Input is `fit.eb`. Makes list of tables `venn.diet` (for each tissue and age), `venn.tissue` for each diet and age.

07.1_venn_diagrams.R:
: Input is `fit.eb`. Makes Venn diagrams (in png) for each diet given tissue and age, and for each tissue given diet and age.

08_ebic_glasso.R
: Input is `mouse.Data`, `gnames.unique`.
: Uses glasso to find precision matrix and its graph, and leiden clustering of graph
: Out: glasso graph in `consensus`, igraph object in `g`, and clustering in `lei`

08_joint-glasso.R
: Input is `mouse.Data`, `gnames.unique`.
: Makes differential graphs using Fused Graphical Lasso (FGL), then calculates several centrality indices.
: Differential graphs comparing CR diets vs AL for each tissue
: Output is a list of lists of FGL objects named `d.fgl`

08.5_graphs.R
: Input is `d.fgl`
: Makes network plots and saves png image.
: Writes Png plots. Saves `communities`

08.75_distribution.R
: Input is `d.fgl`.
: Makes plots of degree distrubiton of differential graphs and saves PNG images.

09_find_targets.R
: Input is `gnames` and `gnames.unique`.
: Finds targets of DE miRNAs.
Output is `targets.unique` and `targets_list` (predicted+validated target list)

10_network-corporation.R
: gathers various sources of interactions.
: miRNA - miRNA connection -> output of `joint_glasso\ebic_glasso`.
: miRNA - gene connection ->  output of `find_targets`.
: miRNA - TF connection -> from public databases `Transmir`
: Gene - gene connection -> from public databases `MINT` and `BioGRID`.
: Gene - TF connection -> from public database `TRRUST`

11_miRNA_enrichment.R
: input is `gnames`. makes enrichment analysis for every contrast.

auth.R
: authorization for gdrive. Deprecated as locale is turkish

chart.Correlation.R
: modified function of correlation

cluster.R
: Input is `mouse.Data`.  ICA clustering on normalized expression.

contrasts.R
: all posible contrasts

covariance.R
covariance_plots_20211006.R
glasso vs huge.R
glasso.R
go enrichment.R
go.R
go_tissue.R
heatmap.R
ica.R
network analysis example.R
network.R
Old Examples
pcs.R
README.md
target gene.R
venn.R
