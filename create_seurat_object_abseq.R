library(Seurat)

# This is an example for merging two separate cartridges. 
# It creates a Seurat object for each cartridge and adds the necessary library information (like SMK).
# Seurat also handles suffixes to the cell labels to prevent collisions.


# First read in the UMI count data from the RSEC_MolsPerCell.csv files (or DBEC_MolsPerCell if preferred)
first <- read.table('cart1_RSEC_MolsPerCell.csv', sep = ',', header = TRUE, row.names = 1)
second <- read.table('cart2_RSEC_MolsPerCell.csv', sep = ',', header = TRUE, row.names = 1)


# Transpose the dataframes using the `t()` function so that the column names are individual cell IDs and the row names are gene names or AbSeqs.
first.tp <- data.frame(t(first), check.names = FALSE)
second.tp <- data.frame(t(second), check.names = FALSE)


# Separate AbSeqs and RNA into separate dataframes
num.abseqs <- 30
first.adt <- first.tp[c(1:num.abseqs),]
first.rna <- first.tp[c((num.abseqs+1):length(rownames(first.tp))),]

second.adt <- second.tp[c(1:num.abseqs),]
second.rna <- second.tp[c((num.abseqs+1):length(rownames(second.tp))),]


# Create the Seurat object using the transposed RNA table 
first.data <- CreateSeuratObject(counts = first.rna, project='cartridge1')
second.data <- CreateSeuratObject(counts = second.rna, project='cartridge2')


# And then add the AbSeq data as a second assay
first.abs.assay <- CreateAssayObject(counts = first.adt)
first.data[["ADT"]] <- first.abs.assay

second.abs.assay <- CreateAssayObject(counts = second.adt)
second.data[["ADT"]] <- second.abs.assay


# Import and add the SMK data table to the Seurat metadata table if applicable
first.smk <- read.table('cart1_Sample_Tag_Calls.csv', sep = ',', header = TRUE, row.names = 1)
second.smk <- read.table('cart2_Sample_Tag_Calls.csv', sep = ',', header = TRUE, row.names = 1)

first.data <- AddMetaData(first.data, metadata = first.smk)
second.data <- AddMetaData(second.data, metadata = second.smk)


# Import and add experimental cell type classifications if applicable
first.cell <- read.table('cart1_cell_type_experimental.csv', sep = ',', header = TRUE, row.names = 1)
second.cell <- read.table('cart2_cell_type_experimental.csv', sep = ',', header = TRUE, row.names = 1)

first.data <- AddMetaData(first.data, metadata = first.cell)
second.data <- AddMetaData(second.data, metadata = second.cell)


# Import and add VDJ data if applicable
first.vdj <- read.table('cart1_VDJ_perCell.csv', sep=',', header = TRUE, row.names = 1)
second.vdj <- read.table('cart1_VDJ_perCell.csv', sep=',', header = TRUE, row.names = 1)

first.data <- AddMetaData(first.data, metadata = first.vdj)
second.data <- AddMetaData(second.data, metadata = second.vdj)


# Combine the objects into one
combined.data <- merge(first.data, y=second.data, add.cell.ids=c('cart1', 'cart2'))


# If you have more than two cartridges, the rest can be added as a list in the y parameter
# combined.data <- merge(first.data, y=c(second.data, third.data), add.cell.ids=c('cart1', 'cart2', 'cart3'))

# More info on merging seurat objects: https://satijalab.org/seurat/articles/merge_vignette.html


# save object as RDS file
saveRDS(combined.data, file='projectname_abseq_seurat.rds')
