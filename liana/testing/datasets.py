import scanpy as sc

def olmer_2022():
    """
    Utility function to load the data HuBMAP Processed Data.

    scRNA-seq of Human Knee Articular Cartilage and Meniscus 
    scRNAseq (10x Genomics v3)
    
    HBM2333.CCCX.767 Data Products Secondary Analysis
    Normalized gene expression with additional metadata, in HDF5 format, readable with the 'anndata' Python package (Format: EDAM_1.24.format_3590)
    
    Returns
    -------

    Returns a largely pre-processed AnnData object with the following attributes:
    Raw counts for  cells;  genes; x samples; y conditions.
    """
    adata = sc.read("secondary_analysis.h5ad", backup_url="https://assets.hubmapconsortium.org/c3e09a5698adbd3308d75c91eb7ff802/secondary_analysis.h5ad")

    # Store the counts for later use
    adata.layers["counts"] = adata.X.copy()
    # Rename label to condition, replicate to patient
    adata.obs = adata.obs.rename({"label": "condition", "replicate": "donor"}, axis=1)

    # assign sample
    adata.obs["sample"] = (adata.obs["condition"].astype("str") + "&" + adata.obs["donor"].str.slice(8, 13))

    # set cell_types abbreviations (recommended given MOFA appends names)
    abbreviations = {'CD4 T cells':'CD4T',
                     'B cells':'B',
                     'NK cells':'NK',
                     'CD8 T cells':'CD8T',
                     'FCGR3A+ Monocytes':'FGR3',
                     'CD14+ Monocytes':'CD14',
                     'Dendritic cells':'DCs',
                     'Megakaryocytes':'Mega',
                     'Chondrocytes':'Chondrocytes'}
    adata.obs["cell_abbr"] = adata.obs["cell_type"].replace(abbreviations)

    return adata
""" 
def kang_2018():

    Utility function to load the data from Kang et al., 2018; GSE96583.
    The data contains ~25k PBMCs cells from 8 pooled patient lupus samples, each before and after IFN-beta stimulation.
    Kang, H., Subramaniam, M., Targ, S. et al. Multiplexed droplet single-cell RNA-sequencing using natural genetic variation. Nat Biotechnol 36, 89–94 (2018). https://doi.org/10.1038/nbt.4042

    The dataset was preprocessed for and is available via pertpy (https://github.com/theislab/pertpy; Heumos et al., In prep.).

    Returns
    -------

    Returns a largely pre-processed AnnData object with the following attributes:
    Raw counts for ~25k cells; ~15k genes; 16 samples; 2 conditions.

    adata = sc.read("kang_counts_25k.h5ad", backup_url="https://figshare.com/ndownloader/files/34464122")

    # Store the counts for later use
    adata.layers["counts"] = adata.X.copy()
    # Rename label to condition, replicate to patient
    adata.obs = adata.obs.rename({"label": "condition", "replicate": "patient"}, axis=1)

    # assign sample
    adata.obs["sample"] = (adata.obs["condition"].astype("str") + "&" + adata.obs["patient"].str.slice(8, 13))

    # set cell_types abbreviations (recommended given MOFA appends names)
    abbreviations = {'CD4 T cells':'CD4T',
                     'B cells':'B',
                     'NK cells':'NK',
                     'CD8 T cells':'CD8T',
                     'FCGR3A+ Monocytes':'FGR3',
                     'CD14+ Monocytes':'CD14',
                     'Dendritic cells':'DCs',
                     'Megakaryocytes':'Mega'}
    adata.obs["cell_abbr"] = adata.obs["cell_type"].replace(abbreviations)

    return adata


"""