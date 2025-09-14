scShinyDash: An Interactive Dashboard for Single-Cell RNA-Seq Analysis
scShinyDash is a user-friendly, web-based application built with R and Shiny for the exploratory analysis of single-cell RNA-sequencing (scRNA-Seq) data. This dashboard provides a complete workflow, from data ingestion and quality control to downstream analysis and interactive visualization, all without needing to write any code.

Features
This application is organized into four main analysis modules:

Data Upload & Processing:

Upload your own Seurat object (.rds file).

Alternatively, load a pre-packaged example dataset for a quick demonstration.

Automatic Validation & Processing: The app intelligently checks if the data has been normalized and if UMAP embeddings exist. It automatically performs these steps if they are missing, ensuring the data is ready for analysis.

View and interact with the cell metadata in a searchable, sortable table.

Differential Gene Expression (DGE):

Perform DGE analysis between any two cell groups defined in your metadata (e.g., compare cluster 1 vs. cluster 2).

Visualize results with an interactive volcano plot, allowing you to dynamically set p-value and log-fold-change thresholds.

View detailed DGE results in a data table.

Pathway Enrichment Analysis:

Directly use the significant genes from the DGE analysis to perform Gene Ontology (GO) enrichment.

Identify the biological processes, molecular functions, and cellular components that are over-represented in your gene list.

Visualize enrichment results with a dot plot and browse them in a table.

Interactive UMAP Visualization:

Explore the UMAP embedding of your cells.

Color the plot by any metadata variable (e.g., cluster, cell type, treatment).

Highlight one or more specific cell populations for focused analysis.

Instantly view a table of the top 10 marker genes for any selected cluster.

Installation & Setup
Follow these steps to set up and run the application on your local machine.

1. Prerequisites

R: Version 4.2 or newer. You can download it from the Comprehensive R Archive Network (CRAN).

RStudio: (Recommended) An integrated development environment for R. Download the free desktop version from Posit.

2. Clone the Repository

Open your terminal or Git Bash and clone this GitHub repository:

git clone <URL_to_your_GitHub_repository>
cd scRNA_dashboard

3. Install Required R Packages

Open R or RStudio and run the following script in the console. This will install all the necessary packages from CRAN and Bioconductor.

# --- Install packages from CRAN ---
install.packages(c("shiny", "Seurat", "ggplot2", "shinythemes", "DT", "plotly", "dplyr"))

# --- Install packages from Bioconductor ---
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "enrichplot", "org.Hs.eg.db"))

How to Run the Application
Set Working Directory: In RStudio, open the app.R file. Your working directory should automatically be set to the project's root folder. Alternatively, you can manually set it via Session > Set Working Directory > To Source File Location.

Launch the App: With app.R open, click the "Run App" button in the top-right corner of the editor pane, or run the following command in the R console:

shiny::runApp()

The application will launch in a new window or in your default web browser.

Using the Test Data

To test the application's full functionality, you can use the provided test dataset.

Download Test Data: pbmc_seurat_v5_processed.rds (Google Drive)

Usage: Once the app is running, click the "Choose File" button on the "1. Data Upload & Metadata" tab and select the .rds file you just downloaded.