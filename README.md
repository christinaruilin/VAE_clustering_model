# VAE_clustering_model
This repository hosts a Variational Autoencoder (VAE) pipeline tailored for clustering CUT&Tag data, from preprocessing through visualization.

Project layout:

- `New_preprocessing_method_train01`: curated cells, fragment summaries, and model-ready inputs, together with VAE run outputs for downstream UMAP embeddings.
- `VAEModel`: Jupyter notebooks that cover data preprocessing, model definition, and training.
- `r`: R scripts for visualization plus the initial preprocessing steps that build the Seurat object and export the TF-IDF matrix.

Overall workflow: preprocess CUT&Tag counts in R to generate the Seurat object and TF-IDF matrix, prepare inputs with `preprocessing.ipynb`, then train and visualize the VAE embeddings in the accompanying notebooks.
