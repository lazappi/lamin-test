import lamindb as ln
import lnschema_bionty as lb
import scanpy as sc
from urllib.request import urlopen
import os.path
import os
import ssl
import certifi
import gzip
import shutil

lb.settings.species = "human"

transform = ln.Transform(name="download_Organoid123")
ln.track(transform)

base_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE114nnn/GSE114802/suppl/"
ssl_context = ssl.create_default_context(cafile=certifi.where())

data_dir = "data/00-raw/Organoid123"
print(f"Downloading dataset from '{base_url}'...")
if not os.path.exists(data_dir):
   os.makedirs(data_dir)

print(f"Downloading genes...")
genes_url = base_url + "GSE114802_org_genes.tsv.gz"
genes_path = os.path.join(data_dir, "genes.tsv")
with urlopen(genes_url, context=ssl_context) as url_file:
    with open(genes_path + ".gz", "wb") as out_file:
        out_file.write(url_file.read())
    with gzip.open(genes_path + ".gz", "rb") as gz_file:
        with open(genes_path, 'wb') as unzipped_file:
            shutil.copyfileobj(gz_file, unzipped_file)
    os.remove(genes_path + ".gz")
genes_file = ln.File(genes_path, key=genes_path)

print(f"Downloading barcodes...")
barcodes_url = base_url + "GSE114802_org_barcodes.tsv.gz"
barcodes_path = os.path.join(data_dir, "barcodes.tsv")
with urlopen(barcodes_url, context=ssl_context) as url_file:
    with open(barcodes_path + ".gz", "wb") as out_file:
        out_file.write(url_file.read())
    with gzip.open(barcodes_path + ".gz", "rb") as gz_file:
        with open(barcodes_path, 'wb') as unzipped_file:
            shutil.copyfileobj(gz_file, unzipped_file)
    os.remove(barcodes_path + ".gz")
barcodes_file = ln.File(barcodes_path, key=barcodes_path)

print(f"Downloading matrix...")
matrix_url = base_url + "GSE114802_org_matrix.mtx.gz"
matrix_path = os.path.join(data_dir, "matrix.mtx")
with urlopen(matrix_url, context=ssl_context) as url_file:
    with open(matrix_path + ".gz", "wb") as out_file:
        out_file.write(url_file.read())
    with gzip.open(matrix_path + ".gz", "rb") as gz_file:
        with open(matrix_path, 'wb') as unzipped_file:
            shutil.copyfileobj(gz_file, unzipped_file)
    os.remove(matrix_path + ".gz")
matrix_file = ln.File(matrix_path, key=matrix_path)

print(f"Creating AnnData object...")
adata = sc.read_10x_mtx(data_dir, var_names="gene_ids")
adata_file = ln.File.from_anndata(
    adata,
    field=lb.Gene.ensembl_gene_id,
    description="Organoid123 raw"
)

print(f"Saving AnnData object...")
adata_file.save()

print("Done!")
