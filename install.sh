#sudo apt-get update

echo "Adding conda channels..."
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

echo "Installing python..."
conda install -y python

echo "Installing python modules..."
conda install -y biopython numpy pandas seaborn bokeh

echo "Installing other tools..."
conda install -y bwa bedtools samtools mafft

# http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/ for gff3togenepred
# http://www.openbioinformatics.org/annovar/annovar_download_form.php annovar