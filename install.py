# install sript for all dependencies except for annovar
import platform
import sys
import subprocess
import os

if __name__ == '__main__':
	# TODO: have a flag for install and a flag for checking dependencies in the same way that lava does 

	print('Please note that this scirpt assumes you have python installed on your system and that you are executing this file from the main lava folder')
	print('If running on Mac OS you also need to have brew properly installed')
	print('currently you are running this install script from:')
	subprocess.call('pwd', shell=True)

	# for super fresh OS systems]
	if platform.system == 'Darwin':
		print('Installing wget')
		subprocess.call('brew install wget', shell = True)

	# Either install pip or upgrade it to the latest version 
	print('Installing pip')
	subprocess.call('wget https://bootstrap.pypa.io/get-pip.py', shell=True)
	subprocess.call('python get-pip.py', shell=True)

	print('Installing python modules...')
	print('First we\'re going to update pip and some other setuptools')
	# Just in case 
	subprocess.call('python -m pip install --upgrade pip setuptools wheel', shell=True)
	print('Installing biopython...')
	subprocess.call('python -m pip install biopython', shell=True)
	print('Installing numpy...')
	subprocess.call('python -m pip install numpy',shell =True)
	print('Installing pandas...')
	subprocess.call('python -m pip install pandas', shell=True)
	print('Installing seaborn...')
	subprocess.call('python -m pip install seaborn', shell=True)
	print('Installing bokeh..')
	subprocess.call('python -m pip install bokeh', shell=True)
	print('Python modules installed')

	
	# TODO: add detection to only do this if picard.jar isn't found in the main lava folder 
	if not os.path.isfile('picard.jar'):
		print('Downloading Picard')
		subprocess.call('wget https://github.com/broadinstitute/picard/releases/download/2.18.15/picard.jar', shell=True)
	if not os.path.isfile('gatk-4.0.11.0.zip'):
		print('Downloading GATK and installing GATK')
		subprocess.call('wget https://github.com/broadinstitute/gatk/releases/download/4.0.11.0/gatk-4.0.11.0.zip', shell=True)
		subprocess.call('tar -xvzf gatk-4.0.11.0.zip', shell=True)
		# clean up
		subprocess.call('rm gatk-4.0.11.0.zip', shell=True)
	if not os.path.isfile('VarScan')
		print('Downloading VarScan')
		subprocess.call('wget https://sourceforge.net/projects/varscan/files/latest/download', shell=True)
		subprocess.call('mv download VarScan', shell=True)

	if platform.system() == 'Linux':
		subprocess.call('wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gff3ToGenePred', shell=True)
		subprocess.call('apt-get install bedtools', shell=True)
		subprocess.call('apt-get install samtools', shell=True)
		subprocess.call('apt-get install bwa', shell=True)
		subprocess.call('apt-get install mafft', shell=True)
	elif platform.system() == 'Darwin':
		subprocess.call('wget http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/gff3ToGenePred', shell = True)
		subprocess.call('brew install bedtools', shell=True)
		subprocess.call('brew install samtools', shell=True)
		subprocess.call('brew install bwa', shell=True)
		subprocess.call('brew install mafft', shell= True)
	subprocess.call('chmod +x gff3ToGenePred', shell =True)

