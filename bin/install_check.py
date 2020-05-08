import sys
import os.path
gatk_version = 'gatk-4.0.11.0'

def check_picard(dir_path):
	dir_path = dir_path
	if os.path.isfile(dir_path + '/picard.jar'):
		return dir_path + '/picard.jar'
	else:
		print('Picard not found - lava is being executed from : ')
		subprocess.call('pwd', shell=True)
		print('LAVA checked for picard in the above folder and the main lava folder.')
		print('To fix this error download picard and unzip it into the main lava directory - for more in-depth help check out the readme.')
		sys.exit(1)

def check_gatk(dir_path):
	dir_path = dir_path
	if os.path.isfile(dir_path + '/' + gatk_version + '/gatk'):
		return dir_path + '/' + gatk_version + '/gatk'
	else:
		print('GATK not found - lava is being executed from : ')
		subprocess.call('pwd', shell=True)
		print('LAVA checked for GATK in the above folder and the main lava folder.')
		print('To fix this error download GATK and unzip it into the main lava directory - for more in-depth help check out the readme.')
		sys.exit(1)

def check_varscan(dir_path):
	dir_path = dir_path
	if os.path.isfile(dir_path + '/VarScan'):
		return dir_path + '/VarScan'
	else:
		print('VarScan not found alias lava.py="python `pwd`/lava.py" - lava is being executed from : ')
		subprocess.call('pwd', shell=True)
		print('LAVA checked for VarScan in the above folder and the main lava folder.')
		print('To fix this error download VarScan and unzip it into the main lava directory. NOTE: the jar file needs to be named VarScan - for more in-depth help check out the readme.')
		sys.exit(1)
