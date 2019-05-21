import os
import subprocess

folders = ['push-new-orig', 'push-pull', 'degree-sorting', 'hubdup', 'RADAR']

for folder in folders:
    origDir = os.getcwd()
    os.chdir(origDir + '/' + folder + '/apps')
    subprocess.call('make clean; make -j 28', shell=True)
    os.chdir(origDir)
    print('~~~~ Compiled ' + folder + ' ~~~~')
