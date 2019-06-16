import os
import subprocess

folders = ['push-new-orig', 'push-pull', 'degree-sorting', 'hubdup', 'RADAR']

for folder in folders:
    origDir = os.getcwd()
    os.chdir(origDir + '/' + folder + '/Ligra/apps')
    subprocess.call('make clean; make -j 4', shell=True)
    os.chdir(origDir)

    if 'push-pull' in folder:
        print('~~~~ Compiled ' + folder + ' ~~~~')
        continue

    os.chdir(origDir + '/' + folder + '/GAP/')
    subprocess.call('make clean; make -j 4', shell=True)
    os.chdir(origDir)
    print('~~~~ Compiled ' + folder + ' ~~~~')
