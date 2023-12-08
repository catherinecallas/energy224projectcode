import os
import subprocess
import shutil


def run(path_curr, run_mode, file_dir):
    ecl_file = os.path.join(path_curr, file_dir)
    os.chdir(ecl_file)
    
    if file_dir == 'PSS_Upsc':
        eclipse_version = 'eclipse'
    if file_dir == 'Coarse_LGR':
        eclipse_version = 'e300'
    
    if run_mode == 'local':
        cmd = 'eclrun ' + eclipse_version + ' ' + 'CO2_ECLIPSE.DATA'
        os.system(cmd)
    if run_mode == 'sherlock':
        cmd = eclipse_version + '.exe' + ' ' + 'CO2_ECLIPSE'
        os.system(cmd)
        
    if file_dir == 'Coarse_LGR':
        source = 'CO2_ECLIPSE.RSM'
        destination = path_curr
        shutil.copy(source, destination)
    os.chdir(path_curr)

    