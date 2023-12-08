import os
import glob

def purge(path_curr):
    path_pss_model = os.path.join(path_curr, 'PSS_Upsc') # path to the psuedo steady state upscaling model
    path_coarse_model = os.path.join(path_curr, 'Coarse_LGR')
    
    for file in glob.glob(os.path.join(path_pss_model,'base.sched')):
        os.remove(file)
    for file in glob.glob(os.path.join(path_coarse_model,'*.ECLIN')):
        os.remove(file)
    for file in glob.glob(os.path.join(path_coarse_model,'base_C.sched')):
        os.remove(file)
    for file in glob.glob(os.path.join(path_curr,'*.RSM')):
        os.remove(file)
    for file in glob.glob(os.path.join(path_curr,'Coarse_well_blocks.txt')):
        os.remove(file)
    for file in glob.glob(os.path.join(path_curr,'LGR_names.txt')):
        os.remove(file)