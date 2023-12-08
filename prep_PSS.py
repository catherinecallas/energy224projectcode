import os
import subprocess
import prep_PSS_helper as helper
import shutil
import glob
import time
import datetime
import sys
import re
import pathlib
import os.path

def createPSSFile(main_path, upsc_fact):
    fine_well_file = os.path.join(main_path,'base.sched')
    Upsc_well_file = os.path.join(main_path,'PSS_Upsc','base.sched')
    
    f= open(Upsc_well_file,"w+")
    f.close()
    
    # number of injectors and their names 
    numOfInj, nameOfWells = helper.readWellsFromFine(fine_well_file)
    
    # given the above read well types and extract compdat from fine
    wellType, f_compdat = helper.readWellTypeFromFine(fine_well_file, nameOfWells)
    
    # read from fine: WELSPECS, WCONINJE
    f_welspecs = helper.readWellInfoFromFine(fine_well_file, 'WELSPECS')
    f_wconinje = helper.readWellInfoFromFine(fine_well_file, 'WCONINJE')
    
    # Create PPS WELSPECS
    pps_welspecs = helper.generateWelspecs(f_welspecs, upsc_fact)
    
    # Create PPS COMPDAT
    pps_compdat = helper.generateCompdat(f_compdat, upsc_fact)
    
    # Create PPS WCONINJE
    pps_wconinj = helper.generateWconinje(f_wconinje, nameOfWells)
    
    remain_base_sched = ['TSTEP','5*400 /',' ','END']
    
    output = open(Upsc_well_file, 'a')
        # WELSPECS
    output.write('WELSPECS'+'\n')
    output.write('\n'.join(map(str, pps_welspecs))+'\n')
    output.write('/'+'\n\n')
    # COMPDAT
    output.write('COMPDAT'+'\n')
    output.write('\n'.join(map(str, pps_compdat))+'\n')
    output.write('/'+'\n\n')
    # WCONINJE
    output.write('WCONINJE'+'\n')
    output.write('\n'.join(map(str, pps_wconinj))+'\n')
    output.write('/'+'\n\n')
    # END
    output.write('\n'.join(map(str, remain_base_sched))+'\n')
    output.close()
    
    
    
    