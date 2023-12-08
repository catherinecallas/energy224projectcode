import os
import numpy as np
import prep_Fine as prep_Fine
import prep_PSS as prep_PSS
import prep_Coarse as prep_Coarse
import run_ECL as run_ECL
import purge as purge
import time

def main():
    upsc_fact = [5,5,5]
    run_mode = 'local' # 'local' or 'sherlock'
    
    path_curr = os.getcwd()
    
    time_1 = time.time()
    # purge last run file
    purge.purge(path_curr)
    
    # prepare fine scale schedule reference depth
    prep_Fine.set_ref_depth(path_curr)
    
    # prepare PSS files
    prep_PSS.createPSSFile(path_curr, upsc_fact)
    print('PSS file done')
    #time_2 = time.time()
    #delta_time = time_2 - time_1
    #print(delta_time)
    
    # run PSS
    print('Running PSS...')
    run_ECL.run(path_curr, run_mode, 'PSS_Upsc')
    print('PSS run done')
    #time_3 = time.time()
    #delta_time = time_3 - time_2
    #print(delta_time)
    
    # prepare Coarse LGR files
    prep_Coarse.createCoarseFile(path_curr, upsc_fact)
    print('Coarse file done')
    #time_4 = time.time()
    #delta_time = time_4 - time_3
    #print(delta_time)
    
    # run Coarse LGR
    run_ECL.run(path_curr, run_mode, 'Coarse_LGR')
    print('Coarse run done')
    #time_5 = time.time()
    #delta_time = time_5 - time_4
    #print(delta_time)
    print('All done')
    
if __name__ == "__main__":
    main()