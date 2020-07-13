"""
run_sim.py
Runs simulation from Shah lab: https://github.com/shahpr/SMoPT
"""

import os
import shlex
import shutil
import subprocess as sp
import time

PATH_TO_SHAH = '../SMoPT'
def run_sim(output_path, linux, sim_time=1500, equil_time=1000, harr_add_time=1500, cyclo_add_time=1500,
             num_rib=200000, num_trna=3300000, num_genes=4839, cyclo_prob=0, cyclo_dissoc_rate=0, harr_rate=0,
             random_seed=1, input_genom_path=f'{PATH_TO_SHAH}/example/input/S.cer.flash-freeze.genom',
             input_trna_path=f'{PATH_TO_SHAH}/example/input/S.cer.tRNA', outputs=[1], logging=False):
    """
    Runs the Shah lab simulation with the given parameters
    Shah simulation parameters shown in parantheses next to variable descriptions
    
    :param output_path: Path to directory to save output files to
    :param linux: Boolean to indicate whether system is linux-based
    :param sim_time: (-Tt) Total simulation time in seconds
    :param equil_time: (-Tb) Burn-in/threshold time. Time spent by the cell to reach equilibrium.
    :param harr_add_time: (-Th) Time at which harringtonine is added to the cell
    :param cyclo_add_time: (-Tc) Time at which cycloheximide is added to the cell
    :param num_rib: (-R) Total number of ribosomes in the cell
    :param num_trna: (-t) Total number of tRNAs in the cell
    :param num_genes: (-N) Total number of genes
    :param cyclo_prob: (-x1) Probability of cycloheximide binding ribosomes
    :param cyclo_dissoc_rate: (-x2) Rate of cycloheximide dissociation from bound ribosomes
    :param harr_rate: (-y) Harringtonine rate for free ribosomes
    :param random_seed: (-s) Random number seed
    
    """
    sim_path = PATH_TO_SHAH + '/bin/SMoPT_v2'
    if linux:
        sim_path += '_linux'
    # Create output directory
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    
    # Create output prefix using -x1, -x2, -y values
    output_name = f'cp{float(cyclo_prob)}_cdr{float(cyclo_dissoc_rate)}_hr{float(harr_rate)}'
    full_path = output_path + output_name + '/'
    
    # Make subdirectory inside output directory for current values
    if not os.path.exists(full_path):
        os.mkdir(full_path)
    
    # Create command
    cmd = shlex.split(f'{sim_path} -Tt {sim_time} -Tb {equil_time} -Th {harr_add_time}\
                        -Tc {cyclo_add_time} -R {num_rib} -t {num_trna} -N {num_genes}\
                        -F {input_genom_path} -C {input_trna_path} -x1 {cyclo_prob}\
                        -x2 {cyclo_dissoc_rate} -y {harr_rate} -s {random_seed}\
                        -O {full_path}{output_name} '\
                        + ' '.join([f'-p{i}' for i in outputs]))
    
    # Run simulation, wait for it to finish
    proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
    if logging:
        print(' '.join(cmd))
    runtime = 0
    while proc.poll() == None:
        time.sleep(1)
        if logging:
            runtime += 1
            print('Running Simulation... Runtime: %d secs' % (runtime), end = '\r')
    if logging:
        print()
        print('Done')
    
    # Return list of saved output files 
    return [full_path + outfile for outfile in os.listdir(full_path)]
