import schimpy
import subprocess
import os
import re
import math

import shutil

import schimpy.param

def tail(f, n, offset=0):
    if isinstance(n, int):
        n = str(n)
    if isinstance(offset, int):
        offset = str(offset)
    proc = subprocess.Popen(['tail', '-n', n + offset, f], stdout=subprocess.PIPE)
    lines = proc.stdout.readlines()
    return b''.join(lines[int(offset):])

class rst_fm_hotstart(object):

    def __init__(self, mod_dir, baro):
        self.mod_dir = mod_dir
        self.baro = baro
        self.param_fn = os.path.join(self.mod_dir,f'param.nml.{self.baro}')
        with open(self.param_fn,"r") as fin:
            content=fin.read()
        self.param = schimpy.param.Params(content)

        self.get_ts_info()
        self.get_mirror_crash()
        self.get_last_hotstart()
    
    def get_ts_info(self):

        self.start_time = self.param.get_run_start()
        self.interval = self.param['dt'] # get computational interval in seconds
        # self.hotstart_freq = self.param.get_hotstart_freq()
        self.nspool = self.param['nspool'] # output is done every nspool steps
        self.ihfskip = self.param['ihfskip'] # new output stack is opened every ihfskip steps
        self.nhot_write = self.param['nhot_write'] # hotstart file is written every nhot_write time steps

    def get_mirror_crash(self):

        mirrortail = tail(os.path.join(self.mod_dir, 'outputs/mirror.out'), 40)
        mmatches = re.findall("TIME STEP=\s+(\d+);", str(mirrortail))

        self.crash_timestep = int(max(mmatches))

    def get_last_hotstart(self):

        num_timesteps = math.floor(self.crash_timestep / self.nhot_write) * self.nhot_write

        # hotstart files are written as hotstart_[process_id]_[time_step].nc
        self.last_hotstart = str(int(num_timesteps))

    def combine_hotstart(self, iteration, machine='hpc5'):
        
        os.chdir('outputs')
        if machine.lower()=='hpc5':
            modld = f'module purge; module load intel/2024.0 hmpt/2.29 hdf5/1.14.3 netcdf-c/4.9.2 netcdf-fortran/4.6.1 schism/5.11.1; '
        else:
            modld = ''
        command = f'{modld}combine_hotstart7 --iteration {iteration}'

        ret = subprocess.run(command, capture_output=True, shell=True)
        os.chdir('../')
        
        print(ret.stdout.decode())

    def param_mod(self, iteration, param_in="param.nml.clinic"):

        # os.chdir(os.path.join(self.mod_dir))
        
        # replace hotstart variables ihot 
        with open(param_in, 'r') as file:
            param_text = file.read()
        param_out = re.sub('ihot = \d+','ihot = 2', param_text)

        param_fn_out = f'{param_in}.hot{iteration}'
        with open(param_fn_out, 'w') as file:
            file.write(param_out)

        # copy the hotstart file
        hotstart_fn = f'outputs/hotstart_it={iteration}.nc'
        hotstart_link_fn = 'hotstart.nc'
        shutil.copyfile(hotstart_fn, hotstart_fn.replace('outputs/',''))

        # copy the mirror.out file
        mirr_in = "outputs/mirror.out"
        if os.path.exists(mirr_in):
            shutil.copyfile(mirr_in, os.path.join(self.mod_dir, f'{mirr_in}.pre-{iteration}'))
        

def create_arg_parser():
    import argparse
    parser = argparse.ArgumentParser(
        description="Combine hotstarts from SCHISM run")
    parser.add_argument('--mod_dir', type=str,
                        help='directory for ', required=True)
    parser.add_argument('--baro', type=str,                        
                        help='"clinic" or "tropic" mode. Looks at relevant param.nml', default=None)
    parser.add_argument('--machine', type=str,                        
                        help='hpc5 or azure, decides which command to submit to cli', default=None)
    return parser


def main():
    # User inputs override the yaml file inputs.
    parser = create_arg_parser()
    args = parser.parse_args()

    rfh = rst_fm_hotstart(args.mod_dir, args.baro)
    rfh.combine_hotstart(rfh.last_hotstart, machine=args.machine)
    rfh.param_mod(rfh.last_hotstart)

    print('Combined Hotstart files')


# if __name__ == '__main__':

#     os.chdir(os.path.dirname(os.path.abspath(__file__)))
#     mod_dir = '.'
#     baro = 'clinic'

#     rfh = rst_fm_hotstart(mod_dir, baro)
#     print(rfh.last_hotstart)
#     rfh.combine_hotstart(rfh.last_hotstart, machine='hpc5')
#     rfh.param_mod(rfh.last_hotstart)

#     print('Combined Hotstart files')