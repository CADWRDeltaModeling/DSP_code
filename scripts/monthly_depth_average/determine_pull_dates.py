import schimpy
import datetime as dt
import os
import re

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

class schism_out_per(object):

    def __init__(self, mod_dir, date_range):
        self.mod_dir = mod_dir
        self.date_range = date_range

        self.param_fn = os.path.join(self.mod_dir,'param.nml.clinic')
        with open(self.param_fn,"r") as fin:
            content=fin.read()
        self.param = schimpy.param.Params(content)

        self.get_ts_info()

        self.det_fns()
    
    def get_ts_info(self):

        self.start_time = self.param.get_run_start()
        self.interval = self.param['dt'] # get computational interval in seconds
        # self.hotstart_freq = self.param.get_hotstart_freq()
        self.nspool = self.param['nspool'] # output is done every nspool steps
        self.ihfskip = self.param['ihfskip'] # new output stack is opened every ihfskip steps
        self.nhot_write = self.param['nhot_write'] # hotstart file is written every nhot_write time steps
        self.start_date = dt.datetime(year=self.param['start_year'], month=self.param['start_month'], day=self.param['start_day'])

    def det_fns(self):

        print('hi')

        self.fns = [(self.date_range[0] - self.start_date).days, (self.date_range[1] - self.start_date).days] # file numbers to pull outputs


def main():
    print('main')



if __name__ == '__main__':

    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    mod_dir = '/scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/baseline_lhc_4'
    date_range = [dt.datetime(2008,8,1), dt.datetime(2008,8,31)]

    sop = schism_out_per(mod_dir, date_range)

    print(f'Collect output from file # {sop.fns[0]} to {sop.fns[1]}')