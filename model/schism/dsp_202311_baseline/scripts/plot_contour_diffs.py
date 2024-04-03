#!/usr/bin/env python
# -*- coding: utf-8 -*-import pandas as pd
import sys

def diff_congruent_gr3(fname0,fname1,fout):
    print(f"Differenceing {fname0} and {fname1}")
    with open(fout,'w') as out:
        with open(fname0,'r') as first:
            with open(fname1,'r') as second:
                first_lines = list(first.readlines())
                nelem,nnode = map(int,first_lines[1].strip().split()[0:2])
                for i,(line0,line1) in enumerate(zip(first_lines,second.readlines())):
                    if i>1 and i<nnode+2:
                        try:
                            n,x,y,z = line0.strip().split()
                        except:
                            print(line1)
                            raise
                        n1,x1,y1,z1 = line1.strip().split()
                        diff = float(z1)-float(z)
                        outline = f"{n} {x} {y} {diff}\n"

                        out.write(outline)
                    else:
                        out.write(line0)
                    
                    
if __name__ == "__main__":
    # fname0,fname1,out = sys.argv[1:4]
    import os

    os.chdir('/scratch/tomkovic/DSP_code/model/schism/dsp_202311_baseline/scripts/')
    fname0 = 'baseline/14dsalt_200708.gr3' 
    fname1 = 'suisun/14dsalt_200708.gr3' 
    out = 'base_suis_diff_200708.gr3'
    diff_congruent_gr3(fname0,fname1,out)

    # conda env w suxarray conda activate them?

    # go to this folder then:
    # python ../calculate_depth_average.py --varname salinity --date_base 2006-11-14 --date_start 2007-07-01 --date_end 2007-09-01 --path_study ../../
    # rename those
    # python ../calculate_14d_salt.py --path_out2d ../../suisun/outputs/out2d_100.nc --path_salinity_nc depth_averaged_salinity.nc --path_prefix_gr3 salt14d
    # python plot_congruent_gr3.py baseline/14dsalt_200708.gr3 suisun/14dsalt_200708.gr3 base_suis_diff_200708.gr3