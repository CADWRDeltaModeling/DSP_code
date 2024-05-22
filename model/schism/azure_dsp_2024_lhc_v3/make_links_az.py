import datetime as dt
import os

def make_links(start_date, end_date, src_dir, src_dir_narr, link_dir):
    # Creates sflux links from start to end _dates using the src and narr dirs
    # the link_dir is where the links will be created
    delt = dt.timedelta(days=1)
    current = start_date
    if (current >= end_date):
        print(f'ERROR: Start date {start_date.strftime("%b %d, %Y")} is after end date {end_date.strftime("%b %d, %Y")}')
    nfile = 0
    print(f'\t\t current: {current}')
    print(f'\t\t end_date: {end_date}')
    while (current <= end_date):
        # Air data
        # Ours
        src_str_air = os.path.join(src_dir,
                                   "baydelta_schism_air_%s%02d%02d.nc" % (current.year, current.month, current.day))
        # NARR
        src_str_rad = os.path.join(src_dir_narr, 
                                   "%4d_%02d/narr_rad.%4d_%02d_%02d.nc" % (current.year, current.month, current.year, current.month, current.day))
        src_str_prc = os.path.join(src_dir_narr, 
                                   "%4d_%02d/narr_prc.%4d_%02d_%02d.nc" % (current.year, current.month, current.year, current.month, current.day))
        nfile += 1
        link_str_air = os.path.join(link_dir, "sflux_air_1.%04d.nc" % (nfile))
        link_str_rad = os.path.join(link_dir, "sflux_rad_1.%04d.nc" % (nfile))
        link_str_prc = os.path.join(link_dir, "sflux_prc_1.%04d.nc" % (nfile))
        if not os.path.islink(link_str_air):
            os.symlink(src_str_air, link_str_air)
        else:
            os.remove(link_str_air)
            os.symlink(src_str_air, link_str_air)
        if not os.path.islink(link_str_rad):
            os.symlink(src_str_rad, link_str_rad)
        else:
            os.remove(link_str_rad)
            os.symlink(src_str_rad, link_str_rad)
        if not os.path.islink(link_str_prc):
            os.symlink(src_str_prc, link_str_prc)
        else:
            os.remove(link_str_prc)
            os.symlink(src_str_prc, link_str_prc)
        current += delt

os.chdir(os.path.dirname(os.path.abspath(__file__)))

print(f'\t Making links from {src_dir} and {src_dir_narr} to {link_dir} relative to {os.getcwd()}')

make_links(dt.datetime({start_date}), dt.datetime({end_date}), '{src_dir}', '{src_dir_narr}', '{link_dir}')

print('\t sflux links have been created. use ls -la sflux to see what links are pointing to.')