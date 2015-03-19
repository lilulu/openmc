#!/usr/bin/env python

import os
import sys
import h5py
import numpy as np

def main():

  # begin loop
  i = 6
  while i < 401:

    # create file name
    filename = 'statepoint.'+str(i)+'.h5'
    filenamebk = 'statepoint.'+str(i-1)+'.h5'

    # Open both HDF5 files
    print '../feed/'+filename
    print '../nofeed/'+filename
    fh_feed = h5py.File('../feed/'+filename,'r')
    fh_nofeed = h5py.File('../nofeed/'+filename,'r')
    fh_feed_bk = h5py.File('../feed/'+filenamebk,'r')
    fh_nofeed_bk = h5py.File('../nofeed/'+filenamebk,'r')


    # Extract numpy datasets
    cmfd_feed = fh_feed['cmfd/cmfd_src']
    openmc_feed = fh_feed['cmfd/openmc_src']
    openmc_nofeed = fh_nofeed['cmfd/openmc_src']
    cmfd_feed_bk = fh_feed_bk['cmfd/cmfd_src']
    openmc_feed_bk = fh_feed_bk['cmfd/openmc_src']
    openmc_nofeed_bk = fh_nofeed_bk['cmfd/openmc_src']

    # subtract off mean
    c_f = cmfd_feed[0,0,0:99,0]*i#- cmfd_feed_bk[0,0,0:99,0]*(i-1))*4e6
    o_f = (openmc_feed[0,0,0:99,0]*i - openmc_feed_bk[0,0,0:99,0]*(i-1))*4e6
    o_nf = (openmc_nofeed[0,0,0:99,0]*i - openmc_nofeed_bk[0,0,0:99,0]*(i-1))*4e6

    # renormalize
    c_f = c_f/np.sum(c_f)*c_f.shape[0]
    o_f = o_f/np.sum(o_f)*o_f.shape[0]
    o_nf = o_nf/np.sum(o_nf)*o_nf.shape[0]

    # Get file name
    filename = filename.split('.')
    name = filename[0]+filename[1]

    # Batch number
    batch = filename[1]

    # Write out gnuplot file
    write_gnuplot(name, c_f, o_f, o_nf, batch)

    i += 1

def write_gnuplot(name, c_f, o_f, o_nf, bat):

  # Header String for GNUPLOT
  headerstr = """#!/usr/bin/env gnuplot

set terminal pdf enhanced dashed
set output '{output}'
set key bottom center
set xlabel 'Slab Position [cm]' 
set ylabel 'Flux [-]'
set yrange [0.0:1.6]
set title  'Batch number {batch}'
set style line 1 lt 1 lc rgb "red" lw 3
set style line 2 lt 2 lc rgb "blue" lw 3
set style line 3 lt 1 lc rgb "orange" lw 3""".format(output=name+'.pdf',batch=bat)

  # Write out the plot string
  pltstr = \
"plot '-' with lines ls 3 title 'CMFD Fission Source - CMFD Case',\\\n\
     '-' with lines ls 1 title 'OpenMC Fission Source - Base Case',\\\n\
     '-' with lines ls 2 title 'OpenMC Fission Source - CMFD Case'"

  # Write out the data string
  i = 0
  datastr = ''
  while i < c_f.shape[0]:
    datastr = datastr + '{0}\n'.format(c_f[i])
    i += 1
  datastr = datastr + 'e\n'
  i = 0
  while i < o_nf.shape[0]:
    datastr = datastr + '{0}\n'.format(o_nf[i])
    i += 1
  datastr = datastr + 'e\n'
  i = 0
  while i < o_f.shape[0]:
    datastr = datastr + '{0}\n'.format(o_f[i])
    i += 1

  # Concatenate all
  outstr = headerstr + '\n' + pltstr + '\n' + datastr
  outstr = outstr.replace('nan','0.0')

  # Write File
  with open(name+".plot",'w') as f:
    f.write(outstr)

  # Run GNUPLOT
  os.system('gnuplot ' + name+".plot")
 
if __name__ == "__main__":
  main()
