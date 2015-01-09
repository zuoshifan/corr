#!/bin/sh
# @ job_name         = mpidipole
# @ job_type         = bluegene
# @ comment          = "BGQ cubep3m dipole"
# @ error            = $(job_name).$(Host).$(jobid).err
# @ output           = $(job_name).$(Host).$(jobid).out
# @ bg_size          = 64
# @ wall_clock_limit = 1:30:00
# @ bg_connectivity  = Torus
# @ queue


export RUN_DIRECTORY='/home/p/pen/zhm/020_3fast'

export RUN_SUFFIX='.log'

runjob --np 64 --ranks-per-node=1 --cwd=$RUN_DIRECTORY : $RUN_DIRECTORY/lp >& logfile
