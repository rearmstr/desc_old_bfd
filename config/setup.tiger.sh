source /tigress/HSC/stack/hsc3-perseus-20171206/loadLSST.bash
module load rh/devtoolset/6
setup hscPipe 6.7
setup old_bfd
setup -jr ./
setup -jr obs_subaru
export TMV_DIR=/tigress/rea3/bfd_update/my_tmv

