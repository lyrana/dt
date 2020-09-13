#! /bin/bash
make -f Makefile.test
make -f Makefile.df
make -f Makefile.mea
make -f Makefile.types
make -f Makefile.sap
#make -f Makefile.upbf
make -f Makefile_x.upbf
make -f Makefile_xy.upbf
make -f Makefile_xyz.upbf
make -f Makefile_2D_e.upbf
make -f Makefile_spherical_r.upbf
#make -f Makefile.part
make -f Makefile_x.part
make -f Makefile_r.part
make -f Makefile_xy.part
make -f Makefile_xyz.part
