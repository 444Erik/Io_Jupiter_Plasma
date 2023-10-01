KPL/MK

   The names and contents of the kernels referenced by this
   meta-kernel are as follows:

   File name                   Contents
   --------------------------  -----------------------------
   naif0012.tls.pc             Generic LSK
   jup365.bsp                  Jupiter and Jupiter moons 1600 to 2200
   de441_part-2.bsp            Planetary positions 1965 to 17191 (not necessary?)
   pck00010.tpc                Rotation of Jupiter

   \begindata
   KERNELS_TO_LOAD = ( 'kernels/lsk/naif0012.tls.pc',
                       'kernels/spk/jup365.bsp',
                       'kernels/spk/de441_part-2.bsp',
                       'kernels/pck/pck00010.tpc' )
   \begintext