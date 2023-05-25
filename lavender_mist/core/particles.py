#/usr/bin/env python
"""
    Build Air Particals Object
    Note this class is an array of air particals

"""

import numpy as np
from ..lib import utils, const
print_prefix='core.particals>>'

class Particles:

    '''
    Construct air parcel array

    Attributes

    part_id:        int, partical id
    itramem:        int, memorized release times (s) of the particles
    xlon, xlat:     float, longitude/latitude of the particles
    xmeter, ymeter: float, x, y displacement (m) of the particles
    ztra1:          float, height of the particles
    topo:           float, topography
    pvi:            float, pressure
    qvi:            float, specific humidity
    rhoi:           float, density
    hmixi:          float, mixing ratio
    tri:            float, temperature
    tti:            float, potential temperature
    xmass:          float*nspec, mass of each specie

    Methods
    -----------


    '''
    
    def __init__(self, nptcls=1, nspec=1):
        """ construct air parcel obj """  
        # initiate partical templates
        __ARFP32=np.zeros(nptcls,dtype=np.float32)
        __ARIP32=np.zeros(nptcls,dtype=np.int32)
        
        self.nptcls,self.nspec=nptcls,nspec
        self.part_id=np.arange(nptcls,dtype=np.int32)
        self.itramem=__ARFP32.copy()
        self.xlon,self.xlat=__ARFP32.copy(),__ARFP32.copy()
        self.ix, self.iy, self.iz=\
            __ARIP32.copy(),__ARIP32.copy(),__ARIP32.copy()
        self.dx, self.dy, self.dz=\
            __ARFP32.copy(),__ARFP32.copy(),__ARFP32.copy()
        self.ztra1,self.topo=__ARFP32.copy(),__ARFP32.copy()
        self.pvi,self.qvi,self.rhoi,self.hmixi=\
            __ARFP32.copy(),__ARFP32.copy(),__ARFP32.copy(),__ARFP32.copy()
        self.tri,self.tti=__ARFP32.copy(),__ARFP32.copy()
        self.xmass=np.zeros((nptcls,nspec),dtype=np.float32)
    
        utils.write_log(
            print_prefix+'array with %d particals initiated!' % nptcls)

    def update(self, emis):
        """ update particle array by emission """ 
        self.itramem=np.linspace(
            -emis.emis_span,const.FP32_NEGZERO,self.nptcls,dtype=np.float32)  
        self.emis_span=emis.emis_span
        self.z0=emis.height
        self.ix, self.iy, self.iz=\
            self.ix+emis.ipos[0],\
            self.iy+emis.ipos[1],\
            self.iz+emis.ipos[2]
        self.ix0, self.iy0,self.iz0=emis.ipos[0],emis.ipos[1],emis.ipos[2]
        self.EOdx, self.EOdy=emis.EOdx, emis.EOdy
        # active particles number
        self.np=0


if __name__ == "__main__":
    pass
