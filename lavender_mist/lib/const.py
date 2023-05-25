#/usr/bin/env python3
"""Commonly used Constants 


"""
import sys
import numpy as np
CWD=sys.path[0]


DEG2RAD=np.pi/180.0
RAD2DEG=180.0/np.pi

# Physical Constants
k=0.4 # Von Karman Constant
r_air=287.05   # individual gas constant for dry air, J/kg/K
cp=1004.67     # specific heat of dry air at constant pressure, J/kg/K
G=9.81 # m/s2
R_EARTH=6370000
T0=273.15 # K

# Diffusivity above PBL
diffz=0.1 # m^2s^-1
diffh=50.0 # m^2s^-1

SF_HEIGHT=100.0   # surface layer height, m
SCALE_VEL=5.0 # m/s for velocity scale
LATDIS=111000.0 # m for 1 deg latitude 

FP32_ISIM=np.float32(1e-30)
FP32_NEGZERO=np.float32(-0.0)

MAX_DT=3600 # s

# Calculate Mesh
L53_LOG={'layers':np.array([10*np.power(1.16, i)-10 for i in range(0,53)]),
         'c0':10.0, 'c1':1.16}

TEST_Z={'layers':np.arange(0, 10000, 200).astype(np.float32),
        'magic_idz':np.array([0]).astype(np.int32),
        'sep_z':np.array([200]).astype(np.float32)}

LAY2_Z={'layers':np.concatenate((
            np.arange(0, 3000, 200),
            np.arange(3000, 15000, 1000))).astype(np.float32),
        'magic_idz':np.array([0, 15]).astype(np.int32),
        'sep_z':np.array([200, 1000]).astype(np.float32)}


DENSE_Z={'layers':np.concatenate((
            np.arange(0, 300, 20), 
            np.arange(300, 1000, 50),
            np.arange(1000, 2000, 100),
            np.arange(2000,6000, 500),
            np.arange(6000, 10000, 1000), 
            np.arange(10000, 20000, 2000))).astype(np.float32),
        'magic_idz':np.array([0, 15]).astype(np.int32),
        'sep_z':np.array([20, 50]).astype(np.float32)}

# Render
## font
BIGFONT=22
MIDFONT=18
SMFONT=14
## fig
FIG_WIDTH=10
FIG_HEIGHT=10
FRM_MARGIN=[0.0, 0.0, 1.0, 1.0]

# Miscellaneous
YMDHM='%Y%m%d%H%M'
Y_M_D_H_M='%Y-%m-%d %H:%M'
YMDHMS='%Y%m%d%H%M%S'


