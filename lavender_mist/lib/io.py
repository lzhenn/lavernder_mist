#!/usr/bin/env python3
"""specific module for IO"""
import datetime, struct, os
import pandas as pd
from scipy.io import FortranFile, FortranEOFError

from ..core import particles
from . import utils, const
import xarray as xr

print_prefix='lib.io>>'

class FileHandler():
    """ file handler class """
    def __init__(self, strt_time_str, end_time_str):
        self.strt_time=datetime.datetime.strptime(
            strt_time_str, const.YMDHM)
        
        self.end_time=datetime.datetime.strptime(
            end_time_str, const.YMDHM)
    
    def construct_file_list(self):
        '''
        construct file list according to file wildcard
        '''
        try:
            time_frms= pd.date_range(
                start=self.strt_time, end=self.end_time, freq=self.frq)
        except:
            utils.throw_error(
                print_prefix+'''cannot generate time frames,
                check init_t, end_t, and output_frq in config file''')
        self.tfrms=time_frms 
        self.file_list=[]
        for ts in time_frms:
            self.file_list.append(
                utils.parse_tswildcard(ts, self.file_wildcard))

class PtclHandler(FileHandler):

    '''
    Construct Particle File Handler for dumped particle data 
    in unformatted fortran binary file
    
    Methods
    -----------
    __init__:   initialize FLEXPHandler and loading data
    write_frame: write flexpart format partical dump file
    '''
    
    def __init__(self, cfg):
        '''
        Initialize OutHandler with config and loading data
        '''
        FileHandler.__init__(
            self, cfg['INPUT']['init_t'], cfg['INPUT']['end_t'])

        utils.write_log(print_prefix+'construct OutFileHandler')
        self.fmt=cfg['INPUT']['file_fmt']
        self.path=cfg['INPUT']['file_root']
        self.file_wildcard=cfg['INPUT']['file_wildcard']
        self.frq=cfg['INPUT']['feed_frq']
        
        self.fig_wildcard=cfg['VISUAL']['fig_wildcard']
        self.fig_fmt=cfg['VISUAL']['fig_fmt']
        
        #specs=cfgparser.cfg_get_varlist(cfg, 'EMISSION', 'specs')
        # May 25, 2023: only one spec is considered
        self.nspec=1

        self.construct_file_list()

    def load_frame(self, fn):
        '''
        load particle dump data according to single file 
        '''
        fn_full=os.path.join(self.path, fn)
        
        if not(os.path.exists(fn_full)):
            utils.write_log(
                print_prefix+fn_full+' not exist, skip...', lvl=30)
            return None 

        utils.write_log(print_prefix+'loading file: '+fn_full)
        
        if self.fmt=='flexpart':
            prtarray=self.__load_flexpart(fn_full)
        elif self.fmt=='nc':
            prtarray=self.__load_nc(fn_full)
        else:
            utils.throw_error(
                print_prefix+'unknown file format: '+self.fmt)
        utils.write_log(
            print_prefix+'%s loaded successfully!' % fn_full)
        return prtarray

    def __load_nc(self, fn): 
        '''
        load particle dump data according to single file 
        '''
       
        ds=xr.open_dataset(fn)
        nptcls=len(ds['parcel_id'])
        
        # construct partical array
        prtarray=particles.Particles(nptcls=nptcls, nspec=self.nspec)

        prtarray.xlon=ds['xlon'].values
        prtarray.xlat=ds['xlat'].values
        prtarray.ztra1=ds['xh'].values
        prtarray.itramem=ds['xtime'].values

        ds.close()

        return prtarray


    def __load_flexpart(self, fn):
        '''
        load flexpart format partical dump file
        '''

        rec_slab_temp=(
            'xlon', 'xlat', 'ztra1', 
            'itramem', 'topo', 'pvi', 'qvi', 'rhoi',
            'hmixi', 'tri', 'tti')
        dump_file = FortranFile(fn, 'r')        
        
        # header
        rec=dump_file.read_record('12c')
        (itime, nparts, iomode_xycoord)=struct.unpack('III', rec)
        
        # construct partical array
        prtarray=particles.Particles(nptcls=nparts, nspec=self.nspec)

        # id, rec_slab, and species
        len_slab_flag=str(
            (1+len(rec_slab_temp)+self.nspec)*4)+'c'
        decode_fmt='IfffIfffffff'+'f'*self.nspec
        for prtid in range(nparts):    
            try: 
                rec=dump_file.read_record(len_slab_flag)
            except FortranEOFError:
                break
            decode_rec=struct.unpack(decode_fmt, rec)
            for itmid, itm in enumerate(rec_slab_temp):
                getattr(prtarray, itm)[prtid]=decode_rec[itmid+1]
        
        return prtarray


