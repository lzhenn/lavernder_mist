#!/usr/bin/env python3
'''
Module Painter


History:
Nov 06, 2022 --- render3d serial version completed!
Nov 06, 2022 --- FLEXPART output support for rendering

'''

import os, gc, subprocess
import numpy as np
import xarray as xr 
from multiprocessing import Pool
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import cmaps
import numba as nb

from . import io, utils, const, mtask_tool, cfgparser

print_prefix='lib.painter>>'

ncmap = ListedColormap(cmaps.GMT_globe(range(0,209,2))) 

def render3d(cfg, show_case):
    '''
    render 3d plot of the dumped partical file 
    '''
    utils.write_log(print_prefix+'render 3d plot...')
    # get topo file path
    topo_path=cfg['VISUAL']['topo_file_path']
    max_ny=float(cfg['VISUAL']['max_ny'])

    # vis domain bdy    
    latN, latS=\
        float(cfg['DOMAIN']['latN']), float(cfg['DOMAIN']['latS'])
    lonW, lonE=\
        float(cfg['DOMAIN']['lonW']), float(cfg['DOMAIN']['lonE'])
    zbot, ztop=\
        float(cfg['DOMAIN']['zbot']), float(cfg['DOMAIN']['ztop'])
    bnd_lim=[lonW,lonE,latS,latN,zbot,ztop]
    
    # cam pos
    cam_pos_lst=cfgparser.cfg_get_varlist(cfg, 'VISUAL', 'cam_pos')
    cam_pos_lst=[float(ele) for ele in cam_pos_lst]
    # get topo data
    topo_lat,topo_lon,topo_z = __get_topo_data(
        cfg,topo_path, latN, latS, lonW, lonE)
    # topo min lim
    topo_min=float(cfg['VISUAL']['topo_min'])
    topo_zoom=float(cfg['VISUAL']['topo_zoom'])
    color_zoom=float(cfg['VISUAL']['color_scale'])
    topo_z=topo_z/topo_zoom
    topo_z=topo_z.where(topo_z>0,topo_min)
    
    ny=topo_lat.shape[0]
    
    if ny>max_ny:
        zoom_r=int(np.floor(ny/max_ny))+1
        topo_lat=topo_lat[::zoom_r]
        topo_lon=topo_lon[::zoom_r]
        topo_z=topo_z[::zoom_r,::zoom_r] 
    topo_lon, topo_lat = np.meshgrid(topo_lon.values, topo_lat.values)
    
    dx = np.diff(topo_lon[1,1:3])
    dy = np.diff(topo_lat[1:3,1])
    # get partical dump data
    file_hdler=io.PtclHandler(cfg)
    
    # mtsks_lst[iproc(strt_subidx,end_subidx)]
    ntasks=int(cfg['VISUAL']['ntasks'])
    mtsks_lst=mtask_tool.alloc_mtasks_lst(
        ntasks, file_hdler.file_list)
    ntasks=len(mtsks_lst) 
    run_plot=cfg['VISUAL'].getboolean('run_plot')
    if run_plot: 
        # MULTIPROCESSING: start process pool
        process_pool = Pool(processes=ntasks)
        results=[]

        # open tasks ID 0 to ntasks-1
        for itsk in range(0, ntasks):  
            # start and end file index for each task
            result=process_pool.apply_async(
                __mtsk_render3d, 
                args=(
                    itsk, mtsks_lst, file_hdler, 
                    topo_lat, topo_lon, topo_z,
                    dx, dy, bnd_lim, 
                    show_case,
                    color_zoom, cam_pos_lst))
            results.append(result)
        print(results[0].get())
        process_pool.close()
        process_pool.join() 
    else:
        utils.write_log(print_prefix+'skip plot')
    
    form_anim=cfg['VISUAL'].getboolean('form_anim')
    if form_anim:
        _form_anim(cfg)

def  _form_anim(cfg):
        '''form animation from pngs'''
        utils.write_log('form animation from pngs')
        anim_name=cfg['VISUAL']['anim_outname']
        anim_fps=int(cfg['VISUAL']['anim_fps'])
        fig_path=os.path.join(const.CWD,'fig','*.png')
        anim_path=os.path.join(const.CWD,'fig',anim_name)
        cmd="ffmpeg -r %d -pattern_type glob -i '%s' -c:v libx264 -vf \"fps=8,format=yuv420p\" %s" % (
            anim_fps, fig_path, anim_path)
        subprocess.call(cmd, shell=True) 


def __mtsk_render3d(
        itsk, seg_list, fh, 
        topo_lat, topo_lon, topo_z, 
        dx, dy, bnd_lim,
        show_case, 
        color_zoom=1.0, cam_pos_lst=[60.0,0,-80.0,0]):
    '''
    multiprocessing render 3d plot of the dumped partical file
    '''
    
    cam_elev0, cam_elevr, cam_azim0, cam_azimr = cam_pos_lst 
    len_seg=len(seg_list[itsk])    
    for subid, ts in enumerate(seg_list[itsk]):
        # get global idx
        gidx=fh.file_list.index(ts)

        # init particle array in current frame
        prtarray=fh.load_frame(ts)
        
        if prtarray is None:
            continue
        
        utils.write_log(
            '%sTASK[%04d]: subtask (%04d/%04d) rendering %s with nptcl=%d' % (
                print_prefix, itsk, subid+1, len_seg, ts, prtarray.nptcls))
        
        plon,plat,pz=prtarray.xlon, prtarray.xlat, prtarray.ztra1
        # Create a mask that selects the points within the specified boundary
        lonW,lonE,latS,latN,zbot,ztop=bnd_lim
        mask = (plat >= latS) & (plat <= latN) & (
            plon >= lonW) & (plon <= lonE) & (pz >= zbot) & (pz <= ztop)

        # Apply the mask to the point positions to extract the points within the boundary
        plon = plon[mask]
        plat = plat[mask]
        pz = pz[mask] 

        fig = plt.figure(
            figsize=[const.FIG_WIDTH, const.FIG_HEIGHT],dpi=150)

        ax = fig.add_axes(const.FRM_MARGIN, projection='3d')
   
        res = len(topo_lon[0,:])
        # camera position in xyz
        xyz = np.array(__sph2cart(*__sphview(ax)), ndmin=3).T   
        # "distance" of bars from camera
        zo = np.multiply(
            [topo_lon, topo_lat, np.zeros_like(topo_z)], xyz).sum(0)  
        bars = np.empty(topo_lon.shape, dtype=object)

        #ncmap = ListedColormap(["red","blue","green"])
        #cnlevels = np.concatenate(
        #    (np.arange(-3150,0,50),np.arange(0,6150,150)))
        cnlevels = np.concatenate(
            (np.arange(-315,0,5),np.arange(0,615*color_zoom,15*color_zoom)))


        for i,(x,y,dz,o) in enumerate(__ravzip(
            topo_lon, topo_lat,topo_z, zo)):
            for nll in range(0,len(cnlevels),1):
                if nll==0 and dz<cnlevels[nll]:
                    color0 = ncmap([nll])
                    break
                elif (nll>0 and nll<(len(cnlevels)-1)) and (
                    dz>=cnlevels[nll-1] and dz<cnlevels[nll]):
                    color0 = ncmap([nll])
                    break
                elif nll==(len(cnlevels)-1) and dz>=cnlevels[nll]:
                    color0 = ncmap([nll+1])

            j, k = divmod(i, res)
            bars[j, k] = pl = ax.bar3d(x, y, 0, dx, dy, dz, color0)
            pl._sort_zpos = o
        if show_case=='':
            ax.scatter(
                plon, plat,pz, marker='.', color='white',
                zorder=999999, s=1, alpha=0.1)   
        elif show_case =='voc+no2':
               
            plon1,plat1,pz1=plon-0.08,plat-0.12,pz+300.0
               
            collision_mask=check_collision(plon,plat,pz,plon1,plat1,pz1,0.001,0.001,100)
            collision_mask2=check_collision(plon1,plat1,pz1,plon,plat,pz,0.001,0.001,100)
            utils.write_log(
                '%srendering %s@%s with ncollision=%d (%d) %6.2f%% (%6.2f%%)' % (
                print_prefix, show_case, ts, 
                np.sum(collision_mask),np.sum(collision_mask2), 
                np.sum(collision_mask)/len(collision_mask)*100,
                np.sum(collision_mask2)/len(collision_mask2)*100))
            plon2,plat2,pz2=plon1[collision_mask],plat1[collision_mask],pz1[collision_mask]
            
            # collision of ozone formation
            ax.scatter(
                plon2, plat2, pz2, marker='.', color='green',
                zorder=99, s=2, alpha=1.0)
            
            # voc source from vegetation
            plon2=plon1[np.logical_not(collision_mask)]
            plat2=plat1[np.logical_not(collision_mask)]
            pz2=pz1[np.logical_not(collision_mask)]
            ax.scatter(
                plon2, plat2,pz2, marker='.', color='blue',
                zorder=10, s=1, alpha=0.1)
            
            # no2 source from ship
            plon2=plon[np.logical_not(collision_mask2)]
            plat2=plat[np.logical_not(collision_mask2)]
            pz2=pz[np.logical_not(collision_mask2)]
            
            ax.scatter(
                plon2, plat2, pz2, marker='.', color='red',
                zorder=5, s=1, alpha=0.1)
            del plon1,plat1,pz1,plon2,plat2,pz2, collision_mask, collision_mask2
        ax.set_facecolor('k')
        ax.set_zscale('log')
        ax.set_xlim(lonW,lonE)
        ax.set_ylim(latS,latN)
        ax.set_zlim(zbot, ztop)
        ax.view_init(
            elev=cam_elev0+cam_elevr*gidx, 
            azim=cam_azim0+cam_azimr*gidx)
        ax.grid(False)
        ax.annotate(
            'AirTracers@%s' % fh.tfrms[gidx].strftime(const.Y_M_D_H_M),
            xy=(0.5, 0.95), xycoords='axes fraction', ha='center', 
            va='top', fontsize=const.MIDFONT, color='white')        
        
        plt.axis('off')
                
        figpath=os.path.join(
            const.CWD,'fig', 
            utils.parse_tswildcard(fh.tfrms[gidx], fh.fig_wildcard))
        plt.savefig(figpath, dpi=150)
        plt.close('all')
        utils.write_log(
            '%sTASK[%04d]: subtask (%04d/%04d) render3D done!' % (
                print_prefix, itsk, subid+1, len_seg))
        gc.collect()
        # break for test
        # break
    return 0


@nb.njit(parallel=True)
def check_collision(
        plat, plon, pz, plat1, plon1, pz1,
        latmin, lonmin, zmin):
    
    # initialize the collision array
    collision = np.zeros_like(plat1, dtype=np.bool_)

    # loop over particles in group2
    for i in nb.prange(plat1.shape[0]):
        # compute the distances between the i-th particle in group2
        # and all particles in group1
        dlat = np.abs(plat1[i] - plat)
        dlon = np.abs(plon1[i] - plon)
        dz = np.abs(pz1[i] - pz)

        # check if any particle in group1 collides with the i-th particle in group2
        if np.any((dlat < latmin) & (dlon < lonmin) & (dz < zmin)): 
            collision[i] = True
    return collision

def __sph2cart(r, theta, phi):
    '''spherical to cartesian transformation.'''
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z

def __sphview(ax):
    '''returns the camera position for 3D axes in spherical coordinates'''
    r = np.square(
        np.max([ax.get_xlim(), ax.get_ylim()], 1)).sum()
    theta, phi = np.radians((90-ax.elev, ax.azim))
    return r, theta, phi

def __ravzip(*itr):
    '''flatten and zip arrays'''
    return zip(*map(np.ravel, itr))

 
    
        
def __get_topo_data(cfg, topo_file, latN, latS, lonW, lonE):
    '''
    get topo data from topo file
    '''
    #topo_file=os.path.join(topo_path,'db/', fn)
    if not(os.path.exists(topo_file)):
        utils.write_log(
            print_prefix+topo_file+' not exist, skip...', lvl=30)
        return None
    else:
        utils.write_log(print_prefix+'load topo data from '+topo_file)
        topo=xr.open_dataset(topo_file)
    topo=topo.sel(y=slice(latS,latN), x=slice(lonW,lonE))
    topo_lat,topo_lon,topo_z=topo['y'], topo['x'], topo['z']
    
    return topo_lat,topo_lon,topo_z