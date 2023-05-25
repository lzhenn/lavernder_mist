#/usr/bin/env python3
"""Commonly used utilities for parallization

    Function    
    ---------------
   
    alloc_mtasks_lst():

"""
from multiprocessing import cpu_count
from . import utils

def alloc_mtasks_lst(ntasks, task_list):
    '''
    allocated tasks to multiple processes
    '''
    len_lst=len(task_list)
    
    # determine number of tasks
    if ntasks ==0:
        ntasks=cpu_count()
    if ntasks > len_lst:
        ntasks=len_lst
        utils.write_log('ntasks reduced to  %4d according list length.' % (ntasks))
    if ntasks > cpu_count():
        ntasks=cpu_count()
        utils.write_log('ntasks reduced to  %4d according cpu_count' % (ntasks))    
    utils.write_log('ntasks set to: %4d' % (ntasks))

    avglen, res=divmod(len(task_list), ntasks)
    seg_lst=[]
    for i in range(0, ntasks):
        if i < res:
            seg_lst.append(task_list[i*avglen+i: (i+1)*avglen+i+1])
        else:
            seg_lst.append(task_list[i*avglen+res: (i+1)*avglen+res])
    return seg_lst