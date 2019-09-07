# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2019-01-22 11:24:06
# @Last modified by:   jsgounot
# @Last Modified time: 2019-02-11 13:47:18

import os, sys
from multiprocessing import Pool, cpu_count
from subprocess import check_output

# ------ Multiprocessing ------
 
class FunArgs() :
 
    def __init__(self, fun, * args, ** kwargs) :
        self.fun = fun
        self.args = args
        self.kwargs = kwargs
 
    def launch_fun(self) :
        return self.fun(* self.args, ** self.kwargs)
 
    def __str__(self) :
        return "%s : (%s %s)" %(str(self.fun), str(self.args), str(self.kwargs))

    def __repr__(self) :
        return str(self)

class Multiprocess() :
 
    @staticmethod
    def lambda_fun(farg) :
        return farg.launch_fun()
 
    def run(self, fargs, ncore=1) :
        print ("ncore : %i - available : %i" %(ncore, cpu_count()))
        if ncore == 1 : return [farg.launch_fun() for farg in fargs]
 
        pool = Pool(ncore)
        func = Multiprocess.lambda_fun
        fargs = ((farg, ) for farg in fargs)
 
        try :
            data = pool.starmap(func, fargs)
            pool.close()
            pool.join()
            return data
 
        except KeyboardInterrupt :
            print("Interrupt childs process")
            pool.terminate()
            sys.exit()
 
def mproc(fargs, ncore=1) :
    mp = Multiprocess()
    return mp.run(fargs, ncore=ncore)

