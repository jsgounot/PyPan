# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2019-01-22 10:42:57
# @Last modified by:   jsgounot
# @Last Modified time: 2019-02-01 16:33:04

import os
import json
import subprocess

dname = os.path.dirname
bname = os.path.basename

CONFIGFILE = None

class PyPanUtilsException(Exception) :
    pass

def set_config(fname) :
    global CONFIGFILE
    with open(fname) as f :
        CONFIGFILE = json.load(f)

def get_config() :
    global CONFIGFILE
    if CONFIGFILE is None :
        raise PyPanUtilsException("Configuration file not set")

    return CONFIGFILE

def get_field(* args) :
    cf = get_config()
    for arg in args :
        try : cf = cf[arg]
        except KeyError : raise PyPanUtilsException("Field not found in configuration file : %s" 
            %(" - ".join(str(field) for field in args)))

    return cf

def get_soft(name) :
    return get_field("SOFTPATH", name)

def check_install(soft) :
    rc = subprocess.call(['which', soft])
    if rc != 0 : raise PyPanUtilsException("Look like software '%s' is not installed" %(soft))

def check_soft(ignore_check_install=False) :
    names = ["snap", "augustus", "seg"]
    for name in names : 
        soft = get_field("SOFTPATH", name)
        if not ignore_check_install : check_install(soft)

def get_outdir(species, * args, mkdir=True) :
    outdir = get_field("PROJECT", species, "outdir")
    outdir = os.path.join(outdir, * args) if args else outdir
    if mkdir : os.makedirs(outdir, exist_ok=True)
    return outdir

def get_parameter(kind) :
    return get_field("PARAMETERS", kind)

def get_prot_db(name) :
    return get_field("BLASTDBPRED", name)

def get_dbs_used(name) :
    return get_field("BLASTDBUSED").get(name, [])