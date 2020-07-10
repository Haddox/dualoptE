
from __future__ import print_function
import numpy as np
import subprocess
from distributed import Client, as_completed, wait
from dask_jobqueue import SLURMCluster

import uuid
import sys,os
import re
from scipy.optimize import minimize
import json
import time
import random
import subprocess
from scipy.stats import entropy
from copy import deepcopy

jobname="HH_run18D"

def run_job(cmd):
    ret = subprocess.call(cmd, shell=True,timeout=2400)
    ret = int(ret)
    if ret != 0:
        raise RuntimeError("command failed", cmd)
    return cmd, ret


# parse results files
def parse_result(file, scalescore, fieldscore, fieldcount):
    score, validate = 0,-1
    if os.path.exists(file):
        with open(file, 'r') as f:
            line = f.readline().strip()
            resfields = re.split("[ |,]+", line)
            if (len(resfields)>fieldcount):
                score = scalescore * float(resfields[fieldscore])
                validate = float(resfields[fieldcount])
            #else:
            #    raise IOError("error parsing %s!"%file)
    #else:
    #    raise IOError("%s doesn't exist!"%file)
    return score,validate



class ParamsOpt:
    def __init__(self):
        self.STRUCTCOUNTER = 1
        self.EVALONLY = 0
        self.RANDOMIZE = 0.0

        ## score weights
        self.wt_rr_int = 0.1;
        self.wt_rr_core = 0.1;
        self.wt_decoy = 1.0;
        self.wt_min_decoy = 1.0;
        self.wt_docking = 1.0;
        self.wt_grad = 1.0;
        self.wt_seqrecov = 1.0;
        self.wt_stab = 0.0000001;
        self.wt_distr = 10.0;
        self.wt_expsol = 1.0;

        ## tether weights
        self.wt_hbdel = 0.02;
        self.wt_lkoffdel = 0.01;

        ## special mins << FD UNUSED FOR NOW
        self.DISTR_MIN = 0.00;
        self.EXPSOL_MIN = 0.00;

        self.joblist = []

        ## PARAMETERS
        self.scorefxn="wts_beta_cart.wts"
        self.joblistfn = "alljobs"

        #jobrun parameters
        self.cluster_obj = SLURMCluster(cores=1, processes=1, memory="4GB",
                queue='long', name=jobname, extra=["--no-nanny", "--no-bokeh"],
                walltime="160:00:00")
        self.cluster_obj.scale(200)
        self.client = Client(self.cluster_obj, timeout=600)
        print(self.client.scheduler_info())

        ## Simplex params
        self.maxiter = 5000
        self.ftol = 1e-6
        self.sttime = time.time()
        self.initialized = False


    def initialize(self):
        params_init = self.read_init_flags("./flags_0")
        #print (params_init)

        values_full, scales, names, active = [],[],[],[]
        for elt in (params_init):
            names.append( elt[0] )
            values_full.append( elt[1] )
            scales.append( elt[2] )
            active.append( True )

        self.names = names
        self.values_full = np.array(values_full)
        active = np.array(active)
        self.x0 = self.values_full[active]
        self.scales = scales
        self.active = active

        self.read_joblist()

        print("#parameters to optimize:", np.sum(active))
        print("number of jobs:", len(self.joblist))

        self.initialized = True

    def read_init_flags(self, flagfn):
        readmode = 0

        params_init = []

        #self.initialize_elec(params_init)
        self.base_charges = {}

        with open(flagfn,'r') as f:
            lines = f.readlines()

        for line in lines:
            line = line.strip()
            if line == "" :
                continue

            if line.startswith("#") and not line.startswith("# rama") :
                continue

            if "LK_DGFREE" in line:
                readmode = 1
            elif "LK_LAMBDA" in line:
                readmode = 11
            elif "LJ_RADIUS" in line:
                readmode = 21
            elif "LK_RADIUS_OFFSET" in line:
                readmode = 31
            elif "LJ_WDEPTH" in line:
                readmode = 41
            elif "chemical:set_atomic_charge" in line:
                readmode = 2
            elif "chemical:set_patch_atomic_charge" in line:
                readmode = 12
            elif "-elec_sigmoidal_die" in line:
                readmode = 22
            elif "hb_don_strength" in line:
                readmode = 3
            elif "hb_acc_strength" in line:
                readmode = 4
            elif "-lk_ball" in line:
                readmode = 5
            elif "-set_weights" in line:
                readmode = 6
            elif line.startswith("# rama"):
                readmode = 7
            elif "hb_env_dep_exposed_scale" in line:
                readmode = 8
            elif "ldsrbb_" in line or "hbond_AHdist_scale" in line:
                readmode = 9
            else:
                sys.exit("unknown: "+line)

            if readmode == 1:
                fields = line.split(':')
                params_init.append( ["LKDGFREE:"+fields[-3], float(fields[-1]), 6.0] )
            elif readmode == 11:
                fields = line.split(':')
                params_init.append( ["LKLAMBDA:"+fields[-3], float(fields[-1]), 1.0] )
            elif readmode == 21:
                fields = line.split(':')
                params_init.append( ["LJRADIUS:"+fields[-3], float(fields[-1]), 0.1] )
            elif readmode == 31:
                fields = line.split(':')
                params_init.append( ["LKRADIUSOFFSET:"+fields[-3], float(fields[-1]), 1.0] )
            elif readmode == 41:
                fields = line.split(':')
                params_init.append( ["LJWDEPTH:"+fields[-3], float(fields[-1]), 0.1] )
            elif readmode == 2:
                fields = line.split(':')
                self.base_charges[ fields[-3]+":"+fields[-2] ] = float(fields[-1])
            elif readmode == 12:
                fields = line.split(':')
                self.base_charges[ fields[-4]+":"+fields[-3]+":"+fields[-2] ] = float(fields[-1])
            elif readmode == 22:
                fields = line.split(' ')
                params_init.append( ["sigmoidal_die:"+fields[0][1:], float(fields[1]), 0.5*float(fields[1])] )
            elif readmode == 3:
                fields = re.split(":|\s",line)
                params_init.append( ["hb_don_strength:"+fields[-2],float(fields[-1]), 0.5] )
            elif readmode == 4:
                fields = re.split(":|\s",line)
                params_init.append( ["hb_acc_strength:"+fields[-2],float(fields[-1]), 0.5] )
            elif readmode == 5:
                fields = re.split(" +",line)
                if (fields[0] == "-lk_ball_waters_sp2"):
                    params_init.append( ["lkball_sp2_ang", float(fields[1]), 24.0]  )
                elif (fields[0] == "-lk_ball_waters_sp3"):
                    params_init.append( ["lkball_sp3_ang", float(fields[1]), 24.0]  )
                elif (fields[0] == "-lk_ball_waters_donor"):
                    params_init.append( ["lkball_len", float(fields[1]), 1.0]  )
                elif (fields[0] == "-lk_ball_ramp_width_A2"):
                    params_init.append( ["lkball_ramp_width", float(fields[1]), 1.0]  )
                elif (fields[0] == "-lk_ball_overlap_width_A2"):
                    params_init.append( ["lkball_overlap_width", float(fields[1]), 1.0]  )
                elif (fields[0] == "-lk_ball_overlap_gap"):
                    params_init.append( ["lkball_overlap_gap", float(fields[1]), 0.2]  )
            elif readmode == 6:
                fields = re.split(" +",line)
                for i in range(1,len(fields),2):
                    params_init.append( ["weights:"+fields[i], float(fields[i+1]), 0.1] )
            elif readmode == 7:
                fields = re.split(" +",line)
                params_init.append( [fields[1], float(fields[2]), 0.5*float(fields[2])] )
            elif readmode == 8:
                fields = re.split(" +",line)
                params_init.append( ["hb_env_dep_exposed_scale",float(fields[1]), 0.2] )
            elif readmode == 9:
                fields = re.split(" +",line)
                params_init.append( [fields[0][1:], float(fields[1]), 0.5*float(fields[1])] )

        return params_init

    def initialize_elec(self, elec_scale_init):
        self.charge_group_net = {
            "NTERM_NH":1.0,
            "NTERM_NH_GLY":1.0,
            "NTERM_NH_PRO":1.0,
            "CTERM_CO":-1.0,
            "ALA_NH":0.0,
            "ALA_CO":0.0,
            "ALA_SC":0.0,
            "ARG_NH":0.0,
            "ARG_CO":0.0,
            "ARG_SC":1.0,
            "ASN_NH":0.0,
            "ASN_CO":0.0,
            "ASN_SC":0.0,
            "ASP_NH":0.0,
            "ASP_CO":0.0,
            "ASP_SC":-1.0,
            "CYD_NH":0.0,
            "CYD_CO":0.0,
            "CYD_SC":0.0,
            "CYS_NH":0.0,
            "CYS_CO":0.0,
            "CYS_SC":0.0,
            "GLN_NH":0.0,
            "GLN_CO":0.0,
            "GLN_SC":0.0,
            "GLU_NH":0.0,
            "GLU_CO":0.0,
            "GLU_SC":-1.0,
            "GLY_NH":0.0,
            "GLY_CO":0.0,
            "HIS_NH":0.0,
            "HIS_CO":0.0,
            "HIS_SC":0.0,
            "HIS_D_NH":0.0,
            "HIS_D_CO":0.0,
            "HIS_D_SC":0.0,
            "ILE_NH":0.0,
            "ILE_CO":0.0,
            "ILE_SC":0.0,
            "LEU_NH":0.0,
            "LEU_CO":0.0,
            "LEU_SC":0.0,
            "LYS_NH":0.0,
            "LYS_CO":0.0,
            "LYS_SC":1.0,
            "MET_NH":0.0,
            "MET_CO":0.0,
            "MET_SC":0.0,
            "PHE_NH":0.0,
            "PHE_CO":0.0,
            "PHE_SC":0.0,
            "PRO_NH":0.0,
            "PRO_CO":0.0,
            "SER_NH":0.0,
            "SER_CO":0.0,
            "SER_SC":0.0,
            "THR_NH":0.0,
            "THR_CO":0.0,
            "THR_SC":0.0,
            "TRP_NH":0.0,
            "TRP_CO":0.0,
            "TRP_SC":0.0,
            "TYR_NH":0.0,
            "TYR_CO":0.0,
            "TYR_SC":0.0,
            "VAL_NH":0.0,
            "VAL_CO":0.0,
            "VAL_SC":0.0,
        }
        self.charge_groups = {
            "NTERM_NH":["ALA:NtermProteinFull:N", "ALA:NtermProteinFull:CA", "ALA:NtermProteinFull:1H", "ALA:NtermProteinFull:2H", "ALA:NtermProteinFull:3H", "ALA:NtermProteinFull:HA"],
            "NTERM_NH_GLY":["GLY:NtermProteinFull:N", "GLY:NtermProteinFull:CA", "GLY:NtermProteinFull:1H", "GLY:NtermProteinFull:2H", "GLY:NtermProteinFull:3H", "GLY:NtermProteinFull:1HA", "GLY:NtermProteinFull:2HA"],
            "NTERM_NH_PRO":["PRO:NtermProteinFull:N", "PRO:NtermProteinFull:CA", "PRO:NtermProteinFull:1H", "PRO:NtermProteinFull:2H", "PRO:NtermProteinFull:HA", "PRO:NtermProteinFull:CB", "PRO:NtermProteinFull:CG", "PRO:NtermProteinFull:1HD",    "PRO:NtermProteinFull:2HD", "PRO:NtermProteinFull:1HG", "PRO:NtermProteinFull:2HG", "PRO:NtermProteinFull:1HB", "PRO:NtermProteinFull:2HB"],
            "CTERM_CO":["ALA:CtermProteinFull:C", "ALA:CtermProteinFull:OXT", "ALA:CtermProteinFull:O"],
            "ALA_NH":["ALA:N", "ALA:CA", "ALA:H", "ALA:HA"],
            "ALA_CO":["ALA:C", "ALA:O"],
            "ALA_SC":["ALA:CB", "ALA:1HB", "ALA:2HB", "ALA:3HB"],
            "ARG_NH":["ARG:N", "ARG:CA", "ARG:H", "ARG:HA"],
            "ARG_CO":["ARG:C", "ARG:O"],
            "ARG_SC":["ARG:CB", "ARG:CG", "ARG:CD", "ARG:NE", "ARG:CZ", "ARG:NH1", "ARG:NH2", "ARG:1HH1", "ARG:2HH1", "ARG:1HH2", "ARG:2HH2", "ARG:HE", "ARG:1HB", "ARG:2HB", "ARG:1HG", "ARG:2HG", "ARG:1HD", "ARG:2HD"],
            "ASN_NH":["ASN:N", "ASN:CA", "ASN:H", "ASN:HA"],
            "ASN_CO":["ASN:C", "ASN:O"],
            "ASN_SC":["ASN:CB", "ASN:CG", "ASN:OD1", "ASN:ND2", "ASN:1HD2", "ASN:2HD2", "ASN:1HB", "ASN:2HB"],
            "ASP_NH":["ASP:N", "ASP:CA", "ASP:H", "ASP:HA"],
            "ASP_CO":["ASP:C", "ASP:O"],
            "ASP_SC":["ASP:CB", "ASP:CG", "ASP:OD1", "ASP:OD2", "ASP:1HB", "ASP:2HB"],
            "CYD_NH":["CYD:N", "CYD:CA", "CYD:H", "CYD:HA"],
            "CYD_CO":["CYD:C", "CYD:O"],
            "CYD_SC":["CYD:CB", "CYD:SG", "CYD:1HB", "CYD:2HB"],
            "CYS_NH":["CYS:N", "CYS:CA", "CYS:H", "CYS:HA"],
            "CYS_CO":["CYS:C", "CYS:O"],
            "CYS_SC":["CYS:CB", "CYS:SG", "CYS:1HB", "CYS:2HB", "CYS:HG"],
            "GLN_NH":["GLN:N", "GLN:CA", "GLN:H", "GLN:HA"],
            "GLN_CO":["GLN:C", "GLN:O"],
            "GLN_SC":["GLN:CB", "GLN:CG", "GLN:CD", "GLN:OE1", "GLN:NE2", "GLN:1HE2", "GLN:2HE2", "GLN:1HB", "GLN:2HB", "GLN:1HG", "GLN:2HG"],
            "GLU_NH":["GLU:N", "GLU:CA", "GLU:H", "GLU:HA"],
            "GLU_CO":["GLU:C", "GLU:O"],
            "GLU_SC":["GLU:CB", "GLU:CG", "GLU:CD", "GLU:OE1", "GLU:OE2", "GLU:1HB", "GLU:2HB", "GLU:1HG", "GLU:2HG"],
            "GLY_NH":["GLY:N", "GLY:CA", "GLY:H", "GLY:1HA", "GLY:2HA"],
            "GLY_CO":["GLY:C", "GLY:O"],
            "HIS_NH":["HIS:N", "HIS:CA", "HIS:H", "HIS:HA"],
            "HIS_CO":["HIS:C", "HIS:O"],
            "HIS_SC":["HIS:CB", "HIS:CG", "HIS:ND1", "HIS:CD2", "HIS:CE1", "HIS:NE2", "HIS:HE2", "HIS:1HB", "HIS:2HB", "HIS:HE1", "HIS:HD2"],
            "HIS_D_NH":["HIS_D:N", "HIS_D:CA", "HIS_D:H", "HIS_D:HA"],
            "HIS_D_CO":["HIS_D:C", "HIS_D:O"],
            "HIS_D_SC":["HIS_D:CB", "HIS_D:CG", "HIS_D:ND1", "HIS_D:CD2", "HIS_D:CE1", "HIS_D:NE2", "HIS_D:HD1", "HIS_D:1HB", "HIS_D:2HB", "HIS_D:HE1", "HIS_D:HD2"],
            "ILE_NH":["ILE:N", "ILE:CA", "ILE:H", "ILE:HA"],
            "ILE_CO":["ILE:C", "ILE:O"],
            "ILE_SC":["ILE:CB", "ILE:CG1", "ILE:CG2", "ILE:CD1", "ILE:HB", "ILE:1HG2", "ILE:2HG2", "ILE:3HG2", "ILE:1HG1", "ILE:2HG1", "ILE:1HD1", "ILE:2HD1", "ILE:3HD1"],
            "LEU_NH":["LEU:N", "LEU:CA", "LEU:H", "LEU:HA"],
            "LEU_CO":["LEU:C", "LEU:O"],
            "LEU_SC":["LEU:CB", "LEU:CG", "LEU:CD1", "LEU:CD2", "LEU:1HB", "LEU:2HB", "LEU:HG", "LEU:1HD1", "LEU:2HD1", "LEU:3HD1", "LEU:1HD2", "LEU:2HD2", "LEU:3HD2"],
            "LYS_NH":["LYS:N", "LYS:CA", "LYS:H", "LYS:HA"],
            "LYS_CO":["LYS:C", "LYS:O"],
            "LYS_SC":["LYS:CB", "LYS:CG", "LYS:CD", "LYS:CE", "LYS:NZ", "LYS:1HZ", "LYS:2HZ", "LYS:3HZ", "LYS:1HB", "LYS:2HB", "LYS:1HG", "LYS:2HG", "LYS:1HD", "LYS:2HD", "LYS:1HE", "LYS:2HE"],
            "MET_NH":["MET:N", "MET:CA", "MET:H", "MET:HA"],
            "MET_CO":["MET:C", "MET:O"],
            "MET_SC":["MET:CB", "MET:CG", "MET:SD", "MET:CE", "MET:1HB", "MET:2HB", "MET:1HG", "MET:2HG", "MET:1HE", "MET:2HE", "MET:3HE"],
            "PHE_NH":["PHE:N", "PHE:CA", "PHE:H", "PHE:HA"],
            "PHE_CO":["PHE:C", "PHE:O"],
            "PHE_SC":["PHE:CB", "PHE:CG", "PHE:CD1", "PHE:CD2", "PHE:CE1", "PHE:CE2", "PHE:CZ", "PHE:HD1", "PHE:HE1", "PHE:HZ", "PHE:HE2", "PHE:HD2", "PHE:1HB", "PHE:2HB"],
            "PRO_NH":["PRO:N", "PRO:CA", "PRO:HA", "PRO:CB", "PRO:CG", "PRO:1HD", "PRO:2HD", "PRO:1HG", "PRO:2HG", "PRO:1HB", "PRO:2HB"],
            "PRO_CO":["PRO:C", "PRO:O"],
            "SER_NH":["SER:N", "SER:CA", "SER:H", "SER:HA"],
            "SER_CO":["SER:C", "SER:O"],
            "SER_SC":["SER:CB", "SER:OG", "SER:HG", "SER:1HB", "SER:2HB"],
            "THR_NH":["THR:N", "THR:CA", "THR:H", "THR:HA"],
            "THR_CO":["THR:C", "THR:O"],
            "THR_SC":["THR:CB", "THR:OG1", "THR:CG2", "THR:HG1", "THR:HB", "THR:1HG2", "THR:2HG2", "THR:3HG2"],
            "TRP_NH":["TRP:N", "TRP:CA", "TRP:H", "TRP:HA"],
            "TRP_CO":["TRP:C", "TRP:O"],
            "TRP_SC":["TRP:CB", "TRP:CG", "TRP:CD1", "TRP:CD2", "TRP:NE1", "TRP:CE2", "TRP:CE3", "TRP:CZ2", "TRP:CZ3", "TRP:CH2", "TRP:HE1", "TRP:HD1", "TRP:HZ2", "TRP:HH2", "TRP:HZ3", "TRP:HE3", "TRP:1HB", "TRP:2HB"],
            "TYR_NH":["TYR:N", "TYR:CA", "TYR:H", "TYR:HA"],
            "TYR_CO":["TYR:C", "TYR:O"],
            "TYR_SC":["TYR:CB", "TYR:CG", "TYR:CD1", "TYR:CD2", "TYR:CE1", "TYR:CE2", "TYR:CZ", "TYR:OH", "TYR:HH", "TYR:HD1", "TYR:HE1", "TYR:HE2", "TYR:HD2", "TYR:1HB", "TYR:2HB"],
            "VAL_NH":["VAL:N", "VAL:CA", "VAL:H", "VAL:HA"],
            "VAL_CO":["VAL:C", "VAL:O"],
            "VAL_SC":["VAL:CB", "VAL:CG1", "VAL:CG2", "VAL:HB", "VAL:1HG1", "VAL:2HG1", "VAL:3HG1", "VAL:1HG2", "VAL:2HG2", "VAL:3HG2"]
        }
        self.coupled_charge_groups = {
            "allaa_CO":["ALA_CO", "ASN_CO", "ASP_CO", "ARG_CO", "GLY_CO", "CYS_CO", "CYD_CO", "ILE_CO", "GLN_CO", "GLU_CO", "HIS_CO", "HIS_D_CO", "LEU_CO", "MET_CO", "PHE_CO", "PRO_CO", "TYR_CO", "TRP_CO", "VAL_CO", "LYS_CO", "SER_CO", "THR_CO"],   ## 1
            "allaa_NH":[ "ALA_NH", "ASN_NH", "ASP_NH", "ARG_NH", "GLY_NH", "CYS_NH", "CYD_NH", "ILE_NH", "GLN_NH", "GLU_NH", "HIS_NH", "HIS_D_NH", "LEU_NH", "MET_NH", "PHE_NH", "PRO_NH", "TYR_NH", "TRP_NH", "VAL_NH", "LYS_NH", "SER_NH", "THR_NH"],   ## 2
            "nonpolar_SC":["ALA_SC", "CYS_SC", "CYD_SC", "ILE_SC", "LEU_SC", "MET_SC", "VAL_SC"],   ## 3
            "ASN_SC":["ASN_SC"],  ## 4
            "ASP_SC":["ASP_SC"],  ## 5
            "ARG_SC":["ARG_SC"],  ## 6
            "GLN_SC":["GLN_SC"],  ## 7
            "GLU_SC":["GLU_SC"],  ## 8
            "HIS_SC":["HIS_SC", "HIS_D_SC"],  ## 9
            "LYS_SC":["LYS_SC"],  ## 10
            "SER_SC":["SER_SC"],  ## 11
            "THR_SC":["THR_SC"],  ## 12
            "PHE_SC":["PHE_SC"],  ## 13
            "TRP_SC":["TRP_SC"],  ## 14
            "TYR_SC":["TYR_SC"],  ## 15
            "nterm_NH":["NTERM_NH", "NTERM_NH_GLY", "NTERM_NH_PRO"], ## THE ORDER HERE IS IMPORTANT... GENERAL CASE MUST COME 1ST!!!!
            "cterm_CO":["CTERM_CO"] ## 17
        }

        for i in self.coupled_charge_groups.keys():
            elec_scale_init.append(["elec:"+i,1.0,0.6])

    def read_joblist(self):
        with open(self.joblistfn, 'r') as infn:
            lines = infn.readlines()
        for line in lines:
            self.joblist.append( line.strip() )

    def setup_params_flag(self, outfn):
        flags_content = ""

        lkb_sp2, lkb_sp3, lkb_len = 0,0,0
        for i,name in enumerate(self.names):
            fields = re.split(':', name)
            param = self.values_full[i]

            if (fields[0] == "elec"):
                for group in self.coupled_charge_groups[fields[1]]:
                    netcharge = self.charge_group_net[group]
                    natoms = len(self.charge_groups[group])
                    for atom in self.charge_groups[group]:
                        newcharge = param * (
                            self.base_charges[atom] - (netcharge/natoms) ) + (netcharge/natoms)
                        atomfields = re.split(':', atom)
                        if (len(atomfields) == 2):
                            flags_content+="""-chemical:set_atomic_charge fa_standard:%s:%.5f\n""" % (atom, newcharge)
                        else:
                            flags_content+="""-chemical:set_patch_atomic_charge fa_standard:%s:%.5f\n""" % (atom, newcharge)
            elif (fields[0] == "LKDGFREE"):
                flags_content+="""-chemical:set_atom_properties fa_standard:%s:LK_DGFREE:%.5f\n""" % (fields[1], param)
            elif (fields[0] == "LKLAMBDA"):
                flags_content+="""-chemical:set_atom_properties fa_standard:%s:LK_LAMBDA:%.5f\n""" % (fields[1], param)
            elif (fields[0] == "LJRADIUS"):
                flags_content+="""-chemical:set_atom_properties fa_standard:%s:LJ_RADIUS:%.5f\n""" % (fields[1], param)
            elif (fields[0] == "LKRADIUSOFFSET"):
                flags_content+="""-chemical:set_atom_properties fa_standard:%s:LK_RADIUS_OFFSET:%.5f\n""" % (fields[1], param)
            elif (fields[0] == "LJWDEPTH"):
                flags_content+="""-chemical:set_atom_properties fa_standard:%s:LJ_WDEPTH:%.5f\n""" % (fields[1], max(0.001, param))
            elif (fields[0] == "hb_don_strength"):
                flags_content+="""-score::hb_don_strength %s:%.3f\n""" % (fields[1], param)
            elif (fields[0] == "hb_acc_strength"):
                flags_content+="""-score::hb_acc_strength %s:%.3f\n""" % (fields[1], param)
            elif (fields[0] == "weights"):
                flags_content+="""-set_weights %s %.4f\n""" % (fields[1], param)
            elif (fields[0] == "sigmoidal_die"):
                flags_content+="""-%s %.4f\n""" % (fields[1], param)
            elif (fields[0] == "lkball_sp2_ang"):
                lkb_sp2 = param
            elif (fields[0] == "lkball_sp3_ang"):
                lkb_sp3 = param
            elif (fields[0] == "lkball_len"):
                lkb_len = param
            elif (fields[0] == "lkball_ramp_width"):
                flags_content+="""-lk_ball_ramp_width_A2 %.4f\n""" % param
            elif (fields[0] == "lkball_overlap_width"):
                flags_content+="""-lk_ball_overlap_width_A2 %.4f\n""" % param
            elif (fields[0] == "lkball_overlap_gap"):
                flags_content+="""-lk_ball_overlap_gap %.4f\n""" % param
            elif (fields[0] == "rama"):
                flags_content+="""# rama:%s %.4f\n""" % (fields[1], param)
            elif (fields[0] == "hb_env_dep_exposed_scale"):
                flags_content+="""-hb_env_dep_exposed_scale %.4f\n""" % (param)
            elif (fields[0] == "ldsrbb_low_scale" or fields[0] == "ldsrbb_high_scale"):
                flags_content+="""-%s %.3f\n""" % (fields[0], param)
            elif (fields[0] == "ldsrbb_minlength" or fields[0] == "ldsrbb_maxlength"):
                flags_content+="""-%s %i\n""" % (fields[0], param)
            elif (fields[0] == "hbond_AHdist_scale"):
                flags_content+="""-%s %.3f\n""" % (fields[0], param)

        # to do: optional?
        if (lkb_len != 0):
            flags_content+="""-lk_ball_waters_sp2 %.4f %.1f %.1f %.4f %.1f %.1f\n""" % (lkb_len, lkb_sp2, 0.0, lkb_len, lkb_sp2, 180.0)
            flags_content+="""-lk_ball_waters_sp3 %.4f %.1f %.1f %.4f %.1f %.1f\n""" % (lkb_len, lkb_sp3, 120.0, lkb_len, lkb_sp3, 240.0)
            flags_content+="""-lk_ball_waters_ring %.4f %.1f %.1f\n""" % (lkb_len, 180.0, 0.0)
            flags_content+="""-lk_ball_waters_donor %.4f\n""" % lkb_len

        with open(outfn, 'w') as flagsout:
            flagsout.write(flags_content)

    def run_jobs(self, params, cmdlist):
        output_dir = "./opt_%d"%self.STRUCTCOUNTER
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        flagsout_fn = "./flags_%d"%self.STRUCTCOUNTER

        self.setup_params_flag(flagsout_fn)
        results = []


        for cmd in cmdlist:
            results.append( self.client.submit(run_job, cmd) )
        #print(self.client.gather(results))
        wait(results)

        #cleanup_cmd = "./run_cleanup.sh ./opt_%d &> /dev/null"%self.STRUCTCOUNTER
        #os.system(cleanup_cmd)
        time.sleep(60)
        cleanup_cmd = "./run_cleanup.sh ./opt_%d"%self.STRUCTCOUNTER
        p = subprocess.Popen(
            cleanup_cmd, shell=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        out, err = p.communicate()
        #print("out:", out)
        #print("err:", err)
        time.sleep(30)

    def validate_results(self,ct_rr_int,ct_rr_core,ct_decoy,ct_min_decoy,ct_docking,ct_grad,ct_seqrecov,ct_distr,ct_expsol):
        if (ct_rr_int < 3380): 
            return False
        if (ct_rr_core < 18850): 
            return False
        if (ct_decoy < 100112): 
            return False
        if (ct_min_decoy < 5340): 
            return False
        if (ct_docking < 8997): 
            return False
        if (ct_seqrecov < 15796): 
            return False
        if (ct_grad < 47): 
            return False
        if (ct_distr < 90):
            return False
        if (ct_expsol < 2000): 
            return False

        return True


    def compute_total_score(self):
        sc_rr_int,sc_rr_core,sc_decoy,sc_min_decoy,sc_docking,sc_grad,sc_seqrecov,sc_distr,sc_expsol,sc_seqstab = 0,0,0,0,0,0,0,0,0,0
        ct_rr_int,ct_rr_core,ct_decoy,ct_min_decoy,ct_docking,ct_grad,ct_seqrecov,ct_distr,ct_expsol = 0,0,0,0,0,0,0,0,0

        # hbond offsets
        hbdon, hbdon2, hbacc, hbacc2 = 0.0,0.0,0.0,0.0
        Ndon, Nacc = 0,0
        for i,name in enumerate(self.names):
            fields = re.split(':', name)
            param = self.values_full[i]
            if (fields[0] == "hb_don_strength"):
                hbdon += param
                hbdon2 += param*param
                Ndon += 1
            elif (fields[0] == "hb_acc_strength"):
                hbacc += param
                hbacc2 += param*param
                Nacc += 1
        if (Ndon == 0) or (Nacc) == 0:
            sc_hbdel = 0.0
        else:
            sc_hbdel = (
                np.sqrt( Ndon*hbdon2 - hbdon*hbdon ) / Ndon
                + np.sqrt( Nacc*hbacc2 - hbacc*hbacc ) / Nacc )

        # lk offsets
        lkro, Nro = 0.0, 0
        for i,name in enumerate(self.names):
            fields = re.split(':', name)
            param = self.values_full[i]
            if (fields[0] == "LKRADIUSOFFSET"):
                lkro += param*param
                Nro += 1
        if Nro == 0:
            sc_lkoffdel = 0
        else:
            sc_lkoffdel = np.sqrt( lkro / Nro )

        sc_rr_int,ct_rr_int = parse_result(
            "opt_%d/rr_interface_result"%self.STRUCTCOUNTER, -1, 0, 2)
        sc_rr_core,ct_rr_core = parse_result(
            "opt_%d/rr_core_result"%self.STRUCTCOUNTER, -1, 0, 2)
        sc_decoy,ct_decoy = parse_result(
            "opt_%d/score_decoy_result"%self.STRUCTCOUNTER, -1, 0, 2)
        sc_min_decoy,ct_min_decoy = parse_result(
            "opt_%d/score_min_result"%self.STRUCTCOUNTER, -1, 0, 2)
        sc_docking,ct_docking = parse_result(
            "opt_%d/score_docking_result"%self.STRUCTCOUNTER, -1, 0, 2)
        sc_grad,ct_grad = parse_result(
            "opt_%d/xtal_grad_result"%self.STRUCTCOUNTER, 1, 0, 1)

        sc_seqrecov,ct_seqrecov = parse_result(
            "opt_%d/seqrecov_result"%self.STRUCTCOUNTER, 1, 2, 0)
        sc_seqstab,_ = parse_result(
            "opt_%d/seqrecov_result"%self.STRUCTCOUNTER, 1, 3, 0)

        sc_distr,ct_distr = parse_result(
            "opt_%d/distr_result"%self.STRUCTCOUNTER, 1, 0, 1)
        sc_expsol,ct_expsol = parse_result(
            "opt_%d/expsol_result"%self.STRUCTCOUNTER, 1, 0, 2)

        ### SPECIAL CASE FOR DISTR!
        if (sc_distr < self.DISTR_MIN):
            sc_distr = self.DISTR_MIN

        ### SPECIAL CASE FOR expsol!
        if (sc_expsol < self.EXPSOL_MIN):
            sc_expsol = self.EXPSOL_MIN

        elapsed_time = time.time() - self.sttime

        if self.validate_results(ct_rr_int,ct_rr_core,ct_decoy,ct_min_decoy,ct_docking,ct_grad,ct_seqrecov,ct_distr,ct_expsol)==0:
            print( "[%4d : %6d s]  Error! rerunning!  %d, %d, %d, %d, %d, %d, %d, %d, %d"%(self.STRUCTCOUNTER,
                elapsed_time, ct_rr_int,ct_rr_core,ct_decoy,ct_min_decoy,ct_docking,ct_grad,ct_seqrecov,ct_distr,ct_expsol) )
            #return 0
            sys.exit()
        else:
            score_tot = (
                self.wt_rr_int * sc_rr_int +
                self.wt_rr_core * sc_rr_core +
                self.wt_decoy * sc_decoy +
                self.wt_min_decoy * sc_min_decoy +
                self.wt_docking * sc_docking +
                self.wt_grad * sc_grad +
                self.wt_seqrecov * sc_seqrecov +
                self.wt_stab * sc_seqstab +
                self.wt_distr * sc_distr +
                self.wt_expsol * sc_expsol +
                self.wt_hbdel * sc_hbdel +
                self.wt_lkoffdel * sc_lkoffdel
            )
            if self.STRUCTCOUNTER==1:
                print( "                   %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s"%(
                    "TOTAL", "RR_INT", "RR_CORE", "DECOY", "MINDECOY", "DOCK", "XTALGRAD", "SEQREC", "STAB", "DISTR", "EXPSOL","HBDEL","LKDEL") )
            print( "[ %4d : %6d s] %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f"%(
                    self.STRUCTCOUNTER, elapsed_time,
                    score_tot,sc_rr_int,sc_rr_core,sc_decoy,sc_min_decoy,sc_docking,sc_grad,sc_seqrecov,sc_seqstab,sc_distr,sc_expsol,sc_hbdel,sc_lkoffdel ))

            self.STRUCTCOUNTER += 1
            if self.EVALONLY != 0:
                sys.exit()

            return score_tot

    def construct_init_simplex(self, x0, scales):
        assert(len(x0)==len(scales))
        simplex = np.tile(x0,(len(x0)+1,1))
        for i in range(len(x0)):
            simplex[i+1][i] += scales[i]
        return simplex

    def func(self, x):
        self.values_full[self.active] = x

        score = 0
        firstoutfile = "opt_%d/rr_interface_result"%self.STRUCTCOUNTER

        if os.path.exists(firstoutfile):
            score = self.compute_total_score()
            if score != 0:
                return score
        while score == 0:
            cmdlist = []
            for cmd in self.joblist:
                cmd = cmd+ " {} {}".format(self.STRUCTCOUNTER, self.scorefxn)
                cmdlist.append(cmd)
            self.run_jobs(x, cmdlist)
            score = self.compute_total_score()

        return score

    def run_opt(self):
        if not self.initialized:
            self.initialize()
        init_simplex = self.construct_init_simplex(self.x0, self.scales)
        res = minimize(self.func, self.x0, method='Nelder-Mead',
                       options={'initial_simplex':init_simplex,
                       'maxiter': self.maxiter,
                       'fatol': self.ftol })
        print( res )


if __name__ == "__main__":
    opt = ParamsOpt()
    opt.run_opt()
