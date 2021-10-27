#!/usr/bin/env python
#### <license> ####
# Copyright (c) 2016-2017, Lawrence Livermore National Security,
# LLC. Produced at the Lawrence Livermore National Laboratory. Written
# by Robert Blake <blake14@llnl.gov>.
# 
# LLNL-CODE-720003.
# All rights reserved.
#
# This file is part of MELODEE. For details, see
# http://github.com/llnl/melodee.
#
# Licensed under the Apache License, Version 2.0 (the "Licensee"); you
# may not use this file except in compliance with the License.  You may
# obtain a copy of the License at:
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
# implied. See the License for the specific language governing
# permissions and limitations under the license.
#### </license> #### 

import sys
import re
import sympy
from sympy.core import S

from melodee.cardioidGenerator import MyCCodeSympyPrinter,CPrintVisitor

from melodee.parser import MelodeeParser,Differentiator
from melodee import utility
from melodee.utility import order


def pretty(symbol):
    return str(symbol)

def generateContinuity(model, targetName):
    template = {}
    template["target"] = targetName

    diffvars = model.diffvars()
    diffvarUpdate = {var : model.diffvarUpdate(var) for var in diffvars}
    gates = model.varsWithAttribute("gate") & diffvars
    params = model.varsWithAttribute("param")
    coupled = model.varsWithAttribute("coupled") & diffvars
    
    inputs=model.inputs()
    
    differ = Differentiator(model, diffvars | params | inputs)
    gateJacobians = {}
    for gate in diffvars:
        (gateJacobians[gate],dontcare) = differ.diff(diffvarUpdate[gate],gate)
    differ.augmentInstructions()

    dt = model.addSymbol("_dt")
    approxvars = set([dt])
    if model.time != None:
        approxvars.add(model.time)

    expensiveVars = model.extractExpensiveFunctions()
    model.makeNamesUnique()

    computeTargets = set()
    computeTargets |= set([diffvarUpdate[var] for var in diffvars])
    computeTargets |= set([gateJacobians[var] for var in diffvars])

    statevars = model.inputs()|diffvars
    computeAllDepend = model.allDependencies(approxvars|statevars, computeTargets)
    constants = model.allExcluding(approxvars, statevars) & computeAllDepend

    varCount = {}
    ii=0
    for var in order(coupled):
        varCount[var] = ii
        ii += 1
    for var in order(diffvars-coupled):
        varCount[var] = ii
        ii += 1

    ii=0
    for var in order(inputs):
        varCount[var] = ii
        ii += 1
    for var in order(params):
        varCount[var] = ii
        ii += 1

    fieldNumbers = {}
    ii=1
    for var in order(inputs):
        fieldNumbers[var] = ii
        ii += 1
    for var in order(coupled):
        fieldNumbers[var] = ii
        ii += 1
    out = utility.Indenter(open(targetName+".txt","w"))
    out("	component			%(target)s			StateVar", template)
    for var in order(coupled):
        out("voltage	spatially coupled	[['Field %d', '']]	d%s_dt = 0	StateVar			%s", fieldNumbers[var], pretty(var), pretty(var))
    for var in order(diffvars-coupled):
        out("voltage	state	[['Value', '1']]	d%s_dt = 0	StateVar			%s", pretty(var),pretty(var))
    out("param root	component			__root__			Parameters")
    out.inc()
    for var in order(inputs):
        out("parameter	[['Field %s', '']]		Parameters		%s", fieldNumbers[var], pretty(var))
        
    for var in order(params):
        out("parameter	[['Value', '%s']]		Parameters			%s", str(model.ssa[var].sympy), pretty(var))
    out.dec()
    

    out = utility.Indenter(open(targetName+".continuity.c","w"))
    out('''
void cpu_advance_be1(REAL _t, REAL _t_end, REAL *_y_global, REAL *_y_global_temp, REAL *_rpar_global, int num_gauss_pts, int num_gpus, int gp_offset)
{
    //REAL bi;
    REAL _dt = _t_end-_t;


''' % template)
    out.inc()
    if model.time:
        out("REAL %s = _t;", pretty(model.time))
    good = set()
    good |= approxvars
    for var in order(diffvars):
        out("REAL %s=_y_global[%d];", pretty(var), varCount[var]);
    good |= statevars
    for var in order(params|inputs):
        out("REAL %s=_rpar_global[%d];", pretty(var), varCount[var]);
    good |= params

    cprinter = CPrintVisitor(out, model.ssa, good, "REAL")
    model.printTarget(good,computeTargets,cprinter)
    

    for var in order(diffvars):
        out("_y_global_temp[%d] = _y_global[%d] + %s*_dt/(1-_dt*%s);", varCount[var], varCount[var], diffvarUpdate[var], gateJacobians[var])
    
    out.dec()
    out('''

//assigned_vars: 10

}
''' % template)

generators = {
    frozenset(["continuity"]): generateContinuity,
}
