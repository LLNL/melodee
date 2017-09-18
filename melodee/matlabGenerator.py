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
from sympy.printing.octave import octave_code


from melodee.parser import MelodeeParser
from melodee import utility
from melodee.utility import order

def pretty(symbol):
    return str(symbol)


class MatlabPrintVisitor:
    def __init__(self, out, ssa, params):
        self.out = out
        self.ssa = ssa
        self.params = params
    def ifPrint(self,printer,ifSymbol,thenList,elseList,choiceList):
        self.out("if (%s)",pretty(ifSymbol))
        self.out.inc()
        printer(thenList)
        for choiceVar in choiceList:
            choice = self.ssa[choiceVar]
            lhs = pretty(choiceVar)
            rhs = pretty(choice.thenVar)
            if lhs != rhs:
                self.out("%s = %s;",lhs,rhs)
        self.out.dec()
        self.out("else")
        self.out.inc()
        printer(elseList)
        for choiceVar in choiceList:
            choice = self.ssa[choiceVar]
            lhs = pretty(choiceVar)
            rhs = pretty(choice.elseVar)
            if lhs != rhs:
                self.out("%s = %s;",lhs,rhs)
        self.out.dec()
        self.out("end")
    def equationPrint(self,lhs,rhs):
        rhsText = octave_code(rhs.sympy)
        if lhs in self.params:
            self.out("""
if (isfield(__params, '%(name)s'))
   %(name)s = __params.%(name)s;
else
   %(name)s = %(rhs)s;
   __params.%(name)s = %(name)s;
end""", name=pretty(lhs), rhs=rhsText)
        else:
            self.out("%s = %s;", pretty(lhs), rhsText)

def numberVars(diffvars):
    ret = {}
    count = 1
    for var in order(diffvars):
        ret[var] = count
        count += 1
    return ret

def generateMatlab(model, targetName, initFile, diffFile, traceFile):
    template = {}
    template["target"] = targetName

    inputs = model.inputs()
    if inputs:
        template["arglower"] = 1
    else:
        template["arglower"] = 0

    out = initFile
    params = model.varsWithAttribute("param")
    printer = MatlabPrintVisitor(out,model.ssa,params)
    out("""
function [__y_init, __ordering, __params] = %(target)s_init(varargin)
   narginchk(0+%(arglower)d,1);
   if (nargin >= 1)
      __params = varargin{1};
   else
      __params = struct();
   end
""" % template)
    out.inc()

    
    #define time
    good = set()
    if model.time != None:
        timename = pretty(model.time)
        out("% define time")
        out("%s = 0;",timename)
        good.add(model.time)
    else:
        timename = "__current_time"
    template["timename"] = timename

    if inputs:
        out("%define the inputs")
        good |= inputs
    for symbol in order(inputs):
        out("%s = __params.%s(%s);",pretty(symbol),pretty(symbol),timename)

    out("\n\n")
    out("%define the initial conditions")
    diffvars = model.diffvars()
    model.printTarget(good,params|diffvars,printer)

    diffvarNumbering = numberVars(diffvars)
    out("__y_init = zeros(%d, 1);", len(diffvarNumbering))
    for var in order(diffvars):
        out("__y_init(%d) = %s;",diffvarNumbering[var],pretty(var))

    out("__ordering = struct();")
    for var in order(diffvars):
        out("__ordering.%s = %d;", pretty(var),diffvarNumbering[var])
        
    out.dec()
    out("""
end
""" % template)

    out = diffFile
    printer = MatlabPrintVisitor(out,model.ssa,params)
    out("""
function __dydt = %(target)s(%(timename)s,__diffvars,varargin)
""" % template)
    out.inc()
    out("__dydt = zeros(%d,1);" % len(diffvarNumbering))
    
    good = set()
    if model.time != None:
        good.add(model.time)

    out("% make copies of the differential vars")
    for var in order(diffvars):
        out("%s = __diffvars(%d);", pretty(var), diffvarNumbering[var])
    good |= diffvars

    out("""
narginchk(2+%(arglower)d,3);
if (nargin >= 3)
    __params = varargin{1};
else
   __params = struct();
end
""" % template)
        
    if inputs:
        out("% define all inputs")
        good |= inputs
        for symbol in order(inputs):
            out("%s = __params.%s(%s);",pretty(symbol),pretty(symbol),timename)

    out("% define the differential update")
    diffvarUpdate = {var : model.diffvarUpdate(var) for var in diffvars}
    model.printTarget(good,set(diffvarUpdate.values()),printer)

    out("% stuff the differential update into an array")
    for var in order(diffvars):
        out("__dydt(%d) = %s;", diffvarNumbering[var], pretty(diffvarUpdate[var]))

    out.dec()
    out("""
end""" % template)

    out = traceFile
    printer = MatlabPrintVisitor(out, model.ssa, params)
    out("""
function __trace = %(target)s_trace(%(timename)s,__diffvars,varargin)
""" % template)
    out.inc()
    out("% make copies of the differential vars")
    for var in order(diffvars):
        out("%s = __diffvars(:,%d);", pretty(var), diffvarNumbering[var])
    out("""
narginchk(2+%(arglower)d,3);
if (nargin >= 3)
    __params = varargin{1};
else
   __params = struct();
end
""" % template)
    if inputs:
        out("% define all inputs")
        good |= inputs
        for symbol in order(inputs):
            out("%s = __params.%s(%s);",pretty(symbol),pretty(symbol),timename)

    tracevars = model.varsWithAttribute("trace")
    out("%calculate the tracing vars we need")
    model.printTarget(good,tracevars-good,printer)
    
    out("%save the tracing vars")
    out("__trace = struct();")
    for var in order(tracevars):
        out('__trace.%s = %s;', var, var)
    out.dec()
    out("""
end""" % template)

def main():
    import sys
    import argparse
    ap=argparse.ArgumentParser(description="Convert Melodee into Matlab functions")
    ap.add_argument("target",
                    help="Comma separated list of subsystems to generate",
                    )
    ap.add_argument("filenames", nargs=argparse.REMAINDER,
                    help="Filenames containing all the targets",
                    )
    options=ap.parse_args()
    target = options.target
    filenames = options.filenames
    
    p = MelodeeParser()

    models = {}
    for filename in filenames:
        p.parse(open(filename,"r").read())
    originalTargets = target.split(",")
    if len(originalTargets) > 1:
        target = "_".join(originalTargets)
        targetModel = ""
        targetModel += "subsystem %s {\n" % target
        for model in originalTargets:
            targetModel += "   use %s;\n" % model
        targetModel += "}\n"
        p.parse(targetModel)
        
    generateMatlab(p.getModel(target), target,
                   #utility.Indenter(sys.stdout),
                   #utility.Indenter(sys.stdout),
                   utility.Indenter(open(target+"_init.m","w")),
                   utility.Indenter(open(target+".m","w")),
                   utility.Indenter(open(target+"_trace.m","w")),
    )
