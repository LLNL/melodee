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
from sympy.printing.octave import OctaveCodePrinter


from melodee.parser import MelodeeParser
from melodee import utility
from melodee.utility import order

def pretty(symbol):
    if str(symbol)[0] == '_':
        return "U"+str(symbol)
    else:
        return str(symbol)
    return str(symbol)

class MyOctaveCodeSympyPrinter(OctaveCodePrinter):
    def __init__(self,*args,**kwargs):
        OctaveCodePrinter.__init__(self, *args, **kwargs)
    def _print_Relational(self,expr):
        if expr.rel_op == "==" or expr.rel_op == "!=":
            PREC = sympy.printing.precedence.precedence(expr)
            return "%s %s %s" % (self.parenthesize(expr.lhs, PREC),
                                 expr.rel_op,
                                 self.parenthesize(expr.rhs, PREC))
        else:
            return super(MyOctaveCodeSympyPrinter, self)._print_Relational(expr)
    def _print_Symbol(self,symbol):
        return pretty(symbol)

class MatlabPrintVisitor:
    def __init__(self, out, ssa, params):
        self.out = out
        self.ssa = ssa
        self.params = params
        self.oprinter = MyOctaveCodeSympyPrinter()
        
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
        rhsText = self.oprinter.doprint(rhs.sympy)
        lhsText = self.oprinter.doprint(lhs)
        if lhs in self.params:
            self.out("""
if (isfield(U_params, '%(name)s'))
   %(name)s = U_params.%(name)s;
else
   %(name)s = %(rhs)s;
   U_params.%(name)s = %(name)s;
end""", name=lhsText, rhs=rhsText)
        else:
            self.out("%s = %s;", lhsText, rhsText)

def numberVars(diffvars):
    ret = {}
    count = 1
    for var in order(diffvars):
        ret[var] = count
        count += 1
    return ret

def generateMatlab(model, targetName):
    template = {}
    template["target"] = targetName

    inputs = model.inputs()
    if inputs:
        template["arglower"] = 1
    else:
        template["arglower"] = 0

    good = set()
    if model.time != None:
        good.add(model.time)
        timename = str(model.time)
    else:
        timename = "U_current_time"
    template["timename"] = timename
        
    out = utility.Indenter(open(targetName+"_init.m","w"))
    params = model.varsWithAttribute("param")
    printer = MatlabPrintVisitor(out,model.ssa,params)
    out("""
function [U_y_init, U_ordering, U_params] = %(target)s_init(%(timename)s,varargin)
   narginchk(1+%(arglower)d,2);
   if (nargin >= 1)
      U_params = varargin{1};
   else
      U_params = struct();
   end
""" % template)
    out.inc()

    if inputs:
        out("%define the inputs")
        good |= inputs
    for symbol in order(inputs):
        out("%s = U_params.%s(%s);",pretty(symbol),pretty(symbol),timename)

    out("\n\n")
    out("%define the initial conditions")
    diffvars = model.diffvars()
    model.printTarget(good,params|diffvars,printer)

    diffvarNumbering = numberVars(diffvars)
    out("U_y_init = zeros(%d, 1);", len(diffvarNumbering))
    for var in order(diffvars):
        out("U_y_init(%d) = %s;",diffvarNumbering[var],pretty(var))

    out("U_ordering = struct();")
    for var in order(diffvars):
        out("U_ordering.%s = %d;", pretty(var),diffvarNumbering[var])
        
    out.dec()
    out("""
end
""" % template)

    out = utility.Indenter(open(targetName+".m","w"))
    printer = MatlabPrintVisitor(out,model.ssa,params)
    out("""
function U_dydt = %(target)s(%(timename)s,U_diffvars,varargin)
""" % template)
    out.inc()
    out("U_dydt = zeros(%d,1);" % len(diffvarNumbering))
    
    good = set()
    if model.time != None:
        good.add(model.time)

    out("% make copies of the differential vars")
    for var in order(diffvars):
        out("%s = U_diffvars(%d);", pretty(var), diffvarNumbering[var])
    good |= diffvars

    out("""
narginchk(2+%(arglower)d,3);
if (nargin >= 3)
    U_params = varargin{1};
else
   U_params = struct();
end
""" % template)
        
    if inputs:
        out("% define all inputs")
        good |= inputs
        for symbol in order(inputs):
            out("%s = U_params.%s(%s);",pretty(symbol),pretty(symbol),timename)

    out("% define the differential update")
    diffvarUpdate = {var : model.diffvarUpdate(var) for var in diffvars}
    model.printTarget(good,set(diffvarUpdate.values()),printer)

    out("% stuff the differential update into an array")
    for var in order(diffvars):
        out("U_dydt(%d) = %s;", diffvarNumbering[var], pretty(diffvarUpdate[var]))

    out.dec()
    out("""
end""" % template)

    out = utility.Indenter(open(targetName+"_trace.m","w"))
    printer = MatlabPrintVisitor(out, model.ssa, params)
    out("""
function U_trace = %(target)s_trace(%(timename)s,U_diffvars,varargin)
""" % template)
    out.inc()
    out("% make copies of the differential vars")
    for var in order(diffvars):
        out("%s = U_diffvars(:,%d);", pretty(var), diffvarNumbering[var])
    out("""
narginchk(2+%(arglower)d,3);
if (nargin >= 3)
    U_params = varargin{1};
else
   U_params = struct();
end
""" % template)
    if inputs:
        out("% define all inputs")
        good |= inputs
        for symbol in order(inputs):
            out("%s = U_params.%s(%s);",pretty(symbol),pretty(symbol),timename)

    tracevars = model.varsWithAttribute("trace")
    out("%calculate the tracing vars we need")
    model.printTarget(good,tracevars-good,printer)
    
    out("%save the tracing vars")
    out("U_trace = struct();")
    for var in order(tracevars):
        out('U_trace.%s = %s;', var, var)
    out.dec()
    out("""
end""" % template)

generators = {
    frozenset(["matlab"]): generateMatlab,
}
