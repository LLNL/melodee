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


from melodee.parser import MelodeeParser,Differentiator
from melodee import utility
from melodee.utility import order

from melodee.cardioidGenerator import MyCCodeSympyPrinter, repeat, CPrintVisitor

def printParamConst(out,thesevars):
    out.inc()
    for var in order(thesevars):
        out("GlobalData_t %s = p->%s;",var,var)
    out.dec()

def stateName(var,statevars):
    if var in statevars:
        return "sv->%s" % var
    else:
        return str(var)

class LookupPrintVisitor(CPrintVisitor):
    def __init__(self,out, ssa, declared, interpTargets, decltype="double"):
        self.interpTargets = interpTargets
        self.fitFromTarget = {}
        for (fit,targets) in interpTargets.items():
            for target in targets:
                self.fitFromTarget[target] = fit
        super(LookupPrintVisitor,self).__init__(out,ssa,declared,decltype)
    
    def printLookup(self,lhs):
        if lhs in self.interpTargets:
            self.out("LUT_data_t %s_row[NROWS_%s];",lhs,lhs)
            self.out("LUT_interpRow(&IF->tables[_%s_TAB], %s, __i, %s_row);",lhs,lhs,lhs) 

    def equationPrint(self,lhs,rhs):
        if lhs in self.fitFromTarget:
            rhsText = "%s_row[%s_idx]" % (self.fitFromTarget[lhs],lhs)
            self.equationPrintWithRhs(lhs,rhsText)
        else:
            super(LookupPrintVisitor,self).equationPrint(lhs,rhs)
        if lhs in self.interpTargets:
            self.printLookup(lhs)
    def printChoice(self, lhs):
        if lhs in self.fitFromTarget:
            self.equationPrint(lhs,None)

def generateCarp(model, targetName):
    si = model.si
    unitFromExternalName = {"Vm" : si.get("mV"),
                            "Lambda" : si.get("1"),
                            "delLambda" : si.get("1")/si.get("ms"),
                            "Tension" : si.get("kPa"),
                            "K_e" : si.get("mM"),
                            "Na_e" : si.get("mM"),
                            "Ca_e" : si.get("uM"),
                            "Iion" : si.get("uA")/si.get("uF"),
                            "tension_component" : si.get("1"),
                            "K_i" : si.get("mM"),
                            "Na_i" : si.get("mM"),
                            "Ca_i" : si.get("uM"),
                            "Ca_b" : si.get("uM"),
                            "CaRel" : si.get("uM"),
                            "CaUp" : si.get("uM"),
                            "ATP" : si.get("mM"),
                            "ADP" : si.get("uM"),
                            }
    publicNames = set(["Vm",
                       "Lambda",
                       "delLambda",
                       "Tension",
                       "K_e",
                       "Na_e",
                       "Ca_e",
                       "Iion",
                       "tension_component",
                       ])
    privateNames = set(["K_i",
                        "Na_i",
                        "Ca_i",
                        "Ca_b",
                        "CaRel",
                        "CaUp",
                        "ATP",
                        "ADP",
                        ])

    template = {}
    template["target"] = targetName
    template["uppercaseModel"] = targetName.upper()
    template["timeFactor"] = "1" #FIXME
    
    
    reqvars=None # anything with a req
    modvars=None # anything with a mod
    statevars=None # anything with a s->
    paramvars=None # anything with a p->
    flagvars=None # anything that can be used as a flag
    nodevars=None
    diffvars=None
    simvars=None
    interpvars=None # interp vars.
    interpTargets=None
    fevars = None
    tracevars = None
    
    ## Vars for printing
    externalNameFromVar = {}
    for var in model.varsWithAttribute("external"):
        externalName = model.info("external",var)
        if not externalName:
            externalName = str(var)
        externalNameFromVar[var] = externalName

    diffvars = model.diffvars()
    exvars = set(externalNameFromVar.keys())
    reqvars = exvars & set(model.inputs())
    modvars = exvars & set(model.ssa.keys())

    paramvars = model.varsWithAttribute("param")
    flagvars = model.varsWithAttribute("flag")
    nodevars = model.varsWithAttribute("nodal")
    tracevars = model.varsWithAttribute("trace")
    
    statevars = (diffvars | nodevars) - exvars

    fevars = set(diffvars)
    
    simvars = set()
    dt = model.addSymbol("_dt")
    if model.time:
        simvars.add(model.time)
    simvars.add(dt)

    diffvarUpdate = {var : model.diffvarUpdate(var) for var in diffvars}
    V = model.input("V")
    V_init = model.output("V_init")
    Iion = model.output("Iion")

    modvars.add(Iion); externalNameFromVar.setdefault(Iion,"Iion")
    reqvars.add(V); externalNameFromVar.setdefault(V,"Vm")
    
    rushvars = model.varsWithAttribute("gate") & diffvars
    differ = Differentiator(model, diffvars | paramvars | reqvars | modvars)
    rushJacobians = {}
    for rush in order(rushvars):
        (rushJacobians[rush],dontcare) = differ.diff(diffvarUpdate[rush],rush)
    differ.augmentInstructions()
    rushTargets = {}
    for rush in order(rushvars):
        F = model.ssa[diffvarUpdate[rush]].sympy
        L = model.ssa[rushJacobians[rush]].sympy
        M = (F-L*rush).simplify()
        
        RLA = model.addInstruction("_%s_RLA" % rush, sympy.exp(dt*L))
        RLB = model.addInstruction("_%s_RLB" % rush, M/L*(RLA-1))
        rushTargets[rush] = model.addInstruction("_%s_update" % rush,rush*RLA+RLB)

    fevars -= rushvars
        
    expensiveVars = model.extractExpensiveFunctions()
    model.makeNamesUnique()

    computeTargets = set()
    computeTargets.add(Iion)
    computeTargets |= set([diffvarUpdate[var] for var in diffvars-rushvars])
    for gate in rushvars:
        computeTargets.add(rushTargets[gate])

    approxvars = set([dt])
    computeAllDepend = model.allDependencies(approxvars|statevars|reqvars, computeTargets)

    constants = model.allExcluding(approxvars, statevars|reqvars) & computeAllDepend

    interpvars = model.varsWithAttribute("interp")
    interpTargets = {}
    allinterps = set()
    for fit in interpvars:
        good = approxvars | set([fit])
        dependsOnlyOnFitvar = model.allExcluding(good, (statevars|reqvars)-good)
        interpCandidates = (dependsOnlyOnFitvar & computeAllDepend) - constants

        expensiveInterps = expensiveVars & interpCandidates
        inexpensiveInterps = model.allExcluding(good, ((statevars|reqvars)-good)|expensiveInterps)
        interpCandidates -= inexpensiveInterps
        
        externallyUsedInterps = (
            model.allDependencies(approxvars|statevars|reqvars|interpCandidates, computeTargets)
            &
            interpCandidates
            )
        interpTargets[fit] = externallyUsedInterps
        allinterps |= externallyUsedInterps
        
    computeAllDepend = model.allDependencies(
        approxvars|statevars|reqvars|allinterps,
        computeTargets)
    constants = model.allExcluding(approxvars,statevars|reqvars) & computeAllDepend

    out = utility.Indenter(open(targetName+".h","w"), "  ")
    cprinter=CPrintVisitor(out, model.ssa, set())
    out('''
//// HEADER GUARD ///////////////////////////
// If automatically generated, keep above
// comment as first line in file.  
#ifndef __%(uppercaseModel)s_H__
#define __%(uppercaseModel)s_H__
//// HEADER GUARD ///////////////////////////
// DO NOT EDIT THIS SOURCE CODE FILE
// ANY CHANGES TO THIS FILE WILL BE OVERWRITTEN!!!!
//   If you need to make a change, do what you have to do and send Rob
//   an email with what you had to change and why.  He'll fix the translator
//   so your fix will be recorded the next time the translator runs.

''', template)
    out('#define %(target)s_REQDAT 0|%(req)s',
        target=targetName,
        req="|".join(["%s_DATA_FLAG"% externalNameFromVar[var] for var in reqvars])
        )
    out('#define %(target)s_MODDAT 0|%(mod)s',
        target=targetName,
        mod="|".join(["%s_DATA_FLAG"% externalNameFromVar[var] for var in modvars])
        )
    out('''

typedef struct %(target)s_params {
''', template)
    out.inc()
    for var in order(paramvars):
        out("GlobalData_t %s;", var)
    out.dec()
    out('''
    char* flags;
} %(target)s_Params;''', template)
    
    out('static char* %(target)s_flags = "%(flags)s";',
        target=targetName,
        flags="|".join([str(var) for var in flagvars])
        )
    out('''

typedef struct %(target)s_state {''', template)
    out.inc()
    for var in order(statevars):
        out("GlobalData_t %s;", var)
    out.dec()
    out('''
} %(target)s_state;

#ifdef __CUDACC__
extern "C" {
#endif

void initialize_params_%(target)s(ION_IF *);
void construct_tables_%(target)s(ION_IF *);
void destroy_%(target)s(ION_IF *);
void initialize_sv_%(target)s(ION_IF *, GlobalData_t**);
GLOBAL void compute_%(target)s(int, int, ION_IF *, GlobalData_t**);

#ifdef __CUDACC__
};
#endif

//// HEADER GUARD ///////////////////////////
#endif
//// HEADER GUARD ///////////////////////////

''', template)

    out = utility.Indenter(open(targetName+".c","w"), "  ")
    out('''
// DO NOT EDIT THIS SOURCE CODE FILE
// ANY CHANGES TO THIS FILE WILL BE OVERWRITTEN!!!!
//   If you need to make a change, do what you have to do and send Rob
//   an email with what you had to change and why.  He'll fix the translator
//   so your fix will be recorded the next time the translator runs.

#include "ION_IF.h"
#include "NumComp.h"
#include "%(target)s.h"


#ifdef __CUDACC__
#define pow powf
#define log logf
#endif

void trace_%(target)s(ION_IF* IF, int node, FILE* file, GlobalData_t** impdata);

void destroy_%(target)s( ION_IF *IF )
{
  destroy_luts( IF );
  SV_free( &IF->sv_tab );
  // rarely need to do anything else
}

void initialize_params_%(target)s( ION_IF *IF )
{
  cell_geom *region = &IF->cgeom;
  %(target)s_Params *p = (%(target)s_Params *)IF->params;

''', template)
    out.inc()

    good = set()
    
    out("//Compute the paramaters")
    cprinter=CPrintVisitor(out, model.ssa, set())
    model.printTarget(good, paramvars, cprinter)
    for var in order(paramvars):
        out("p->%s = %s;", var, var)

    out.dec()
    out('''
    //p->Bufc = 0.2;
    //p->cell_type = EPI;
    //if (0) ;
    //else if (flag_set(p->flags, "EPI")) p->cell_type = EPI;
    //else if (flag_set(p->flags, "MCELL")) p->cell_type = MCELL;
    //else if (flag_set(p->flags, "ENDO")) p->cell_type = ENDO;

  IF->type->trace = trace_%(target)s;

}


// Define the parameters for the lookup tables
enum Tables {
''', template)
    out.inc()
    for var in order(interpvars):
        out("_%s_TAB,",var)
    out.dec()
    out('''

  N_TABS
};

// Define the indices into the lookup tables.
''', template)
    for var in order(interpvars):
        out("enum %s_TableIndex {" % var)
        out.inc()
        for target in order(interpTargets[var]):
            out("%s_idx,", target)
        out("NROWS_%s", var)
        out.dec()
        out("};")
    out('''

void construct_tables_%(target)s( ION_IF *IF )
{
  GlobalData_t _dt = IF->dt*%(timeFactor)s;
  cell_geom *region = &IF->cgeom;
  %(target)s_Params *p = (%(target)s_Params *)IF->params;

  IF->numLUT = N_TABS;
  IF->tables = (LUT *)IMP_malloc( N_TABS, sizeof(LUT) );''', template)
    printParamConst(out,paramvars-nodevars)
    good = set()
    good |= paramvars-nodevars
    good.add(dt)
    out('''
  
''', template)
    for var in order(interpvars):
        (lb,ub,inc) = model.info("interp",var).split(",")
        out('''
  // Create the %(var)s lookup table
  LUT* %(var)s_tab = &IF->tables[_%(var)s_TAB];
  {
     double epsilon=fabs(%(lb)s)*1e-7;
     LUT_alloc(%(var)s_tab, NROWS_%(var)s, %(lb)s+epsilon, %(ub)s+epsilon, %(inc)s, "%(target)s %(var)s");
  }
  for (int __i=%(var)s_tab->mn_ind; __i<=%(var)s_tab->mx_ind; __i++) {
    double %(var)s = %(var)s_tab->res*__i;
    LUT_data_t* %(var)s_row = %(var)s_tab->tab[__i];
''',
            target=targetName,
            var=var,
            lb=lb,
            ub=ub,
            inc=inc,
            )
        out.inc(2)
        cprinter=CPrintVisitor(out, model.ssa, set())
        model.printTarget(good | set([var]), interpTargets[var], cprinter)
        for target in order(interpTargets[var]):
            out('%(var)s_row[%(target)s_idx] = %(target)s;',
                var=var,
                target=target)
        out.dec()
        out('}')
        out('check_LUT(%(var)s_tab);', var=var)
        out.dec()

    out('''
}



void    initialize_sv_%(target)s( ION_IF *IF, GlobalData_t **impdata )
{
  GlobalData_t _dt = IF->dt*%(timeFactor)s;
  cell_geom *region = &IF->cgeom;
  %(target)s_Params *p = (%(target)s_Params *)IF->params;

  SV_alloc( &IF->sv_tab, IF->numNode, sizeof(%(target)s_state) );
  %(target)s_state *sv_base = (%(target)s_state *)IF->sv_tab.y;
  GlobalData_t t = 0;''', template)
    printParamConst(out,paramvars-nodevars)
    good = set()
    good |= paramvars-nodevars
    good.add(dt)
    out('''

  //Prepare all the public arrays.''', template)
    out.inc()
    for var in order(reqvars|modvars):
        out("GlobalData_t* %s_ext = impdata[%s];",var,externalNameFromVar[var])
    out.dec()
    out('''
  //Prepare all the private functions.

  //set the initial values
  for(int __i=0; __i<IF->sv_tab.numSeg; __i++ ){
    %(target)s_state *sv = sv_base+__i;''', template)
    out.inc(2)

    out("//Get the definition of nodal vars")
    for var in order(paramvars & nodevars):
        out("GlobalData_t %s = p->%s;", var, var)
    good |= paramvars & nodevars

    out("//Read in the input vars")
    for var in order(reqvars-set([V])):
        out("GlobalData_t %s = %s_ext[__i];",var, var)
    good |= reqvars-set([V])

    out("//Define voltage")
    valid = paramvars|(reqvars-set([V]))
    cprinter=CPrintVisitor(out, model.ssa, set())
    model.printTarget(valid,set([V_init]),cprinter)
    out("double %s = %s;",V,V_init)
    good |= set([V])
    
    out("//Define statevars and modvars")
    valid |= set([V])
    model.printTarget(valid,(statevars|modvars)-valid,cprinter)

    out("//Save all the state variables")
    for var in order(statevars):
        out("sv->%s = %s;",var,var)
    out("//Save all external vars")
    for var in order(modvars):
        out("%s_ext[__i] = %s;",var,var)
    out("//Save voltage")
    out("%s_ext[__i] = %s;",V,V);

    out.dec(2)
    out('''
    //Change the units of external variables as appropriate.
    //sv->Ca_i *= 1e-3;
  }

}

/** compute the  current
 *
 * param start   index of first node
 * param end     index of last node
 * param IF      IMP
 * param plgdata external data needed by IMP
 */
GLOBAL void compute_%(target)s(int start, int end, ION_IF *IF, GlobalData_t **impdata )
{
  GlobalData_t _dt = IF->dt*%(timeFactor)s;
  cell_geom *region = &IF->cgeom;
  %(target)s_Params *p  = (%(target)s_Params *)IF->params;
  %(target)s_state *sv_base = (%(target)s_state *)IF->sv_tab.y;

  GlobalData_t t = IF->tstp.cnt*_dt;

''', template)
    printParamConst(out,paramvars-nodevars)
    good = set()
    good |= paramvars-nodevars
    good.add(dt)
    out('''

  //Prepare all the public arrays.''', template)
    out.inc()
    for var in order(reqvars|modvars):
        out("GlobalData_t* %s_ext = impdata[%s];",var,externalNameFromVar[var])
    out.dec()
    out('''

#ifdef __CUDACC__
  int __i = blockDim.x * blockIdx.x + threadIdx.x;
  if ( __i>=end ) { return; }
  {
#else
#pragma omp parallel for schedule(static)
  for (int __i=start; __i<end; __i++) {
#endif
    %(target)s_state *sv = sv_base+__i;

''', template)
    out.inc(2)
    
    out("//Read in the input vars")
    for var in order(reqvars):
        out("GlobalData_t %s = %s_ext[__i];",var, var)
    good |= reqvars
        
    out("//Get the definition of statevars")
    for var in order(statevars):
        out("GlobalData_t %s = sv->%s;", var, var)
    good |= statevars
        
    out('//Change the units of external variables as appropriate.')
    out('//sv->Ca_i *= 1e-3;')
    
    cprinter=LookupPrintVisitor(out, model.ssa, set(),interpTargets)
    out('//Compute lookup tables for things that have already been defined.')
    for var in order(good & set(interpTargets.keys())):
        cprinter.printLookup(var)

    out('//Compute external modvars')
    justDefined = model.allDependencies(good|allinterps,modvars-diffvars)-good
    model.printSet(justDefined,cprinter)
    good |= justDefined
    
    out('//Complete Forward Euler Update', template)
    justDefined = model.allDependencies(good|allinterps,set([diffvarUpdate[var] for var in fevars]))-good
    model.printSet(justDefined,cprinter)
    good |= justDefined
    
    out('//Complete Rush Larsen Update', template)
    justDefined = model.allDependencies(good|allinterps,set(rushTargets.values()))-good
    model.printSet(justDefined,cprinter)
    good |= justDefined

    out('//Finish the update', template)
    for var in order(fevars):
        out("%s += %s*_dt;",stateName(var,statevars),diffvarUpdate[var])
    for var in order(rushvars):
        out("%s = %s;",stateName(var,statevars),rushTargets[var])
    out('//Change the units of external variables as appropriate.', template)    
    out('//sv->Ca_i *= 1e3;',template)

    for var in order(modvars):
        out("%s_ext[__i] = %s;", var,var)
    out.dec(2)
    out('''
  }    


}


void trace_%(target)s(ION_IF* IF, int node, FILE* file, GlobalData_t** impdata) 
{
  static bool first = true;
  if (first) {
    first = false;
    FILE_SPEC theader = f_open("%(target)s_trace_header.txt","wt");
    f_printf( theader, 
''', template)
    out.inc(2)
    for var in order(tracevars):
        out('"%s\\n"',var)
    out.dec(2)
    out('''
      );

    f_close(theader);
  }
  
  CUDA_TRACE( %(target)s, IF, node, file )
#ifdef HAVE_CUDA_LIMPET
#define fprintf(A,F, ...)  printf( F, __VA_ARGS__ )
#endif

  GlobalData_t _dt = IF->dt*%(timeFactor)s;
  cell_geom *region = &IF->cgeom;
  %(target)s_Params *p  = (%(target)s_Params *)IF->params;

  %(target)s_state *sv_base = (%(target)s_state *)IF->sv_tab.y;
  %(target)s_state *sv = sv_base+node;
  int __i = node;

  GlobalData_t t = IF->tstp.cnt*_dt;

''', template)
    out.inc()
    printParamConst(out,paramvars-nodevars)
    good = set()
    good |= paramvars-nodevars
    good.add(dt)
    out('//Prepare all the public arrays.', template)

    for var in order(reqvars|modvars):
        out("GlobalData_t* %s_ext = impdata[%s];",var,externalNameFromVar[var])

    out("//Read in the input vars")
    for var in order(reqvars):
        out("GlobalData_t %s = %s_ext[__i];",var, var)
    good |= reqvars
        
    out("//Get the definition of statevars")
    for var in order(statevars):
        out("GlobalData_t %s = sv->%s;", var, var)
    good |= statevars
        
    out('//Change the units of external variables as appropriate.')
    out('//sv->Ca_i *= 1e-3;')

    cprinter=CPrintVisitor(out, model.ssa, set())
    model.printTarget(good,tracevars-good,cprinter)

    out('//Output the desired variables')
    for var in order(tracevars):
        out('fprintf(file, "%%4.12f\\t", %s);',var)

    out.dec()
    out('''
  //Change the units of external variables as appropriate.
  //sv->Ca_i *= 1e3;

#ifdef HAVE_CUDA_LIMPET
#undef fprintf
#endif

}''', template)

generators = {
    frozenset(["carp"]) : generateCarp,
}
