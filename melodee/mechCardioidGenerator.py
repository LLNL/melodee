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

from melodee.parser import MelodeeParser,Differentiator
from melodee import utility
from melodee.utility import order,itemsOrderedByKey

from melodee.cardioidGenerator import CPrintVisitor

def generateCardioidMech(model, targetName):
    template = {}
    template["target"] = targetName

    diffvars = model.diffvars()
    diffvarUpdate = {var : model.diffvarUpdate(var) for var in diffvars}
    flagvars = model.varsWithAttribute("flag")
    paramvars = model.varsWithAttribute("param")
    tracevars = model.varsWithAttribute("trace")
    simvars = set()
    
    stretch = model.input("stretch")
    tension = model.output("tension")
    if "stretchVel" not in model._inputs:
        stretchVel = model.addSymbol("stretchVel")
    else:
        stretchVel = model.input("stretchVel")
    inputvars = set([stretch,stretchVel])|model.inputs()
    
    partialFromDiff = {}
    diffPartialFromDiff = {}

    differ = Differentiator(model, inputvars-set([stretch]))
    partialvars = set()
    for diffvar in diffvars:
        partial = model.addSymbol("_partial_"+str(diffvar))
        partialvars.add(partial)
        differ.cacheResult(diffvar, stretchVel, partial, None)
        partialFromDiff[diffvar] = partial
    partialStretch = model.addSymbol("_partial_stretch")
    differ.cacheResult(stretch, stretchVel, partialStretch, None)

    for diffvar in diffvars:
        (diffvarUpdate[partialFromDiff[diffvar]],dontcare) = differ.diff(diffvarUpdate[diffvar],stretchVel)
    (dtension_temp,dontcare) = differ.diff(tension, stretchVel)
    differ.augmentInstructions()
    dtension_dstretchVel = model.addInstruction("dtension_dstretchVel",dtension_temp)
    
    outputvars = set([tension,dtension_dstretchVel])    

    out = utility.Indenter(open(targetName+".hh","w"))
    out(r'''
#include "ExcitationContraction.hh"
#include <vector>
#include <string>

namespace %(target)s
{

struct State
{''',template)
    out.inc()
    for var in order(diffvars|partialvars):
        out("double %s;", var)
    out.dec()
    out('''
};

class ThisModel : public ExcitationContraction
{
 public:
   ThisModel(const int numPoints);
   virtual std::string name() const;
   virtual int getHandle(const std::string& varName) const;
   virtual std::string getUnit(const int varHandle) const;

   virtual std::vector<std::string> getVarGroup(const std::string type) const;
   
   virtual double get(const int varHandle) const;
   virtual void   set(const int varHandle, const double value);

   virtual double get(const int varHandle, const int iCell) const;
   virtual void   set(const int varHandle, const int iCell, const double value);

   virtual double get(const int varHandle, const int iCell, const double* inputs[]) const;

   virtual void resetDefaultParameters();
   virtual void initialize(const double *const inputs[]);
   virtual void tryTimestep(const double dt, const double *const inputs[]);
   virtual void outputForNextTimestep(const double dt, const double* const inputs[], double* const outputs[]);
   virtual void commitTimestep();

 private:
   int _numPoints;
   //double _dt;

   //FLAGS
''', template)
    out.inc()
    for var in order(flagvars):
        out("double %s;", var)
    out.dec()
    out(r'''
   //PARAMETERS
''', template)
    out.inc()
    for var in order(paramvars):
        out("double %s;", var)
    out.dec()
    out(r'''

   //STATE
   std::vector<State> _state[2];
   std::vector<State> _partial;
   int _oldIdx;

   inline int _old() const { return _oldIdx; }
   inline int _new() const { return 1 ^ _oldIdx; }
   inline void swapOldAndNew() { _oldIdx = _new(); }
  
};

}
   
''', template)

    out = utility.Indenter(open(targetName+".cpp","w"))
    out(r'''
#include "%(target)s.hh"
#include <cassert>
#include <cmath>

using namespace std;

namespace %(target)s
{

string ThisModel::name() const
{
   return "%(target)s";
}

enum varHandles
{
''', template)
    handlevars = diffvars|inputvars|outputvars|flagvars|paramvars
    out.inc()
    for var in order(handlevars):
        out('_handle_%s,',var)
    out.dec()
    out(r'''
   NUM_HANDLES
};

const string varNames[] = {
''', template)
    out.inc()
    for var in order(handlevars):
        out('"%s",',var)
    out.dec()    
    out(r'''
   ""
};

const string varUnits[] = {
''', template)
    out.inc()
    for var in order(handlevars):
        out('"not_implemented",')
    out.dec()    
    out(r'''
   ""
};

enum inputOrder
{
''', template)
    out.inc()
    for var in order(inputvars):
        out('_inIdx_%s,',var)
    out.dec()
    out(r'''
   NUM_INPUTS
};

enum outputOrder
{
''', template)
    out.inc()
    for var in order(outputvars):
        out('_outIdx_%s,',var)
    out.dec()
    out(r'''
   NUM_OUTPUTS
};

int ThisModel::getHandle(const string& varName) const
{
   for (int ii=0; ii<NUM_HANDLES; ii++)
   {
      if (varName == varNames[ii])
      {
         return ii;
      }
   }
   return -1;
}

string ThisModel::getUnit(const int varHandle) const
{
   assert(varHandle > 0 && varHandle < NUM_HANDLES);
   return varUnits[varHandle];
}

vector<string> ThisModel::getVarGroup(const std::string type) const
{
   vector<string> names;
   if (0) {}
 ''', template)
    out.inc()
    varGroups = {
        "input" : inputvars,
        "output" : outputvars,
        "checkpoint" : diffvars,
        "flag" : flagvars,
        "param" : paramvars,
        }
    for (group,groupvars) in itemsOrderedByKey(varGroups):
        out(r'else if (type == "%s")' % group)
        out(r'{')
        out.inc()
        out(r'names.reserve(%d);', len(groupvars))
        for var in order(groupvars):
            out(r'names.push_back("%s");', var)
        out.dec()
        out(r'}')
    out.dec()
    out(r'''
   return names;
}


double ThisModel::get(const int varHandle) const
{
   if (0) {}
''', template)
    out.inc()
    for var in order(flagvars|paramvars):
        out('else if (varHandle == _handle_%s)',var)
        out('{')
        out('   return %s;',var)
        out('}')
    out.dec()
    out(r'''
   return NAN;
}
void ThisModel::set(const int varHandle, const double value)
{
   if (0) {}
''', template)
    out.inc()
    for var in order(flagvars|paramvars):
        out('else if (varHandle == _handle_%s)',var)
        out('{')
        out('   %s = value;',var)
        out('}')
    out.dec()
    out(r'''
   assert(0 && "Can't set a value for parameter that doesn't exist");
}

double ThisModel::get(const int varHandle, const int iCell) const
{
   if (0) {}
''', template)
    out.inc()
    for var in order(diffvars):
        out('else if (varHandle == _handle_%s)',var)
        out('{')
        out('   return _state[_old()][iCell].%s;',var)
        out('}')
    out.dec()
    out(r'''
   return get(varHandle);
}

void ThisModel::set(const int varHandle, const int iCell, const double value)
{
   if (0) {}
''', template)
    out.inc()
    for var in order(diffvars):
        out('else if (varHandle == _handle_%s)',var)
        out('{')
        out('   _state[_old()][iCell].%s = value;',var)
        out('}')
    out.dec()
    out(r'''
   set(varHandle, value);
}

double ThisModel::get(const int varHandle, const int iCell, const double* inputs[]) const
{
''', template)
    out.inc()
    out('if (0) {}')
    for var in order(tracevars):
        out('else if (varHandle == _handle_%s)',var)
        out('{')
        out('}')
    out.dec()        
    out(r'''
   return get(varHandle, iCell);
}

ThisModel::ThisModel(const int numPoints)
{
   _numPoints = numPoints;
   _oldIdx = 0;
   _state[0].resize(numPoints);
   _state[1].resize(numPoints);

   //set the flags
''',template)
    out.inc()
    good = set(simvars)
    cprinter = CPrintVisitor(out, model.ssa, paramvars|flagvars|simvars)
    model.printTarget(good,flagvars,cprinter)
    good |= flagvars
    out.dec()
    out(r'''
   //set the default parameters
   resetDefaultParameters();
}

void ThisModel::resetDefaultParameters()
{
''', template)
    out.inc()
    cprinter = CPrintVisitor(out, model.ssa, paramvars|flagvars|simvars)
    model.printTarget(flagvars|simvars,paramvars,cprinter)
    good |= paramvars
    out.dec()
    out(r'''
}

void ThisModel::initialize(const double* const _inputs[])
{
   for (int _icell=0; _icell < _numPoints; _icell++)
   {''', template)
    out.inc(2)
    for var in order(inputvars):
        out("const double %s = _inputs[_inIdx_%s][_icell];", var, var)
    good |= inputvars
    cprinter = CPrintVisitor(out, model.ssa, paramvars|flagvars)
    model.printTarget(flagvars,diffvars,cprinter)
    for var in order(diffvars):
        out("_state[_old()][_icell].%s = %s;",var,var)
    good |= diffvars
    out.dec(2)
    out(r'''
   }
}

const int internalTimestep = 20;

void ThisModel::tryTimestep(const double dt, const double* const _inputs[])
{
   for (int _icell=0; _icell < _numPoints; _icell++)
   {''', template)
    out.inc(2)
    for var in order(inputvars-set([stretch])):
        out("const double %s = _inputs[_inIdx_%s][_icell];", var, var)
    out("double %s = _inputs[_inIdx_%s][_icell];",stretch,stretch)
    good |= inputvars
    for var in order(partialvars):
        out("double %s = 0;",var)
    good |= partialvars
    for var in order(diffvars):
        out("double %s = _state[_old()][_icell].%s;", var, var)
    out("for (int itime=0; itime<internalTimestep; itime++)")
    out("{")
    out.inc()
    out("double %s = itime*(dt/internalTimestep);",partialStretch)
    good |= set([partialStretch])
    cprinter = CPrintVisitor(out, model.ssa, paramvars|flagvars)
    model.printTarget(good,set(diffvarUpdate.values()),cprinter)
    for var in order(diffvars|partialvars):
        out("%s += %s*(dt/internalTimestep);",var,diffvarUpdate[var])
    out("%s += %s*(dt/internalTimestep);",stretch,stretchVel)
    out.dec()
    out("}")
    for var in order(diffvars|partialvars):
        out("_state[_new()][_icell].%s = %s;",var,var)
    good |= diffvars
    out.dec(2)
    out(r'''
  }
}

void ThisModel::outputForNextTimestep(const double dt, const double* const _inputs[], double* const outputs[])
{
   for (int _icell=0; _icell < _numPoints; _icell++)
   {''', template)
    out.inc(2)
    for var in order(inputvars-set([stretch])):
        out("const double %s = _inputs[_inIdx_%s][_icell];", var, var)
    for var in order(diffvars|partialvars):
        out("const double %s = _state[_new()][_icell].%s;",var,var)
    out("const double %s = _inputs[_inIdx_%s][_icell]+dt*%s;",stretch,stretch,stretchVel)
    out("const double %s = dt;",partialStretch)
    cprinter = CPrintVisitor(out, model.ssa, paramvars|flagvars)
    model.printTarget(good,outputvars,cprinter)
    for var in order(outputvars):
        out("outputs[_outIdx_%s][_icell] = %s;",var,var)
    out.dec(2)
    out(r'''
   }
}

void ThisModel::commitTimestep()
{
   swapOldAndNew();
}

}
''', template)

generators = {
    frozenset(["cardioid","mech"]): generateCardioidMech,
}
