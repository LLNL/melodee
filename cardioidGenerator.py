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
from sympy.printing.ccode import ccode

from parser import MelodeeParser,Differentiator
import utility
from utility import order

def pretty(symbol):
    return str(symbol)

class CPrintVisitor:
    def __init__(self, out, ssa, declared, decltype="double"):
        self.out = out
        self.ssa = ssa
        self.declaredStack = [set(declared)]
        self.decltype = decltype
    def pushStack(self):
        self.declaredStack.append(set(self.declaredStack[-1]))
    def popStack(self):
        return self.declaredStack.pop()
    def ifPrint(self,printer,ifSymbol,thenList,elseList,choiceList):
        for choice in choiceList:
            if pretty(choice) not in self.declaredStack[-1]:
                self.out("%s %s;",self.decltype,choice)
            self.declaredStack[-1].add(pretty(choice))
        self.out("if (%s)",pretty(ifSymbol))
        self.out("{")
        self.out.inc()
        self.pushStack()
        printer(thenList)
        for choiceVar in choiceList:
            choice = self.ssa[choiceVar]
            lhs = pretty(choiceVar)
            rhs = pretty(choice.thenVar)
            if lhs != rhs:
                self.out("%s = %s;",lhs,rhs)
        self.popStack()
        self.out.dec()
        self.out("}")
        self.out("else")
        self.out("{")
        self.out.inc()
        self.pushStack()
        printer(elseList)
        for choiceVar in choiceList:
            choice = self.ssa[choiceVar]
            lhs = pretty(choiceVar)
            rhs = pretty(choice.elseVar)
            if lhs != rhs:
                self.out("%s = %s;",lhs,rhs)
        self.popStack()
        self.out.dec()
        self.out("}")
    def equationPrint(self,lhs,rhs):
        rhsText = ccode(rhs.sympy)
        if pretty(lhs) not in self.declaredStack[-1]:
            self.out("%s %s = %s;",self.decltype,pretty(lhs),rhsText)
            self.declaredStack[-1].add(pretty(lhs))
        else:
            self.out("%s = %s;",pretty(lhs),rhsText)

class ParamPrintVisitor:
    def __init__(self, out, otherPrintVisitor, params):
        self.out  = out
        self.other = otherPrintVisitor
        self.params = params
    def ifPrint(self,printer,ifSymbol,thenList,elseList,choiceList):
        self.other.ifPrint(printer,ifSymbol,thenList,elseList,choiceList)
    def equationPrint(self,lhs,rhs):
        self.other.declaredStack[-1].add(pretty(lhs))
        rhsText = ccode(rhs.sympy)
        if lhs in self.params:
            self.out("setDefault(%s, %s);", pretty(lhs),rhsText)
        else:
            self.other.equationPrint(lhs,rhs)

def generateCardioid(model, targetName, headerFile, sourceFile):
    template = {}
    template["target"] = targetName

    params = model.varsWithAttribute("param")
    diffvars = model.diffvars()
    gates = model.varsWithAttribute("gate") & diffvars
    diffvarUpdate = {var : model.diffvarUpdate(var) for var in diffvars}
    #gatevars = model.getVars("gate")
    #tracevars = model.getVars("trace")
    
    V = model.input("V")
    V_init = model.output("V_init")
    Iion = model.output("Iion")

    differ = Differentiator(model, diffvars | params | set([V]))
    gateJacobians = {}
    for gate in gates:
        (gateJacobians[gate],dontcare) = differ.diff(diffvarUpdate[gate],gate)
    differ.augmentInstructions()
    
    out = headerFile
    out('''
#include "Reaction.hh"
#include "object.h"
#include <vector>

namespace %(target)s
{

   struct State
   {''', template)

    out.inc(2)
    for var in order(diffvars):
        out("double %s;", pretty(var))
    out.dec(2)
    out('''
   };

   struct PerCellFlags
   { };
   struct PerCellParameters
   { };

   class ThisReaction : public Reaction
   {
    public:
      ThisReaction(const int numPoints);
      std::string methodName() const;
      
      void calc(double dt,
                const VectorDouble32& Vm,
                const std::vector<double>& iStim,
                VectorDouble32& dVm);
      void initializeMembraneVoltage(VectorDouble32& Vm);
      virtual void getCheckpointInfo(std::vector<std::string>& fieldNames,
                                     std::vector<std::string>& fieldUnits) const;
      virtual int getVarHandle(const std::string& varName) const;
      virtual void setValue(int iCell, int varHandle, double value);
      virtual double getValue(int iCell, int varHandle) const;
      virtual const std::string getUnit(const std::string& varName) const;

    public:
      //PARAMETERS''', template)
    out.inc(2)
    for var in order(params):
        out("double %s;", pretty(var))
    out.dec(2)
    out('''
      //per-cell flags
      std::vector<PerCellFlags> perCellFlags_;
      std::vector<PerCellParameters> perCellParameters_;

    private:
      unsigned nCells_;
      std::vector<State> state_;
   };
}

namespace scanReaction 
{
   Reaction* scan%(target)s(OBJECT* obj, const int numPoints);
}
''', template)

    out = sourceFile
    
    out('''
/**

   How to convert this code to work for any other model:

   - Search/Replace the model name with your own specific string in the header and source files
   - Add your own code to EDIT_FLAGS and EDIT_PARAMETERS
   - Add your own code to EDIT_PERCELL_FLAGS and EDIT_PERCELL_PARAMETERS
   - Add your own states to EDIT_STATE
   - Add your computation code to the main calc routine, copy pasting frmo matlab.
   
 */


#include "%(target)s.hh"
#include "object_cc.hh"
#include <cmath>

using namespace std;

#define Eq(x,y) ((x)==(y))


namespace scanReaction 
{

#define setDefault(name, value) objectGet(obj, #name, reaction->name, #value)
   
   Reaction* scan%(target)s(OBJECT* obj, const int numPoints)
   {
      %(target)s::ThisReaction* reaction = new %(target)s::ThisReaction(numPoints);

      //override the defaults
      //EDIT_PARAMETERS''', template)
    out.inc(2)
    good = set()
    cprinter = CPrintVisitor(out, model.ssa, params)
    paramPrinter = ParamPrintVisitor(out, cprinter, params)
    model.printTarget(good, params, paramPrinter)
    good |= params
    out.dec(2)
    out('''
      
      return reaction;
   }
#undef setDefault

}

namespace %(target)s 
{
   
string ThisReaction::methodName() const
{
   return "%(target)s";
}
const char* varNames[] = 
{''', template)
    out.inc()
    varNames = ['"'+pretty(var)+'"' for var in order(diffvars)]
    out(",\n".join(varNames))
    out.dec()
    out('''
};
const char* varUnits[] =
{''', template)
    out.inc()
    #varUnits = ['"'+units[var]+'"' for var in order(diffvars)]
    varUnits = ['"%s"'%str(model.ssa[var].astUnit.rawUnit) for var in order(diffvars)]
    out(",\n".join(varUnits))
    out.dec()
    out('''
};

#define NUMVARS (sizeof(varNames)/sizeof(char*))

int getVarOffset(const std::string& varName)
{
   for (int ivar=0; ivar<NUMVARS; ivar++) 
   {
      if (varNames[ivar] == varName) 
      {
         return ivar;
      }
   }
   assert(0 && "Control should never get here.");
   return -1;
}

void assertStateOrderAndVarNamesAgree(void)
{
   State s;
#define checkVarOrder(x) assert(reinterpret_cast<double*>(&s)+getVarOffset(#x) == &s . x)

   int STATIC_ASSERT_checkAllDouble[(NUMVARS == sizeof(s)/sizeof(double))? 1: 0];

   //EDIT_STATE''', template)
    out.inc()
    for var in order(diffvars):
        out('checkVarOrder(%s);',pretty(var))
    out.dec()
    out('''
}
   
ThisReaction::ThisReaction(const int numPoints)
: nCells_(numPoints)
{
   assertStateOrderAndVarNamesAgree();
   state_.resize(nCells_);
   perCellFlags_.resize(nCells_);
   perCellParameters_.resize(nCells_);
}

void ThisReaction::calc(double dt, const VectorDouble32& __Vm,
                       const vector<double>& __iStim , VectorDouble32& __dVm)
{
   for (unsigned __ii=0; __ii<nCells_; ++__ii)
   {
      //set Vm
      const double V = __Vm[__ii];
      const double iStim = __iStim[__ii];

      //set all state variables''', template)
    out.inc(2)
    for var in order(diffvars):
        out('const double %s=state_[__ii].%s;',pretty(var),pretty(var))
    good |= diffvars
    good.add(V)

    target = set(diffvarUpdate.values()) | set(gateJacobians.values())
    target.add(Iion)
    cprinter = CPrintVisitor(out, model.ssa, params)
    model.printTarget(good,target,cprinter)
    
    out.dec(2)
    out(r'''


      if (1) 
      {
         bool foundError=false;
#define CHECK_BLOWUP(x) do { if (!isfinite(x)) { fprintf(stderr, "Error in node %%d, variable " #x " = %%g\n", ii, (x)); foundError=true; } } while(0)
         CHECK_BLOWUP(Iion);
            
         //EDIT_STATE''', template)
    out.inc(3)
    for var in order(diffvars):
        out('CHECK_BLOWUP(%s);',pretty(diffvarUpdate[var]))
    out.dec(3)
    out(r'''

#undef CHECK_BLOWUP
         
         if (foundError) 
         {
#define PRINT_STATE(x) do { fprintf(stderr, "node %%d: " #x " = %%g\n", ii, (x)); } while(0)
            //EDIT_STATE''', template)
    out.inc(4)
    for var in order(diffvars):
        out('PRINT_STATE(%s);',pretty(var))
    out.dec(4)
    out('''            
#undef PRINT_STATE
            
            exit(255);
         }
      }
      
      
      //EDIT_STATE''', template)
    out.inc(2)
    for var in order(diffvars-gates):
        out('state_[__ii].%s += dt*%s;',pretty(var),pretty(diffvarUpdate[var]))
    for var in order(gates):
        out('state_[__ii].%(v)s += %(d)s/%(j)s*expm1(dt*%(j)s);',
            v=pretty(var),
            d=pretty(diffvarUpdate[var]),
            j=gateJacobians[var],
        )
                     
    out.dec(2)
    out('''      
   }
}

void ThisReaction::initializeMembraneVoltage(VectorDouble32& __Vm)
{
   assert(Vm.size() >= nCells_);

''', template)
    out.inc()

    cprinter = CPrintVisitor(out, model.ssa, params)
    good = set()
    good |= params
    target = set([V_init])
    model.printTarget(good,target,cprinter)

    good.add(V)
    out("double V = V_init;") #FIXME, possibly make more general?
    model.printTarget(good,diffvars,cprinter)
    
    out.dec()
    out('''      

   State __initState;
''', template)
    out.inc()
    for var in order(diffvars):
        out('__initState.%s = %s;',pretty(var),pretty(var))
    out.dec()
    out('''      
   

   __Vm.assign(__Vm.size(), V_init);

   state_.resize(nCells_);
   state_.assign(state_.size(), __initState);
}

const string ThisReaction::getUnit(const std::string& varName) const
{
   return varUnits[getVarOffset(varName)];
}

#define HANDLE_OFFSET 1000
int ThisReaction::getVarHandle(const std::string& varName) const
{
   return getVarOffset(varName)+HANDLE_OFFSET;
}

void ThisReaction::setValue(int iCell, int varHandle, double value) 
{
   reinterpret_cast<double*>(&state_[iCell])[varHandle-HANDLE_OFFSET] = value;
}


double ThisReaction::getValue(int iCell, int varHandle) const
{
   return reinterpret_cast<const double*>(&state_[iCell])[varHandle-HANDLE_OFFSET];
}

void ThisReaction::getCheckpointInfo(vector<string>& fieldNames,
                                     vector<string>& fieldUnits) const
{
   fieldNames.resize(NUMVARS);
   fieldUnits.resize(NUMVARS);
   for (int ivar=0; ivar<NUMVARS; ivar++) 
   {
      fieldNames[ivar] = varNames[ivar];
      fieldUnits[ivar] = varUnits[ivar];
   }
}

}''', template)

if __name__=="__main__":
    import sys
    
    p = MelodeeParser()

    sys.argv.pop(0)
    target = sys.argv.pop(0)
    models = {}
    for filename in sys.argv:
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
        
    generateCardioid(p.getModel(target), target,
                     #utility.Indenter(sys.stdout),
                     #utility.Indenter(sys.stdout),
                     utility.Indenter(open(target+".hh","w")),
                     utility.Indenter(open(target+".cc","w")),
    )
