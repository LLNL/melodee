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
from sympy.printing.ccode import CCodePrinter
from sympy.core import S


from parser import MelodeeParser,Differentiator
import utility
from utility import order

def repeat(thing, repetitions):
    return (thing,) * repetitions

class MyCCodeSympyPrinter(CCodePrinter):
    def __init__(self,*args,**kwargs):
        CCodePrinter.__init__(self,*args,**kwargs)
    def _print_Pow(self, expr):
        PREC = sympy.printing.precedence.precedence(expr)
        if expr.exp == 0:
            return 1
        elif expr.exp == 0.5:
            return 'sqrt(%s)' % self._print(expr.base)
        elif expr.exp.is_constant and int(expr.exp) == expr.exp:
            result = self.parenthesize(expr.base,PREC)
            if expr.exp > 0:
                return "*".join(repeat(result, int(expr.exp)))
            else:
                return "1.0/"+"/".join(repeat(result, -int(expr.exp)))
        return 'pow(%s, %s)' % (self._print(expr.base),
                                self._print(expr.exp))
    def _print_Relational(self,expr):
        if expr.rel_op == "==" or expr.rel_op == "!=":
            PREC = sympy.printing.precedence.precedence(expr)
            return "%s %s %s" % (self.parenthesize(expr.lhs, PREC),
                                 expr.rel_op,
                                 self.parenthesize(expr.rhs, PREC))
        else:
            return super(MyCCodeSympyPrinter, self)._print_Relational(expr)

def pretty(symbol):
    return str(symbol)

class CPrintVisitor:
    def __init__(self, out, ssa, declared, decltype="double"):
        self.out = out
        self.ssa = ssa
        self.declaredStack = [set(declared)]
        self.decltype = decltype
        self.cprinter = MyCCodeSympyPrinter()
                      
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
        rhsText = self.cprinter.doprint(rhs.sympy, None)
        self.equationPrintWithRhs(lhs,rhsText)
    def equationPrintWithRhs(self,lhs,rhsText):
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
        self.cprinter = MyCCodeSympyPrinter()
    def ifPrint(self,printer,ifSymbol,thenList,elseList,choiceList):
        self.other.ifPrint(printer,ifSymbol,thenList,elseList,choiceList)
    def equationPrint(self,lhs,rhs):
        self.other.declaredStack[-1].add(pretty(lhs))
        rhsText = self.cprinter.doprint(rhs.sympy, None)
        if lhs in self.params:
            self.out("setDefault(%s, %s);", pretty(lhs),rhsText)
        else:
            self.other.equationPrint(lhs,rhs)

class InterpolatePrintVisitor:
    def __init__(self, otherPrintVisitor, interps):
        self.other = otherPrintVisitor
        self.interps = interps
    def ifPrint(self,printer,ifSymbol,thenList,elseList,choiceList):
        self.other.ifPrint(printer,ifSymbol,thenList,elseList,choiceList)
    def equationPrint(self,lhs,rhs):
        if lhs in self.interps:
            (interpVar, count) = self.interps[lhs]
            rhsText = "_interpolant[%d].eval(%s)" % (count, pretty(interpVar))
            self.other.equationPrintWithRhs(lhs,rhsText)
        else:
            self.other.equationPrint(lhs,rhs)

def generateCardioid(model, targetName, headerFile, sourceFile):
    template = {}
    template["target"] = targetName

    diffvars = model.diffvars()
    diffvarUpdate = {var : model.diffvarUpdate(var) for var in diffvars}
    params = model.varsWithAttribute("param")
    gates = model.varsWithAttribute("gate") & diffvars
    polyfits = model.varsWithAttribute("interp")
    
    V = model.input("V")
    V_init = model.output("V_init")
    Iion = model.output("Iion")

    differ = Differentiator(model, diffvars | params | set([V]))
    gateJacobians = {}
    for gate in gates:
        (gateJacobians[gate],dontcare) = differ.diff(diffvarUpdate[gate],gate)
    differ.augmentInstructions()

    dt = model.addSymbol("_dt")
    gateTargets = {}
    for gate in gates:
        F = model.ssa[diffvarUpdate[gate]].sympy
        L = model.ssa[gateJacobians[gate]].sympy
        M = (F-L*gate).simplify()
        
        RLA = model.addInstruction("_%s_RLA" % gate, sympy.exp(dt*L))
        RLB = model.addInstruction("_%s_RLB" % gate, M/L*(RLA-1))
        gateTargets[gate] = (RLA,RLB)

    expensiveVars = model.extractExpensiveFunctions()
    model.makeNamesUnique()

    computeTargets = set()
    computeTargets.add(Iion)
    computeTargets |= set([diffvarUpdate[var] for var in diffvars-gates])
    for gate in gates:
        (RLA,RLB) = gateTargets[gate]
        computeTargets.add(RLA)
        computeTargets.add(RLB)

    approxvars = set([dt])
    statevars = model.inputs()|diffvars
    computeAllDepend = model.allDependencies(approxvars|statevars, computeTargets)

    constants = model.allExcluding(approxvars, statevars) & computeAllDepend

    polyfitTargets = {}
    allfits = set()
    for fit in polyfits:
        good = approxvars | set([fit])
        dependsOnlyOnFit = model.allExcluding(good, statevars-good)
        polyfitCandidates = (dependsOnlyOnFit & computeAllDepend) - constants

        expensiveFit = expensiveVars & polyfitCandidates
        inexpensiveFit = model.allExcluding(good, (statevars-good)|expensiveFit)
        polyfitCandidates -= inexpensiveFit
        
        externallyUsedFits = (
            model.allDependencies(approxvars|statevars|polyfitCandidates, computeTargets)
            &
            polyfitCandidates
            )
        polyfitTargets[fit] = externallyUsedFits
        allfits |= externallyUsedFits

    fitCount = 0
    interps = {}
    for fit in order(polyfits):
        for target in order(polyfitTargets[fit]):
            interps[target] = (fit, fitCount)
            fitCount += 1

    computeAllDepend = model.allDependencies(
        approxvars|statevars|allfits,
        computeTargets)
    constants = model.allExcluding(approxvars,statevars) & computeAllDepend

    out = headerFile
    out('''
#include "Reaction.hh"
#include "Interpolation.hh"
#include "object.h"
#include <vector>
#include <sstream>

namespace scanReaction 
{
    Reaction* scan%(target)s(OBJECT* obj, const int numPoints, const double __dt);
}

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
      ThisReaction(const int numPoints, const double __dt);
      std::string methodName() const;
      
      void createInterpolants(const double _dt);
      void calc(double dt,
                const VectorDouble32& Vm,
                const std::vector<double>& iStim,
                VectorDouble32& dVm);
      //void updateNonGate(double dt, const VectorDouble32&Vm, VectorDouble32&dVR);
      //void updateGate   (double dt, const VectorDouble32&Vm) ;
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
      double __cachedDt;''', template)
    out.inc(2)
    out("Interpolation _interpolant[%d];" % fitCount)
    out.dec(2)
    out('''
      friend Reaction* scanReaction::scan%(target)s(OBJECT* obj, const int numPoints, const double __dt);
   };
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
#include "mpiUtils.h"
#include <cmath>
#include <cassert>
#include <fstream>
#include <iostream>

using namespace std;

namespace scanReaction 
{

#define setDefault(name, value) objectGet(obj, #name, reaction->name, #value)
   
   Reaction* scan%(target)s(OBJECT* obj, const int numPoints, const double _dt)
   {
      %(target)s::ThisReaction* reaction = new %(target)s::ThisReaction(numPoints, _dt);

      //override the defaults
      //EDIT_PARAMETERS''', template)
    out.inc(2)
    good = set([dt])
    cprinter = CPrintVisitor(out, model.ssa, params)
    paramPrinter = ParamPrintVisitor(out, cprinter, params)
    model.printTarget(good, params, paramPrinter)
    good |= params
    out.dec(2)
    out(r'''
      string fitName;
      objectGet(obj, "fit", fitName, "");
      int funcCount = sizeof(reaction->_interpolant)/sizeof(reaction->_interpolant[0]);
      if (fitName != "")
      {
         OBJECT* fitObj = objectFind(fitName, "FIT");
         double _fit_dt; objectGet(fitObj, "dt", _fit_dt, "nan");''' % template)
    out.inc(3)
    for var in order(params):
        out('double _fit_%s; objectGet(fitObj, "%s", _fit_%s, "nan");''',var,var,var)
    out.dec(3)
    out(r'''
         if (_dt == _fit_dt''' % template)
    out.inc(4)
    for var in order(params):
        out("&& _fit_%s == reaction->%s",var,var)
    out.dec(4)
    out(r'''
         )
         {
            vector<string> functions;
            objectGet(fitObj, "functions", functions);
            OBJECT* funcObj;
            if (functions.size() == funcCount)
            {
               for (int _ii=0; _ii<functions.size(); _ii++)
               {
                  OBJECT* funcObj = objectFind(functions[_ii], "FUNCTION");
                  objectGet(funcObj, "numer", reaction->_interpolant[_ii].numNumer_, "-1");
                  objectGet(funcObj, "denom", reaction->_interpolant[_ii].numDenom_, "-1");
                  objectGet(funcObj, "coeff", reaction->_interpolant[_ii].coeff_);
               }
               return reaction;
            }
         }
      }

      reaction->createInterpolants(_dt);

      //save the interpolants
      if (funcCount > 0 && getRank(0) == 0)
      {
         ofstream outfile((string(obj->name) +".fit.data").c_str());
         outfile.precision(16);
         fitName = string(obj->name) + "_fit";
         outfile << obj->name << " REACTION { fit=" << fitName << "; }\n";
         outfile << fitName << " FIT {\n";
         outfile << "   dt = " << _dt << ";\n";''' % template)
    out.inc(3)
    for var in order(params):
        out(r'outfile << "   %s = " << reaction->%s << ";\n";', var, var)
    out.dec(3)
    out(r'''
         outfile << "   functions = ";
         for (int _ii=0; _ii<funcCount; _ii++) {
            outfile << obj->name << "_interpFunc" << _ii << " ";
         }
         outfile << ";\n";
         outfile << "}\n";

         for (int _ii=0; _ii<funcCount; _ii++)
         {
            outfile << obj->name << "_interpFunc" << _ii << " FUNCTION { "
                    << "numer=" << reaction->_interpolant[_ii].numNumer_ << "; "
                    << "denom=" << reaction->_interpolant[_ii].numDenom_ << "; "
                    << "coeff=";
            for (int _jj=0; _jj<reaction->_interpolant[_ii].coeff_.size(); _jj++)
            {
               outfile << reaction->_interpolant[_ii].coeff_[_jj] << " ";
            }
            outfile << "; }\n";
         }
         outfile.close();
      }
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
   
ThisReaction::ThisReaction(const int numPoints, const double __dt)
: nCells_(numPoints)
{
   assertStateOrderAndVarNamesAgree();
   state_.resize(nCells_);
   perCellFlags_.resize(nCells_);
   perCellParameters_.resize(nCells_);
   __cachedDt = __dt;
}

void ThisReaction::createInterpolants(const double _dt) {
''', template)

    out.inc()
    fitMap = []
    for fit in range(0,len(interps)):
        fitMap.append(0)
    for target,(fit,fitCount) in interps.items():
        fitMap[fitCount] = (target,fit)
    for fitCount in range(0,len(interps)):
        target,fit = fitMap[fitCount]
        (lb,ub,inc) = model.info("interp",fit).split(",")
        lookup = {
            "target" : pretty(target),
            "ub" : ub,
            "lb" : lb,
            "inc" : inc,
            "fit" : pretty(fit),
            "fitCount" : fitCount,
        }
        out("{")
        out.inc()
        out("int _numPoints = (%(ub)s - %(lb)s)/%(inc)s;", lookup)
        out("vector<double> _inputs(_numPoints);")
        out("vector<double> _outputs(_numPoints);")
        out("for (int _ii=0; _ii<_numPoints; _ii++)", lookup)
        out("{")
        out.inc()
        out("double %(fit)s = %(lb)s + (%(ub)s - %(lb)s)*(_ii+0.5)/_numPoints;", lookup)
        out("_inputs[_ii] = %(fit)s;", lookup)
        cprinter = CPrintVisitor(out, model.ssa, params)
        model.printTarget(good|set([fit]),set([target]),cprinter)
        out("_outputs[_ii] = %(target)s;", lookup)
        out.dec()
        out("}")
        out(r'''
double relError = 1e-4;
        double actualTolerance = _interpolant[%(fitCount)d].create(_inputs,_outputs, relError);
if (actualTolerance > relError  && getRank(0) == 0)
{
   cerr << "Warning: Could not meet tolerance for %(target)s: " 
        << actualTolerance " > " << relError
        << " target" << endl;
}''', lookup)
        out.dec()
        out("}")
    out.dec()
    out('''
}

void ThisReaction::calc(double _dt, const VectorDouble32& __Vm,
                       const vector<double>& __iStim , VectorDouble32& __dVm)
{
   //define the constants''', template)
    out.inc()

    cprinter = CPrintVisitor(out, model.ssa, params)
    model.printTarget(good,computeAllDepend&constants,cprinter)

    good |= constants
    
    out.dec()
    out('''
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

    iprinter = InterpolatePrintVisitor(cprinter, interps)
    model.printSet(model.allDependencies(good|allfits,computeTargets)-good, iprinter)
    
    out.dec(2)
    out(r'''
      
      //EDIT_STATE''', template)
    out.inc(2)
    for var in order(diffvars-gates):
        out('state_[__ii].%s += _dt*%s;',pretty(var),pretty(diffvarUpdate[var]))
    for var in order(gates):
        (RLA,RLB) = gateTargets[var]
        out('state_[__ii].%(v)s = %(a)s*%(v)s + %(b)s;',
            v=pretty(var),
            a=pretty(RLA),
            b=pretty(RLB),
        )
                     
    out.dec(2)
    out('''      
      __dVm[__ii] = -Iion;
   }
}

void ThisReaction::initializeMembraneVoltage(VectorDouble32& __Vm)
{
   assert(__Vm.size() >= nCells_);

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
