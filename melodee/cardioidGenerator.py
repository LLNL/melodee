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
import io
from sympy.printing.ccode import C99CodePrinter
from sympy.core import S


from melodee.parser import MelodeeParser,Differentiator
from melodee import utility
from melodee.utility import order

def repeat(thing, repetitions):
    return (thing,) * repetitions

class MyCCodeSympyPrinter(C99CodePrinter):
    def __init__(self,*args,**kwargs):
        C99CodePrinter.__init__(self,*args,**kwargs)
    def _print_Pow(self, expr):
        PREC = sympy.printing.precedence.precedence(expr)
        if expr.exp == 0:
            return "1.0"
        elif expr.exp == 0.5:
            return 'sqrt(%s)' % self._print(expr.base)
        elif expr.exp == 1:
            return self.parenthesize(expr.base,PREC)
        elif expr.exp.is_constant and int(expr.exp) == expr.exp:
            base = self.parenthesize(expr.base,PREC)
            if expr.exp > 0:
                return "(" + "*".join(repeat(base,int(expr.exp))) + ")"
            else:
                return "(1.0/" + "/".join(repeat(base,-int(expr.exp))) + ")"
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

class CPrintVisitor(object):
    def __init__(self, out, ssa, declared, decltype="double"):
        self.out = out
        self.ssa = ssa
        self.declaredStack = [set([str(x) for x in declared])]
        self.decltype = decltype
        self.cprinter = MyCCodeSympyPrinter()

    def currentlyDeclared(self):
        return self.declaredStack[-1]
    def pushStack(self):
        self.declaredStack.append(set(self.currentlyDeclared()))
    def popStack(self):
        return self.declaredStack.pop()
    def ifPrint(self,printer,ifSymbol,thenList,elseList,choiceList):
        for choice in choiceList:
            if str(choice) not in self.currentlyDeclared():
                self.out("%s %s;",self.decltype,choice)
            self.currentlyDeclared().add(str(choice))
        if ifSymbol != None:
            self.out("if (%s)",ifSymbol)
            self.out("{")
            self.out.inc()
            self.pushStack()
            printer(thenList)
            for choiceVar in choiceList:
                choice = self.ssa[choiceVar]
                lhs = str(choiceVar)
                rhs = str(choice.thenVar)
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
                lhs = str(choiceVar)
                rhs = str(choice.elseVar)
                if lhs != rhs:
                    self.out("%s = %s;",lhs,rhs)
            self.popStack()
            self.out.dec()
            self.out("}")
        for choice in choiceList:
            self.printChoice(choice)
    def equationPrint(self,lhs,rhs):
        rhsText = self.cprinter.doprint(rhs.sympy, None)
        self.equationPrintWithRhs(lhs,rhsText)
    def equationPrintWithRhs(self,lhs,rhsText):
        if str(lhs) not in self.currentlyDeclared():
            self.out("%s %s = %s;",self.decltype,lhs,rhsText)
            self.currentlyDeclared().add(str(lhs))
        else:
            self.out("%s = %s;",lhs,rhsText)
    def printChoice(self,lhs):
        pass

class ParamPrintVisitor(CPrintVisitor):
    def __init__(self, out, ssa, declared, params, decltype="double"):
        self.params = params
        super(ParamPrintVisitor,self).__init__(out, ssa, declared, decltype)
    def equationPrint(self,lhs,rhs):
        rhsText = self.cprinter.doprint(rhs.sympy, None)
        if lhs in self.params:
            self.currentlyDeclared().add(str(lhs))
            self.out("setDefault(%s, %s);", lhs,rhsText)
        else:
            super(ParamPrintVisitor,self).equationPrint(lhs,rhs)
    def printChoice(self, lhs):
        if lhs in self.params:
            self.out("setDefault(%s, %s);", lhs,lhs)
class InterpolatePrintCPUVisitor(CPrintVisitor):
    def __init__(self, out, ssa, declared, interps, decltype="double"):
        self.interps = interps
        super(InterpolatePrintCPUVisitor, self).__init__(out,ssa,declared,decltype)
    def equationPrint(self,lhs,rhs):
        if lhs in self.interps:
            (interpVar, count) = self.interps[lhs]
            rhsText = "_interpolant[%d].eval(%s)" % (count, interpVar)
            self.equationPrintWithRhs(lhs,rhsText)
        else:
            super(InterpolatePrintCPUVisitor,self).equationPrint(lhs,rhs)
    def printChoice(self, lhs):
        if lhs in self.interps:
            self.equationPrint(lhs,None)

class InterpolatePrintNvidiaVisitor(CPrintVisitor):
    def __init__(self, out, ssa, declared, interps, decltype="double"):
        self.interps = interps
        super(InterpolatePrintNvidiaVisitor, self).__init__(out,ssa,declared,decltype)
    def equationPrint(self,lhs,rhs):
        if lhs in self.interps:
            (interpVar, count) = self.interps[lhs]
            self.out('"; generateInterpString(ss,_interpolant[%d], "%s"); ss << "',count,interpVar)
            self.equationPrintWithRhs(lhs,"_ratPoly")
        else:
            super(InterpolatePrintNvidiaVisitor,self).equationPrint(lhs,rhs)
    def printChoice(self, lhs):
        if lhs in self.interps:
            self.equationPrint(lhs,None)

def generateCardioid(model, targetName, arch="cpu"):
    template = {}
    template["target"] = targetName

    diffvars = model.diffvars()
    diffvarUpdate = {var : model.diffvarUpdate(var) for var in diffvars}
    params = model.varsWithAttribute("param")
    markovs = model.varsWithAttribute("markov") & diffvars
    gates = model.varsWithAttribute("gate") & diffvars
    polyfits = model.varsWithAttribute("interp")
    tracevars = model.varsWithAttribute("trace")
    nointerps = model.varsWithAttribute("nointerp")
    
    V = model.input("V")
    V_init = model.output("V_init")
    Iion = model.output("Iion")

    differ = Differentiator(model, diffvars | params | set([V]))
    gateJacobians = {}
    for gate in order(gates):
        (gateJacobians[gate],dontcare) = differ.diff(diffvarUpdate[gate],gate)
    markovJacobians = { var : {} for var in markovs}
    for imarkov in order(markovs):
        for jmarkov in order(markovs):
            (dontcare,markovJacobians[imarkov][jmarkov]) = differ.diff(diffvarUpdate[imarkov],jmarkov)
    differ.augmentInstructions()

    dt = model.addSymbol("_dt")
    gateTargets = {}
    for gate in order(gates):
        F = model.ssa[diffvarUpdate[gate]].sympy
        L = model.ssa[gateJacobians[gate]].sympy
        M = (F-L*gate).simplify()
        
        RLA = model.addInstruction("_%s_RLA" % gate, sympy.exp(dt*L)-1)
        RLB = model.addInstruction("_%s_RLB" % gate, M/L)
        gateTargets[gate] = (RLA,RLB)

    markovOld = {}
    for markov in order(markovs):
        markovOld[markov] = model.addSymbol("_mi_old_%s" % markov)
    markovTargets = {}
    for imarkov in order(markovs):
        summation = 0
        for jmarkov in order(markovs):
            if imarkov == jmarkov:
                continue
            if jmarkov in markovTargets:
                thisSym = markovTargets[jmarkov]
            else:
                thisSym = markovOld[jmarkov]
            summation += markovJacobians[imarkov][jmarkov].sympy*thisSym
        sss = (diffvarUpdate[imarkov]+dt*summation)/(1-dt*markovJacobians[imarkov][imarkov].sympy)
        markovTargets[imarkov] = model.addInstruction("_mi_new_%s" % imarkov,sss)

    expensiveVars = model.extractExpensiveFunctions()
    model.makeNamesUnique()

    computeTargets = set()
    computeTargets.add(Iion)
    computeTargets |= set([diffvarUpdate[var] for var in diffvars-gates-markovs])
    for gate in gates:
        (RLA,RLB) = gateTargets[gate]
        computeTargets.add(RLA)
        computeTargets.add(RLB)
    for markov in markovs:
        computeTargets.add(markovTargets[markov])

    approxvars = set([dt])
    statevars = model.inputs()|diffvars
    computeAllDepend = model.allDependencies(approxvars|statevars|params, computeTargets)

    constants = model.allExcluding(approxvars, statevars) & computeAllDepend

    polyfitTargets = {}
    allfits = set()
    for fit in polyfits:
        good = approxvars | params | set([fit])
        dependsOnlyOnFit = model.allExcluding(good, (statevars-good)|nointerps)
        polyfitCandidates = (dependsOnlyOnFit & computeAllDepend) - constants

        expensiveFit = expensiveVars & polyfitCandidates
        inexpensiveFit = model.allExcluding(good, (statevars-good)|expensiveFit)
        polyfitCandidates -= inexpensiveFit
        
        externallyUsedFits = (
            model.allDependencies(approxvars|statevars|params|polyfitCandidates, computeTargets)
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
        approxvars|statevars|params|allfits,
        computeTargets)
    constants = model.allExcluding(approxvars,statevars) & computeAllDepend

    out = utility.Indenter(open(targetName+".hh","w"))
    out('''
#include "Reaction.hh"
#include "Interpolation.hh"
#include "object.h"
#include <vector>
#include <sstream>
''', template)
    if arch=="nvidia":
        out(r'''
#include "TransportCoordinator.hh"
#include <nvrtc.h>
#include <cuda.h>
''', template)
    elif arch=="cpu":
        pass
    else:
        assert(False)
    out(r'''    
namespace scanReaction 
{
    Reaction* scan%(target)s(OBJECT* obj, const int numPoints, const double __dt);
}

namespace %(target)s
{

''', template)
    if arch=="nvidia":
        pass
    elif arch=="cpu":
        out.inc()
        out(r'struct State')
        out("{")
        out.inc()
        for var in order(diffvars):
            out("double %s;", var)
        out.dec()
        out("};")
        out.dec()
    else:
        assert(False)
    out(r'''

   class ThisReaction : public Reaction
   {
    public:
      ThisReaction(const int numPoints, const double __dt);
      std::string methodName() const;
      
      void createInterpolants(const double _dt);
      //void updateNonGate(double dt, const VectorDouble32&Vm, VectorDouble32&dVR);
      //void updateGate   (double dt, const VectorDouble32&Vm) ;
      virtual void getCheckpointInfo(std::vector<std::string>& fieldNames,
                                     std::vector<std::string>& fieldUnits) const;
      virtual int getVarHandle(const std::string& varName) const;
      virtual void setValue(int iCell, int varHandle, double value);
      virtual double getValue(int iCell, int varHandle) const;
      virtual double getValue(int iCell, int varHandle, double V) const;
      virtual const std::string getUnit(const std::string& varName) const;

    private:
      unsigned nCells_;
      double __cachedDt;

    public:
      //PARAMETERS''', template)
    out.inc(2)
    for var in order(params):
        out("double %s;", var)
    out.dec(2)
    out('''
''',template)
    if arch=="nvidia":
        out(r'''
    public:
      void calc(double dt,
                const Managed<ArrayView<double>> Vm_m,
                const Managed<ArrayView<double>> iStim_m,
                Managed<ArrayView<double>> dVm_m);
      void initializeMembraneVoltage(ArrayView<double> Vm);
      virtual ~ThisReaction();
      void constructKernel();
    private:
      TransportCoordinator<PinnedVector<double> > stateTransport_;
      std::string _program_code;
      nvrtcProgram _program;
      std::vector<char> _ptx;
      CUmodule _module;
      CUfunction _kernel;
      int blockSize_;
''', template)
    elif arch=="cpu":
        out(r'''
    public:
      void calc(double dt,
                const VectorDouble32& Vm,
                const std::vector<double>& iStim,
                VectorDouble32& dVm);
      void initializeMembraneVoltage(VectorDouble32& Vm);
    private:
      std::vector<State> state_;
''',template)
    else:
        assert(False)
    out.inc(2)
    out("Interpolation _interpolant[%d];" % fitCount)
    out.dec(2)
    out('''
      friend Reaction* scanReaction::scan%(target)s(OBJECT* obj, const int numPoints, const double __dt);
   };
}

''', template)

    out = utility.Indenter(open(targetName+".cc","w"))
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
    for var in order(params):
        out("double %s;", var)
    good = set([dt])
    paramPrinter = ParamPrintVisitor(out, model.ssa, params, params)
    model.printTarget(good, params, paramPrinter)
    good |= params

    for var in order(params):
        out("reaction->%s = %s;",var,var)
    
    out.dec(2)

    fitDependencies = model.allDependencies(good|polyfits,allfits) & good
    out(r'''
      bool reusingInterpolants = false;
      string fitName;
      objectGet(obj, "fit", fitName, "");
      int funcCount = sizeof(reaction->_interpolant)/sizeof(reaction->_interpolant[0]);
      if (fitName != "")
      {
         OBJECT* fitObj = objectFind(fitName, "FIT");''' % template)
    out.inc(3)
    if dt in fitDependencies:
        out('double _fit_dt; objectGet(fitObj, "dt", _fit_dt, "nan");')
    for var in order(fitDependencies-set([dt])):
        out('double _fit_%s; objectGet(fitObj, "%s", _fit_%s, "nan");''',var,var,var)
    out.dec(3)
    out(r'''
         if (1''' % template)
    out.inc(4)
    if dt in fitDependencies:
        out("&& _fit_dt == _dt")
    for var in order(fitDependencies-set([dt])):
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
               reusingInterpolants = true;
            }
         }
      }

      if (!reusingInterpolants)
      {
      reaction->createInterpolants(_dt);

      //save the interpolants
      if (funcCount > 0 && getRank(0) == 0)
      {
         ofstream outfile((string(obj->name) +".fit.data").c_str());
         outfile.precision(16);
         fitName = string(obj->name) + "_fit";
         outfile << obj->name << " REACTION { fit=" << fitName << "; }\n";
         outfile << fitName << " FIT {\n";''' % template)

    out.inc(4)
    if dt in fitDependencies:
        out(r'outfile << "   dt = " << _dt << ";\n";')
    for var in order(fitDependencies-set([dt])):
        out(r'outfile << "   %s = " << reaction->%s << ";\n";', var, var)
    out.dec(4)
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
      }''',template)
    if arch=="nvidia":
        out.inc(2)
        out("reaction->constructKernel();")
        out.dec(2)
    elif arch=="cpu":
        pass
    else:
        assert(False)
    out('''
      return reaction;
   }
#undef setDefault

}

namespace %(target)s 
{

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
            "target" : target,
            "ub" : ub,
            "lb" : lb,
            "inc" : inc,
            "fit" : fit,
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
        << actualTolerance << " > " << relError
        << " target" << endl;
}''', lookup)
        out.dec()
        out("}")
    out.dec()
    out('''
}''',template)
    if arch=="nvidia":
        out(r'''
void generateInterpString(stringstream& ss, const Interpolation& interp, const char* interpVar)
{
   ss <<
   "{\n"
   "   const double _numerCoeff[]={";
    for (int _ii=interp.numNumer_-1; _ii>=0; _ii--)
   {
      if (_ii != interp.numNumer_-1) { ss << ", "; }
      ss << interp.coeff_[_ii];
   }
   ss<< "};\n"
   "   const double _denomCoeff[]={";
   for (int _ii=interp.numDenom_+interp.numNumer_-2; _ii>=interp.numNumer_; _ii--)
   {
      ss << interp.coeff_[_ii] << ", ";
   }
   ss<< "1};\n"
   "   double _inVal = " << interpVar << ";\n"
   "   double _numerator=_numerCoeff[0];\n"
   "   for (int _jj=1; _jj<sizeof(_numerCoeff)/sizeof(_numerCoeff[0]); _jj++)\n"
   "   {\n"
   "      _numerator = _numerCoeff[_jj] + _inVal*_numerator;\n"
   "   }\n"
   "   if (sizeof(_denomCoeff)/sizeof(_denomCoeff[0]) == 1)\n"
   "   {\n"
   "      _ratPoly = _numerator;\n"
   "   }\n"
   "   else\n"
   "   {\n"
   "      double _denominator=_denomCoeff[0];\n"
   "      for (int _jj=1; _jj<sizeof(_denomCoeff)/sizeof(_denomCoeff[0]); _jj++)\n"
   "      {\n"
   "         _denominator = _denomCoeff[_jj] + _inVal*_denominator;\n"
   "      }\n"
   "      _ratPoly = _numerator/_denominator;\n"
   "   }\n"
   "}"
      ;
}

void ThisReaction::constructKernel()
{

   stringstream ss;
   ss.precision(16);
   ss <<
   "enum StateOffset {\n"
''',template)
        out.inc()
        for var in order(diffvars):
            out(r'"   %s_off,\n"',var)
        out.dec()
        out(r'''
   "   NUMSTATES\n"
   "};\n"
   "extern \"C\"\n"
   "__global__ void %(target)s_kernel(const double* _Vm, const double* _iStim, double* _dVm, double* _state) {\n"
   "const double _dt = " << __cachedDt << ";\n"
   "const int _nCells = " << nCells_ << ";\n"
''',template)


        out.inc(1)
        for var in order(params):
            out(r'"const double %s = " << %s << ";\n"',var,var)
        out.dec(1)
        out(r'''
   "const int _ii = threadIdx.x + blockIdx.x*blockDim.x;\n"
   "if (_ii >= _nCells) { return; }\n"
   "const double V = _Vm[_ii];\n"
   "double _ratPoly;\n"
''',template)
        out.inc(1)
        for var in order(diffvars):
            out(r'"const double %s = _state[_ii+%s_off*_nCells];\n"',var,var)

        good |= diffvars
        good.add(V)

        calcCodeBuf = io.BytesIO()
        calcOut = utility.Indenter(calcCodeBuf)
        iprinter = InterpolatePrintNvidiaVisitor(calcOut, model.ssa, params, interps)
        model.printSet(model.allDependencies(good|allfits,computeTargets)-good, iprinter)
        for line in calcCodeBuf.getvalue().splitlines():
            out('"'+line+r'\n"')

        out(r'"\n\n//EDIT STATE\n"')
        for var in order(diffvars-gates):
            out(r'"_state[_ii+%s_off*_nCells] += _dt*%s;\n"',var,diffvarUpdate[var])
        for var in order(gates):
            (RLA,RLB) = gateTargets[var]
            out(r'"_state[_ii+%(v)s_off*_nCells] += %(a)s*(%(v)s+%(b)s);\n"',
            v=var,
            a=RLA,
            b=RLB,
        )
        out(r'"_dVm[_ii] = -%s;\n"', Iion)
        out.dec()
        out(r'''
   "}\n";

   _program_code = ss.str();
   //cout << ss.str();
   nvrtcCreateProgram(&_program,
                      _program_code.c_str(),
                      "%(target)s_program",
                      0,
                      NULL,
                      NULL);
   nvrtcCompileProgram(_program,
                       0,
                       NULL);
   std::size_t size;
   nvrtcGetPTXSize(_program, &size);
   _ptx.resize(size);
   nvrtcGetPTX(_program, &_ptx[0]);

   cuModuleLoadDataEx(&_module, &_ptx[0], 0, 0, 0);
   cuModuleGetFunction(&_kernel, _module, "%(target)s_kernel");
}

void ThisReaction::calc(double dt,
                const Managed<ArrayView<double>> Vm_m,
                const Managed<ArrayView<double>> iStim_m,
                Managed<ArrayView<double>> dVm_m)
{
   ArrayView<double> state = stateTransport_.modifyOnDevice();

   {
      int errorCode=-1;
      if (blockSize_ == -1) { blockSize_ = 1024; }
      while(1)
      {
         ConstArrayView<double> Vm = Vm_m.readOnDevice();
         ConstArrayView<double> iStim = iStim_m.readOnDevice();
         ArrayView<double> dVm = dVm_m.modifyOnDevice();
         double* VmRaw = const_cast<double*>(&Vm[0]);
         double* iStimRaw = const_cast<double*>(&iStim[0]);
         double* dVmRaw = &dVm[0];
         double* stateRaw= &state[0];
         void* args[] = { &VmRaw,
                          &iStimRaw,
                          &dVmRaw,
                          &stateRaw};
         int errorCode = cuLaunchKernel(_kernel,
                                        (nCells_+blockSize_-1)/blockSize_, 1, 1,
                                        blockSize_,1,1,
                                        0, NULL,
                                        args, 0);
         if (errorCode == CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES && blockSize_ > 0)
         {
            blockSize_ /= 2;
            continue;
         }
         else if (errorCode != CUDA_SUCCESS)
         {
            printf("Cuda return %%d;\n", errorCode);
            assert(0 && "Launch failed");
         }
         else
         {
            break;
         }
         //cuCtxSynchronize();
      }
   }
}

enum StateOffset {''', template)
        out.inc()
        for var in order(diffvars):
            out("_%s_off,",var)
        out.dec()
        out('''
   NUMSTATES
};

ThisReaction::ThisReaction(const int numPoints, const double __dt)
: nCells_(numPoints)
{
   stateTransport_.setup(PinnedVector<double>(nCells_*NUMSTATES));
   __cachedDt = __dt;
   blockSize_ = -1;
   _program = NULL;
   _module = NULL;
}

ThisReaction::~ThisReaction() {
   if (_program)
   {
      nvrtcDestroyProgram(&_program);
      _program = NULL;
   }
   if (_module)
   {
      cuModuleUnload(_module);
      _module = NULL;
   }
}

''',template)
    elif arch=="cpu":
        out('''
ThisReaction::ThisReaction(const int numPoints, const double __dt)
: nCells_(numPoints)
{
   state_.resize(nCells_);
   __cachedDt = __dt;
}

void ThisReaction::calc(double _dt, const VectorDouble32& __Vm,
                       const vector<double>& __iStim , VectorDouble32& __dVm)
{
   //define the constants''', template)
        out.inc()

        iprinter = InterpolatePrintCPUVisitor(out, model.ssa, params, interps)
        model.printTarget(good,computeAllDepend&constants,iprinter)

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
            out('const double %s=state_[__ii].%s;',var,var)
        good |= diffvars
        good.add(V)

        out("//get the gate updates (diagonalized exponential integrator)")
        gateGoals = set()
        for RLA,RLB in gateTargets.values():
            gateGoals.add(RLA)
            gateGoals.add(RLB)
        good.add(dt)
        gateSet = model.allDependencies(good,gateGoals)-good
        model.printSet(gateSet,iprinter)
        good |= gateSet

        out("//get the other differential updates")
        diffGoals = set([diffvarUpdate[var] for var in diffvars-gates])
        diffSet = model.allDependencies(good,diffGoals)-good
        model.printSet(diffSet, iprinter)
        good |= diffSet

        out("//get Iion")
        IionSet = model.allDependencies(good|allfits,set([Iion]))-good
        model.printSet(IionSet, iprinter)
        good |= IionSet
        
        out("//Do the markov update (1 step rosenbrock with gauss siedel)")
        for var in order(markovs):
            out("double %s = %s;",markovTargets[var],diffvarUpdate[var])
        out("int _count=0;")
        out("double _error;")
        out("do")
        out("{")
        out.inc()
        iprinter.pushStack()

        for var in order(markovs):
            out("double %s = %s;",markovOld[var], markovTargets[var])
        markovGoals = set(markovOld.values())
        good |= markovGoals
        cprinter.declaredStack[-1] |= set(["%s" % var for var in markovTargets.values()])
        markovSet = model.allDependencies(good,set(markovTargets.values()))
        model.printSet(markovSet,cprinter)

        out("_error = 0;")
        for var in order(markovs):
            old = markovOld[var]
            new = markovTargets[var]
            out("_error += (%s-%s)*(%s-%s);",old,new,old,new)
        out("_count++;")

        iprinter.popStack()
        out.dec()
        out("} while (_error > 1e-100 && _count<50);")

        out("//EDIT_STATE")
        for var in order(diffvars-gates-markovs):
            out('state_[__ii].%s += _dt*%s;',var,diffvarUpdate[var])
        for var in order(gates):
            (RLA,RLB) = gateTargets[var]
            out('state_[__ii].%(v)s += %(a)s*(%(v)s+%(b)s);',
                v=var,
                a=RLA,
                b=RLB,
            )
        for var in order(markovs):
            out("state_[__jj].%s += _dt*%s;",var,markovTargets[var])
                     
        out("__dVm[__ii] = -%s;", Iion)
        out.dec(2)

        out('''
   }
}''',template)
    else:
        assert(False)
    out('''

   
string ThisReaction::methodName() const
{
   return "%(target)s";
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

    if arch=="nvidia":
        def stateName(var,index):
            return "stateData[_%s_off*nCells_+%s]" % (var,index)
    elif arch=="cpu":
        def stateName(var,index):
            return "state_[%s].%s" % (index,var)
    else:
        assert(False)
        
    if arch=="nvidia":
        out("ArrayView<double> stateData = stateTransport_;")
    elif arch=="cpu":
        out("state_.resize(nCells_);")
    else:
        assert(False)
    out(r'for (int iCell=0; iCell<nCells_; iCell++)')
    out('{')
    out.inc()
    for var in order(diffvars):
        out("%s = %s;",stateName(var,"iCell"),var)
    out.dec(2)
    out('''
   }

   __Vm.assign(__Vm.size(), V_init);
}

enum varHandles
{''', template)
    out.inc()
    for var in order(diffvars)+order(tracevars-diffvars):
        out('%s_handle,',var)
    out.dec()
    out('''
   NUMHANDLES
};

const string ThisReaction::getUnit(const std::string& varName) const
{
   if(0) {}''', template)
    out.inc()
    for var in order(diffvars|tracevars):
        varUnit = model.ssa[var].astUnit
        if varUnit.isNull():
            unitText = "INVALID"
        else:
            unitText = str(varUnit.rawUnit)
        out('else if (varName == "%s") { return "%s"; }',var,unitText)
    out.dec()
    out('''
   return "INVALID";
}

int ThisReaction::getVarHandle(const std::string& varName) const
{
   if (0) {}''', template)
    out.inc()
    for var in order(diffvars|tracevars):
        out('else if (varName == "%s") { return %s_handle; }', var, var)
    out.dec()
    out('''

   return -1;
}

void ThisReaction::setValue(int iCell, int varHandle, double value) 
{''', template)
    out.inc()
    if arch=="nvidia":
        out("ArrayView<double> stateData = stateTransport_;")
    elif arch=="cpu":
        pass
    else:
        assert(False)
    out("if (0) {}")
    for var in order(diffvars):
        out('else if (varHandle == %s_handle) { %s = value; }', var, stateName(var,"iCell"))
    out.dec()
    out('''
}


double ThisReaction::getValue(int iCell, int varHandle) const
{''', template)
    out.inc()
    if arch=="nvidia":
        out("ConstArrayView<double> stateData = stateTransport_;")
    elif arch=="cpu":
        pass
    else:
        assert(False)
    out("if (0) {}")
    for var in order(diffvars):
        out('else if (varHandle == %s_handle) { return %s; }', var, stateName(var,"iCell"))
    out.dec()
    out('''

   return NAN;
}

double ThisReaction::getValue(int iCell, int varHandle, double V) const
{
''', template)
    out.inc()
    if arch=="nvidia":
        out("ConstArrayView<double> stateData = stateTransport_;")
    elif arch=="cpu":
        pass
    else:
        assert(False)
    for var in order(diffvars):
        out('const double %s=%s;',var,stateName(var,"iCell"))
    out('if (0) {}')
    for var in order(diffvars|tracevars):
        out('else if (varHandle == %s_handle)', var)
        out('{')
        out.inc()
        good=diffvars|set([V])|params
        cprinter = CPrintVisitor(out, model.ssa, good)
        model.printTarget(good,set([var])-params,cprinter)
        out('return %s;',var)
        out.dec()
        out('}')
    out.dec()
    out('''

   return NAN;
}

void ThisReaction::getCheckpointInfo(vector<string>& fieldNames,
                                     vector<string>& fieldUnits) const
{
   fieldNames.clear();
   fieldUnits.clear();''', template)
    out.inc()
    for var in order(diffvars):
        out('fieldNames.push_back("%s");',var)
        out('fieldUnits.push_back(getUnit("%s"));',var)
    out.dec()
    out('''
}

}''', template)


generators = {
    frozenset(["cardioid"]) : lambda model, targetName : generateCardioid(model,targetName,"cpu"),
    frozenset(["cardioid", "cpu"]) : lambda model, targetName : generateCardioid(model,targetName,"cpu"),
    frozenset(["cardioid","nvidia"]) : lambda model, targetName : generateCardioid(model,targetName,"nvidia"), 
}
