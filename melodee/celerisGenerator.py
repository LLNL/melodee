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
from sympy.printing.fortran import FCodePrinter
from sympy.core import S


from melodee.parser import MelodeeParser,Differentiator
from melodee import utility
from melodee.utility import order

def repeat(thing, repetitions):
    return (thing,) * repetitions

class MyFCodeSympyPrinter(FCodePrinter):
    def __init__(self,*args,**kwargs):
        FCodePrinter.__init__(self,*args,**kwargs)
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
        return '%s**%s' % (self._print(expr.base),
                                self._print(expr.exp))
    def _print_Relational(self,expr):
        if expr.rel_op == "==" or expr.rel_op == "!=":
            PREC = sympy.printing.precedence.precedence(expr)
            return "%s %s %s" % (self.parenthesize(expr.lhs, PREC),
                                 expr.rel_op,
                                 self.parenthesize(expr.rhs, PREC))
        else:
            return super(MyFCodeSympyPrinter, self)._print_Relational(expr)
    def fortranifyCode(self, code):
        codeLines = code.splitlines()
        paddedLines = self._pad_leading_columns(codeLines)
        wrappedLines = self._wrap_fortran(paddedLines)
        return "\n"+"\n".join(wrappedLines)
        
def pretty(symbol):
    return str(symbol)

class FPrintVisitor:
    def __init__(self, ssa, decltype="double"):
        self.newBuffer()
        self.ssa = ssa
        self.decltype = decltype
        self.cprinter = MyFCodeSympyPrinter()

    def newBuffer(self):
        self.ioBuf = io.BytesIO()
        self.out = utility.Indenter(self.ioBuf)

    def getBuffer(self):
        return self.ioBuf.getvalue()
        
    def ifPrint(self,printer,ifSymbol,thenList,elseList,choiceList):
        self.out("if (%s)",pretty(ifSymbol))
        self.out("then")
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
        self.out("end if")
    def equationPrint(self,lhs,rhs):
        rhsText = self.cprinter._print(rhs.sympy)
        self.equationPrintWithRhs(lhs,rhsText)
    def equationPrintWithRhs(self,lhs,rhsText):
        self.out("%s = %s;",pretty(lhs),rhsText)
    def finalize(self):
        retval = self.fortranifyCode(self.getBuffer())
        self.newBuffer()
        return retval
    def fortranifyCode(self, code):
        return self.cprinter.fortranifyCode(code)
        
def generateCeleris(model, targetName):
    template = {}
    template["target"] = targetName

    diffvars = model.diffvars()
    diffvarUpdate = {var : model.diffvarUpdate(var) for var in diffvars}

    llambda = model.input("lambda")
    actTime = model.input("actTime")
    if "stretchVel" not in model._inputs:
        stretchVel = model.addSymbol("_stretchVel")
    else:
        stretchVel = model.input("stretchVel")
    tension = model.output("tension")

    inputs = set([llambda,actTime,stretchVel])
    ############

    partialFromDiff = {}
    diffPartialFromDiff = {}

    differ = Differentiator(model, set([actTime,stretchVel]))
    partialvars = set()
    for diffvar in diffvars:
        partial = model.addSymbol("_partial_"+str(diffvar))
        differ.cacheResult(diffvar, stretchVel, partial, None)
        partialFromDiff[diffvar] = partial
    partialLambda = model.addSymbol("_partial_lambda")
    differ.cacheResult(llambda, stretchVel, partialLambda, None)
        
    for diffvar in diffvars:
        (diffvarUpdate[partialFromDiff[diffvar]],dontcare) = differ.diff(diffvarUpdate[diffvar],stretchVel)
    (dtension_dstretchVel,dontcare) = differ.diff(tension, stretchVel)
    differ.augmentInstructions()
    ###############

    partialvars = set(partialFromDiff.values()+[partialLambda])
    good=inputs|diffvars|partialvars
        
    indexFromVar = {}
    numRhsVars=0
    for var in order(diffvars):
        numRhsVars += 1
        indexFromVar[var] = numRhsVars        
    for var in order(good-diffvars):
        numRhsVars += 1
        indexFromVar[var] = numRhsVars        

    
    diffTargets = set()
    diffTargets |= set(diffvarUpdate.values())

    tensionTargets = set()
    tensionTargets |= set([tension])
    tensionTargets |= set([dtension_dstretchVel])

    fprinter = FPrintVisitor(model.ssa)    
    out = utility.Indenter(open(targetName+".f","w"))
    out('''
      SUBROUTINE %(target)sRHS (M_NEQ, M_T, M_SV, M_DSV)
!     Declare inputs and outputs
      REAL (KRSIM) :: M_T, M_SV(*), M_DSV(*)

''', template)
    
    out('''
!     Declare the local variables
''', template)

    good=inputs|diffvars|partialvars
    diffDepend = model.allDependencies(good,diffTargets)|good
    varnames = set([pretty(var) for var in diffDepend])
    code = ""
    for name in order(varnames):
        code += "REAL (KRSIM) :: %s\n" % name 
    out(fprinter.fortranifyCode(code))
    
    out('''
!     Set up the initial conditions from SV
''', template)
    code = ""
    for var in order(good):
        code += "%s = M_SV[%d];\n" % (pretty(var), indexFromVar[var])
    out(fprinter.fortranifyCode(code))
    out('''
!     Compute the derivatives
''', template)

    model.printTarget(diffvars|partialvars|inputs, diffTargets, fprinter)
    out(fprinter.finalize())
    out('''
!     Store the derivatives
''', template)
    code = ""
    for diffvar in order(diffvars):
        code += "M_DSV[%d] = %s;\n" %(indexFromVar[diffvar], pretty(diffvarUpdate[diffvar]))
        partial = partialFromDiff[diffvar]
        code += "M_DSV[%d] = %s;\n" %(indexFromVar[partial], pretty(diffvarUpdate[partial]))
    code += "M_DSV[%d] = 1;\n" % indexFromVar[actTime]
    code += "M_DSV[%d] = %s;\n" % (indexFromVar[llambda], pretty(stretchVel))
    code += "M_DSV[%d] = 1;\n" % indexFromVar[partialLambda]
    code += "M_DSV[%d] = 0;\n" % indexFromVar[stretchVel]
    out(fprinter.fortranifyCode(code))

    template["numRhsVars"] = numRhsVars
    template["numDiffVars"] = len(diffvars)
    out('''

      RETURN
      END SUBROUTINE %(target)sRHS


      SUBROUTINE %(target)s (SV_0, PREV_LAMBDA, NEXT_LAMBDA, DT, 
     @  PREV_ACTTIME, SV_1, TENSION, DTENSION)

      REAL (KRSIM), INTENT (IN) :: SV_0(%(numDiffVars)d), THIS_LAMBDA
      REAL (KRSIM), INTENT (OUT) :: SV_1(%(numDiffVars)d), TENSION, DTENSION
      REAL (KRSIM), INTENT (IN) :: PREV_LAMBDA, NEXT_LAMBDA, DT, 
     @ PREV_ACTTIME
      REAL (KRSIM), INTENT (OUT) :: TENSION, DTENSION

!     Declare the local variables
      REAL (KRSIM) :: ODEPACK_VARS(%(numRhsVars)d)
''', template)
    
    varnames = set([pretty(var) for var in model.allDependencies(good,tensionTargets)|good])
    code = ""
    for name in order(varnames):
        code += "REAL (KRSIM) :: %s\n" % name 
    out(fprinter.fortranifyCode(code))
    out('''
!    Setup the local variables correctly
''', template)
    code = ""
    for diffvar in order(diffvars):
        code += "ODEPACK_VARS[%d] = SV_0[%d];\n" %(indexFromVar[diffvar],indexFromVar[diffvar])
        partial = partialFromDiff[diffvar]
        code += "ODEPACK_VARS[%d] = 0;\n" %(indexFromVar[partial])
    code += "ODEPACK_VARS[%d] = PREV_ACTIME;\n" % indexFromVar[actTime]
    code += "ODEPACK_VARS[%d] = PREV_LAMBDA;\n" % indexFromVar[llambda]
    code += "ODEPACK_VARS[%d] = 0;\n" % indexFromVar[partialLambda]
    code += "ODEPACK_VARS[%d] = (NEXT_LAMBDA-PREV_LAMBDA)/DT;\n" % indexFromVar[stretchVel]    
    out(fprinter.fortranifyCode(code))
    out('''
    
!     Integrate the equations
''', template)
    out('''

!     Evaluate tension
''', template)
    code = ""
    for diffvar in order(diffvars):
        code += "%s = ODEPACK_VARS[%d];\n" %(pretty(diffvar) , indexFromVar[diffvar])
        partial = partialFromDiff[diffvar]
        code += "%s = ODEPACK_VARS[%d];\n" %(pretty(partial), indexFromVar[partial])
    code += "%s = PREV_ACTTIME+DT;\n" % pretty(actTime)
    code += "%s = NEXT_LAMBDA;\n" % pretty(llambda)
    code += "%s = DT;\n" % pretty(partialLambda)
    code += "%s = (NEXT_LAMBDA-PREV_LAMBDA)/DT;\n" % pretty(stretchVel)
    out(fprinter.fortranifyCode(code))
    

    model.printTarget(diffvars|partialvars|inputs, tensionTargets, fprinter)
    out(fprinter.finalize())
    out('''
!     Set up the final outputs
''', template)
    code = ""
    for diffvar in order(diffvars):
        code += "SV_1[%d] = %s;\n" %(indexFromVar[diffvar], pretty(diffvar))
    code += "TENSION = %s;\n" % pretty(tension)
    code += "DTENSION = %s/DT;\n" % pretty(dtension_dstretchVel)
    out(fprinter.fortranifyCode(code))
    out('''

      RETURN
      END SUBROUTINE %(target)s

''', template)

generators = {
    frozenset(["celeris"]) : generateCeleris,
}
