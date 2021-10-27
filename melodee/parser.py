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

import re
import sympy


import ply.lex as lex
import ply.yacc as yacc

from melodee import units
from melodee.utility import order,itemsOrderedByKey,itemsOrderedByValue

###################################################################
#these parts make up the external API

class ConsolidatedSystem:
    '''Output a system that reduces an encap to a useable core'''
    def __init__(self):
        self._namedVariables = {}
        self._inputs = {}
        self._outputs = {}
        self._diffvars = {}
        self.ssa = SSA()
        self.time = None
        self.instructions = []
    
    def dependencies(self, var):
        return self.ssa.dependencies(var)
    
    def allDependencies(self, good, target):
        depend = set()
        front = set(target)
        while front:
            newFront = set()
            for var in front:
                depend.add(var)
                if var not in good:
                    newFront |= self.dependencies(var)
            front = newFront
        return depend

    def allExcluding(self,origReachable,bad,instructions=None):
        if instructions == None:
            instructions = self.instructions
        reachable = set(origReachable)
        for inst in instructions:
            if isinstance(inst, IfInstruction):
                if inst.ifVar not in reachable:
                    continue
                reachable |= self.allExcluding(reachable,bad,inst.thenInstructions)
                reachable |= self.allExcluding(reachable,bad,inst.elseInstructions)
                reachable |= self.allExcluding(reachable,bad,inst.choiceInstructions)
            else:
                if (inst not in reachable and
                    inst not in bad and
                    self.dependencies(inst)<=reachable):
                    reachable.add(inst)
        return reachable-origReachable
    
    
    def printTarget(self, good, target, printVisitor):
        allUpdate = self.allDependencies(good, target)-good
        self.printSet(allUpdate, printVisitor)
        
    def printSet(self, allUpdate, printVisitor):
        def printer(instructions, allUpdate=allUpdate, printVisitor=printVisitor):
            for inst in instructions:
                if isinstance(inst, IfInstruction) and instructionsInList([inst]) & allUpdate:
                    subInstructions = (instructionsInList(inst.thenInstructions) | instructionsInList(inst.elseInstructions)) & allUpdate
                    if not subInstructions:
                        printVisitor.ifPrint(printer, None, None, None, inst.choiceInstructions)
                    else:
                        printVisitor.ifPrint(printer, inst.ifVar, inst.thenInstructions, inst.elseInstructions, inst.choiceInstructions)
                else:
                    if inst in allUpdate:
                        printVisitor.equationPrint(inst, self.ssa[inst])
        printer(self.instructions)

    def varsWithAttribute(self, name):
        return set(self._namedVariables.get(name, {}).keys())

    def input(self, name):
        return self._inputs[name]
    def inputs(self):
        return set(self._inputs.values())
    def output(self,name): 
        return self._outputs[name]
    def outputs(self):
        return set(self._outputs.values())
    def diffvars(self):
        return set(self._diffvars.keys())
    def diffvarUpdate(self, diffvar):
        return self._diffvars[diffvar]
        
    def make_a(self, name, symbol, info=None):
        self._namedVariables.setdefault(name, {})[symbol] = info

    def info(self, name, symbol):
        return self._namedVariables[name][symbol]

    def addSymbol(self, name):
        if isinstance(name, str):
            return Symbol(name)
        else:
            return name

    def addSSA(self, name, sympy):
        var = self.addSymbol(name)
        self.ssa[var] = AST(sympy,ASTUnit.null())
        return var        
    
    def addInstruction(self, name, sympy):
        var = self.addSymbol(name)
        self.instructions.append(var)
        self.ssa[var] = AST(sympy,ASTUnit.null())
        return var

    def makeNamesUnique(self):
        def uniqueNames(instructions, currentNames=set()):
            for inst in instructions:
                if isinstance(inst, IfInstruction):
                    uniqueNames(inst.thenInstructions,set(currentNames))
                    uniqueNames(inst.elseInstructions,set(currentNames))
                    uniqueNames(inst.choiceInstructions,currentNames)
                else:
                    if inst.name in currentNames:
                        count=1
                        base=inst.name
                        while True:
                            inst.name = "%s_%03d" % (base, count)
                            if inst.name not in currentNames:
                                break
                            else:
                                count += 1
                    currentNames.add(inst.name)
        uniqueNames(self.instructions)

    def extractExpensiveFunctions(self):
        (newInstructions,retval) = self._extractExpensiveFunctions(self.instructions)
        self.instructions = newInstructions
        return retval
    
    def _extractExpensiveFunctions(self, instructions):
        newInstructions = []
        expensiveInstructions = set()
        for inst in instructions:
            if isinstance(inst, IfInstruction):
                newThenList, thenExpense = self._extractExpensiveFunctions(inst.thenInstructions)
                newElseList, elseExpense = self._extractExpensiveFunctions(inst.elseInstructions)
                newInstructions.append(IfInstruction(inst.ifVar,newThenList,newElseList,inst.choiceInstructions))
                expensiveInstructions |= thenExpense
                expensiveInstructions |= elseExpense
            else:
                (thisVars, thisExpr, isExpensive) = self._extractFunctionsFromSympy(
                    self.ssa[inst].sympy, True)
                self.ssa[inst].sympy = thisExpr
                newInstructions += thisVars
                newInstructions.append(inst)
                expensiveInstructions |= set(thisVars)
                if isExpensive:
                    expensiveInstructions.add(inst)
        return (newInstructions, expensiveInstructions)

    def _extractFunctionsFromSympy(self, expr, top=False):
        newVars = []
        newArgs = []
        isExpensive = False
        for arg in expr.args:
            (newVarsFromArg, newArg,dontcare) = self._extractFunctionsFromSympy(arg)
            newVars += newVarsFromArg
            newArgs.append(newArg)
        if newVars:
            expr = expr.func(*newArgs)
        if isinstance(expr.func, type(sympy.Function)) or (
                expr.func == sympy.Pow and (
                    not expr.exp.is_constant() or
                    int(expr.exp) != expr.exp)):
            if top:
                isExpensive = True
            else:
                newSym = self.addSSA("_expensive_functions", expr)
                expr = newSym
                newVars.append(newSym)
        return (newVars, expr, isExpensive)
        
        
class MelodeeParser:
    def __init__(self, *args, **kwargs):
        self._parser = InternalMelodeeParser(*args, **kwargs)
    def parse(self, text):
        return self._parser.parse(text)
    def getModel(self, target):
        return self._parser.getModel(target)
    def getUnits(self):
        return self._parser.si

class Differentiator:
    '''This class differentiates a ConsolidatedSystem, caching results as
    needed.  It can also commit the derivatives into the target
    Consolidated system.
    '''
    def __init__(self, cs, independentVars):
        self.cs = cs
        self.independentVars = independentVars
        self.cache = {}
        self.zero = AST(sympy.sympify("0"),ASTUnit.null())
        self.one = AST(sympy.sympify("1"),ASTUnit.null())
        self.epsilon = Symbol("epsilon")

    def cacheResult(self, var, wrt, resultVar, expr):
        self.cache.setdefault(var, {})[wrt] = (resultVar, expr)
        return (resultVar, expr)
        
    def diff(self, var, wrt):
        if self.cache.get(var, {}).get(wrt, None) != None:
            return self.cache[var][wrt]

        resultVar = Symbol("_d_%s_wrt_%s" %(var,wrt))
        
        if var == wrt:
            return self.cacheResult(var,wrt,resultVar,self.one)
        if var in self.independentVars:
            return self.cacheResult(var,wrt,resultVar,self.zero)

        ast = self.cs.ssa[var]
        if isinstance(ast, Choice):
            (thenVar, thenExpr) = self.diff(ast.thenVar,wrt)
            (elseVar, elseExpr) = self.diff(ast.elseVar,wrt)
            if (not isinstance(thenExpr, Choice) and
                not isinstance(elseExpr, Choice) and
                thenExpr.sympy == elseExpr.sympy):
                return self.cacheResult(var,wrt,resultVar,thenExpr)
            else:
                resultAst = Choice(ast.ifVar,thenVar,elseVar,
                                   ASTUnit.null()
                                   )
                return self.cacheResult(var,wrt,resultVar,resultAst)
        else:
            deps = ast.dependencies()
            subMap = {}
            for dep in deps:
                (diffVar,diffExpr) = self.diff(dep,wrt)
                result=None
                if not (isinstance(diffExpr, Choice) or diffExpr == None):
                    if diffExpr.sympy==self.zero.sympy:
                        result = dep
                    elif diffExpr.sympy==self.one.sympy:
                        result = dep + self.epsilon
                if result==None:
                    result = dep + self.epsilon * diffVar
                subMap[dep] = result
            diffSympy = ast.sympy.subs(subMap)
            diffSympy = diffSympy.diff(self.epsilon).subs(self.epsilon,0)
            diffAst = AST(diffSympy,ASTUnit.null())
            return self.cacheResult(var,wrt,resultVar,diffAst)

    def _augmentInstructions(self,instructions):
        newInstructions = []
        for inst in instructions:
            if isinstance(inst, IfInstruction):
                newThenList = self._augmentInstructions(inst.thenInstructions)
                newElseList = self._augmentInstructions(inst.elseInstructions)
                newChoiceList = []
                for choice in inst.choiceInstructions:
                    newChoiceList.append(choice)
                    if choice in self.cache:
                        for (dontcare, (diffvar,diffexpr)) in self.cache[choice].items():
                            if isinstance(diffexpr, Choice):
                                newChoiceList.append(diffvar)
                            else:
                                newInstructions.append(diffvar)
                            self.cs.ssa[diffvar] = diffexpr
                newInstructions.append(IfInstruction(inst.ifVar,newThenList,newElseList,newChoiceList))
            else:
                newInstructions.append(inst)
                if inst in self.cache:
                    for (dontcare, (diffvar, diffexpr)) in self.cache[inst].items():
                        newInstructions.append(diffvar)
                        self.cs.ssa[diffvar] = diffexpr
        return newInstructions

    def augmentInstructions(self):
        self.cs.instructions = self._augmentInstructions(self.cs.instructions)



#############################################################################

class XXXSyntaxError(SyntaxError):
    def __init__(self, text):
        self.text = text
        print(text)
    def __str__(self):
        return self.text

class Symbol(sympy.Dummy):
    def __new__(cls, *args, **kwargs):
        return super(Symbol,cls).__new__(cls,*args,**kwargs)
    def __str__(self):
        ret = sympy.Dummy.__str__(self)
        #return ret[1:]
        return ret

class SSA(dict):
    def __setitem__(self, key, value):
        if key in self:
            print(key, value)
            raise MultipleSymbolAssignment(key)
        dict.__setitem__(self,key,value)
    def dependencies(self, symbolMaybeList):
        # get dependecies for a list or a single symbol
        try:
            iter(symbolMaybeList)
        except TypeError:
            #just a regular symbol
            return self[symbolMaybeList].dependencies()
        else:
            ret = set()
            for sym in symbolMaybeList:
                ret |= self[sym].dependencies()
            return ret


class IfInstruction:
    def __init__(self, ifVar, thenInstructions, elseInstructions, choiceInstructions):
        self.ifVar = ifVar
        self.thenInstructions = thenInstructions
        self.elseInstructions = elseInstructions
        self.choiceInstructions = choiceInstructions

def instructionsInList(instructions):
    retval = set()
    for inst in instructions:
        if isinstance(inst, IfInstruction):
            retval |= instructionsInList(inst.thenInstructions)
            retval |= instructionsInList(inst.elseInstructions)
            retval |= instructionsInList(inst.choiceInstructions)
        else:
            retval.add(inst)
    return retval
        
class Choice:
    def __init__(self, ifVar, thenVar, elseVar, unit):
        self.ifVar = ifVar
        self.thenVar = thenVar
        self.elseVar = elseVar
        self.astUnit = unit
    def __repr__(self):
        return 'Choice(%s,%s,%s,%s)' % (repr(self.ifVar),
                                        repr(self.thenVar),
                                        repr(self.elseVar),
                                        repr(self.astUnit),
        )
    def dependencies(self):
        return set([self.ifVar, self.thenVar, self.elseVar])
    def subs(self, subMap):
        newIfVar = subMap.get(self.ifVar, self.ifVar)
        newThenVar = subMap.get(self.thenVar, self.thenVar)
        newElseVar = subMap.get(self.elseVar, self.elseVar)
        return Choice(newIfVar,newThenVar,newElseVar,self.astUnit)

class AST:
    def __init__(self, sympyExpr, unit=None):
        self.sympy = sympyExpr
        self.astUnit = unit
    def __repr__(self):
        return 'AST(sympify(' + repr(self.sympy) +'), '+ repr(self.astUnit) +')'
    def dependencies(self):
        return set(self.sympy.atoms(Symbol))
    def subs(self, subMap):
        return AST(self.sympy.subs(subMap),self.astUnit)

def textToAST(text, unit):
    return AST(sympy.sympify(text), unit)

class ASTUnit:
    def __init__(self, unit, explicit):
        self.rawUnit = unit
        self.explicit = explicit
    def isNull(self):
        return self.rawUnit is None
    def __mul__(self, other):
        if other.isNull() or self.isNull():
            return ASTUnit.null()
        else:
            return ASTUnit(self.rawUnit*other.rawUnit, self.explicit or other.explicit)
    def __div__(self, other):
        if other.isNull() or self.isNull():
            return ASTUnit.null()
        else:
            return ASTUnit(self.rawUnit/other.rawUnit, self.explicit or other.explicit)
    def __pow__(self, number):
        if self.isNull():
            return ASTUnit.null()
        else:
            return ASTUnit(self.rawUnit ** number, self.explicit)
    def __repr__(self):
        return 'ASTUnit('+repr(self.rawUnit)+','+repr(self.explicit)+')'
    @staticmethod
    def null():
        return ASTUnit(None, False)

class Scope:
    def __init__(self, parent=None):
        self.parent = parent
        self.symbols = {}
        self.units = {}
        self.instructions = []
    def makeChild(self):
        return Scope(self)
    def getParent(self):
        return self.parent
    def getSymbol(self,name):
        if name in self.symbols:
            return self.symbols[name]
        elif self.parent != None:
            return self.parent.getSymbol(name)
        else:
            raise KeyError(name)
    def hasSymbol(self, name):
        try:
            self.getSymbol(name)
        except KeyError as e:
            return False
        return True
    def getUnit(self, name):
        if name in self.units:
            return self.units[name]
        elif self.parent != None:
            return self.parent.getUnit(name)
        else:
            raise KeyError(name)
    def hasUnit(self,name):
        try:
            self.getUnit(name)
        except KeyError as e:
            return False
        return True
    def setSymbol(self, name, symbol):
        self.symbols[name] = symbol
    def setUnit(self, name, unit):
        self.units[name] = unit
    def addInstruction(self, inst):
        self.instructions.append(inst)

class Port:
    def __init__(self, unit):
        self.rawUnit = unit
        
class Junction:
    def __init__(self,unit):
        self.rawUnit = unit

class Connections:
    def __init__(self):
        self._junctionFromPort = {}
        self._portsFromJunction = {}
    def makeJunction(self,unit):
        junction = Junction(unit)
        self._portsFromJunction[junction] = set()
        return junction
    def makePort(self, junction):
        port = Port(junction.rawUnit)
        self._junctionFromPort[port] = junction
        self._portsFromJunction[junction].add(port)
        return port
    def junctionFromPort(self, port):
        return self._junctionFromPort[port]
    def portsFromJunction(self, junction):
        return set(self._portsFromJunction[junction])
    def removePort(self, port):
        self._portsFromJunction[self._junctionFromPort[port]].remove(port)
        del self._junctionFromPort[port]
    def removeJunction(self, junction):
        assert(not self._portsFromJunction[junction])
        del self._portsFromJunction[junction]
    def substituteJunction(self, matchJunction, replacementJunction):
        if matchJunction != replacementJunction:
            ports = self._portsFromJunction[matchJunction]
            for port in ports:
                self._junctionFromPort[port] = replacementJunction
            self._portsFromJunction[replacementJunction] |= ports
            self._portsFromJunction[matchJunction] = set()
    
class Encapsulation:
    def __init__(self, subsystem=None):
        self.subsystem = subsystem
        self.children = {}
        self.junctions = {}
        self.ports = {}
        self.externalJunctions = set()
    def hasJunction(self, name):
        return name in self.junctions
    def getJunction(self, name):
        return self.junctions[name]
    def setJunction(self, name, junction, external=False):
        self.junctions[name] = junction
        if external:
            self.externalJunctions.add(junction)
    def clearJunctions(self):
        self.junctions = {}
    def isExternal(self, junction):
        return junction in self.externalJunctions
    def addChild(self, name, encapsulation):
        if name in self.children:
            raise MultipleAssignmentDisallowed(name)
        self.children[name] = encapsulation
    def copy(self, connections, newJuncFromOld=None):
        if newJuncFromOld == None:
            newJuncFromOld = {}
        for (name,oldJunction) in self.junctions.items():
            if oldJunction not in newJuncFromOld:
                newJuncFromOld[oldJunction] = connections.makeJunction(oldJunction.rawUnit)

        other = Encapsulation()
        other.subsystem = self.subsystem
        for (name,encap) in self.children.items():
            other.children[name] = encap.copy(connections, newJuncFromOld)
        oldPorts = self.allLocalPorts()
        newPortFromOld = {}
        for oldPort in oldPorts:
            oldJunction = connections.junctionFromPort(oldPort)
            if oldJunction not in newJuncFromOld:
                newJuncFromOld[oldJunction] = connections.makeJunction(oldJunction.rawUnit)
            newJunction = newJuncFromOld[oldJunction]
            newPortFromOld[oldPort] = connections.makePort(newJunction)
        
        for (name,oldJunction) in self.junctions.items():            
            other.junctions[name] = newJuncFromOld[oldJunction]
        for (name,oldPort) in self.ports.items():
            other.ports[name] = newPortFromOld[oldPort]
        for oldJunction in self.externalJunctions:
            other.externalJunctions.add(newJuncFromOld[oldJunction])
        return other
    def allLocalPorts(self):
        return set(self.ports.values())
    def allLocalInternalJunctions(self):
        return set(self.junctions.values()) - self.externalJunctions
    def allEncap(self, ret=None):
        if ret==None:
            ret = []
        ret.append(self)
        for child in self.children.values():
            child.allEncap(ret)
        return ret
    def nameFromJunction(self, lookupJunc, connections):
        #print "---", lookupJunc
        #print "---", self.junctions
        for (name,junction) in self.junctions.items():
            if junction==lookupJunc:
                return name
        #print "---", self.ports    
        for (name,port) in self.ports.items():
            junction = connections.junctionFromPort(port)
            if junction==lookupJunc:
                return name
        possibleNames = set()
        retName = None
        for child in self.children.values():
            possibleName = child.nameFromJunction(lookupJunc,connections)
            if possibleName != None:
                possibleNames.add(possibleName)
                retName = possibleName
        if len(possibleNames) == 1:
            return retName
        else:
            return None
    def nameEncapFromPort(self, ret=None):
        if ret==None:
            ret={}
        for (var, port) in self.ports.items():
            ret[port] = (var, self)
        for child in self.children.values():
            child.nameEncapFromPort(ret)
        return ret
    def allInternalJunctions(self):
        ret = self.allLocalInternalJunctions()
        for child in self.children.values():
            ret |= child.allInternalJunctions()
        return ret
    def allPorts(self):
        ret = self.allLocalPorts()
        for child in self.children.values():
            ret |= child.allPorts()
        return ret
    def removeChild(self, name, connections):
        child = self.children[name]
        child.removeConnections(connections)
        del self.children[name]
    def removeConnections(self, connections):
        for child in self.children.values():
            child.removeConnections(connections)
        for port in self.allLocalPorts():
            connections.removePort(port)
        for junction in self.junctions.values():
            connections.removeJunction(junction)
    def isInput(self, var):
        return var in self.subsystem.inputs
    def isOutput(self, var):
        return var in self.subsystem.outputs
    def isAccum(self, var):
        return var in self.subsystem.accums
    def instructions(self):
        return self.subsystem.scope.instructions
    
def fullName(tupleName):
    return ":".join(tupleName)

def consolidateSystem(timeUnit, rootEncap, connections,si):
        self = ConsolidatedSystem()
        self.si = si
        nameFromEncap={}
        def nameAllEncaps(encap, myName=(), retval=None):
            retval[encap] = myName
            for (postfix, child) in list(encap.children.items()):
                nameAllEncaps(child, myName+(postfix,), retval)
            return retval
        nameFromEncap = nameAllEncaps(rootEncap, (), {})

        allEncap = list(nameFromEncap.keys())

        #go ahead and name all the junctions, get their symbols
        symbolFromJunction = {}
        for (encap, tupName) in nameFromEncap.items():
            for (name,junction) in encap.junctions.items():
                symbolFromJunction[junction] = Symbol(fullName(tupName + (name,)))
                
        # deal with the assignment issues.
        isAccum = set()
        isAssign = set()
        for encap in allEncap:
            for (name,port) in encap.ports.items():
                junction = connections.junctionFromPort(port)
                symbol = symbolFromJunction[junction]
                if encap.isAccum(name):
                    assert(encap.isOutput(name))
                    if junction in isAssign:
                        raise XXXSyntaxError("Can't both assign and accumulate variable %s." % str(symbol))
                    isAccum.add(junction)
                elif encap.isOutput(name):
                    if junction in isAssign:
                        raise XXXSyntaxError("Can't assign to variable %s more than once." % str(symbol))
                    if junction in isAccum:
                        raise XXXSyntaxError("Can't both assign and accumulate variable %s." % str(symbol))
                    isAssign.add(junction)

        # see if we ever mess with time
        timeName = None
        for encap in allEncap:
            if encap.subsystem.time != None:
                timeName = encap.subsystem.time
                if timeUnit == None:
                    timeUnit = encap.subsystem.getUnit(timeName)
                break
        if timeName != None:
            timeSym = Symbol(timeName)
        else:
            timeSym = None

        #build the symbol correspondence
        class SymbolMapping:
            def __init__(self):
                self._map = {}
            def get(self, encap, symbol,tupName):
                if (encap, symbol) not in self._map:
                    self._map[(encap,symbol)] = Symbol(fullName(tupName + (str(symbol),)))
                return self._map[(encap,symbol)]
            def set(self, encap, oldSym, newsymbol):
                self._map[(encap,oldSym)] = newsymbol
            def getList(self, encap, symbolList, tupName):
                return [self.get(encap, oldSym, tupName) for oldSym in symbolList]

        #build the SSA
        convertedInstructions = []
        symbolMap = SymbolMapping()
        accumsFromJunction = {}
        for (encap,tupName) in itemsOrderedByValue(nameFromEncap):
            subsystem = encap.subsystem
            for (name,port) in encap.ports.items():
                junction = connections.junctionFromPort(port)
                newsymbol = symbolFromJunction[junction]
                if encap.isInput(name):
                    symbolMap.set(encap,subsystem.getVar(name),newsymbol)
                elif encap.isAccum(name):
                    newsymbol = symbolMap.get(encap,subsystem.getVar(name),tupName)
                    symbolMap.set(encap,subsystem.getVar(name),newsymbol)
                    accumsFromJunction.setdefault(junction, []).append(newsymbol)
                else:
                    assert(encap.isOutput(name))
                    symbolMap.set(encap,subsystem.getVar(name),newsymbol)

            #process time vars here.
            if subsystem.time != None:
                assert(timeSym != None)
                oldSym = subsystem.getVar(subsystem.time)
                oldUnit = subsystem.getUnit(subsystem.time)
                if oldUnit == timeUnit:
                    symbolMap.set(encap, oldSym, timeSym)
                else:
                    newSym = symbolMap.get(encap, oldSym,tupName)
                    convertedInstructions.append(newSym)
                    self.ssa[newSym] = AST(sympy.Mul(timeSym,oldUnit.convertTo(timeUnit)),ASTUnit(timeUnit,False))
            #process other named variables
            #for oldSym in subsystem.allSymbols():
            #    symbolMap.get(encap,oldSym,tupName)
            #process the remaining symbols
            for (oldSym, oldAst) in subsystem.ssa.items():
                newSym = symbolMap.get(encap,oldSym,tupName)
                localSubs = {}
                for depSym in oldAst.dependencies():
                    localSubs[depSym] = symbolMap.get(encap,depSym,tupName)
                self.ssa[newSym] = oldAst.subs(localSubs)
        
        # Add SSA symbols for accumulations
        for (junction,accums) in accumsFromJunction.items():
            assert(len(accums) >= 1)
            newSymbol = symbolFromJunction[junction]
            convertedInstructions.append(newSymbol)
            self.ssa[newSymbol] = AST(sympy.Add(*accums), ASTUnit(junction.rawUnit,False))

        #get the inputs/outputs
        for (name,junction) in itemsOrderedByKey(rootEncap.junctions):
            if rootEncap.isExternal(junction):
                symbol = symbolFromJunction[junction]
                if junction in isAssign or junction in isAccum:
                    self._outputs[str(symbol)] = symbol
                else:
                    self._inputs[str(symbol)] = symbol

        # fix time
        self.time = timeSym

        #get the diffvars
        for (encap,tupName) in itemsOrderedByValue(nameFromEncap):
            subsystem = encap.subsystem
            appendAfter = {}
            for name in subsystem.diffvars:
                oldSym = subsystem.getVar(name)
                newSym = symbolMap.get(encap,oldSym,tupName)

                diffUnit = subsystem.getUnit(name)
                initSym = symbolMap.get(encap,subsystem.getVar(name+".init"),tupName)
                self.ssa[newSym] = AST(initSym,ASTUnit(diffUnit,False))
                appendAfter[initSym] = newSym

                updateUnit = subsystem.getUnit(name+".diff")
                updateSym = symbolMap.get(encap,subsystem.getVar(name+".diff"),tupName)
                if timeUnit == None:
                    timeUnit = diffUnit/updateUnit
                if diffUnit/updateUnit == timeUnit:
                    diffvarUpdate = updateSym
                else:
                    newName = fullName(tupName + (name+".diff_convertUnit",))
                    newDiff = Symbol(newName)
                    appendAfter[updateSym] = newDiff
                    self.ssa[newDiff] = AST(sympy.Mul(updateSym,updateUnit.convertTo(diffUnit/timeUnit)),
                                           ASTUnit(diffUnit/timeUnit,False));
                    diffvarUpdate = newDiff
                self._diffvars[newSym] = diffvarUpdate
            def convert(inst, symbolMap=symbolMap, encap=encap, tupName=tupName):
                if isinstance(inst,IfInstruction):
                    return IfInstruction(
                        symbolMap.get(encap,inst.ifVar,tupName),
                        [convert(sub) for sub in inst.thenInstructions],
                        [convert(sub) for sub in inst.elseInstructions],
                        [convert(sub) for sub in inst.choiceInstructions],
                        )
                else:
                    return symbolMap.get(encap,inst,tupName)
            for inst in encap.instructions():
                converted = convert(inst)
                convertedInstructions.append(converted)
                if converted in appendAfter:
                    convertedInstructions.append(appendAfter[converted])

        for (encap,tupName) in itemsOrderedByValue(nameFromEncap):
            subsystem = encap.subsystem
            for (name,thisAttrMap) in subsystem.attributeMap.items():
                oldSym = subsystem.getVar(name)
                for (attr,info) in thisAttrMap.items():
                    self.make_a(attr, symbolMap.get(encap,oldSym,tupName), info)

        #order the converted instructions
        def symsFromInstructionGenerator(root):
            if isinstance(root, IfInstruction):
                ret = set()
                for inst in root.thenInstructions:
                    ret |= symsFromInstructionGenerator(inst)
                for inst in root.elseInstructions:
                    ret |= symsFromInstructionGenerator(inst)
                ret |= set(root.choiceInstructions)
                return ret
            else:
                return set([root])
        symsFromInstruction= {}
        for inst in convertedInstructions:
            symsFromInstruction[inst] = symsFromInstructionGenerator(inst)
        dependFromInstruction={}
        for inst in convertedInstructions:
            syms = symsFromInstruction[inst]
            dependFromInstruction[inst] = self.ssa.dependencies(syms)-syms

        self.instructions = []
        undefinedInstructions = set(convertedInstructions)
        definedSymbols = set()
        if self.time != None:
            definedSymbols.add(self.time)
        definedSymbols |= set(self._inputs.values())

        while undefinedInstructions:
            removedThisIter = False
            for inst in convertedInstructions:
                if inst not in undefinedInstructions:
                    continue
                syms = symsFromInstruction[inst]
                depend = dependFromInstruction[inst]
                #print "===",inst
                #print "+++",symsFromInstruction[inst]
                #print "+++",depend
                if depend <= definedSymbols:
                    #print "DEFINING", syms
                    removedThisIter = True
                    undefinedInstructions.remove(inst)
                    definedSymbols |= syms
                    self.instructions.append(inst)

            if not removedThisIter:
                print(undefinedInstructions)
                print(len(undefinedInstructions))
                print(definedSymbols)
                assert(removedThisIter)
            #print "---------------------------------------------------------------"

        #now we need to fix variable names.  ideally, we'd like to
        #use the shortest name everywhere if at all possible
            
        #this is just a first pass.  Ideally, we'd like to
        #recognize when a name is no longer used but I'm not
        #requiring that here.
        def uniqueNames(instructions, nameMap={}, currentNames=set()):
            for inst in instructions:
                if isinstance(inst,IfInstruction):
                    uniqueNames(inst.thenInstructions,nameMap,set(currentNames))
                    uniqueNames(inst.elseInstructions,nameMap,set(currentNames))
                    uniqueNames(inst.choiceInstructions,nameMap,currentNames)
                else:
                    simpleName = re.sub(r'\.','_',str(inst))
                    nameArray = simpleName.split(r':')
                    lb=len(nameArray)-1
                    while lb >= 0 and "_".join(nameArray[lb:]) in currentNames:
                        lb += -1
                    if lb == -1:
                        lb = 0
                    name = "_".join(nameArray[lb:])
                    nameMap[inst] = name
                    currentNames.add(name)
            return nameMap
        nameMap = uniqueNames(self.instructions)
        #now rename all the variables.
        for (var,name) in nameMap.items():
            var.name = name

        self.makeNamesUnique()
        return self

class Subsystem:
    def __init__(self):
        self.ssa = SSA()
        self.scope = Scope()
        self.inputs = set()
        self.outputs = set()
        self.diffvars = set()
        self.accums = set()
        self.attributeMap = {}
        self.time = None
        self.frozen = set()

    def getAllDependencies(self, target):
        try:
            iter(target)
            iterableTarget = target
        except TypeError:
            iterableTarget = [target]
        front = set(iterableTarget)
        dependencies = set()
        while front:
            dependencies |= front
            newFront = set()
            for symbol in front:
                if symbol in self.ssa:
                    newFront |= self.ssa[symbol].dependencies()
            front = newFront-dependencies
        return dependencies

    def getVar(self, name):
        return self.scope.getSymbol(name)
    def getUnit(self,name):
        return self.scope.getUnit(name)
            
def strifyInstructions(ilist, ssa, indent=0):
    myIndent = "   "*indent
    ret = ""
    for instruction in ilist:
        if isinstance(instruction,IfInstruction):
            ret += myIndent + "if (%s) {\n" % instruction.ifVar
            ret += strifyInstructions(instruction.thenInstructions,ssa,indent+1)
            ret += myIndent + "} else {\n"
            ret += strifyInstructions(instruction.elseInstructions,ssa,indent+1)
            ret += myIndent + "}\n"
            for choice in instruction.choiceInstructions:
                ret += myIndent + "%s = %s;\n" % (choice, ssa[choice])
        else:
            ret += myIndent + "%s = %s;\n" % (instruction, ssa[instruction])
    return ret

        
class InternalMelodeeParser:
    def __init__(self, **kw):
        self.debug = kw.get('debug', 0),
        self.start = kw.get('start', 'topLevelStatementsOpt')
        self.lexer = lex.lex(module=self, debug=self.debug)
        self.parser = yacc.yacc(module=self,
                                debug=self.debug,
                                write_tables=0,
                                start=self.start,
        )
        self.si = units.Si()
        self.connections = Connections()
        self.scopeStack = []
        self.tempCount = 0
        self.enumerations = {}
        self.encapsulationStack = [Encapsulation()]
        self.timeUnitFromEncapName = {}
        self.clearEnvironment()

    def parse(self, text):
        self.clearEnvironment()
        self.parser.parse(text)
    def getModel(self, name):
        rootEncap = self.encapsulationStack[0].children[name]
        timeUnit = self.timeUnitFromEncapName[name]
        return consolidateSystem(timeUnit,rootEncap,self.connections,self.si)
        
    def clearEnvironment(self):
        self.timeVar = None
        self.timeUnit = None
        self.encapsulationStack[0].clearJunctions()
        self.lexer.lineno = 1
        
    def currentSubsystem(self):
        return self.currentEncapsulation().subsystem
    def currentScope(self):
        return self.scopeStack[-1]
    def currentEncapsulation(self):
        return self.encapsulationStack[-1]
    def pushScope(self):
        self.scopeStack[-1] = self.scopeStack[-1].makeChild()
        return self.scopeStack[-1]
    def popScope(self):
        retval = self.scopeStack[-1]
        self.scopeStack[-1] = self.scopeStack[-1].getParent()
        return retval
    def readAccessVar(self, var):
        encapsulation = self.currentEncapsulation()
        subsystem = self.currentSubsystem()
        scope = self.currentScope()
        #if the var is an accum
        if var in subsystem.accums:
            raise XXXSyntaxError("Can't read from accumulation variable '%s'." % var)
        #if the var has a symbol already
        if scope.hasSymbol(var):
            #get the symbol
            symbol = scope.getSymbol(var)
            #if this symbol has a unit defined
            if scope.hasUnit(var):
                rawUnit = scope.getUnit(var)
            elif symbol in subsystem.ssa:
                rawUnit = subsystem.ssa[symbol].astUnit.rawUnit
            else:
                rawUnit = None
            return AST(symbol, ASTUnit(rawUnit, explicit=False))
        else:
            #FIXME the approach here breaks in if conditions if we allow vars to be overwritten (like enumerations)
            #see if this is a shared var
            (exists, junction) = self.searchForJunction(var)
            if exists:
                port = self.connections.makePort(junction)
                encapsulation.ports[var] = port
                #make a symbol for this guy, mark it as an input
                symbol = Symbol(var)
                subsystem.scope.setSymbol(var, symbol)
                subsystem.scope.setUnit(var, junction.rawUnit)
                subsystem.inputs.add(var)
                #repeat the lookup command now that we have a symbol
                return self.readAccessVar(var)
            #check for a timevar
            elif var == self.timeVar:
                #make a symbol for this guy, mark it as time
                symbol = Symbol(var)
                subsystem.scope.setSymbol(var, symbol)
                subsystem.scope.setUnit(var, self.timeUnit)
                subsystem.time = var
                #repeat the lookup for this command now that we have a symbol
                return self.readAccessVar(var)
            #check for an enumeration
            elif var in self.enumerations:
                return self.enumerations[var]
            else:
                raise XXXSyntaxError("Variable '%s' used but not defined." % var)

    def processDeclaration(self, var, structure, attributes, unitOpt):
        if "provides" in structure:
            junction = self.markProvides(var,unitOpt)
            self.currentEncapsulation().ports[var] = self.connections.makePort(junction)
        else:
            self.checkDeclarable(var)
            self.processUnitAssignment(var, unitOpt)
            
        if "accum" in structure:
            self.currentSubsystem().accums.add(var)
        elif "diffvar" in structure:
            if self.timeUnit == None:
                raise XXXSyntaxError("Must include an 'integrate' statement before declaring subsystems.")
            if not self.currentScope().hasUnit(var):
                raise XXXSyntaxError("Differential variable %s must have a unit." % var)

            self.currentSubsystem().diffvars.add(var)
            self.currentScope().setSymbol(var, Symbol(var))

            rawUnit = self.currentScope().getUnit(var)
            self.currentScope().setUnit(var+".init", rawUnit)
            self.currentScope().setUnit(var+".diff", rawUnit/self.timeUnit)
        
        for (name, info) in attributes:
            self.currentSubsystem().attributeMap.setdefault(var, {})[name] = info

    def processUnitAssignment(self, var, unitOpt):
        if unitOpt != None:
            if self.currentScope().hasUnit(var):
                if self.currentScope().getUnit(var) != unitOpt:
                    raise XXXSyntaxError("Declared units for '%s' don't match previous declaration." % var)
            else:
                self.currentScope().setUnit(var, unitOpt)

    def processAssignment(self, var, operand, rhs):
        subsystem = self.currentSubsystem()
        scope = self.currentScope()
        #error if this is a diffvar
        if var in subsystem.diffvars:
            raise XXXSyntaxError("Can't assign to differential variable '%s'." % var)
        if var in subsystem.inputs:
            raise XXXSyntaxError("Can't assign to shared variable '%s'." % var)
        if var == self.timeVar:
            raise XXXSyntaxError("Can't assign to integration variable '%s'." % var)
        if var in subsystem.accums:
            if operand != '+=' and operand != '-=':
                raise XXXSyntaxError("Can only += or -= accumulation variable '%s'."%var)
        if var not in subsystem.outputs and not scope.hasSymbol(var):
            (exists, dontcare) = self.searchForJunction(var)
            if exists:
                raise XXXSyntaxError("Can't assign to shared variable '%s'." % var)

        if operand != '=':
            if not scope.hasSymbol(var):
                if var in subsystem.accums:
                    lhs = textToAST("0",ASTUnit(scope.getUnit(var),False))
                else:
                    raise XXXSyntaxError("'%s' used before assignment."%var)
            else:
                lhs = self.readAccessVar(var)

            if operand == '+=' or operand == '-=':
                if operand == '-=':
                    rhs = AST(sympy.Mul(sympy.Integer(-1),rhs.sympy), rhs.astUnit)
                rhs = AST(sympy.Add(lhs.sympy,rhs.sympy), self.checkExactUnits(lhs.astUnit,rhs.astUnit))
            elif operand == '*=' or operand == '/=':
                if operand == "/=":
                    rhs = AST(sympy.Pow(rhs.sympy,sympy.Integer(-1)), rhs.astUnit ** -1)
                rhs = AST(sympy.Mul(lhs.sympy,rhs.sympy), lhs.astUnit*rhs.astUnit)
            elif operand == '^=':
                rhs = self.powerProcess(lhs, rhs)
            else:
                assert(0)

        #Time to do final unit checks
        if scope.hasUnit(var):
            #print var
            #print scope.getUnit(var)
            #print rhs
            rhs = self.checkExplicitCast(scope.getUnit(var), rhs)
        
        #ok, we're ready to do the assignment!
        if var in subsystem.frozen:
            XXXSyntaxError('"%s" cannot be assigned once it has been read.' % var)
        symbol = Symbol(var)
        subsystem.ssa[symbol] = rhs
        scope.addInstruction(symbol)
        scope.setSymbol(var, symbol)

        #parameter checking code.
        #FIXME
        nameDepend = set([symbol.name for symbol in rhs.dependencies()])
        newlyReadParams = (set(subsystem.attributeMap.keys()) & nameDepend) - subsystem.frozen
        newlyReadParams -= set([var])
        subsystem.frozen |= newlyReadParams
        
    def checkDeclarable(self, var):
        #if the variable is already type defined in this subsystem
        ss = self.currentSubsystem()
        if var in ss.inputs:
            raise XXXSyntaxError("'%s' already used as a shared variable in this subsystem." % var)
        if var in ss.diffvars:
            raise XXXSyntaxError("'%s' already defined as a differential variable in this subsystem." % var)
        if var in ss.accums:
            raise XXXSyntaxError("'%s' already defined as an accumulation output for this subsystem." % var)
        if var in ss.outputs:
            raise XXXSyntaxError("'%s' already an output for this subsystem." % var)
        if var == self.timeVar:
            raise XXXSyntaxError("'%s' already used as the integration variable." % var) 
        #FIXME
        #if the variable has already been assigned in this subsystem
        #if ss.scope.hasSymbol(var) and ss.scope.getSymbol(var) in ss.ssa:
        #    raise XXXSyntaxError("'%s' has already been assigned to, it can't now be declared." % var)
        #if the variable is already unit defined in this subsystem
        #if ss.scope.hasUnit(var):
        #    raise XXXSyntaxError("Units previously defined for '%s'" % var) 
    
    def markProvides(self, varname, unitOpt):
        self.checkDeclarable(varname)
        #if this variable connects to another variable
        (exists, junction) = self.searchForJunction(varname)
        if exists:
            #get the unit for that junction
            junctionUnit = junction.rawUnit
            #if a unit was defined and that unit does not match the junction unit
            if unitOpt != None and unitOpt != junctionUnit:
                raise XXXSyntaxError("Provided units for '%s' don't match the shared units")
            unit = junctionUnit
        else:
            #if a unit was not provided
            if unitOpt == None:
                raise XXXSyntaxError("Must have a unit for provides variable '%s'." % varname)
            unit = unitOpt
            junction = self.connections.makeJunction(unit)
            #make a new junction with a unit
            self.currentEncapsulation().setJunction(varname, junction)
        self.currentScope().setUnit(varname, unit)
        self.currentSubsystem().outputs.add(varname)
        #return the junction
        return junction

    def searchForJunction(self, name):
        assert(len(self.encapsulationStack) >= 1)
        for ii in range(len(self.encapsulationStack)-1,-1,-1):
            if self.encapsulationStack[ii].hasJunction(name):
                return (True, self.encapsulationStack[ii].getJunction(name))
        return (False, None)

    def newTempVar(self):
        current = self.tempCount
        self.tempCount += 1
        return "__melodee_temp_%03d" % current
    
    def processIfCondition(self, ifExpr, thenScope, elseScope):
        choiceInstructions = []

        self.checkExactUnits(ifExpr.astUnit, self.boolean())
        var = self.newTempVar()
        symbol = Symbol(var)
        self.currentSubsystem().ssa[symbol] = ifExpr
        self.currentScope().addInstruction(symbol)
        ifSymbol = symbol
        
        #iterate over local symbols
        for var in order(set(thenScope.symbols.keys()) | set(elseScope.symbols.keys())):
            if not thenScope.hasSymbol(var) or not elseScope.hasSymbol(var):
                continue
            #make sure the units match
            if thenScope.hasUnit(var) and elseScope.hasUnit(var):
                unit = thenScope.getUnit(var)
                if unit != elseScope.getUnit(var):
                    raise XXXSyntaxError("if condition for '%s' declares different units depending on branch." % var)
                self.currentScope().setUnit(var, unit)
            elif not thenScope.hasUnit(var) and not elseScope.hasUnit(var):
                unit = None
            else:
                raise XXXSyntaxError("if condition for '%s' declares different units depending on branch." % var)
            
            if var in self.currentSubsystem().frozen:
                XXXSyntaxError("Due to '%s' being used in an if condition, I don't know what it's default value should be." % var)
            choice = Choice(ifSymbol, thenScope.getSymbol(var), elseScope.getSymbol(var),
                            ASTUnit(unit, False))
            symbol = Symbol(var)
            self.currentSubsystem().ssa[symbol] = choice
            choiceInstructions.append(symbol)
            self.currentScope().setSymbol(var, symbol)

        iif = IfInstruction(ifSymbol,thenScope.instructions,elseScope.instructions,choiceInstructions)
        self.currentScope().addInstruction(iif)

    def lookupEncapsulation(self, name):
        assert(len(self.encapsulationStack) >= 1)
        for ii in range(len(self.encapsulationStack)-1,-1,-1):
            if name in self.encapsulationStack[ii].children:
                return self.encapsulationStack[ii].children[name]
        raise XXXSyntaxError("Can't find encapsulation '%s'." % name)
    
    def processUseStatement(self, childName, encap, removeList, exportList):
        newEncap = encap.copy(self.connections)
        # delete from the encap as necessary
        for removeItem in removeList:
            removeEncap = newEncap
            while len(removeItem) > 1:
                removeEncap = removeEncap.children[removeItem[0]]
                removeItem = removeItem[1:]
            removeEncap.removeChild(removeItem[0], self.connections)
        # do the exporting that has been manually specified.
        junctionFromExportedName = {}
        for item in exportList:
            (speccedName, exportName) = item
            exportEncap = newEncap
            while len(speccedName) > 1:
                exportEncap = exportEncap.children[speccedName[0]]
                speccedName = speccedName[1:]
            name = speccedName[0]
            if name in exportEncap.junctions:
                newJunction = exportEncap.junctions[name]
            elif name in exportEncap.ports:
                newJunction = self.connections.junctionFromPort(exportEncap.ports[name])
            else:
                raise XXXSyntaxError("Can't find variable '%s' for exporting." % name)
            if exportName in junctionFromExportedName:
                raise XXXSyntaxError('Cannot export "%s" from "%s": name is already exported.' % (exportName,childName))
            junctionFromExportedName[exportName] = newJunction
        # find all the ports
        allPorts = newEncap.allPorts()
        # split those ports up between input, output, and accum ports
        # compute the junctions for each of the port lists
        nameEncapFromPort = newEncap.nameEncapFromPort()
        inputJunctions = set()
        outputJunctions = set()
        accumJunctions = set()
        for port in allPorts:
            junction = self.connections.junctionFromPort(port)
            (name, subEncap) = nameEncapFromPort[port]
            if subEncap.isInput(name):
                inputJunctions.add(junction)
            elif subEncap.isOutput(name):
                outputJunctions.add(junction)
            elif subEncap.isAccum(name):
                accumJunctions.add(junction)
            else:
                assert(False)

        # get a list of junctions declared externally to the encapsulation.
        internalJunctions = newEncap.allInternalJunctions()
        allJunctions = inputJunctions|outputJunctions|accumJunctions
        externalJunctions = allJunctions - internalJunctions
        # figure out the implicitly exported junctions
        danglingInputJunctions = inputJunctions - outputJunctions - accumJunctions
        exportedJunctions = set(junctionFromExportedName.values())
        implicitlyExportedJunctions = (danglingInputJunctions | externalJunctions) - exportedJunctions
        #print "+++++", childName
        #print danglingInputJunctions
        #print externalJunctions
        #print junctionFromExportedName
        #print newEncap.junctions
        #print newEncap.externalJunctions
        
        # get names for each of the implicitly exported junctions
        for junction in implicitlyExportedJunctions:
            #find the name for this junction
            name = newEncap.nameFromJunction(junction, self.connections)
            if name == None:
                allNames = set()
                for port in self.connections.portsFromJunction(junction):
                    (name,encap) = nameEncapFromPort[port]
                    allNames.add(name)
                errorMessage = "Can't get a unique name for an implicily exported junction in child \"%s\".\n" % childName
                errorMessage += "Please export manually. Candidate names are:\n"
                for name in allNames:
                    errorMessage += name + "\n"
                raise XXXSyntaxError(errorMessage)
            #check for possible name conflicts
            if name in junctionFromExportedName:
                raise XXXSyntaxError('Cannot implicitly export "%s" from "%s": another variable has already been exported with this name.' % (name,childName))
            # add the export
            junctionFromExportedName[name] = junction

        # export the junctions
        for (exportName,newJunction) in junctionFromExportedName.items():
            (exists, oldJunction) = self.searchForJunction(exportName)
            if not exists:
                self.currentEncapsulation().setJunction(exportName, newJunction,
                                                        external=(newJunction in externalJunctions))
            else:
                self.connections.substituteJunction(newJunction,oldJunction)

        self.currentEncapsulation().addChild(childName, newEncap)

    def checkConnections(self, rootEncap):
        #Get a list of all internal junctions
        internalJunctions = rootEncap.allInternalJunctions()
        #make a mapping between ports and encapsulation/name pairs
        nameEncapFromPort = rootEncap.nameEncapFromPort()
        #for each internal junction
        errors = ""
        for thisJunction in internalJunctions:
            #get the ports for this junction
            thesePorts = self.connections.portsFromJunction(thisJunction)
            #if any of the ports are input ports
            foundInputPorts = False
            for (var, encap) in [nameEncapFromPort[port] for port in thesePorts]:
                if encap.isInput(var):
                    foundInputPorts = True
                    break;
            if foundInputPorts:
                #make sure the port is assigned a value by someone
                outputs = []
                accums = []
                for (var, encap) in [nameEncapFromPort[port] for port in thesePorts]:
                    if encap.isOutput(var):
                        outputs.append((var,encap))
                    if encap.isAccum(var):
                        accums.append((var,encap))
                assert(len(accums) <= len(outputs))
                if len(outputs) < 1:
                    errors += "No definition for %s in %s\n"  %(var, rootEncap)
                if len(outputs) > 1 and len(outputs) != len(accums):
                    errors += "Can't mix output and accum for %s" % var
        if errors:
            raise XXXSyntaxError(errors)
        
    def nodim(self):
        return ASTUnit(self.si.get("1"), False)
    def boolean(self):
        return ASTUnit(self.si.get("bool"), False)
    
    def checkExplicitCast(self, castRaw, ast):
        retval=AST(ast.sympy, ASTUnit(castRaw, explicit=True))
        if ast.astUnit.isNull():
            return retval
        elif ast.astUnit.rawUnit == castRaw:
            return retval
        else:
            raise XXXSyntaxError("Can't cast "+str(ast.astUnit)+" to unit of "+str(castRaw)+")")

    def checkExactUnits(self, lunit, runit):
        if lunit.isNull() or runit.isNull():
            if lunit.explicit or runit.explicit:
                raise XXXSyntaxError("Explicit units cannot combine with null units!")
            else:
                return ASTUnit.null()
        else:
            # now we have something with units!  Check if they match.
            resultIsExplicit = lunit.explicit or runit.explicit
            if lunit.rawUnit == runit.rawUnit:
                return ASTUnit(lunit.rawUnit, resultIsExplicit)
            else:
                if resultIsExplicit:
                    raise XXXSyntaxError("Units don't match, maybe you need a conversion?")
                else:
                    return ASTUnit.null()
    def convertUnitTo(self, expr, newUnit):
        if expr.astUnit.isNull():
            raise XXXSyntaxError("Can't convert a null unit!")
        elif not expr.astUnit.rawUnit.isCompatibleWith(newUnit):
            raise XXXSyntaxError("Incompatible unit conversion requested.")
        else:
            factor = expr.astUnit.rawUnit.convertTo(newUnit)
            return AST(sympy.Mul(factor,expr.sympy), ASTUnit(newUnit, explicit=False))
    
    #def astToVar(self, var, ast):
    #    self.currentSubsystem().ssa[var] = ast
    #    self.currentScope().addInstruction(var)
    #    return (var, ast.astUnit)
        
    #def astToTemp(self, ast):
    #    return self.astToVar(self.newTempVar(),ast)
    #def astToSymbol(self, name, ast):
    #    (var, astUnit) = self.astToVar(Symbol(name),ast)
    #    self.currentScope().setSymbol(name, var)
    #    return (var, astUnit)
        
    
    t_ignore = " \t\r"
                                    
    def t_newline(self, t):
        r'\n+'
        t.lexer.lineno += len(t.value)

    def t_ignore_C_COMMENT(self, t):
        r'/\*(?:.|\n)*?\*/'
        #count the number of new lines
        t.lexer.lineno += t.value.count("\n")
        
    def t_error(self, t):
        print("Illegal character '%s'" % t.value[0])
        t.lexer.skip(1)

    reserved = {
        "if" : "IF",
        "else" : "ELSE",
        "bool" : "BOOL",
        "true" : "TRUE",
        "false" : "FALSE",
        "accum" : "ACCUM",
        "diffvar" : "DIFFVAR",
        "diff" : "DIFF",
        "init" : "INIT",
        "shared" : "SHARED",
        "provides" : "PROVIDES",
        "subsystem" : "SUBSYSTEM",
        "pow" : "POW",
        "convert" : "CONVERT",
        "enum" : "ENUM",
        
        "integrate" : "INTEGRATE",

        "use" : "USE",
        "export" : "EXPORT",
        "as" : "AS",
    }

    tokens = (
        "LEQ",
        "GEQ",
        "NEQ",
        "BOOLEQ",
        "AND",
        "OR",
        "NOT",
        "PLUSEQ",
        "MINUSEQ",
        "TIMESEQ",
        "DIVIDEEQ",
        "EXPONEQ",
        "NUMBER",
        "NAME",
        "ONE",
    ) + tuple(reserved.values())

    literals = "(){}<>!;?:+-/*^=.,@"

    def t_AND(self, t):
        r'(?:and)|(?:&&)'
        return t
    def t_OR(self, t):
        r'(?:or)|(?:\|\|)'
        return t
    def t_NOT(self, t):
        r'!|(?:not)'
        return t

    def t_NAME(self, t):
        r'[_a-zA-Z][_a-zA-Z0-9?]*'
        t.type = self.reserved.get(t.value,"NAME")
        return t

    def t_NUMBER(self, t):
        r'(?:(?:[0-9]+(?:\.[0-9]*)?)|(?:\.[0-9]+))(?:[eE][-\+]?[0-9]+)?'
        if t.value == "1":
            t.type = "ONE"
        return t

    t_ignore_CPP_COMMENT = r'//.*?(?=\n)'
    t_ignore_PYTHON_COMMENT = r'\#.*?(?=\n)'
    t_ignore_MATLAB_COMMENT = r'%.*?(?=\n)'

    t_LEQ = r'<='
    t_GEQ = r'>='
    t_NEQ = r'!='
    t_BOOLEQ = r'=='
    
    t_PLUSEQ = r'\+='
    t_MINUSEQ = r'-='
    t_TIMESEQ = r'\*='
    t_DIVIDEEQ = r'/='
    t_EXPONEQ = r'^='
    

    precedence = (
        #("nonassoc", "UNITRESOLVE"),
        #("nonassoc", "BOOLRESOLVE"),
        #("left", '?', ':', 'TERNARY'),
        #("left", "OR"),
        #("left", "AND"),
        #("nonassoc", "<", ">", "BOOLEQ", "LEQ", "GEQ", "NEQ"),
        ('left', '+', '-'),
        ('left', '*', '/'),
        #('right', 'ANNOTATE_UNIT'),
        #('right', 'UMINUS'),
        ('left', '^'),
       )

    def p_topLevelStatementsOpt(self, p):
        '''topLevelStatementsOpt : topLevelStatement topLevelStatementsOpt'''
        pass

    def p_topLevelStatementsOpt_term(self, p):
        '''topLevelStatementsOpt : empty'''
        pass

    def p_topLevelStatement(self, p):
        '''topLevelStatement : sharedStatement
                             | integrateStatement
        '''
        pass
    def p_topLevelSubsystem(self, p):
        '''topLevelStatement : subSystemBegin topLevelEnvironmentFix subSystemStatementsOpt '}' '''
        name = p[1]
        encap = self.closeSubsystem()
        self.timeUnitFromEncapName[name] = self.timeUnit
        self.currentEncapsulation().addChild(name,encap)
        self.checkConnections(encap)
        
    def p_topLevelEnvironmentFix(self, p):
        '''topLevelEnvironmentFix : empty'''
        thisEncap = self.currentEncapsulation()
        topEncap = self.encapsulationStack[-2]
        for (name,oldJunction) in topEncap.junctions.items():
            newJunction = self.connections.makeJunction(oldJunction.rawUnit)
            thisEncap.setJunction(name,newJunction,external=True)
        
    def p_integrateStatement(self, p):
        '''integrateStatement : INTEGRATE var unitDef ';'
        '''
        self.timeVar = p[2]
        self.timeUnit = p[3]
        
    def p_sharedStatement_unit(self, p):
        '''sharedStatement : SHARED var unitDef ';'
        '''
        junction = self.connections.makeJunction(p[3])
        self.currentEncapsulation().setJunction(p[2], junction)

    def p_nameList_term(self, p):
        '''nameList : NAME'''
        p[0] = [p[1]]
    def p_nameList_shift(self, p):
        '''nameList : NAME ',' nameList'''
        p[0] = [p[1]] + p[3]

    def p_subSystemDefinition(self, p):
        '''subSystemDefinition : subSystemBegin subSystemStatementsOpt '}' '''
        p[0] = (p[1],self.closeSubsystem())

    def closeSubsystem(self):
        self.scopeStack.pop()
        thisEncap = self.encapsulationStack.pop()
        return thisEncap    
        
    def p_subSystemBegin(self, p):
        '''subSystemBegin : SUBSYSTEM NAME '{' '''
        thisName = p[2]
        self.encapsulationStack.append(Encapsulation(Subsystem()))
        self.scopeStack.append(self.currentSubsystem().scope)
        p[0] = thisName
        
    def p_subSystemStatementsOpt(self, p):
        '''subSystemStatementsOpt : subSystemStatement subSystemStatementsOpt
                                  | empty 
        '''
        pass

    def p_subSystemStatement_shared(self, p):
        '''subSystemStatement : sharedStatement'''
        pass

    def p_subSystemStatement_subSystem(self, p):
        '''subSystemStatement : subSystemDefinition'''
        (name, encap) = p[1]
        self.currentEncapsulation().addChild(name,encap)

    def p_attributeList_single(self,p):
        '''attributeList : attribute'''
        p[0] = [p[1]]
    def p_attributeList_multiple(self, p):
        '''attributeList : attributeList attribute'''
        p[0] = p[1] + [p[2]]
    def p_attribute_bare(self, p):
        '''attribute : '@' NAME '''
        p[0] = (p[2], "")
    def p_attribute_withCapture(self, p):
        '''attribute : '@' NAME '(' nestedParenCaptured ')' '''
        p[0] = (p[2], p[4])
    def p_nestedParenCaptured_empty(self, p):
        '''nestedParenCaptured : empty'''
        p[0] = ""
    def p_nestedParenCaptured_notParen(self, p):
        '''nestedParenCaptured : nestedParenCaptured notParen'''
        p[0] = p[1] + p[2]
    def p_nestedParenCaptured_paren(self, p):
        '''nestedParenCaptured : nestedParenCaptured '(' nestedParenCaptured ')' '''
        p[0] = p[1] + '(' + p[3] + ')'
    def p_notParen(self, p):
        '''notParen : IF
        | ELSE
        | BOOL
        | TRUE
        | FALSE
        | ACCUM
        | DIFFVAR
        | DIFF
        | INIT
        | SHARED
        | PROVIDES
        | SUBSYSTEM
        | POW
        | CONVERT
        | ENUM
        | INTEGRATE
        | USE
        | EXPORT
        | AS
        | LEQ
        | GEQ
        | NEQ
        | BOOLEQ
        | AND
        | OR
        | NOT
        | PLUSEQ
        | MINUSEQ
        | TIMESEQ
        | DIVIDEEQ
        | EXPONEQ
        | NUMBER
        | ONE
        | NAME
        | '{'
        | '}'
        | '<'
        | '>'
        | '!'
        | ';'
        | '?'
        | ':'
        | '+'
        | '-'
        | '/'
        | '*'
        | '^'
        | '='
        | '.'
        | ','
        | '@'
        '''
        p[0] = p[1]

    def p_subSystemStatement_singleDeclaration(self, p):
        '''subSystemStatement : declaration declList ';' '''
        (structure,attributes) = p[1]
        declList = p[2]
        for (var,unitOpt) in declList:
            self.processUnitAssignment(var, unitOpt)
            self.processDeclaration(var, structure, attributes, unitOpt)
    def p_subSystemStatement_declDefinition(self, p):
        '''subSystemStatement : declaration generalizedAssignment ';' '''
        (structure,attributes) = p[1]
        (var, unitOpt, operand, rhs) = p[2]
        self.processUnitAssignment(var, unitOpt)
        self.processDeclaration(var, structure,attributes, unitOpt)
        self.processAssignment(var, operand, rhs)        

    def p_declaration_both(self, p):
        '''declaration : structureDecl attributeList'''
        p[0] = (p[1],p[2])
    def p_declaration_nostructure(self, p):
        '''declaration : attributeList'''
        p[0] = (set(),p[1])
    def p_declaration_noattribute(self, p):
        '''declaration : structureDecl'''
        p[0] = (p[1],[])

    def p_structureDecl_two(self, p):
        '''structureDecl : PROVIDES ACCUM
                         | PROVIDES DIFFVAR '''
        p[0] = set([p[1], p[2]])
    def p_structureDecl_one(self, p):
        '''structureDecl : PROVIDES
                         | DIFFVAR
        '''
        p[0] = set([p[1]])

    def p_unitDecl(self,p):
        '''unitDecl : varDiffOpt unitOpt '''
        p[0] = (p[1],p[2])
        
    def p_declList_single(self, p):
        '''declList : declItem '''
        p[0] = p[1]
    def p_declList_multiple(self, p):
        '''declList : declList ',' declItem '''
        p[0] = p[1] + p[3]
    def p_declItem_single(self, p):
        '''declItem : unitDecl'''
        p[0] = [p[1]]
    def p_declItem_multiple(self, p):
        '''declItem : var varDiffList unitOpt'''
        varDiffList = [p[1]] + p[2]
        unitOpt = p[3]
        p[0] = [(var,unitOpt) for var in varDiffList]

    def p_varDiffList_single(self, p):
        '''varDiffList : varDiffOpt'''
        p[0] = [p[1]]
    def p_varDiffList_multiple(self, p):
        '''varDiffList : varDiffList varDiffOpt'''
        p[0] = p[1] + [p[2]]

    
    def p_unitOpt_empty(self, p):
        '''unitOpt : empty '''
        p[0] = None
    def p_unitOpt_unit(self, p):
        '''unitOpt : unitDef '''
        p[0] = p[1]


    def p_subSystemStatement_definition(self, p):
        '''subSystemStatement : conditionalStatement '''
        pass

    def p_conditionalStatementsOpt_root(self, p):
        '''conditionalStatementsOpt : empty'''
        pass
    def p_conditionalStatementsOpt_recuse(self, p):
        '''conditionalStatementsOpt : conditionalStatementsOpt conditionalStatement'''
        pass

    def p_conditionalStatement_varDefn(self, p):
        '''conditionalStatement : generalizedAssignment ';' '''
        (var, unitOpt, operand, rhs) = p[1]
        self.processUnitAssignment(var, unitOpt)
        self.processAssignment(var, operand, rhs)        

    def p_generalizedAssignment(self, p):
        '''generalizedAssignment : unitDecl '=' realExpr
                                 | unitDecl PLUSEQ realExpr
                                 | unitDecl MINUSEQ realExpr
                                 | unitDecl TIMESEQ realExpr
                                 | unitDecl DIVIDEEQ realExpr
                                 | unitDecl EXPONEQ realExpr
        '''
        (var,unit) = p[1]
        p[0] = (var, unit, p[2], p[3])
         
        
    ########################################################################
    def p_subSystemStatement_if(self, p):
        '''conditionalStatement : ifStatement'''
        pass

    def p_ifStatement(self, p):
        '''ifStatement : initialIfCond thenBody elseOpt'''
        self.processIfCondition(p[1],p[2],p[3])
    def p_initialIfCond(self,p):
        '''initialIfCond : IF '(' realExpr ')' '''
        p[0] = p[3]
    def p_thenBody(self,p):
        '''thenBody : '{' ifScopeBegin conditionalStatementsOpt '}' '''
        p[0] = self.popScope()
    def p_elseOpt(self,p):
        '''elseOpt : ifScopeBegin 
                   | ELSE '{' ifScopeBegin conditionalStatementsOpt '}'
                   '''
        p[0] = self.popScope()

    def p_elseIfCond(self, p):
        '''elseIfCond : ELSE IF '(' realExpr ')' '''
        p[0] = p[4]
    def p_elseOpt_continue(self, p):
        '''elseOpt : elseIfCond ifScopeBegin thenBody elseOpt'''
        self.processIfCondition(p[1],p[3],p[4])
        p[0] = self.popScope()

    def p_ifScopeBegin(self, p):
        '''ifScopeBegin : empty'''
        self.pushScope()

    ##################################################################


    def p_subSystemStatement_use(self, p):
        '''subSystemStatement : useStatement'''
        pass
    def p_useStatement_noBlock(self, p):
        '''useStatement : useEncapsulationSpec ';' '''
        (name, encap, removeList) = p[1]
        exportList = []
        self.processUseStatement(name, encap, removeList, exportList)
    def p_useStatement_withBlock(self, p):
        '''useStatement : useEncapsulationSpec '{' useBlockStatementListOpt '}' '''
        (name, encap, removeList) = p[1]
        exportList = p[3]
        self.processUseStatement(name, encap, removeList, exportList)
    def p_useEncapsulationSpec_simpleSamename(self, p):
        '''useEncapsulationSpec : useBegin useRemoveListOpt'''
        (lookupName, encap) = p[1]
        useName = lookupName
        p[0] = (useName, encap, p[2])
    def p_useEncapsulationSpec_simpleRenamed(self, p):
        '''useEncapsulationSpec : useBegin useRemoveListOpt AS NAME'''
        (lookupName, encap) = p[1]
        useName = p[4]
        p[0] = (useName, self.lookupEncapsulation(lookupName), p[2])
    def p_useBegin_simple(self, p):
        '''useBegin : USE NAME '''
        lookupName = p[2]
        p[0] = (lookupName, self.lookupEncapsulation(lookupName))
    def p_useBegin_complex(self, p):
        '''useBegin : USE NAME '.' speccedName '''
        lookupName = p[2]
        retEncap = self.lookupEncapsulation(lookupName)
        for name in p[4]:
            lookupName = name
            retEncap = retEncap.children[lookupName]
        p[0] = (lookupName, retEncap)
    def p_useRemoveListOpt_empty(self, p):
        '''useRemoveListOpt : empty'''
        p[0] = []
    def p_useRemoveListOpt_recurse(self, p):
        '''useRemoveListOpt : useRemoveListOpt '-' '.' speccedName '''
        p[0] = p[1] + [p[4]]
    def p_speccedName_default(self, p):
        '''speccedName : NAME'''
        p[0] = [p[1]]
    def p_specceedName_specified(self, p):
        '''speccedName : speccedName '.' NAME'''
        p[0] = p[1] + [p[3]]
    def p_useBlockStatementListOpt_empty(self, p):
        '''useBlockStatementListOpt : empty'''
        p[0] = []
    def p_useBlockStatementListOpt_statement(self, p):
        '''useBlockStatementListOpt : useBlockStatementListOpt useBlockStatement'''
        p[0] = p[1] + [p[2]]
    def p_useBlockStatement_exportSimple(self, p):
        '''useBlockStatement : EXPORT speccedName ';' '''
        p[0] = (p[2], p[2][-1])
    def p_useBlockStatement_exportComplex(self, p):
        '''useBlockStatement : EXPORT speccedName AS NAME ';' '''
        p[0] = (p[2], p[4])

    #def p_unitDefBar(self, p):
    #    '''unitDefBar : unitExpr '|' '''
    #    p[0] = p[1]
    def p_unitDefBracket(self, p):
        '''unitDef : '{' unitExpr '}' '''
        p[0] = p[2]

    def p_unitExpr_literal(self, p):
        '''unitExpr : NAME'''
        p[0] = self.si.get(p[1])
    def p_unitExpr_bool(self, p):
        '''unitExpr : BOOL'''
        p[0] = self.si.get(p[1])
    def p_unitExpr_enum(self, p):
        '''unitExpr : ENUM '(' nameList ')' '''
        unitName = "enum("+(",".join(p[3]))+")"
        self.si.addBase(unitName)
        unit = self.si.get(unitName)
        ii=0
        for name in p[3]:
            self.enumerations[name] = textToAST(str(ii), ASTUnit(unit, False))
            ii += 1
        p[0] = unit
    def p_unitExpr_1(self, p):
        '''unitExpr : ONE'''
        p[0] = self.si.get('1')
    def p_unitExpr_paren(self, p):
        '''unitExpr : '(' unitExpr ')' '''
        p[0] = p[2]
    def p_unitExpr_op(self, p):
        '''unitExpr : unitExpr '*' unitExpr
                    | unitExpr '/' unitExpr
        '''
        if p[2] == "*":
            p[0] = p[1]*p[3]
        else:
            p[0] = p[1]/p[3]
    def p_unitExpr_expon(self, p):
        '''unitExpr : unitExpr '^' numberLiteral'''
        p[0] = p[1] ** float(p[3])
    def p_var(self, p):
        '''var : NAME'''
        p[0] = p[1]
    def p_varDiffOpt_base(self, p):
        '''varDiffOpt : var'''
        p[0] = p[1]
    def p_varDiffOpt_diff(self, p):
        '''varDiffOpt : var '.' DIFF
                      | var '.' INIT
        '''
        p[0] = p[1] +'.'+ p[3]
        
    def powerProcess(self, x, y):
        if (not x.astUnit.isNull()) and y.sympy.is_constant():
            newUnit = x.astUnit ** float(y.sympy)
        else:
            newUnit = ASTUnit.null()
        return AST(sympy.Pow(x.sympy, y.sympy), newUnit)
    
    ############################################

    def p_realExpr_pass(self,p):
        '''realExpr : ternaryExpr'''
        p[0] = p[1]

    def p_ternaryOp_pass(self, p):
        '''ternaryExpr : orExpr'''
        p[0] = p[1]
    def p_ternaryOp_impl(self, p):
        '''ternaryExpr : orExpr '?' realExpr ':' realExpr'''
        retUnit = self.checkExactUnits(p[3].astUnit,p[5].astUnit)

        self.checkExactUnits(p[1].astUnit, self.boolean())
        var = self.newTempVar()
        symbol = Symbol(var)
        self.currentSubsystem().ssa[symbol] = p[1]
        self.currentScope().addInstruction(symbol)
        ifSymbol = symbol
        
        var = self.newTempVar()

        symbol = Symbol(var)
        self.currentSubsystem().ssa[symbol] = p[3]
        thenSymbol = symbol

        symbol = Symbol(var)
        self.currentSubsystem().ssa[symbol] = p[5]
        elseSymbol = symbol

        choice = Choice(ifSymbol,thenSymbol,elseSymbol,retUnit)

        symbol = Symbol(var)
        self.currentSubsystem().ssa[symbol] = choice
        choiceSymbol = symbol

        iif = IfInstruction(ifSymbol,[thenSymbol],[elseSymbol],[choiceSymbol])
        self.currentScope().addInstruction(iif)

        p[0] = AST(choiceSymbol, retUnit)
                
    def p_orExpr_pass(self,p):
        '''orExpr : andExpr'''
        p[0] = p[1]
    def p_orExpr_impl(self, p):
        '''orExpr : orExpr OR andExpr'''
        p[0] = AST(sympy.Or(p[1].sympy,p[3].sympy),self.boolean())

    def p_andExpr_pass(self, p):
        '''andExpr : booleqExpr'''
        p[0] = p[1]
    def p_andExpr_impl(self, p):
        '''andExpr : andExpr AND booleqExpr'''
        p[0] = AST(sympy.And(p[1].sympy,p[3].sympy),self.boolean())

    def p_booleqExpr_pass(self, p):
        '''booleqExpr : relationExpr'''
        p[0] = p[1]
    def p_booleqExpr_impl(self,p):
        '''booleqExpr : relationExpr BOOLEQ relationExpr
                      | relationExpr NEQ    relationExpr
        '''
        if p[2] == "NEQ":
            boolOp = sympy.Ne
        else:
            boolOp = sympy.Eq
        self.checkExactUnits(p[1].astUnit,p[3].astUnit)
        p[0] = AST(boolOp(p[1].sympy,p[3].sympy),self.boolean())

    def p_relationExpr_pass(self, p):
        '''relationExpr : additiveExpr'''
        p[0] = p[1]
    def p_relationExpr_impl(self, p):
        '''relationExpr : additiveExpr "<" additiveExpr
                        | additiveExpr ">" additiveExpr
                        | additiveExpr LEQ additiveExpr
                        | additiveExpr GEQ additiveExpr
        '''
        if p[2] == "<":
            boolOp = sympy.Lt
        elif p[2] == ">":
            boolOp = sympy.Gt
        elif p[2] == "LEQ":
            boolOp = sympy.Le
        else:
            boolOp = sympy.Ge
        self.checkExactUnits(p[1].astUnit,p[3].astUnit)
        p[0] = AST(boolOp(p[1].sympy,p[3].sympy),self.boolean())

    def p_additiveExpr_pass(self, p):
        '''additiveExpr : multiplicitiveExpr'''
        p[0] = p[1]
    def p_additiveExpr_impl(self, p):
        '''additiveExpr : additiveExpr '+' multiplicitiveExpr
                        | additiveExpr '-' multiplicitiveExpr
        '''
        lhs = p[1]
        rhs = p[3]
        if p[2] == '-':
            rhs = AST(sympy.Mul(sympy.Integer(-1),rhs.sympy), rhs.astUnit)
        p[0] = AST(sympy.Add(lhs.sympy,rhs.sympy), self.checkExactUnits(lhs.astUnit,rhs.astUnit))

    def p_multiplicitiveExpr_pass(self, p):
        '''multiplicitiveExpr : unaryExpr'''
        p[0] = p[1]
    def p_multiplicitiveExpr_impl(self, p):
        '''multiplicitiveExpr : multiplicitiveExpr '*' unaryExpr
                              | multiplicitiveExpr '/' unaryExpr
        '''
        lhs = p[1]
        rhs = p[3]
        if p[2] == '/':
            rhs = AST(sympy.Pow(rhs.sympy,sympy.Integer(-1)), rhs.astUnit ** -1)
        ret = AST(sympy.Mul(lhs.sympy,rhs.sympy), lhs.astUnit*rhs.astUnit)
        p[0] = ret

    def p_unaryExpr_pass(self, p):
        '''unaryExpr : unitLabelExpr'''
        p[0] = p[1]
    def p_unaryExpr_uminus(self, p):
        '''unaryExpr : '-' unaryExpr'''
        p[0] = AST(sympy.Mul(sympy.Integer(-1),p[2].sympy), p[2].astUnit)
    def p_unaryExpr_not(self, p):
        '''unaryExpr : NOT unaryExpr
                     | '!' unaryExpr
        '''
        p[0] = AST(sympy.Not(p[2].sympy),self.boolean())

    def p_unitExpr_pass(self, p):
        '''unitLabelExpr : exponentExpr'''
        p[0] = p[1]
    def p_unitExpr_impl(self, p):
        '''unitLabelExpr : exponentExpr unitDef '''
        p[0] = self.checkExplicitCast(p[2], p[1])

    def p_exponentExpr_pass(self, p):
        '''exponentExpr : functionExpr'''
        p[0] = p[1]
    def p_exponentExpr_impl(self, p):
        '''exponentExpr : exponentExpr '^' functionExpr'''
        p[0] = self.powerProcess(p[1], p[3])

    def p_functionExpr_impl(self,p):
        '''functionExpr : parenExpr'''
        p[0] = p[1]
    def p_functionExpr_pow(self,p):
        '''functionExpr : POW '(' realExpr ',' realExpr ')' '''
        p[0] = self.powerProcess(p[3],p[5])
    def p_functionExpr_convert(self,p):
        '''functionExpr : CONVERT '(' realExpr ',' unitExpr ')' '''
        p[0] = self.convertUnitTo(p[3], p[5])
    def p_functionExpr_func(self, p):
        '''functionExpr : NAME '(' funcArgListOpt ')' '''
        functionName = p[1]
        if 0:
            pass
        elif functionName == "fabs":
            functionName = "Abs"
        p[0] = AST(getattr(sympy,functionName)(*[x.sympy for x in p[3]]), self.nodim())
    def p_funcArgListOpt_zero(self, p):
        '''funcArgListOpt : empty'''
        p[0] = []
    def p_funcArgListOpt_nonzero(self, p):
        '''funcArgListOpt : funcArgList'''
        p[0] = p[1]
    def p_funcArgList_term(self, p):
        '''funcArgList : realExpr'''
        self.checkExactUnits(p[1].astUnit,self.nodim())
        p[0] = [p[1]]
    def p_funcArgList_shift(self, p):
        '''funcArgList : realExpr ',' funcArgList'''
        self.checkExactUnits(p[1].astUnit, self.nodim())
        p[0] = [p[1]] + p[3]
    def p_parenExpr_pass(self, p):
        '''parenExpr : primaryExpr'''
        p[0] = p[1]
    def p_parenExpr_impl(self, p):
        '''parenExpr : '(' realExpr ')' '''
        p[0] = p[2]

    def p_primaryExpr_var(self, p):
        '''primaryExpr : varDiffOpt'''
        p[0] = self.readAccessVar(p[1])
    def p_primaryExpr_boolLiteral(self, p):
        '''primaryExpr : boolLiteral
        '''
        p[0] = textToAST(p[1], self.boolean())
    def p_boolLiteral(self, p):
        '''boolLiteral : TRUE
                       | FALSE
        '''
        p[0] = p[1]
    def p_primaryExpr_const(self, p):
        '''primaryExpr : numberLiteralPlusOpt'''
        p[0] = textToAST(p[1], ASTUnit.null())
    def p_numberLiteralPlusOpt_noPlus(self, p):
        '''numberLiteralPlusOpt : numberLiteral
        '''
        p[0] = p[1]
    def p_numberLiteralPlusOpt_plus(self, p):
        '''numberLiteralPlusOpt : '+' numberLiteral
        '''
        p[0] = p[2]
    def p_numberLiteral(self, p):
        '''numberLiteral : NUMBER
                         | ONE
        '''
        p[0] = p[1]


    ####################################

    def p_empty(self, p):
        '''empty :'''
        pass

    def p_statementRecovery(self,p):
        '''subSystemStatement : error ';' '''
        print("moving on from syntax error on line number "+str(p.lineno(2)))
    def p_subsystemRecovery(self,p):
        '''subSystemDefinition : subSystemBegin error '}' '''
        print("moving on from syntax error on line number "+str(p.lineno(3)))
        self.scopeStack.pop()
        thisEncapsulation = self.encapsulationStack.pop()
    def p_useRecovery(self,p):
        '''useBlockStatement : error ';' '''
        print("moving on from syntax error on line number "+str(p.lineno(3)))
    #def p_error(self,p):
    #    if p:
    #        print "SyntaxError on line number "+str(self.lexer.lineno)+" at token", p.type, p.value
    #        self.parser.errok()


    #def p_subSystemStatement_renameStatement(self, p):
    #    '''subSystemStatement : renameStatement'''
    #    pass

    
if __name__=="__main__":
    p = InternalMelodeeParser()
    data = '''
and && or || not ! 0 2.0 .3 40. 5e+6 if myID */* bljsadfj */ */
'''
    p.lexer.input(data)

    while True:
        tok = p.lexer.token()
        if not tok:
            break
        print(tok)

    p = InternalMelodeeParser(start="unitExpr")
    print(p.parse("mV/ms"))
    print(p.parse("uA/uF"))

    p = InternalMelodeeParser(start="realExpr")
    p.p_subSystemBegin([None, 'subsystem', "testing", '{'])
    p.currentScope().setSymbol("a", Symbol("a"))
    p.currentScope().setSymbol("b", Symbol("b"))
    p.currentScope().setSymbol("c", Symbol("c"))
    p.currentScope().setSymbol("d", Symbol("d"))
    p.currentScope().setUnit("a", p.si.get("unitless"))
    p.currentScope().setUnit("b", p.si.get("unitless"))
    p.currentScope().setUnit("c", p.si.get("ms"))
    p.currentScope().setUnit("d", p.si.get("s"))
    print(p.parse("a+b/c+d"))
    print(p.parse("a+(b/c)+d"))
    print(p.parse("a+b"))
    print(p.parse("(a+b)"))
    print(p.parse("a"))
    print(p.parse("(a)"))
    print(p.parse("((a))"))
    print(p.parse("((a+b))/((c+d))"))
    print(p.parse("1 {ms}+ c"))
    print(p.parse("convert(1 {ms}, s)+ d {s}"))
    print(p.parse("a == b"))

    HH = '''
integrate time {ms};

shared V {mV};
subsystem HH {
   shared Iion {uA/cm^2};
   provides E_R {mV} = -75;
   shared E_Na {mV};
   subsystem leakage_current {
      E_L {mV} = (E_R+10.613{mV});
      param g_L {mS/cm^2} = 0.3;
      i_L {uA/cm^2} = g_L*(V-E_L);
      provides accum Iion += i_L;
   }
   subsystem potassium_channel {
      shared n {1};
      subsystem potassium_channel_n_gate {
         provides diffvar n {1};
         alpha_n {1/ms} = -0.01{1/mV/ms}*(V+65{mV})/(exp(-(V+65{mV})/10{mV})-1{1});
         beta_n {1/ms} = 0.125{1/ms}*exp((V+75{mV})/80{mV});
         n.init = 0.325;
         n.diff = (alpha_n*(1-n)-beta_n*n);
      }
      E_K {mV} = (E_R-12{mV});
      g_K {mS/cm^2} = 36;
      i_K {uA/cm^2} = g_K*n^4*(V-E_K);
      provides accum Iion += i_K;
   }
   subsystem sodium_channel {
      shared h {1};
      shared m {1};
      subsystem sodium_channel_h_gate {
         provides diffvar h {1};
         alpha_h {1/ms} = 0.07{1/ms}*exp(-(V+75{mV})/20{mV});
         beta_h {1/ms} = 1{1/ms}/(exp(-(V+45{mV})/10{mV})+1);
         h.init = 0.6;
         h.diff = (alpha_h*(1-h)-beta_h*h);
      }
      subsystem sodium_channel_m_gate {
         provides diffvar m {1};
         alpha_m {1/ms} = -0.1{1/mV/ms}*(V+50{mV})/(exp(-(V+50{mV})/10{mV})-1);
         beta_m {1/ms} = 4{1/ms}*exp(-(V+75{mV})/18{mV});
         m.init = 0.05;
         m.diff = (alpha_m*(1-m)-beta_m*m);
      }
      g_Na {mS/cm^2} = 120;
      i_Na {uA/cm^2} = g_Na*m^3*h*(V-E_Na);
      provides accum Iion += i_Na;
   }
   subsystem stimulus {
      i_Stim {uA/cm^2} = ((time >= 10{ms} && time <= 10.5{ms}) ? 20{uA/cm^2} : 0{uA/cm^2});
      provides accum Iion += i_Stim;
   }
   subsystem membrane {
      provides diffvar V {mV};
      Cm {uF/cm^2} = 1;
      provides E_Na {mV} = (E_R+115{mV});
      V.init = -75;
      V.diff = -Iion/Cm;
   }
}


shared E_Na {mV};
subsystem newINa {
  diffvar m {1};
  diffvar j {1};
  diffvar h {1};

  alpha_h = ((V < -40) ? 0.057*exp(-(V+80)/6.8) : 0);
  beta_h = ((V < -40) ? (2.7*exp(0.079*V)+310000*exp(0.3485*V)) : 0.77/0.13*(1+exp((V+10.66)/-11.1)));
  h_inf = 1/pow((1+exp((V+71.55)/7.43)),2);
  tau_h = 1/(alpha_h+beta_h);
  h.init = 0.7573;
  h.diff = (h_inf-h)/tau_h;

  alpha_j = ((V < -40) ? (-25428*exp(0.2444*V)-6.948e-6*exp(-0.04391*V))*(V+37.78)/1/(1+exp(0.311*(V+79.23))) : 0);
  beta_j = ((V < -40) ? 0.02424*exp(-0.01052*V)/(1+exp(-0.1378*(V+40.14))) : 0.6*exp(0.057*V)/(1+exp(-0.1*(V+32))));
  j_inf = 1/pow((1+exp((V+71.55)/7.43)),2);
  tau_j = 1/(alpha_j+beta_j);
  j.init = 0.7225;
  j.diff = (j_inf-j)/tau_j;

  alpha_m = 1/(1+exp((-60-V)/5));
  beta_m = (0.1/(1+exp((V+35)/5))+0.1/(1+exp((V-50)/200)));
  m_inf = 1/pow((1+exp((-56.86-V)/9.03)),2);
  tau_m = 1*alpha_m*beta_m;
  m.init = 0.00155;
  m.diff = (m_inf-m)/tau_m;

  g_Na = 14.838;
  i_Na = g_Na*pow(m,3)*h*j*(V-E_Na);
  provides accum i_Natot {uA/cm^2} += i_Na;
}

subsystem modifiedModel {
  use HH - .sodium_channel {
    export leakage_current.V;
    export Iion;
    export E_Na;
  }
  use newINa as INa {
    export V;
    export i_Natot as Iion;
    export E_Na;
  } 
}
'''
    p = InternalMelodeeParser(start="topLevelStatementsOpt")
    p.parse(HH)
    model = p.getModel("modifiedModel")
    print(strifyInstructions(model.instructions, model.ssa))
