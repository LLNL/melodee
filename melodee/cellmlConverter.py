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

from xml.etree import cElementTree as ET
import sys
import re
from melodee.utility import Indenter,order,unzip

inlineUnits = True

xmlns = {"http://www.cellml.org/cellml/1.0#" : "cellml",
         "http://www.w3.org/1998/Math/MathML" : "math",
         }

def registerNamespaces():
    for longName, shortName in xmlns.items():
        ET.register_namespace(shortName, longName)

def stripNamespaces(elementRoot):
    for el in elementRoot.getiterator():
        if '}' in el.tag:
            el.tag = el.tag.split('}', 1)[1]  # strip all namespaces
        for at in list(el.attrib.keys()): # strip namespaces of attributes too
            if '}' in at:
                newat = at.split('}', 1)[1]
                el.attrib[newat] = el.attrib[at]
                del el.attrib[at]
    return elementRoot

class Equation:
    def __init__(self, lhs, rhs, unit=None, dependencies=set(), isDiff=False):
        self.lhs = lhs
        self.rhs = rhs
        self.hasUnit = (unit!=None)
        self.unit = unit
        self.dependencies = dependencies
        self.isDiff = isDiff
        self.initialCondition = None

    def addInitialCondition(self, eqn):
        self.initialCondition = eqn.rhs
        self.dependencies |= eqn.dependencies

    def toCode(self, out):
        if self.hasUnit:
            unitText = " {"+printUnit(self.unit)+"}"
        else:
            unitText = ""
        if self.isDiff:
            if self.initialCondition:
                out("%s.init%s = %s;" % (self.lhs, unitText, self.initialCondition))
            out("%s.diff = %s;" % (self.lhs, self.rhs))
        else:
            out("%s%s = %s;" % (self.lhs, unitText, self.rhs))

def parseEquation(eqnElement, varToUnit):
    if eqnElement[0].tag != "eq":
        print(ET.tostring(eqnElement))
    assert(eqnElement[0].tag == "eq")
    lhsElement = eqnElement[1]
    rhsElement = eqnElement[2]

    isDiff = 0
    if lhsElement.tag == "ci":
        lhs = lhsElement.text.strip()
    else:
        assert(lhsElement.tag == "apply")
        assert(lhsElement.find("diff") != None)
        assert(lhsElement.find("bvar/ci") != None)
        assert(lhsElement.find("ci") != None)
        lhs = lhsElement.find("ci").text.strip()
        isDiff=1

    (rhs, depend) = parseRhs(rhsElement)

    unit = varToUnit[lhs]
    return Equation(lhs, rhs, unit, depend, isDiff)

def checkParen(sss):
    count=0
    for char in sss:
        if char == '(':
            count += 1
        if char == ')':
            count -= 1
            if count < 0:
                return False
    return count==0
        

def protectDivide(text):
    ret = text
    if not (ret[0] == '(' and ret[-1] == ')' and checkParen(ret[1:-1])):
        unitless = re.sub(r'\{.*?\}','',ret)
        if unitless.count('*') or unitless.count('/'):
            ret = '(' + ret + ')'
    return ret

def parseRhs(rhsElement):
    multiplicitiveOps = { "times" : "*", "divide": "/" }
    boolOps = { "eq" : "==",
                "neq" : "!=",
                "gt" : ">",
                "lt" : "<",
                "geq" : ">=",
                "leq" : "<=",
                "and" : "&&",
                "or" : "||",
                "xor" : "^(XOR_FIXCONVERTER)",
                }
    functionOps = { "power": "pow",
                    "root": "sqrt",
                    "abs": "fabs",
                    "exp": "exp",
                    "ln": "log",
                    "floor": "floor",
                    "ceiling": "ceil",
                    "factorial": "fact_FIXCONVERTER",
                    "sin": "sin",
                    "cos": "cos",
                    "tan": "tan",
                    "sec": "sec_FIXCONVERTER",
                    "csc": "csc_FIXCONVERTER",
                    "cot": "cot_FIXCONVERTER",
                    "sinh": "sinh",
                    "cosh": "cosh",
                    "tanh": "tanh",
                    "sech": "sech_FIXCONVERTER",
                    "csch": "csch_FIXCONVERTER",
                    "coth": "coth_FIXCONVERTER",
                    "arcsin": "asin",
                    "arccos": "acos",
                    "arctan": "atan_CHECKME",
                    "arcsec": "asec_FIXCONVERTER",
                    "arccsc": "acsc_FIXCONVERTER",
                    "arccot": "acot_FIXCONVERTER",
                    "arcsinh": "asinh",
                    "arccosh": "acosh",
                    "arctanh": "atanh",
                    "arcsech": "asech_FIXCONVERTER",
                    "arccsch": "acsch_FIXCONVERTER",
                    "arccoth": "acoth_FIXCONVERTER",
                    }
                    
                    
    if rhsElement.tag == "ci":
        return (rhsElement.text.strip(), set([rhsElement.text.strip()]))
    elif rhsElement.tag == "cn":
        text = rhsElement.text.strip()
        sep = rhsElement.find("sep")
        if sep != None:
            text += 'e'+sep.tail.strip()
        if rhsElement.get("units",""):
            if inlineUnits:
                text += "{"+printUnit(units[rhsElement.get("units")])+"}"
        return (text, set())
    elif rhsElement.tag == "true":
        return ("true", set())
    elif rhsElement.tag == "false":
        return ("false", set())
    elif rhsElement.tag == "pi":
        return ("3.1415926535897932", set())
    elif rhsElement.tag == "exponentiale":
        return ("2.7182818284590452", set())
    elif rhsElement.tag == "infinity":
        return ("M_INFINITY",set())
    elif rhsElement.tag == "apply":
        op = rhsElement[0].tag
        if op == "log":
            assert(rhsElement[1].tag == "logbase")
            (base, baseDepend) = parseRhs(rhsElement[1][0])
            (arg, argDepend) = parseRhs(rhsElement[2])
            try: 
                if float(base) == 10:
                    return ("log10("+arg+")",baseDepend|argDepend)
            except ValueError:
                pass
            return ("log("+arg+")/log("+base+")", baseDepend|argDepend)

        if op == "diff":
            assert(rhsElement.find("diff") != None)
            assert(rhsElement.find("bvar/ci") != None)
            assert(rhsElement.find("ci") != None)
            var = rhsElement.find("ci").text.strip()+".diff"
            return (var, set([var]))
        (texts, dependencies) = unzip([parseRhs(item) for item in rhsElement[1:]])
        allDependencies = set()
        for depend in dependencies:
            allDependencies |= depend
        if op == "plus":
            if len(texts)==1:
                return (texts[0], allDependencies)
            else:
                return ("("+"+".join(texts)+")",allDependencies)
        elif op == "minus":
            if len(texts)==1:
                return ("-"+texts[0],allDependencies)
            else:
                return ("("+"-".join(texts)+")",allDependencies)
        elif op == "times":
            return ("*".join(texts), allDependencies)
        elif op == "divide":
            subs = []
            subs.append(texts[0])
            subs += [protectDivide(text) for text in texts[1:]]
            return ('/'.join(subs),allDependencies)
        elif op in boolOps:
            return ((" "+ boolOps[op]+" ").join(texts), allDependencies)
        elif op == "not":
            assert(len(texts) == 1)
            return ("!("+texts[0]+")", allDependencies)
        elif op in functionOps:
            return (functionOps[op]+"("+",".join(texts)+")",
                    allDependencies)
        else:
            print(rhsElement)
            print(op)
            assert(False)
    elif rhsElement.tag == "piecewise":
        #grab all the pieces
        ifClauses = []
        elseClause = ("NOELSE_CLAUSE_CHECKME", set())
        for child in rhsElement:
            (valueClause, valueDepend) = parseRhs(child[0])
            if child.tag == "piece":
                (ifValue, ifDepend) = parseRhs(child[1])
                ifClauses.append((ifValue, valueClause, valueDepend|ifDepend))
            else:
                elseClause = (valueClause, valueDepend)
        #now, print it all out.
        allText = elseClause[0]
        allDepend = elseClause[1].copy()
        while ifClauses:
            thisIfClause = ifClauses.pop()
            allText = "(("+thisIfClause[0]+") ? "+thisIfClause[1]+" : "+allText+")"
            allDepend |= thisIfClause[2]
        return (allText, allDepend)
    else:
        print(rhsElement.tag)
        assert(False)


def parseComponentRef(thisElement, parentComponent, componentMap):
    thisComponent = componentMap[thisElement.get("component")]
    parentComponent.addSubComponent(thisComponent)
    for childElement in thisElement.findall("component_ref"):
        parseComponentRef(childElement, thisComponent, componentMap)

class Component:
    def __init__(self, root):
        self.name = root.get("name")
        self.variables = set()
        self.inputs = set()
        self.outputs = set()
        
        self.eqns = {}
        diffeqns = {}

        parseUnits(units, root.findall("units"))

        self.varToUnit = {}
        #loop through all variables
        for var in root.findall("variable"):
            name = var.get("name")
            if var.get("public_interface", "none") == "in":
                self.inputs.add(name)
            elif var.get("public_interface", "none") == "out":
                self.outputs.add(name)

            if var.get("initial_value"):
                self.eqns[name] = Equation(name, var.get("initial_value"), units[var.get("units")])
            self.varToUnit[name] = units[var.get("units")]
        for mathElement in root.findall("math"):
            for eqnElement in mathElement.findall("apply"):
                eqn = parseEquation(eqnElement, self.varToUnit)
                if eqn.isDiff:
                    diffeqns[eqn.lhs] = eqn
                else:
                    self.eqns[eqn.lhs] = eqn

        for var in diffeqns.keys():
            if var in self.eqns:
                diffeqns[var].addInitialCondition(self.eqns[var])
                self.eqns[var] = diffeqns[var]
            else:
                diffeqns[var].addInitialCondition(
                    Equation(var, "NO INITIAL CONDITION", self.varToUnit[var], set(), False))
                self.eqns[var] = diffeqns[var]

        self.subComponents = {}

    def lookupUnit(self, var):
        if var in self.varToUnit:
            return printUnit(self.varToUnit[var])
        for subComponent in self.subComponents.values():
            possibleAnswer = subComponent.lookupUnit(var)
            if possibleAnswer:
                return possibleAnswer
        return None
        
    def addSubComponent(self, subComponent):
        self.subComponents[subComponent.name] = subComponent

    def toCode(self, out, declared=set()):
        declared = declared.copy()

        out("subsystem %s {", self.name)
        out.inc()

        #find all the variables I need to specify for me.
        allDecl = set()
        for component in self.subComponents.values():
            allDecl |= component.outputs
            allDecl |= component.inputs
        allDecl |= self.localDiffvars()
        allDecl |= self.inputs
        allDecl |= self.outputs
        
        allDiffvars = set()
        for componentName in order(self.subComponents.keys()):
            subComponent = self.subComponents[componentName]
            allDiffvars |= subComponent.outputDiffvars()
        allDiffvars |= self.localDiffvars()
        
        allConstants = set()
        for componentName in order(self.subComponents.keys()):
            subComponent = self.subComponents[componentName]
            allConstants |= subComponent.outputConstants()
        allConstants |= self.localConstants()

        #separate those declared variables into computed and state
        if True:
            allDecl -= declared
        localDiffvars = self.localDiffvars() - self.outputs
        for var in order(self.outputs & allDiffvars):
            out("provides diffvar %s {%s};" , var, self.lookupUnit(var))
        for var in order(self.outputs & allConstants):
            out("provides %s {%s};", var, self.lookupUnit(var))
        for var in order(self.outputs - allDiffvars - allConstants):
            out("provides %s {%s};", var, self.lookupUnit(var))
        for var in order((allDecl - allConstants) - self.outputs - self.inputs - localDiffvars):
            out("shared %s {%s};", var, self.lookupUnit(var))
        for var in order((allDecl & allConstants) - self.outputs - self.inputs - localDiffvars):
            out("shared %s {%s};", var, self.lookupUnit(var))
        for var in order(localDiffvars):
            out("diffvar %s {%s};", var, self.lookupUnit(var))
            
        declared |= allDecl

        #find all the differential variables
        invoked = set()
        good = self.inputs.copy() | allDiffvars | (allConstants-self.localConstants())
        defined = set()
        while True:
            newFront = set()
            #go through all the invocations and see if they are good to go
            for componentName in order(set(self.subComponents.keys())-invoked):
                subComponent = self.subComponents[componentName]
                if subComponent.inputs <= good:
                    subComponent.toCode(out,declared)
                    newFront |= subComponent.outputs
                    invoked.add(componentName)
            #go through all the equations looking for things with met dependencies.
            for var in order((set(self.eqns.keys()) - defined)):
                eqn = self.eqns[var]
                #print var
                #print eqn.dependencies
                #print good
                #print defined
                if eqn.dependencies <= good: 
                    newFront.add(var)
                    eqn.toCode(out)
                    if eqn.isDiff:
                        newFront.add(var+".diff")
            if not newFront:
                break
            else:
                good |= newFront
                defined |= newFront

        if not self.outputs <= defined or invoked < set(self.subComponents.keys()):
            #out("//WARNING, circular dependency detected, dumping the rest.  You'll have to fix this manually.")
            for var in order((set(self.eqns.keys()) - defined)):
                self.eqns[var].toCode(out)
            for componentName in order(set(self.subComponents.keys())-invoked):
                subComponent = self.subComponents[componentName]
                subComponent.toCode(out,declared)
            debug = False
            if debug:
                print(list(self.eqns.keys()))
                print(self.outputs)
                print(self.outputs - defined)
                print(good)
                print("-----")
                print(invoked)
                print(set(self.subComponents.keys()) - invoked)
                for componentName in set(self.subComponents.keys()) - invoked:
                    print("=====")
                    component = self.subComponents[componentName]
                    print(component.name)
                    print(component.inputs - good)
                    print(component.outputs - good)
                assert(False)

        out.dec()
        out("}")

    def outputDiffvars(self):
        myvars = set()
        for (name,component) in self.subComponents.items():
            myvars |= component.outputDiffvars()
        myvars |= self.localDiffvars()
        myvars &= self.outputs
        return myvars

    def outputConstants(self):
        myvars = set()
        for (name,component) in self.subComponents.items():
            myvars |= component.outputConstants()
        myvars |= self.localConstants()
        myvars &= self.outputs
        return myvars        
    
    def localDiffvars(self):
        return set([ eqn.lhs for eqn in self.eqns.values() if eqn.isDiff])

    def localConstants(self):
        return set([ eqn.lhs for eqn in self.eqns.values() if not eqn.dependencies])


class RootComponent(Component):
    def __init__(self, name, inputs, outputs, varToUnit):
        self.inputs = inputs;
        self.outputs = outputs;
        self.varToUnit = varToUnit;
        self.eqns = {}
        self.subComponents = {}
        self.name = name

def printUnit(unit):
    numeratorList = []
    denominatorList = []
    for (name,subExp) in unit.items():
        if subExp == 1:
            numeratorList.append(name)
        elif subExp == -1:
            denominatorList.append(name)
        elif subExp >= 2:
            numeratorList.append("%s^%d" % (name, subExp))
        elif subExp <= -2:
            denominatorList.append("%s^%d" % (name, -subExp))
    numerator = "*".join(numeratorList)
    denominator = "/".join(denominatorList)
    result = numerator
    if numerator:
        result = numerator
    else:
        result = "1"
    if denominator:
        result += "/"+denominator
    return result

def parseUnits(units, unitElements):
    cellmlPrefixes = {
        "yotta" : "yotta",
        "zetta" : "zeta",
        "exa" : 'E',
        "peta" : 'P',
        "tera" : 'T',
        "giga" : 'G',
        "mega" : 'M',
        "kilo" : 'k',
        "hecto" : 'h',
        "deka" : 'da',
        "deci" : 'd',
        "centi" : 'c',
        "milli" : 'm',
        "micro" : 'u',
        "nano" : 'n',
        "pico" : 'p',
        "femto" : 'f',
        "atto" : 'a',
        "zepto" : 'zepto',
        "yocto" : 'yocto',
        }

    #parse all the base units
    baseUnitElements = set([unitElement for unitElement in unitElements
                            if unitElement.get("base_unit","no") == "yes"])
    for unitElement in baseUnitElements:
        units[unitElement.get("name")] = { unitElement.get("name"), 1}
        
    moreUnitsToParse = True
    while moreUnitsToParse:
        moreUnitsToParse = False
        for unitElement in set(unitElements)-baseUnitElements:
            if unitElement.get("name") in units:
                continue
            #are all the subunits for this guy defined?
            undefinedSubunits = False
            for subUnitElement in unitElement.findall("unit"):
                if subUnitElement.get("units") not in units:
                    undefinedSubunits = True
                    break
            if undefinedSubunits:
                moreUnitsToParse = True
                continue

            #check for stuff we don't want to handle
            for subUnitElement in unitElement.findall("unit"):
                assert(float(subUnitElement.get("offset","0")) == 0)
                assert(float(subUnitElement.get("scale","1")) == 1)
            
            #all the subunits are defined, so we're good to go.
            thisUnit = {}            
            for subUnitElement in unitElement.findall("unit"):
                exponent = int(subUnitElement.get("exponent", "1"))
                prefix = cellmlPrefixes.get(subUnitElement.get("prefix", ""),"")
                subUnit = units[subUnitElement.get("units")]
                appliedPrefix = False
                for (name,subExp) in subUnit.items():
                    if subExp > 0:
                        thisName = prefix+name
                        appliedPrefix = True
                    else:
                        thisName = name
                    thisUnit[thisName] = thisUnit.get(thisName,0) + subExp*exponent
                assert(not prefix or appliedPrefix)
            units[unitElement.get("name")] = thisUnit

units = {
    "ampere" : {"A" : 1},
    "farad" : {"F" : 1},
    "katal" : {"katal" : 1},
    "lux" : {"lux" : 1},
    "pascal" : {"Pa" : 1},
    "tesla" : {"tesla": 1},
    "becquerel" : {"becquerel": 1},
    "gram" : {"g": 1},
    "kelvin" : {"K": 1},
    "meter" : {"m": 1},
    "radian" : {"radian": 1},
    "volt" : {"V": 1},
    "candela" : {"candela": 1},
    "gray" : {"gray": 1},
    "kilogram" : {"kg": 1},
    "metre" : {"m": 1},
    "second" : {"s": 1},
    "watt" : {"watt": 1},
    "celsius" : {"celsius": 1},
    "henry" : {"henry": 1},
    "liter" : {"L": 1},
    "mole" : {"mol": 1},
    "siemens" : {"S": 1},
    "weber" : {"weber": 1},
    "coulomb" : {"C": 1},
    "hertz" : {"Hz": 1},
    "litre" : {"L": 1},
    "newton" : {"N": 1},
    "sievert" : {"sievert": 1},
    "dimensionless" : {},
    "joule" : {"J": 1},
    "lumen" : {"lumen": 1},
    "ohm" : {"ohm": 1},
    "steradian" : {"steradian": 1}
    }

def main():
    import sys;
    import os;

    import argparse
    parser=argparse.ArgumentParser(description="Convert CellML files to Melodee.")
    parser.add_argument("--inline-units", "-u",
                        action="store_true",
                        help="Print out the units of all constants in the CellML file",
                        dest="inlineUnits")
    parser.add_argument("filename")
    result = parser.parse_args()
    cellmlFilename = result.filename
    global inlineUnits
    inlineUnits=result.inlineUnits

    root = stripNamespaces(ET.parse(cellmlFilename))

    parseUnits(units, root.findall("units"))
    
    #get all the components from memory
    components = {}
    for element in root.findall("component"):
        component = Component(element)
        components[component.name] = component
        
    #get the calling hierarchy.
    for group in root.findall("group"):
        assert(group[0].tag == "relationship_ref")
        rref = group[0];
        if rref.get("relationship", "none") == "encapsulation":
            for parentElement in group.findall("component_ref"):
                parentComponent = components[parentElement.get("component")]
                for childElement in parentElement.findall("component_ref"):
                    parseComponentRef(childElement,parentComponent,components)

    #anything that doesn't have a parent is parented by the root comp.
    #make a list of everyone's children.
    allChildren = set()
    for component in components.values():
        allChildren |= set(component.subComponents.values())

    environmentNames = set(["environment","Environment","Myofilaments"])
    rootInputs = set()
    rootOutputs = set()
    rootInputConstants = {}
    rootVarToUnit = {}
    for environmentName in environmentNames:
        if environmentName in components:
            for var in components[environmentName].outputConstants():
                rootInputConstants[var] = components[environmentName].eqns[var]
            rootInputs |= components[environmentName].outputs-set(rootInputConstants.keys())
            rootOutputs |= components[environmentName].inputs
            for (name, unit) in components[environmentName].varToUnit.items():
                rootVarToUnit[name] = unit
    (rootName, rootExt) = os.path.splitext(os.path.basename(sys.argv[1]))
    root = RootComponent(rootName, rootInputs, rootOutputs, rootVarToUnit)
    root.eqns = rootInputConstants
    for component in components.values():
        if component.name not in environmentNames and component not in allChildren:
            root.addSubComponent(component)
    
    out = Indenter(sys.stdout)
    for var in root.inputs:
        if var == "time" or var == "t":
            out("integrate %s {%s};", var, printUnit(rootVarToUnit[var]))
        else:
            out("shared ephemeral %s {%s};", var, printUnit(rootVarToUnit[var]))
    root.toCode(out)
