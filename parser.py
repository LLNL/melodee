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

import ply.lex as lex
import ply.yacc as yacc
import sympy
import units

class SSA(dict):
    def __setitem__(self, key, value):
        if key in self:
            raise MultipleSymbolAssignment(key)
        self[key]=value

class InstructionList(list):
    def addInstruction(self, instruction):
        self.append(instruction)

class IfInstruction:
    def __init__(self, ifVar, thenInstructions, elseInstructions, choiceInstructions):
        self.ifVar = ifVar
        self.thenInstructions = thenInstructions
        self.elseInstructions = elseInstructions
        self.choiceInstructions = choiceInstructions

class Choice:
    def __init__(self, ifVar, thenVar, elseVar, unit):
        self.ifVar = ifVar
        self.thenVar = thenVar
        self.elseVar = elseVar
        self.astUnit = unit

class AST:
    def __init__(self, sympyExpr, unit=None):
        self.sympy = sympyExpr
        self.astUnit = unit
    def __repr__(self):
        return 'AST(sympify(' + repr(self.sympy) +'), '+ repr(self.astUnit) +')'

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
        except KeyError,e:
            return False
        return True
    def getUnit(self, name):
        if name in self.units:
            return self.units[name]
        elif self.parent != None:
            return self.parents.getUnit(name)
        else:
            raise KeyError(name)
    def hasUnit(self,name):
        try:
            self.getUnit(name)
        except KeyError,e:
            return False
        return True
    def setSymbol(self, name, symbol):
        self.symbols[name] = symbol
    def setUnit(self, name, unit):
        self.units[name] = unit
    def addInstruction(self, inst):
        self.instructions.addInstruction(inst)

class SharedScope(Scope):
    def __init__(self):
        Scope.__init__(self, None)
        self.junctions = {}
    def hasJunction(self, name):
        return name in self.junctions
    def getJunction(self, name):
        return self.junctions[name]
    def setJunction(self, name, junction):
        self.junctions[name] = junction

class Junction:
    def __init__(self,unit):
        self.rawUnit = unit
    
class Subsystem:
    def __init__(self, name):
        self.name = name
        self.ssa = SSA()
        self.scope = SharedScope()
        self.inputs = set()
        self.outputs = set()
        self.diffvars = set()
        self.accums = set()
        self.params = set()
        self.time = None
        self.paramDefault = {}
        

class Parser:
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
        self.subsystemStack = []
        self.scopeStack = [SharedScope()]
        self.timeVar = None
        self.timeUnit = None
        
    def currentSubsystem(self):
        return self.subsystemStack[-1]
    def currentScope(self):
        return self.scopeStack[-1]
    def currentInstructions(self):
        return self.currentScope().instructions
    def pushScope(self):
        self.scopeStack[-1] = self.scopeStack[-1].makeChild()
        return self.scopeStack[-1]
    def popScope(self):
        retval = self.scopeStack[-1]
        self.scopeStack[-1] = self.scopeStack[-1].getParent()
        return retval
    def readAccessVar(self, name):
        #if the var is an accum
        if 0:
            #syntax error: can't read from an accum
            pass
        #if the var has a symbol already
        if 0:
            # if this is a parameter that has never been read before
            if 0:
                #mark this parameter as frozen
                pass
            #get the symbol
            #if this symbol has a unit defined
            if 0:
                #set the unit to this symbol
                pass
            else:
                #lookup the unit for the symbol if it's been assigned to
                pass
        else:
            #see if this is a shared var
            if 0:
                #make a symbol for this guy, mark it as an input
                #repeat the lookup command now that we have a symbol
                pass
            #check for a timevar
            elif 0:
                #make a symbol for this guy, mark it as time
                #repeat the lookup for this command now that we have a symbol
                pass
            #check for an enumeration
            elif 0:
                #make a symbol for this guy, mark it as frozen
                #repeat the lookup for this command now that we have a symbol
                pass
            else:
                #syntax error
                pass
            
        
        #search for a local variable
        if self.currentScope().hasSymbol(name):
            symbol = self.currentScope().getSymbol(name)
            rawUnit = None
            if self.currentScope().hasUnit(name):
                rawUnit = self.currentScope().getUnit(name)
            else:
                rawUnit = self.currentSubsystem().ssa[symbol].astUnit.rawUnit
            return AST(symbol, ASTUnit(rawUnit, explicit=False))
        else:
            #check previous scopes for junctions
            assert(len(self.scopeStack) >= 2)
            (exists, junction) = self.searchForJunction(name)
            if exists:
                symbol = sympy.symbols(name)                
                self.currentSubsystem().scope.setSymbol(name, symbol)
                self.currentSubsystem().scope.setUnit(name, junction.rawUnit)
                return self.readAccessVar(name)
            else:
                raise SyntaxError("Variable '"+name+"' used but not defined.")
    def checkDeclarable(self, var):
        #if the variable is already type defined in this subsystem
        ss = self.currentSubsystem()
        if var in ss.inputs:
            raise SyntaxError("'%s' already used as a shared variable in this subsystem." % var)
        if var in ss.diffvars:
            raise SyntaxError("'%s' already defined as a differential variable in this subsystem." % var)
        if var in ss.accums:
            raise SyntaxError("'%s' already defined as an accumulation output for this subsystem." % var)
        if var in ss.params:
            raise SyntaxError("'%s' already defined as a parameter for this subsystem." % var)
        if var in ss.outputs:
            raise SyntaxError("'%s' already an output for this subsystem." % var)
        if var == ss.time:
            raise SyntaxError("'%s' already used as the integration variable for this subsystem." % var) 
        #if the variable has already been assigned in this subsystem
        if ss.scope.hasSymbol(var) and ss.scope.getSymbol(var) in ss.ssa:
            raise SyntaxError("'%s' has already been assigned to, it can't now be declared." % var)
        #if the variable is already unit defined in this subsystem
        if ss.scope.hasUnit(var):
            raise SyntaxError("Units previously defined for '%s'" % var) 
    
    def markProvides(self, varname, unitOpt):
        self.checkDeclarable(varname)
        #if this variable connects to another variable
        (exists, junction) = self.searchForJunction(varname)
        if exists:
            #get the unit for that junction
            junctionUnit = junction.rawUnit
            #if a unit was defined and that unit does not match the junction unit
            if unitOpt != None and unitOpt != junctionUnit:
                raise SyntaxError("Provided units for '%s' don't match the shared units")
            unit = junctionUnit
        else:
            #if a unit was not provided
            if unitOpt == None:
                raise SyntaxError("Must have a unit for provides variable '%s'." % varname)
            unit = unitOpt
            junction = Junction(unit)
            #make a new junction with a unit
            self.currentSubsystem().scope.junction[varname] = junction
        self.currentScope().setUnit(varname, unit)
        self.currentSubsystem().outputs.add(varname)
        #return the junction
        return junction

    def searchForJunction(self, name):
        assert(len(self.subsystemStack) >= 1)
        for ii in range(len(self.scopeStack)-2,0,-1):
            if self.scopeStack[ii].hasJunction(name):
                return (True, self.scopeStack[ii].getJunction(name))
        return (False, None)

    def processIfCondition(self, ifVar, thenBody, elseBody):

        choiceInstructions = InstructionList()
        thenSymbols = thenBody[1]
        elseSymbols = elseBody[1]

        for (key,(thenSympy, thenUnit)) in thenSymbols.local().items():
            if key in elseSymbols: #no local here, to catch both then and else
                (elseSympy, elseUnit) = elseSymbols[key]
                choice = Choice(ifVar,
                                thenSympy,
                                elseSympy,
                                ASTUnit.null() # FIXME?
                                )
                self.astToSymbol(key,choice)
        
        for (key, (elseSympy,elseUnit)) in elseSymbols.local().items():
            if key in thenSymbols.local(): #local here so we don't violate ssa
                (thenSympy, thenUnit) = thenSymbols[key]
                choice = Choice(ifVar,
                                thenSympy,
                                elseSympy,
                                ASTUnit.null() # FIXME?
                                )
                self.astToSymbol(key,choice)

        thenInstructions = thenBody[0]
        elseInstructions = elseBody[0]
        iif = IfInstruction(ifVar,thenInstructions,elseInstructions,choiceInstructions)
        self.currentInstructions().addInstruction(iif)
        #fixme fill in the current symbol table with the choice stuff?

        
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
            raise SyntaxError("Can't cast "+str(ast.astUnit)+" to unit of "+str(castRaw)+")")

    def checkExactUnits(self, lunit, runit):
        if lunit.isNull() or runit.isNull():
            if lunit.explicit or runit.explicit:
                raise SyntaxError("Explicit units cannot combine with null units!")
            else:
                return ASTUnit.null()
        else:
            # now we have something with units!  Check if they match.
            resultIsExplicit = lunit.explicit or runit.explicit
            if lunit.rawUnit == runit.rawUnit:
                return ASTUnit(lunit.rawUnit, resultIsExplicit)
            else:
                if resultIsExplicit:
                    raise SyntaxError("Units don't match, maybe you need a conversion?")
                else:
                    return ASTUnit.null()
    def convertUnitTo(self, expr, newUnit):
        if expr.astUnit.isNull():
            raise SyntaxError("Can't convert a null unit!")
        elif not expr.astUnit.rawUnit.isCompatibleWith(newUnit):
            raise SyntaxError("Incompatible unit conversion requested.")
        else:
            factor = expr.astUnit.rawUnit.convertTo(newUnit)
            return AST(sympy.Mul(factor,expr.sympy), ASTUnit(newUnit, explicit=False))
    
    def parse(self, s):
        return self.parser.parse(s)

    def astToVar(self, var, ast):
        self.currentSubsystem().ssa[var] = ast
        self.currentInstructions().addInstruction(var)
        return (var, ast.astUnit)
        
    def astToTemp(self, ast):
        return self.astToVar(self.newTempVar(),ast)
    def astToSymbol(self, name, ast):
        (var, astUnit) = self.astToVar(sympy.symbols(name),ast)
        self.currentScope().setSymbol(name, var)
        return (var, astUnit)
        
    
    t_ignore = " \t\r"
                                    
    def t_newline(self, t):
        r'\n+'
        t.lexer.lineno += len(t.value)

    def t_error(self, t):
        print("Illegal character '%s'" % t.value[0])
        t.lexer.skip(1)

    reserved = {
        "if" : "IF",
        "else" : "ELSE",
        "elseif" : "ELSEIF",
        "bool" : "BOOL",
        "true" : "TRUE",
        "false" : "FALSE",
        "param" : "PARAM",
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
        #"from" : "FROM",

        #"unit" : "UNIT",
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
        "PLUSPERCENTEQ",
        "MINUSPERCENTEQ",
    ) + tuple(reserved.values())

    literals = "(){}<>!;?:+-/*^=.,~"

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

    t_LEQ = r'<='
    t_GEQ = r'>='
    t_NEQ = r'!='
    t_BOOLEQ = r'=='
    
    t_PLUSEQ = r'\+='
    t_MINUSEQ = r'-='
    t_TIMESEQ = r'\*='
    t_DIVIDEEQ = r'/='
    t_EXPONEQ = r'^='
    t_PLUSPERCENTEQ = r'\+%='
    t_MINUSPERCENTEQ = r'-%='
    
    t_ignore_CPP_COMMENT = r'//.*?\n'
    t_ignore_C_COMMENT = r'/\*.*?\*/'
    t_ignore_PYTHON_COMMENT = r'\#.*?\n'
    t_ignore_MATLAB_COMMENT = r'%.*?\n'

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
        ('right', 'UMINUS'),
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
                             | subSystemDefinition
                             | integrateStatement
        '''
        pass

    def p_integrateStatement(self, p):
        '''integrateStatement : INTEGRATE var unitDef ';'
        '''
        self.timeVar = p[2]
        self.timeUnit = p[3]
        
    def p_sharedStatement_unit(self, p):
        '''sharedStatement : SHARED var unitDef ';'
        '''
        self.currentScope().setJunction(p[2], Junction(p[3]))

    def p_nameList_term(self, p):
        '''nameList : NAME'''
        pass
    def p_nameList_shift(self, p):
        '''nameList : NAME ',' nameList'''
        pass

    def p_subSystemDefinition(self, p):
        '''subSystemDefinition : subSystemBegin subSystemStatementsOpt '}' '''
        p[0] = self.subsystemStack.pop()
    def p_subSystemBegin(self, p):
        '''subSystemBegin : SUBSYSTEM NAME '{' '''
        thisName = p[2]
        if self.subsystemStack:
            thisName = self.currentSubsystem().name + "." + thisName
        thisSubsystem = Subsystem(thisName)            
        self.subsystemStack.append(Subsystem(thisName))
        self.scopeStack.append(self.currentSubsystem().scope)
        
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
        pass

    def p_subSystemStatement_accum(self, p):
        '''subSystemStatement : PROVIDES ACCUM var unitOpt accumDefOpt ';' '''
        junction = self.markProvides(p[3],p[4])
        #mark that junction as an accum junction
        #FIXME
        #mark the variable as an accum
        self.currentSubsystem().accums.add(p[3])
        #if there's a definition, process it.
        #FIXME
        pass
        
    def p_accumDefOpt_empty(self, p):
        '''accumDefOpt : empty'''
        p[0] = None
    def p_accumDefOpt_def(self, p):
        '''accumDefOpt : PLUSEQ  realExpr
                       | MINUSEQ realExpr
        '''
        rhs = p[2]
        if p[1] == '-=':
            rhs = AST(sympy.Mul(sympy.Integer(-1),rhs.sympy, rhs.astUnit))
        p[0] = rhs

    def p_subSystemStatement_param(self, p):
        '''subSystemStatement : PROVIDES PARAM var unitOpt assignOpt ';' 
                              |          PARAM var unitDef assignOpt ';' 
        '''
        #if this is a provides
        if p[1] == "provides":
            (var, unit, assignOpt) = (p[3],p[4],p[5])
            self.markProvides(var,unit)
        else:
            (var, unit, assignOpt) = (p[2],p[3],p[4])
            self.checkDeclarable(var)
            self.currentScope().setUnit(var,unit)
        #mark the variable as a parameter.
        self.currentSubsystem().params.add(var)
        #if there's a definition, process it.
        #FIXME
        pass
    def p_assignOpt_empty(self, p):
        '''assignOpt : empty'''
        p[0] = None
    def p_assignOpt_def(self, p):
        '''assignOpt : '=' realExpr'''
        p[0] = p[1]

    def p_subSystemStatement_diffvar(self, p):
        '''subSystemStatement : PROVIDES DIFFVAR var unitOpt ';'
                              |          DIFFVAR var unitDef ';' 
        '''
        if self.timeUnit is None:
            raise SyntaxError("Must include an 'integrate' statement before declaring subsystems.")
        #if this is a provides
        if p[1] == "provides":
            (var, unit) = (p[3],p[4])
            self.markProvides(var,unit)
        else:
            (var, unit) = (p[2],p[3])
            self.checkDeclarable(var)
            self.currentScope().setUnit(var,unit)
        #mark the variable as a diffvar
        self.currentSubsystem().diffvars.add(var)
        #mark the units for the initial and diff conditions
        rawUnit = self.currentScope().getUnit(var)
        self.currentScope().setUnit(var+".init", rawUnit)
        self.currentScope().setUnit(var+".diff", rawUnit/self.timeUnit)

    def p_subSystemStatement_output(self, p):
        '''subSystemStatement : PROVIDES var unitOpt assignOpt ';' '''
        self.markProvides(p[2],p[3])
        #if there's a definition, process it.
        #FIXME
        pass

    
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
        '''conditionalStatement : varDiffOpt unitOpt '=' realExpr ';'
                                | varDiffOpt unitOpt PLUSEQ realExpr ';'
                                | varDiffOpt unitOpt MINUSEQ realExpr ';'
                                | varDiffOpt unitOpt TIMESEQ realExpr ';'
                                | varDiffOpt unitOpt DIVIDEEQ realExpr ';'
                                | varDiffOpt unitOpt EXPONEQ realExpr ';'
        '''
        pass

    ########################################################################
    def p_subSystemStatement_if(self, p):
        '''conditionalStatement : ifStatement'''
        pass

    def p_ifStatement(self, p):
        '''ifStatement : initialIfCond thenBody elseOpt'''
        self.processIfCondition(p[1],p[2],p[3])
    def p_initialIfCond(self,p):
        '''initialIfCond : IF '(' realExprToTemp ')' '''
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
        '''elseIfCond : ELSEIF '(' realExprToTemp ')' '''
        p[0] = p[4]
    def p_elseOpt_continue(self, p):
        '''elseOpt : ifScopeBegin elseIfCond thenBody elseOpt'''
        self.processIfCondition(p[2],p[3],p[4])
        p[0] = self.popScope()

    def p_ifScopeBegin(self, p):
        '''ifScopeBegin : empty'''
        self.pushScope()

    ##################################################################


    def p_subSystemStatement_use(self, p):
        '''subSystemStatement : useStatement'''
        pass
    def p_useStatement(self, p):
        '''useStatement : USE useList '{' useBlockStatementList '}' '''
        pass
    def p_useList_default(self, p):
        '''useList : speccedName'''
        p[0] = (p[1],set())
    def p_useList_subtract(self, p):
        '''useList : useList '-' speccedName'''
        p[0] = (p[1][0], p[1][1] | set([p[2]]))
    def p_speccedName_default(self, p):
        '''speccedName : NAME'''
        p[0] = p[1]
    def p_specceedName_specified(self, p):
        '''speccedName : speccedName '.' NAME'''
        p[0] = p[1] + '.' + p[2]
    def p_useBlockStatementList_single(self, p):
        '''useBlockStatementList : useBlockStatement'''
        pass
    def p_useBlockStatementList_mult(self, p):
        '''useBlockStatementList : useBlockStatementList useBlockStatement'''
        pass
    def p_useBlockStatement_explicitBindLeft(self, p):
        '''useBlockStatement : NAME '~' '.' speccedName ';' '''
        pass
    def p_useBlockStatement_explicitBindRight(self,p):
        '''useBlockStatement : '.' speccedName '~' NAME ';' '''
        pass
    def p_useBlockStatement_simpleBind(self,p):
        '''useBlockStatement : '~' '.' NAME ';' '''
        pass
    def p_useBlockStatement_implicitBind(self,p):
        '''useBlockStatement : '~' '.' ';' '''
        pass
    def p_useBlockStatement_parameterSet(self, p):
        '''useBlockStatement : '.' speccedName useOp useExpr ';' '''
        pass
    def p_useOp(self, p):
        '''useOp : '='
                 | PLUSEQ
                 | MINUSEQ
                 | PLUSPERCENTEQ
                 | MINUSPERCENTEQ
                 | TIMESEQ
                 | DIVIDEEQ
        '''
        pass
    def p_useExpr_literal(self,p):
        '''useExpr : numberLiteralPlusOpt '''
        pass
    def p_useExpr_boolop(self,p):
        '''useExpr : useExpr '+' useExpr
                   | useExpr '-' useExpr
                   | useExpr '*' useExpr
                   | useExpr '/' useExpr
                   | useExpr '^' useExpr
        '''
        pass
    def p_useExpr_uminus(self,p):
        '''useExpr : '-' useExpr %prec UMINUS'''
        pass

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
        p[0] = self.si.get(unitName)
        #FIXME, all the unitNames need to be tagged with this unit.
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
        '''ternaryExpr : ternaryIf ternaryThen ternaryElse'''
        resultUnit = self.checkExactUnits(p[1][1],p[2][1])
        ast = Choice(p[0],p[1][0],p[2][0],resultUnit)

        self.instructionStack.append(InstructionList())
        (var,unit) = self.astToTemp(ast)
        choiceInstructions = self.instructionStack.pop()

        self.currentInstructions().addInstruction(IfInstruction(p[0],p[1][2], p[2][2]),choiceInstructions)
        p[0] = AST(var,unit)
        
    def p_ternaryIf(self, p):
        '''ternaryIf : orExpr '?' '''
        (var, unit) = self.astToTemp(p[1])
        p[0] = var

    def p_realExprToTemp(self,p):
        '''realExprToTemp : realExpr'''
        p[0] = self.astToTemp(p[1])
        
    def p_ternaryThen(self, p):
        '''ternaryThen : pushNewInstructionList realExprToTemp'''
        (var,unit) = p[2]
        instructions = self.instructionStack.pop()
        p[0] = (var,unit,instructions)
    def p_ternaryElse(self, p):
        '''ternaryElse : ':' pushNewInstructionList realExprToTemp'''
        (var,unit) = p[3]
        instructions = self.instructionStack.pop()
        p[0] = (var,unit,instructions)
        
    def p_pushNewInstructionList(self, p):
        '''pushNewInstructionList : empty'''
        self.instructionStack.append(InstructionList())
        
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
        p[0] = AST(sympy.Mul(lhs.sympy,rhs.sympy), lhs.astUnit*rhs.astUnit)

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
        p[0] = AST(sympy.getattr(p[1])(*[x.sympy for x in p[3]]), self.nodim())
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

    def p_error(self,p):
        if p:
            print "SyntaxError on line number "+str(self.lexer.lineno)+" at token", p.type, p.value
            self.parser.errok()


    #def p_subSystemStatement_renameStatement(self, p):
    #    '''subSystemStatement : renameStatement'''
    #    pass


if __name__=="__main__":
    p = Parser()
    data = '''
and && or || not ! 0 2.0 .3 40. 5e+6 if myID */* bljsadfj */ */
'''
    p.lexer.input(data)

    while True:
        tok = p.lexer.token()
        if not tok:
            break
        print tok

    p = Parser(start="unitExpr")
    print p.parse("mV/ms")
    print p.parse("uA/uF")

    p = Parser(start="realExpr")
    p.p_subSystemBegin("testing")
    p.currentScope().setSymbol("a", sympy.symbols("a"))
    p.currentScope().setSymbol("b", sympy.symbols("b"))
    p.currentScope().setSymbol("c", sympy.symbols("c"))
    p.currentScope().setSymbol("d", sympy.symbols("d"))
    p.currentScope().setUnit("a", p.si.get("unitless"))
    p.currentScope().setUnit("b", p.si.get("unitless"))
    p.currentScope().setUnit("c", p.si.get("ms"))
    p.currentScope().setUnit("d", p.si.get("s"))
    print p.parse("a+b/c+d")
    print p.parse("a+(b/c)+d")
    print p.parse("a+b")
    print p.parse("(a+b)")
    print p.parse("a")
    print p.parse("(a)")
    print p.parse("((a))")
    print p.parse("((a+b))/((c+d))")
    print p.parse("1 {ms}+ c")
    print p.parse("convert(1 {ms}, s)+ d {s}")
    print p.parse("a == b")

    HH = '''
integrate time {ms};
subsystem hodgkin_huxley_1952 {
   shared V {mV};
   shared Iion {uA/cm^2};
   shared E_R {mV};
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
      E_Na {mV} = (E_R+115{mV});
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
      provides E_R {mV};
      Cm {uF/cm^2} = 1;
      E_R {mV} = -75;
      V.init = -75;
      V.diff = -Iion/Cm;
   }
}
'''
    p = Parser(start="topLevelStatementsOpt")
    p.parse(HH)
