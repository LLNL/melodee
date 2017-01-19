#!/usr/bin/env python

import ply.lex as lex
import ply.yacc as yacc
import sympy
import units

class Var:
    def __init__(self, name, unit):
        self.name = name
        self.sympy = sympy.symbols(name)
        self.unit = unit

class ASTUnit:
    def __init__(self, unit, explicit):
        self.unit = unit
        self.explicit = explicit
    def isNull(self):
        return self.unit is None
    def __mul__(self, other):
        if other.isNull() or self.isNull():
            return ASTUnit.null()
        else:
            return ASTUnit(self.unit*other.unit, self.explicit or other.explicit)
    def __div__(self, other):
        if other.isNull() or self.isNull():
            return ASTUnit.null()
        else:
            return ASTUnit(self.unit/other.unit, self.explicit or other.explicit)
    def __pow__(self, number):
        if self.isNull():
            return ASTUnit.null()
        else:
            return ASTUnit(self.unit ** number, self.explicit)
    def __repr__(self):
        return 'ASTUnit('+repr(self.unit)+','+repr(self.explicit)+')'
    @staticmethod
    def null():
        return ASTUnit(None, False)

class SymbolTable:
    def __init__(self, parent=None):
        self.parent = parent
        self.localNames = {}
    def lookup(self, name):
        if name in self.localNames:
            return self.localNames[name]
        elif self.parent==None:
            return None
        else:
            return self.parent.lookup(name)
    def create(self, name, unit=None):
        v = Var(name, unit)
        self.localNames[name] = v
        return v
    
class Statement:
    def __init__(self, lhs, rhs):
        self.lhs = lhs
        self.rhs = rhs

class BasicBlock:
    def __init__(self):
        self.instructions = []
        self.conditionalExit = {}
        self.unconditionalExit = None

class AST:
    def __init__(self, sympyExpr, unit=None):
        self.sympy = sympyExpr
        self.unit = unit
    def __repr__(self):
        return 'AST(sympify(' + repr(self.sympy) +'), '+ repr(self.unit) +')'

def textToAST(text, unit):
    return AST(sympy.sympify(text), unit)

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
        self.currentTable = SymbolTable()

    def nodim(self):
        return ASTUnit(self.si.get("1"), False)

    def checkExactUnits(self, lunit, runit):
        if lunit.isNull() or runit.isNull():
            if lunit.explicit or runit.explicit:
                raise SyntaxError("Explicit units cannot combine with null units!")
            else:
                return ASTUnit.null()
        else:
            # now we have something with units!  Check if they match.
            resultIsExplicit = lunit.explicit or runit.explicit
            if lunit.unit == runit.unit:
                return ASTUnit(lunit.unit, resultIsExplicit)
            else:
                if resultIsExplicit:
                    raise SyntaxError("Units don't match, maybe you need a conversion?")
                else:
                    return ASTUnit.null()
    def convertUnitTo(self, expr, newUnit):
        if expr.unit.isNull():
            raise SyntaxError("Can't convert a null unit!")
        elif not expr.unit.unit.isCompatibleWith(newUnit):
            raise SyntaxError("Incompatible unit conversion requested.")
        else:
            factor = expr.unit.unit.convertTo(newUnit)
            return AST(sympy.Mul(factor,expr.sympy), ASTUnit(newUnit, explicit=False))
    
    def parse(self, s):
        return self.parser.parse(s)
                                    
    
    t_ignore = " \t"
                                    
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
        "true" : "TRUE",
        "false" : "FALSE",
        "stable" : "STABLE",
        "ephemeral" : "EPHEMERAL",
        "flag" : "FLAG",
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
        #"rename" : "RENAME",
        #"to" : "TO",
        #"from" : "FROM",
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

    literals = "(){}<>!;?:+-/*^=.,"

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
        #('left', '+', '-'),
        ('left', '*', '/'),
        #('right', 'ANNOTATE_UNIT'),
        #('right', 'NOT', 'UMINUS'),
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
        '''
        pass

    def p_sharedStatement_unit(self, p):
        '''sharedStatement : SHARED STABLE var '{' unitExpr '}' ';'
                           | SHARED EPHEMERAL var '{' unitExpr '}' ';'
        '''
        pass
    def p_sharedStatement_flag(self, p):
        '''sharedStatement : SHARED flagDeclBool ';'
                           | SHARED flagDeclEnum ';'
        '''
        pass
    def p_flagDeclBool(self, p):
        '''flagDeclBool : FLAG var'''
        pass
    def p_flagDeclEnum(self, p):
        '''flagDeclEnum : FLAG var '{' nameList '}' '''
        pass

    def p_nameList_term(self, p):
        '''nameList : NAME'''
        pass
    def p_nameList_shift(self, p):
        '''nameList : NAME ',' nameList'''
        pass

    def p_subSystemDefinition(self, p):
        '''subSystemDefinition : SUBSYSTEM '{' subSystemStatementsOpt '}' '''
        pass

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

    def p_subSystemStatement_provide(self, p):
        '''subSystemStatement : providesStatement'''
        pass

    def p_providesStatement_flagDecl(self, p):
        '''providesStatement : PROVIDES flagDeclBool ';'
                             | PROVIDES flagDeclEnum ';'
        '''
        pass
    def p_provdesStatement_flagBoolDefn(self, p):
        '''providesStatement : PROVIDES flagDeclBool '=' boolLiteral ';' '''
        pass
    
    def p_providesStatement_flagEnumDefn(self, p):
        '''providesStatement : PROVIDES flagDeclEnum '=' NAME ';' '''
        pass

    def p_providesStatement_DeclSubNoUnit(self, p):
        '''providesStatement : PROVIDES ACCUM var ';'
                             | PROVIDES DIFFVAR var ';'
                             | PROVIDES PARAM var ';'
                             | PROVIDES STABLE var ';'
                             | PROVIDES EPHEMERAL var ';'
                             | PROVIDES var ';'
        '''
        pass
    def p_providesStatement_DeclSubUnit(self, p):
        '''providesStatement : PROVIDES ACCUM var '{' unitExpr '}' ';'
                             | PROVIDES DIFFVAR var '{' unitExpr '}' ';'
                             | PROVIDES PARAM var '{' unitExpr '}' ';'
                             | PROVIDES STABLE var '{' unitExpr '}' ';'
                             | PROVIDES EPHEMERAL var '{' unitExpr '}' ';'
                             | PROVIDES var '{' unitExpr '}' ';'
        '''
        pass

    def p_providesStatement_Defn(self, p):
        '''providesStatement : PROVIDES ACCUM accumDef
                             | PROVIDES PARAM assignDef
                             | PROVIDES assignDef
        '''
        pass

    def p_subSystemStatement_definition(self, p):
        '''subSystemStatement : varDef ';' '''
        pass

    def p_varDef_assign(self, p):
        '''varDef : assignDef'''
        pass

    def p_varDef_accum(self, p):
        '''varDef : accumDef'''
        pass

    def p_varDef_modAssign(self, p):
        '''varDef : var TIMESEQ realExpr ';'
                  | var DIVIDEEQ realExpr ';'
                  | var EXPONEQ realExpr ';'
        '''
        pass

    def p_assignDef_withoutUnit(self, p):
        '''assignDef : var '=' realExpr ';' '''
        pass
    def p_assignDef(self, p):
        '''assignDef : var '{' unitExpr '}' '=' realExpr ';' '''
        pass
    def p_accumDef(self, p):
        '''accumDef : var PLUSEQ realExpr ';'
                    | var MINUSEQ realExpr ';'
        '''
        pass

    def p_subSystemStatement_diffDef(self, p):
        '''subSystemStatement : diffDef'''
        pass

    def p_diffDef(self, p):
        '''diffDef : var '.' INIT '=' realExpr ';'
                   | var '.' DIFF '=' realExpr ';'
        '''
        pass

    def p_subSystemStatement_if(self, p):
        '''subSystemStatement : ifStatement'''
        pass

    def p_ifStatement_noElse(self, p):
        '''ifStatement : ifClauses'''
        pass
    def p_ifStatement_else(self, p):
        '''ifStatement : ifClauses elseClause'''
        pass
    def p_ifClauses(self, p):
        '''ifClauses : ifClause elseifClausesOpt'''
        pass

    def p_elseifClausesOpt(self, p):
        '''elseifClausesOpt : elseifClause elseifClausesOpt'''
        pass
    def p_elseIfClausesOpt_term(self, p):
        '''elseifClausesOpt : empty'''
        pass

    def p_ifClause(self, p):
        '''ifClause : IF realExpr '{' subSystemStatementsOpt '}'
        '''
        pass
    def p_elseifClause(self, p):
        '''elseifClause : ELSEIF realExpr '{' subSystemStatementsOpt '}'
        '''
        pass
    def p_elseClause(self, p):
        '''elseClause : ELSE '{' subSystemStatementsOpt '}' '''
        pass

    def p_unitExpr_literal(self, p):
        '''unitExpr : NAME'''
        p[0] = self.si.get(p[1])
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

    def powerProcess(self, x, y):
        if (not x.unit.isNull()) and y.sympy.is_constant():
            newUnit = x.unit ** float(y.sympy)
        else:
            newUnit = ASTUnit.null()
        return AST(sympy.Pow(x.sympy, y.sympy), newUnit)
    
    ############################################

    def p_realExpr_pass(self,p):
        '''realExpr : orExpr'''
        p[0] = p[1]

    def p_ternaryOp_pass(self, p):
        '''ternaryExpr : orExpr'''
        p[0] = p[1]
    def p_ternaryOp_impl(self, p):
        '''ternaryExpr : orExpr '?' realExpr ':' ternaryExpr'''
        pass

    def p_orExpr_pass(self,p):
        '''orExpr : andExpr'''
        p[0] = p[1]
    def p_orExpr_impl(self, p):
        '''orExpr : orExpr OR andExpr'''
        p[0] = AST(sympy.Or(p[1].sympy,p[3].sympy),self.nodim())

    def p_andExpr_pass(self, p):
        '''andExpr : booleqExpr'''
        p[0] = p[1]
    def p_andExpr_impl(self, p):
        '''andExpr : andExpr AND booleqExpr'''
        p[0] = AST(sympy.And(p[1].sympy,p[3].sympy),self.nodim())

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
        self.checkExactUnits(p[1].unit,p[3].unit)
        p[0] = AST(boolOp(p[1].sympy,p[3].sympy),self.nodim())

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
        self.checkExactUnits(p[1].unit,p[3].unit)
        p[0] = AST(boolOp(p[1].sympy,p[3].sympy),self.nodim())

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
            rhs = AST(sympy.Mul(sympy.Integer(-1),rhs.sympy), rhs.unit)
        p[0] = AST(sympy.Add(lhs.sympy,rhs.sympy), self.checkExactUnits(lhs.unit,rhs.unit))

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
            rhs = AST(sympy.Pow(rhs.sympy,sympy.Integer(-1)), rhs.unit ** -1)
        p[0] = AST(sympy.Mul(lhs.sympy,rhs.sympy), lhs.unit*rhs.unit)

    def p_unaryExpr_pass(self, p):
        '''unaryExpr : unitLabelExpr'''
        p[0] = p[1]
    def p_unaryExpr_uminus(self, p):
        '''unaryExpr : '-' unaryExpr'''
        p[0] = AST(sympy.Mul(sympy.Integer(-1),p[2].sympy), p[2].unit)
    def p_unaryExpr_not(self, p):
        '''unaryExpr : NOT unaryExpr
                     | '!' unaryExpr
        '''
        p[0] = AST(sympy.Not(p[2].sympy),self.nodim())

    def p_unitExpr_pass(self, p):
        '''unitLabelExpr : exponentExpr'''
        p[0] = p[1]
    def p_unitExpr_impl(self, p):
        '''unitLabelExpr : exponentExpr '{' unitExpr '}' '''
        newUnit = ASTUnit(p[3],explicit=True)
        if not p[1].unit.isNull():
            self.checkExactUnits(p[1].unit, newUnit)
        p[0] = AST(p[1].sympy, newUnit)

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
        self.checkExactUnits(p[1],self.nodim())
        p[0] = [p[1]]
    def p_funcArgList_shift(self, p):
        '''funcArgList : realExpr ',' funcArgList'''
        self.checkExactUnits(p[1].unit, self.nodim())
        p[0] = [p[1]] + p[3]
    def p_parenExpr_pass(self, p):
        '''parenExpr : primaryExpr'''
        p[0] = p[1]
    def p_parenExpr_impl(self, p):
        '''parenExpr : '(' realExpr ')' '''
        p[0] = p[2]

    def p_primaryExpr_var(self, p):
        '''primaryExpr : var'''
        v = self.currentTable.create(p[1])
        p[0] = AST(v.sympy, ASTUnit(v.unit,explicit=False))
    def p_primaryExpr_boolLiteral(self, p):
        '''primaryExpr : boolLiteral
        '''
        p[0] = textToAST(p[1], self.nodim())
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
            print "SyntaxError at token", p.type
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
    print p.parse("a+b/c+d")
    print p.parse("a+(b/c)+d")
    print p.parse("a+b")
    print p.parse("(a+b)")
    print p.parse("a")
    print p.parse("(a)")
    print p.parse("((a))")
    print p.parse("((a+b))/((c+d))")
    print p.parse("1 {ms}+ c {ms}")
    print p.parse("convert(1 {ms}, s)+ c {s}")
    print p.parse("a == b")
