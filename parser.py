#!/usr/bin/env python

import ply.lex as lex
import ply.yacc as yacc

class Parser:
    def __init__(self, **kw):
        self.debug = kw.get('debug', 0)
        self.lexer = lex.lex(module=self, debug=self.debug)
        self.parser = yacc.yacc(module=self,
                                debug=self.debug,
                                write_tables=0,
        )
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
        #"rename" : "RENAME",
        #"to" : "TO",
        #"from" : "FROM",
    }

    tokens = (
        "LEQ",
        "GEQ",
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

    literals = "(){}<>!;?:+-/*^=."

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
        ("left", '?', ':', 'TERNARY'),
        ("left", "OR"),
        ("left", "AND"),
        ("nonassoc", "<", ">", "BOOLEQ", "LEQ", "GEQ"),
        ('left', '+', '-'),
        ('left', '*', '/'),
        ('nonassoc', 'ANNOTATE_UNIT'),
        ('right', 'NOT', 'UMINUS'),
        ('left', '^'),
       )

    def p_topLevelStatementsOpt(self, p):
        '''topLevelStatementsOpt : topLevelStatement topLevelStatementsOpt'''
        pass

    def p_topLevelStatementsOpt_term(self, p):
        '''topLevelStatementsOpt : empty'''
        pass
    def p_topLevelStatementOpt_scope(self, p):
        '''topLevelStatementsOpt : '{' topLevelStatementsOpt '}' '''
        pass

    def p_topLevelStatement(self, p):
        '''topLevelStatement : sharedStatement
                             | subSystemDefinition
        '''
        pass

    def p_sharedStatement_unit(self, p):
        '''sharedStatement : SHARED STABLE varWithUnit ';'
                                 | SHARED EPHEMERAL varWithUnit ';'
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

    def p_subSystemStatement_scope(self, p):
        '''subSystemStatement : scopeDefinition'''
        pass

    def p_scopeDefinition_scope(self, p):
        '''scopeDefinition : '{' subSystemStatementsOpt '}' '''
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
        '''providesStatement : PROVIDES ACCUM varMaybeUnit ';'
                             | PROVIDES DIFFVAR varMaybeUnit ';'
                             | PROVIDES PARAM varMaybeUnit ';'
                             | PROVIDES STABLE varMaybeUnit ';'
                             | PROVIDES EPHEMERAL varMaybeUnit ';'
                             | PROVIDES varMaybeUnit ';'
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
        '''varDef : var TIMESEQ realExprMaybeUnit ';'
                  | var DIVIDEEQ realExprMaybeUnit ';'
                  | var EXPONEQ realExprMaybeUnit ';'
        '''
        pass

    def p_assignDef(self, p):
        '''assignDef : varMaybeUnit '=' realExprMaybeUnit ';' '''
        pass
    def p_accumDef(self, p):
        '''accumDef : var PLUSEQ realExprMaybeUnit ';'
                    | var MINUSEQ realExprMaybeUnit ';'
        '''
        pass

    def p_subSystemStatement_diffDef(self, p):
        '''subSystemStatement : diffDef'''
        pass

    def p_diffDef(self, p):
        '''diffDef : var '.' INIT '=' realExprMaybeUnit ';'
                   | var '.' DIFF '=' realExprMaybeUnit ';'
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
        '''ifClause : IF realExprUnclearUnit '{' subSystemStatementsOpt '}'
                    | IF realExprWithoutUnit '{' subSystemStatementsOpt '}'
        '''
        pass
    def p_elseifClause(self, p):
        '''elseifClause : ELSEIF realExprUnclearUnit '{' subSystemStatementsOpt '}'
                        | ELSEIF realExprWithoutUnit '{' subSystemStatementsOpt '}'
        '''
        pass
    def p_elseClause(self, p):
        '''elseClause : ELSE '{' subSystemStatementsOpt '}' '''
        pass

    def p_unitExpr_literal(self, p):
        '''unitExpr : NAME'''
        pass
    def p_unitExpr_1(self, p):
        '''unitExpr : ONE'''
        pass
    def p_unitExpr_paren(self, p):
        '''unitExpr : '(' unitExpr ')' '''
        pass

    def p_unitExpr_op(self, p):
        '''unitExpr : unitExpr '*' unitExpr
                    | unitExpr '/' unitExpr
        '''
        pass

    def p_unitExpr_expon(self, p):
        '''unitExpr : unitExpr '^' numberLiteral'''
        pass

    def p_var(self, p):
        '''var : NAME'''
        pass

    def p_varMaybeUnit_noUnit(self, p):
        '''varMaybeUnit : var'''
        pass
    def p_varMaybeUnit_withUnit(self, p):
        '''varMaybeUnit : varWithUnit'''
        pass
    def p_varWithUnit(self, p):
        '''varWithUnit : var '{' unitExpr '}' '''
        pass

    ############################################
    def p_realExprUnclearUnit_base(self, p):
        '''realExprUnclearUnit : var'''
        pass
    def p_realExprUnclearUnit_binop(self, p):
        '''realExprUnclearUnit : realExprUnclearUnit '+' realExprUnclearUnit
                               | realExprUnclearUnit '-' realExprUnclearUnit
                               | realExprUnclearUnit '*' realExprUnclearUnit
                               | realExprUnclearUnit '/' realExprUnclearUnit
        '''
        pass
    def p_realExprUnclearUnit_pow(self, p):
        '''realExprUnclearUnit : POW '(' realExprUnclearUnit ',' realExprMaybeUnit ')'
                               | realExprUnclearUnit '^' realExprUnclearUnit
                               | realExprUnclearUnit '^' realExprWithoutUnit
                               | realExprUnclearUnit '^' realExprWithUnit
        '''
        pass
    def p_realExpr_paren(self, p):
        '''realExprUnclearUnit : '(' realExprUnclearUnit ')' '''
        pass
    def p_realExprUnclearUnit_unaryMinus(self, p):
        '''realExprUnclearUnit : '-' realExprUnclearUnit %prec UMINUS'''
        pass
    def p_realExprUnclearUnit_ternary(self, p):
        '''realExprUnclearUnit : realExprUnclearUnit '?' realExprUnclearUnit ':' realExprUnclearUnit %prec TERNARY
                               | realExprWithoutUnit '?' realExprUnclearUnit ':' realExprUnclearUnit %prec TERNARY
        '''
        pass
    def p_realExprUnclearUnit_func(self, p):
        '''realExprUnclearUnit : NAME '(' funcArgListOpt ')' '''
        pass
    def p_funcArgListOpt_zero(self, p):
        '''funcArgListOpt : empty'''
        pass
    def p_funcArgListOpt_nonzero(self, p):
        '''funcArgListOpt : funcArgList'''
        pass
    def p_funcArgList_term(self, p):
        '''funcArgList : realExprMaybeUnit'''
        pass
    def p_funcArgList_shift(self, p):
        '''funcArgList : realExprMaybeUnit ',' funcArgList'''
        pass



    ####################################
    def p_numberLiteralPlusOpt(self, p):
        '''numberLiteralPlusOpt : numberLiteral
                                | '+' numberLiteral
        '''
        pass
    def p_numberLiteral(self, p):
        '''numberLiteral : NUMBER
                         | ONE
        '''
        pass
    def p_realExprWithoutUnit_base(self, p):
        '''realExprWithoutUnit : numberLiteralPlusOpt'''
        pass
    def p_realExprUnclearUnit_boolLiteral(self, p):
        '''realExprWithoutUnit : boolLiteral
        '''
        pass
    def p_boolLiteral(self, p):
        '''boolLiteral : TRUE
                       | FALSE
        '''
        pass
    def p_realExprUnclearUnit_boolop(self, p):
        '''realExprWithoutUnit : realExprUnclearUnit boolOp realExprUnclearUnit %prec BOOLEQ
                               | realExprWithoutUnit boolOp realExprWithoutUnit %prec BOOLEQ
                               | realExprWithoutUnit boolOp realExprUnclearUnit %prec BOOLEQ
                               | realExprUnclearUnit boolOp realExprWithoutUnit %prec BOOLEQ
                               | realExprWithUnit boolOp realExprWithUnit       %prec BOOLEQ
                               | realExprWithUnit boolOp realExprUnclearUnit    %prec BOOLEQ
                               | realExprUnclearUnit boolOp realExprWithUnit    %prec BOOLEQ
        '''
        pass
    def p_boolOp(self, p):
        '''boolOp : BOOLEQ
                  | LEQ
                  | GEQ
                  | '<'
                  | '>'
        '''
        pass
    def p_realExprWithoutUnit_and(self, p):
        '''realExprWithoutUnit : realExprWithoutUnit AND realExprWithoutUnit
                               | realExprWithoutUnit AND realExprUnclearUnit
                               | realExprUnclearUnit AND realExprWithoutUnit
        '''
        pass
    def p_realExprWithoutUnit_or(self, p):
        '''realExprWithoutUnit : realExprWithoutUnit OR realExprWithoutUnit
                               | realExprWithoutUnit OR realExprUnclearUnit
                               | realExprUnclearUnit OR realExprWithoutUnit
        '''
        pass
    def p_realExprWithoutUnit_not(self, p):
        '''realExprWithoutUnit : NOT realExprWithoutUnit
                               | NOT realExprUnclearUnit
        '''
        pass
    #def p_realExprWithoutUnit_resolve(self, p):
    #    '''boolExpr : realExprWithoutUnit
    #                | realExprUnclearUnit
    #    '''

    def p_realExprWithoutUnit_binop(self, p):
        '''realExprWithoutUnit : realExprUnclearUnit '+' realExprWithoutUnit
                               | realExprUnclearUnit '-' realExprWithoutUnit
                               | realExprUnclearUnit '*' realExprWithoutUnit
                               | realExprUnclearUnit '/' realExprWithoutUnit
                               | realExprWithoutUnit '+' realExprUnclearUnit
                               | realExprWithoutUnit '-' realExprUnclearUnit
                               | realExprWithoutUnit '*' realExprUnclearUnit
                               | realExprWithoutUnit '/' realExprUnclearUnit
                               | realExprWithoutUnit '+' realExprWithoutUnit
                               | realExprWithoutUnit '-' realExprWithoutUnit
                               | realExprWithoutUnit '*' realExprWithoutUnit
                               | realExprWithoutUnit '/' realExprWithoutUnit
        '''
        pass
    def p_realExprWithoutUnit_ternary(self, p):
        '''realExprWithoutUnit : realExprUnclearUnit '?' realExprWithoutUnit ':' realExprWithoutUnit %prec TERNARY
                               | realExprUnclearUnit '?' realExprWithoutUnit ':' realExprUnclearUnit %prec TERNARY
                               | realExprUnclearUnit '?' realExprUnclearUnit ':' realExprWithoutUnit %prec TERNARY
                               | realExprWithoutUnit '?' realExprWithoutUnit ':' realExprWithoutUnit %prec TERNARY
                               | realExprWithoutUnit '?' realExprWithoutUnit ':' realExprUnclearUnit %prec TERNARY
                               | realExprWithoutUnit '?' realExprUnclearUnit ':' realExprWithoutUnit %prec TERNARY
        '''
        pass
    def p_realExprWithoutUnit_pow(self, p):
        '''realExprWithoutUnit : POW '(' realExprWithoutUnit ',' realExprMaybeUnit ')'
                               | realExprWithoutUnit '^' realExprUnclearUnit
                               | realExprWithoutUnit '^' realExprWithoutUnit
                               | realExprWithoutUnit '^' realExprWithUnit
        '''
        pass
    def p_realExprWithout_paren(self, p):
        '''realExprWithoutUnit : '(' realExprWithoutUnit ')' '''
        pass
    def p_realExprWithoutUnit_unaryMinus(self, p):
        '''realExprWithoutUnit : '-' realExprWithoutUnit %prec UMINUS'''
        pass



    ####################
    def p_realExprWithUnit_FromWithout(self, p):
        '''realExprWithUnit : realExprWithoutUnit '{' unitExpr '}' %prec ANNOTATE_UNIT'''
        pass
    def p_realExprWithUnit_FromUnclear(self, p):
        '''realExprWithUnit : realExprUnclearUnit '{' unitExpr '}' %prec ANNOTATE_UNIT'''
        pass
    def p_realExprWithUnit_binop(self, p):
        '''realExprWithUnit : realExprUnclearUnit '+' realExprWithUnit
                            | realExprUnclearUnit '-' realExprWithUnit
                            | realExprUnclearUnit '*' realExprWithUnit
                            | realExprUnclearUnit '/' realExprWithUnit
                            | realExprWithUnit    '+' realExprUnclearUnit
                            | realExprWithUnit    '-' realExprUnclearUnit
                            | realExprWithUnit    '*' realExprUnclearUnit
                            | realExprWithUnit    '/' realExprUnclearUnit
                            | realExprWithUnit    '+' realExprWithUnit
                            | realExprWithUnit    '-' realExprWithUnit
                            | realExprWithUnit    '*' realExprWithUnit
                            | realExprWithUnit    '/' realExprWithUnit
        '''
        pass
    def p_realExprWithUnit_ternary(self, p):
        '''realExprWithUnit : realExprUnclearUnit '?' realExprWithUnit    ':' realExprWithUnit %prec TERNARY
                            | realExprUnclearUnit '?' realExprWithUnit    ':' realExprUnclearUnit %prec TERNARY
                            | realExprUnclearUnit '?' realExprUnclearUnit ':' realExprWithUnit %prec TERNARY
                            | realExprWithoutUnit '?' realExprWithUnit    ':' realExprWithUnit %prec TERNARY
                            | realExprWithoutUnit '?' realExprWithUnit    ':' realExprUnclearUnit %prec TERNARY
                            | realExprWithoutUnit '?' realExprUnclearUnit ':' realExprWithUnit %prec TERNARY
        '''
        pass
    def p_realExprWithUnit_pow(self, p):
        '''realExprWithUnit : POW '(' realExprWithUnit ',' realExprMaybeUnit ')'
                               | realExprWithUnit '^' realExprUnclearUnit
                               | realExprWithUnit '^' realExprWithoutUnit
                               | realExprWithUnit '^' realExprWithUnit
        '''
        pass
    def p_realExprWith_paren(self, p):
        '''realExprWithUnit : '(' realExprWithUnit ')' '''
        pass
    def p_realExprWithUnit_unaryMinus(self, p):
        '''realExprWithUnit : '-' realExprWithUnit %prec UMINUS'''
        pass

    ######################

    def p_realExprMaybeUnit_Unclear(self, p):
        '''realExprMaybeUnit : realExprUnclearUnit '''#%prec UNITRESOLVE'''
        pass
    def p_realExprMaybeUnit_Without(self, p):
        '''realExprMaybeUnit : realExprWithoutUnit '''#%prec UNITRESOLVE'''
        pass
    def p_realExprMaybeUnit_With(self, p):
        '''realExprMaybeUnit : realExprWithUnit '''#%prec UNITRESOLVE'''
        pass

    def p_empty(self, p):
        '''empty :'''
        pass



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

    
