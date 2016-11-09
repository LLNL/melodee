#!/usr/bin/env python

t_ignore = " \t"

def t_newline(t):
    r'\n+'
    t.lexer.lineno += len(t.value)

def t_error(t):
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

literals = "(){}<>!;+-/*^=."

def t_AND(t):
    r'(?:and)|(?:&&)'
    return t
def t_OR(t):
    r'(?:or)|(?:\|\|)'
    return t
def t_NOT(t):
    r'!|(?:not)'
    return t

def t_NAME(t):
    r'[_a-zA-Z][_a-zA-Z0-9]*'
    t.type = reserved.get(t.value,"NAME")
    return t

def t_NUMBER(t):
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
    ("nonassoc", "RESOLVE"),
    ("left", "OR"),
    ("left", "AND"),
    ("nonassoc", "<", ">", "BOOLEQ", "LEQ", "GEQ"),
    ('left', '+', '-'),
    ('left', '*', '/'),
    ('nonassoc', 'ANNOTATE_UNIT'),
    ('right', 'NOT', 'UMINUS'),
    ('left', '^'),
)

def p_topLevelStatementsOpt(t):
    '''topLevelStatementsOpt : topLevelStatement topLevelStatementsOpt'''
    pass

def p_topLevelStatementsOpt_term(t):
    '''topLevelStatementsOpt : empty'''
    pass
def p_topLevelStatementOpt_scope(t):
    '''topLevelStatementsOpt : '{' topLevelStatementsOpt '}' '''
    pass

def p_topLevelStatement(t):
    '''topLevelStatement : sharedStatement
                         | subSystemDefinition
    '''
    pass

def p_sharedStatement_unit(t):
    '''sharedStatement : SHARED STABLE varWithUnit ';'
                             | SHARED EPHEMERAL varWithUnit ';'
    '''
    pass
def p_sharedStatement_flag(t):
    '''sharedStatement : SHARED flagDeclBool ';'
                             | SHARED flagDeclEnum ';'
    '''
    pass
def p_flagDeclBool(t):
    '''flagDeclBool : FLAG var'''
    pass
def p_flagDeclEnum(t):
    '''flagDeclEnum : FLAG var '{' nameList '}' '''
    pass

def p_nameList_term(t):
    '''nameList : NAME'''
    pass
def p_nameList_shift(t):
    '''nameList : NAME ',' nameList'''
    pass

def p_subSystemDefinition(t):
    '''subSystemDefinition : SUBSYSTEM '{' subSystemStatementsOpt '}' '''
    pass

def p_subSystemStatementsOpt(t):
    '''subSystemStatementsOpt : subSystemStatement subSystemStatementsOpt
                              | empty 
    '''
    pass

def p_subSystemStatement_scope(t):
    '''subSystemStatement : scopeDefinition'''
    pass

def p_scopeDefinition_scope(t):
    '''scopeDefinition : '{' subSystemStatementsOpt '}' '''
    pass

def p_subSystemStatement_shared(t):
    '''subSystemStatement : sharedStatement'''
    pass

def p_subSystemStatement_subSystem(t):
    '''subSystemStatement : subSystemDefinition'''
    pass

def p_subSystemStatement_provide(t):
    '''subSystemStatement : providesStatement'''
    pass

def p_providesStatement_flagDecl(t):
    '''providesStatement : PROVIDES flagDeclBool ';'
                         | PROVIDES flagDeclEnum ';'
    '''
    pass
def p_provdesStatement_flagBoolDefn(t):
    '''providesStatement : PROVIDES flagDeclBool '=' boolLiteral ';' '''
    
def p_providesStatement_flagEnumDefn(t):
    '''providesStatement : PROVIDES flagDeclEnum '=' NAME ';' '''
    pass

def p_providesStatement_DeclSubNoUnit(t):
    '''providesStatement : PROVIDES ACCUM varMaybeUnit ';'
                         | PROVIDES DIFFVAR varMaybeUnit ';'
                         | PROVIDES PARAM varMaybeUnit ';'
                         | PROVIDES STABLE varMaybeUnit ';'
                         | PROVIDES EPHEMERAL varMaybeUnit ';'
                         | PROVIDES varMaybeUnit ';'
    '''
    pass

def p_providesStatement_Defn(t):
    '''providesStatement : PROVIDES ACCUM accumDef
                         | PROVIDES PARAM assignDef
                         | PROVIDES assignDef
    '''
    pass

def p_subSystemStatement_definition(t):
    '''subSystemStatement : varDef ';' '''
    pass

def p_varDef_assign(t):
    '''varDef : assignDef'''
    pass

def p_varDef_accum(t):
    '''varDef : accumDef'''
    pass

def p_varDef_modAssign(t):
    '''varDef : var TIMESEQ realExprMaybeUnit ';'
              | var DIVIDEEQ realExprMaybeUnit ';'
              | var EXPONEQ realExprMaybeUnit ';'
    '''
    pass

def p_assignDef(t):
    '''assignDef : varMaybeUnit '=' realExprMaybeUnit ';' '''
    pass

def p_accumDef(t):
    '''accumDef : var PLUSEQ realExprMaybeUnit ';'
                | var MINUSEQ realExprMaybeUnit ';'
    '''
    pass

def p_subSystemStatement_diffDef(t):
    '''subSystemStatement : diffDef'''
    pass

def p_diffDef(t):
    '''diffDef : var '.' INIT '=' realExprMaybeUnit ';'
               | var '.' DIFF '=' realExprMaybeUnit ';'
    '''
    pass

def p_subSystemStatement_if(t):
    '''subSystemStatement : ifStatement'''
    pass

def p_ifStatement_noElse(t):
    '''ifStatement : ifClauses'''
    pass
def p_ifStatement_else(t):
    '''ifStatement : ifClauses elseClause'''
def p_ifClauses(t):
    '''ifClauses : ifClause elseifClausesOpt'''
    pass

def p_elseifClausesOpt(t):
    '''elseifClausesOpt : elseifClause elseifClausesOpt'''
    pass
def p_elseIfClausesOpt_term(t):
    '''elseifClausesOpt : empty'''
    pass

def p_ifClause(t):
    '''ifClause : IF '(' boolExpr ')' '{' subSystemStatementsOpt '}' '''
    pass
def p_elseifClause(t):
    '''elseifClause : ELSEIF '(' boolExpr ')' '{' subSystemStatementsOpt '}' '''
    pass
def p_elseClause(t):
    '''elseClause : ELSE '{' subSystemStatementsOpt '}' '''
    pass

def p_boolLiteral(t):
    '''boolLiteral : TRUE
                   | FALSE
    '''
    pass

def p_boolExpr_literal(t):
    '''boolExpr : boolLiteral'''
    pass
def p_boolExpr_binary(t):
    '''boolExpr : boolExpr AND boolExpr
                | boolExpr OR boolExpr'''
    pass
def p_boolExpr_net(t):
    '''boolExpr : NOT boolExpr'''
    pass
def p_boolExpr_paren(t):
    '''boolExpr : '(' boolExpr ')' '''
    pass
def p_boolExpr_boolCompare(t):
    '''boolExpr : boolCompare'''
    pass
def p_boolCompare_UnclearUnit(t):
    '''boolCompare : realExprUnclearUnit boolOp realExprUnclearUnit'''
    pass
def p_boolCompare_WithoutUnit(t):
    '''boolCompare : realExprWithoutUnit boolOp realExprWithoutUnit
                   | realExprWithoutUnit boolOp realExprUnclearUnit
                   | realExprUnclearUnit boolOp realExprWithoutUnit
    '''
    pass
def p_boolCompare_WithUnit(t):
    '''boolCompare : realExprWithUnit boolOp realExprWithUnit
                   | realExprWithUnit boolOp realExprUnclearUnit
                   | realExprUnclearUnit boolOp realExprWithUnit
    '''
    pass
def p_boolOp(t):
    '''boolOp : BOOLEQ
              | LEQ
              | GEQ
              | '<'
              | '>'
    '''
    pass
#def p_boolExpr_real(t):
#    '''boolExpr : realExpr'''
#    pass
#def p_realExpr_neg(t):
#    '''realExpr : NOT realExpr '''
#    pass
#def p_realExpr_compare(t):
#    '''realExpr : boolComparison'''
#    pass
#def p_realExpr_binary(t):
#    '''realExpr : realExpr AND realExpr
#                | realExpr OR realExpr
#    '''
#    pass
#def p_boolCompare(t):
#    '''boolComparison : realExpr BOOLEQ realExpr
#                      | realExpr LEQ realExpr
#                      | realExpr GEQ realExpr
#                      | realExpr '<' realExpr
#                      | realExpr '>' realExpr
#    '''
#    pass

def p_unitExpr_literal(t):
    '''unitExpr : NAME'''
    pass
def p_unitExpr_1(t):
    '''unitExpr : ONE'''
    pass
def p_unitExpr_paren(t):
    '''unitExpr : '(' unitExpr ')' '''
    pass

def p_unitExpr_op(t):
    '''unitExpr : unitExpr '*' unitExpr
                | unitExpr '/' unitExpr
    '''
    pass

def p_unitExpr_expon(t):
    '''unitExpr : unitExpr '^' numberLiteral'''
    pass

def p_var(t):
    '''var : NAME'''
    pass

def p_varMaybeUnit_noUnit(t):
    '''varMaybeUnit : var'''
    pass
def p_varMaybeUnit_withUnit(t):
    '''varMaybeUnit : varWithUnit'''
    pass
def p_varWithUnit(t):
    '''varWithUnit : var '{' unitExpr '}' '''
    pass

############################################
def p_realExprUnclearUnit_base(t):
    '''realExprUnclearUnit : var'''
    pass
def p_realExprUnclearUnit_binop(t):
    '''realExprUnclearUnit : realExprUnclearUnit '+' realExprUnclearUnit
                           | realExprUnclearUnit '-' realExprUnclearUnit
                           | realExprUnclearUnit '*' realExprUnclearUnit
                           | realExprUnclearUnit '/' realExprUnclearUnit
    '''
    pass
def p_realExprUnclearUnit_pow(t):
    '''realExprUnclearUnit : POW '(' realExprUnclearUnit ',' realExprMaybeUnit ')'
                           | realExprUnclearUnit '^' realExprMaybeUnit
    '''
    pass
def p_realExpr_paren(t):
    '''realExprUnclearUnit : '(' realExprUnclearUnit ')' '''
    pass
def p_realExprUnclearUnit_unaryMinus(t):
    '''realExprUnclearUnit : '-' realExprUnclearUnit %prec UMINUS'''
    pass
def p_realExprUnclearUnit_func(t):
    '''realExprUnclearUnit : NAME '(' funcArgListOpt ')' '''
    pass
def p_funcArgListOpt_zero(t):
    '''funcArgListOpt : empty'''
    pass
def p_funcArgListOpt_nonzero(t):
    '''funcArgListOpt : funcArgList'''
    pass
def p_funcArgList_term(t):
    '''funcArgList : realExprMaybeUnit'''
    pass
def p_funcArgList_shift(t):
    '''funcArgList : realExprMaybeUnit ',' funcArgList'''
    pass

####################################
def p_numberLiteralPlusOpt(t):
    '''numberLiteralPlusOpt : numberLiteral
                            | '+' numberLiteral
    '''
    pass
def p_numberLiteral(t):
    '''numberLiteral : NUMBER
                     | ONE
    '''
    pass

def p_realExprWithoutUnit_base(t):
    '''realExprWithoutUnit : numberLiteralPlusOpt'''
    pass
def p_realExprWithoutUnit_binop(t):
    '''realExprWithoutUnit : realExprUnclearUnit '+' realExprWithoutUnit
                           | realExprUnclearUnit '-' realExprWithoutUnit
                           | realExprUnclearUnit '*' realExprWithoutUnit
                           | realExprUnclearUnit '/' realExprWithoutUnit
                           | realExprWithoutUnit '+' realExprUnclearUnit
                           | realExprWithoutUnit '-' realExprUnclearUnit
                           | realExprWithoutUnit '*' realExprUnclearUnit
                           | realExprWithoutUnit '/' realExprUnclearUnit
    '''
    pass
def p_realExprWithoutUnit_pow(t):
    '''realExprWithoutUnit : POW '(' realExprWithoutUnit ',' realExprMaybeUnit ')'
                           | realExprWithoutUnit '^' realExprMaybeUnit
    '''
    pass
def p_realExprWithout_paren(t):
    '''realExprWithoutUnit : '(' realExprWithoutUnit ')' '''
    pass
def p_realExprWithoutUnit_unaryMinus(t):
    '''realExprWithoutUnit : '-' realExprWithoutUnit %prec UMINUS'''
    pass



####################
def p_realExprWithUnit_FromWithout(t):
    '''realExprWithUnit : realExprWithoutUnit '{' unitExpr '}' %prec ANNOTATE_UNIT'''
    pass
def p_realExprWithUnit_FromUnclear(t):
    '''realExprWithUnit : realExprUnclearUnit '{' unitExpr '}' %prec ANNOTATE_UNIT'''
    pass
def p_realExprWithUnit_binop(t):
    '''realExprWithUnit : realExprUnclearUnit '+' realExprWithUnit
                        | realExprUnclearUnit '-' realExprWithUnit
                        | realExprUnclearUnit '*' realExprWithUnit
                        | realExprUnclearUnit '/' realExprWithUnit
                        | realExprWithUnit    '+' realExprUnclearUnit
                        | realExprWithUnit    '-' realExprUnclearUnit
                        | realExprWithUnit    '*' realExprUnclearUnit
                        | realExprWithUnit    '/' realExprUnclearUnit
    '''
    pass
def p_realExprWithUnit_pow(t):
    '''realExprWithUnit : POW '(' realExprWithUnit ',' realExprMaybeUnit ')'
                           | realExprWithUnit '^' realExprMaybeUnit
    '''
    pass
def p_realExprWith_paren(t):
    '''realExprWithUnit : '(' realExprWithUnit ')' '''
    pass
def p_realExprWithUnit_unaryMinus(t):
    '''realExprWithUnit : '-' realExprWithUnit %prec UMINUS'''
    pass

######################

def p_realExprMaybeUnit_Unclear(t):
    '''realExprMaybeUnit : realExprUnclearUnit %prec RESOLVE'''
    pass
def p_realExprMaybeUnit_Without(t):
    '''realExprMaybeUnit : realExprWithoutUnit %prec RESOLVE'''
    pass
def p_realExprMaybeUnit_With(t):
    '''realExprMaybeUnit : realExprWithUnit %prec RESOLVE'''
    pass

def p_empty(t):
    '''empty :'''
    pass



#def p_subSystemStatement_renameStatement(t):
#    '''subSystemStatement : renameStatement'''
#    pass


if __name__=="__main__":
    import ply.lex as lex
    lexer = lex.lex()
    data = '''
and && or || not ! 0 2.0 .3 40. 5e+6 if myID */* bljsadfj */ */
'''
    lexer.input(data)

    while True:
        tok = lexer.token()
        if not tok:
            break
        print tok

    import ply.yacc as yacc
    parser = yacc.yacc()
    
