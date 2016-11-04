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
    "avail" : "AVAIL",
    "provides" : "PROVIDES",
    "subsystem" : "SUBSYSTEM",
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
    "NAME"
) + tuple(reserved.values())

literals = "(){}<>!;:+-/*^=."

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



t_LEQ = r'<='
t_GEQ = r'>='
t_BOOLEQ = r'=='

t_PLUSEQ = r'\+='
t_MINUSEQ = r'-='
t_TIMESEQ = r'\*='
t_DIVIDEEQ = r'/='
t_EXPONEQ = r'^='

t_NUMBER = r'(?:(?:[0-9]+(?:\.[0-9]*)?)|(?:\.[0-9]+))(?:[eE][-\+]?[0-9]+)?'

t_ignore_CPP_COMMENT = r'//.*?\n'
t_ignore_C_COMMENT = r'/\*.*?\*/'
t_ignore_PYTHON_COMMENT = r'\#.*?\n'
t_ignore_MATLAB_COMMENT = r'%.*?\n'

precedence = (
    ("nonassoc", "AND", "OR"),
    ("nonassoc", "<", ">", "BOOLEQ", "LEQ", "GEQ"),
    ('left', '+', '-'),
    ('left', '*', '/'),
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
    '''topLevelStatement : availabilityStatement
                         | subSystemDefinition
    '''
    pass

def p_availabilityStatement_unit(t):
    '''availabilityStatement : AVAIL STABLE varWithUnit ';'
                             | AVAIL EPHEMERAL varWithUnit ';'
    '''
    pass
def p_availabilityStatement_flag(t):
    '''availabilityStatement : AVAIL flagDeclBool ';'
                             | AVAIL flagDeclEnum ';'
    '''
    pass
def p_flagDeclBool(t):
    '''flagDeclBool : FLAG var'''
    pass
def p_flagDeclEnum(t):
    '''flagDeclEnum : FLAG var ':' '{' nameList '}' '''
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

def p_subSystemStatement_avail(t):
    '''subSystemStatement : availabilityStatement'''
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
    '''varDef : varMaybeUnit TIMESEQ realExpr ';'
              | varMaybeUnit DIVIDEEQ realExpr ';'
              | varMaybeUnit EXPONEQ realExpr ';'
    '''
    pass

def p_assignDef(t):
    '''assignDef : varMaybeUnit '=' realExpr ';' '''
    pass

def p_accumDef(t):
    '''accumDef : varMaybeUnit PLUSEQ realExpr ';'
                | varMaybeUnit MINUSEQ realExpr ';'
    '''
    pass

def p_subSystemStatement_diffDef(t):
    '''subSystemStatement : diffDef'''
    pass

def p_diffDef(t):
    '''diffDef : var '.' INIT '=' realExpr ';'
               | var '.' DIFF '=' realExpr ';'
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
    '''ifClause : IF boolExpr '{' subSystemStatementsOpt '}' '''
    pass
def p_elseifClause(t):
    '''elseifClause : ELSEIF boolExpr '{' subSystemStatementsOpt '}' '''
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
def p_boolExpr_varLiteral(t):
    '''boolExpr : var'''
    pass
def p_boolExpr_neg(t):
    '''boolExpr : NOT boolExpr '''
    pass
def p_boolExpr_compare(t):
    '''boolExpr : boolComparison'''
    pass
def p_boolExpr_binary(t):
    '''boolExpr : boolExpr AND boolExpr
                | boolExpr OR boolExpr
    '''
    pass
def p_boolCompare(t):
    '''boolComparison : realExpr BOOLEQ realExpr
                      | realExpr LEQ realExpr
                      | realExpr GEQ realExpr
                      | realExpr '<' realExpr
                      | realExpr '>' realExpr
    '''
    pass
def p_boolExpr_paren(t):
    '''boolExpr : '(' boolExpr ')' '''
    pass

def p_unitExpr_literal(t):
    '''unitExpr : NAME'''
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
    '''unitExpr : unitExpr '^' NUMBER'''
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
    '''varWithUnit : var ':' unitExpr'''
    pass
def p_realExpr_literal(t):
    '''realExpr : numberLiteral'''
    pass
def p_realExpr_numberLiteral(t):
    '''numberLiteral : NUMBER
                     | '+' NUMBER
    '''
    pass
def p_realExpr_unaryMinus(t):
    '''realExpr : '-' realExpr %prec UMINUS'''
    pass
def p_realExpr_binaryOp(t):
    '''realExpr : realExpr '+' realExpr
                | realExpr '-' realExpr
                | realExpr '*' realExpr
                | realExpr '/' realExpr
                | realExpr '^' realExpr
    '''
    pass
def p_realExpr_paren(t):
    '''realExpr : '(' realExpr ')' '''
    pass
def p_realExpr_func(t):
    '''realExpr : NAME '(' funcArgListOpt ')' '''
    pass
def p_funcArgListOpt_zero(t):
    '''funcArgListOpt : empty'''
    pass
def p_funcArgListOpt_nonzero(t):
    '''funcArgListOpt : funcArgList'''
    pass
def p_funcArgList_term(t):
    '''funcArgList : realExpr'''
    pass
def p_funcArgList_shift(t):
    '''funcArgList : realExpr ',' funcArgList'''
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
    
