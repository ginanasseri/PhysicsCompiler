%{
#include "tokens.hpp"
%}

%option yylineno

%%
subject         { return SUBJECT; }
graph           { return GRAPH; }
equation        { return EQUATION; }
data            { return DATA; }
mechanics       { yylval.strval = strdup(yytext); return SUBJECT_TYPE; }
stats           { yylval.strval = strdup(yytext); return SUBJECT_TYPE; }
true            { yylval.strval = strdup(yytext); return GRAPH_VALUE; }
false           { yylval.strval = strdup(yytext); return GRAPH_VALUE; }
:               { return COLON; }
\{              { return LCURL; }
\}              { return RCURL; }
\(              { return LPAR; }
\)              { return RPAR; }
,               { return COMMA; }
=               { return ASSIGN; }
[a-zA-Z0-9_-]+\.csv    { yylval.strval = strdup(yytext); return DATA_FILE; }
theta                  { yylval.strval = strdup(yytext); return VARIABLE_NAME; }
rho                    { yylval.strval = strdup(yytext); return VARIABLE_NAME; }
([A-Z]|[a-z])[a-z]?[0-9]*   { yylval.strval = strdup(yytext); return VARIABLE_NAME; }
-?[0-9]+               { yylval.intval = atoi(yytext); return INT; }
-?[0-9]+\.[0-9]*       { yylval.floatval = atof(yytext); return FLOAT; }
-?[0-9]+(\.[0-9]*)?([eE][-+]?[0-9]+)?  { yylval.floatval = atof(yytext); return FLOAT; }
[ \t\r]+               { /* Skip spaces, tabs and carriage returns */ }
\n                     { return NEWLINE; }
.                      { yylval.strval = strdup(yytext); return UNREC; }
%%

int yywrap(void) {
    return 1;
}
