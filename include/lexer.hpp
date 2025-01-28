#ifndef LEXER_HPP
#define LEXER_HPP

#include "tokens.hpp"

YYSTYPE yylval;

extern int yylex();
extern int yylineno;
extern FILE *yyin;

#endif // LEXER_HPP
