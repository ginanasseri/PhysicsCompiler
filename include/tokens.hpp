#ifndef TOKENS_H
#define TOKENS_H
#undef BEGIN

enum TokenType {
    SUBJECT,
    GRAPH,
    EQUATION,
    DATA,
    SUBJECT_TYPE,
    GRAPH_VALUE,
    COLON,
    LCURL,
    RCURL,
    LPAR,
    RPAR,
    COMMA,
    ASSIGN, 
    VARIABLE_NAME,
    INT,
    FLOAT,
    DATA_FILE,
    NEWLINE,
    UNREC, // invalid token 
};

union YYSTYPE {
    int intval;
    float floatval;
    char *strval;
};

extern YYSTYPE yylval;

#endif // TOKENS_H
