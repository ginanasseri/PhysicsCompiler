#include "parser.hpp"
#include "lexer.hpp"

#define EOF_TOKEN 0

void printToken(std::string token_type); // prints token type, used for debugging

// constructor 
Parser::Parser() : subject(""), graph(""), equation_id(0), data_values(""), token(0) {} 


/* ======================== public section ========================== 
 * ==================================================================
 * Function: parse()
 * ------------------------------------------------------------------
 * Parameters: FILE* - the physics source code file being parsed.
 *  
 * Purpose:  Serves as the entry and exit point of the parser, sets the 
 *           input file for the lexer, starts the parsing process, and 
 *           handles the control flow of the process. If an error occurs,
 *           or a syntax error in the source code is discovered at any 
 *           point, the error is propagated to the parse() function and 
 *           returend to the caller. 
 *
 *           The parsing routine for each keyword is called in order and,
 *           upon successful return, the keyword's value will be saved 
 *           Parser's associated data member, and the current token for
 *           the lexer will be the first token in the next line of the 
 *           source code file. 
 *
 * Returns:  Each parse routine and method, unless otherwise stated, returns
 *           an integer value:
 *
 *           0 if the parse was successful 
 *           1 if a syntax error was discovered or a parser error occurred 
 *
 */
int Parser::parse(FILE* input_file) {

    yyin = input_file; // set input file for lexer

    // SUBJECT, GRAPH
    if (parseImports() != 0) {  
        return 1;
    }
    if (getNextToken(NEWLINE, "expected new line") != 0) {
        return 1;
    }
    // EQUATION
    if (parseEquationID() != 0) {
        return 1;
    }
    if (getNextToken(NEWLINE, "expected new line") != 0) {
        return 1;
    }
    // DATA  
    if (parseData() != 0) {
        return 1;
    }
    // EOF_TOKEN (ignore newlines)
    while (token != EOF_TOKEN) {
        if (token == NEWLINE) {
            token = yylex();
            continue;
        }
        printErrorMessage("Additional text found after data section.");
        return 1;
    } 
    return 0;
}

/* ------------------------------------------------------------
 *           Return Data Members Functions
 * ------------------------------------------------------------*/
std::string Parser::getSubject() const {
    return subject;
}

std::string Parser::getGraphValue() const {
    return graph;
}

int Parser::getEquationID() const {
    return equation_id;
}

std::string Parser::getData() const {
    return data_values;
}

// ==================================================================
/* ==================== private seciton =============================
 * ==================================================================
 *  - Settting data member values 
 */
void Parser::setSubject(const std::string& subject_type) {
    subject = subject_type;
}

void Parser::setGraphValue(const std::string& graph_value) {
    graph = graph_value;
}

void Parser::setEquationID(int i) {
    equation_id = i;
}

/* --------------------------------------
 * update mechanics data: addVariable()
 * --------------------------------------
 * Builds the data_values data member (std::string) to match Python
 * sytax: each variable assignment statement is added to the 
 * data_members string separated by a new line and indented once 
 * to match the indentation of the main function. 
 */
void Parser::addVariable(const std::string& value) {
    if (!data_values.empty()) {
        data_values += "\n";
    }
    data_values += "    ";
    data_values += value;
}                          

/* ==================================================================
 * Token Advancing Method: getNextToken()
 * ------------------------------------------------------------------
 * Parameters: int - expected type of current token 
 *             const char* - error message to print if the current token
 *             is not the expected type.
 * 
 *  Purpose: Calls the lexer to get the next token in the file and checks
 *           if the token is the correct type according to the syntax 
 *           of the line being parsed. If the current token type and the
 *           expected token type are not the same, then an error message 
 *           indicating the syntax error type along with the line number 
 *           where it occurred is printed, and the return value of 1
 *           is propagated back to the parse() function which stops the
 *           parsing process. 
 *         
 *  Example: The correct sequence of tokens for the first line is 
 *                   SUBJECT COLON SUBJECT_TYPE
 *         
 *           The syntax errors below will result in the following 
 *           messages (to the right of the arrows) being printed 
 *           and 1 being returned. 
 *           
 *           "subject mechanics"  --> Error: expected ':' (line 1)
 *           "subject; mechanics" --> Error: unexpected ';' (line 1)
 *
 *  Returns:  0 if the current token type matches the expected type
 *            1 if the current token type does not match the expected type
 *  
 */
int Parser::getNextToken(int expected_type, const std::string& type_error) {
    token = yylex();
    if (token == UNREC) {
        std::cerr << "Error: unexpected '" << yylval.strval << "'"
                  << std::endl;
        return 1;
    }
    if (token != expected_type) {
        std::cerr << "Error: " << type_error << " (line "
                  << lineNumber() << ")" << std::endl;
        return 1;
    }
    return 0;
}

/* ------------------------------------------------------------------- 
 * Print Line Number: lineNumber()
 * -------------------------------------------------------------------
 * Purpose: Called in the case where a syntax error was discovered so that
 *          the correct line number of the error in the source code is 
 *          printed with the error message. 
 *          
 *          If the error was at the end of a line, then the current token 
 *          is a new line character, and the lexer's line number count,
 *          yylineno, will be ahead by one line.
 *
 * Returns: yylineno - 1 if the current token is a newline character 
 *          yylineno otherwise
 */
int Parser::lineNumber() {
    if (token == NEWLINE) {
        return yylineno - 1;
    }
    return yylineno;
}

/* ------------------------------------------------------------------- 
 * Print Error Message: printErrorMessage() 
 * -------------------------------------------------------------------
 * Parameters: error_message - details of the type of syntax error which 
 *             occurred: missing token or unrecognized token
 *
 * Purpose: Prints the details of the error along with the line number.
 */
void Parser::printErrorMessage(const std::string& error_message) {
    std::cerr << "Error: " << error_message << " (line " 
              << lineNumber() << ")" << std::endl;
}

/* ------------------------------------------------------------------- 
 * Check Curl Balance: checkCurlBalance()
 * -------------------------------------------------------------------
 * Parameters: curl_balance - set to 1 indicating we have read 1 LCURL
 *
 * Purpose: Checks that we have exactly one RCURL to close the LCURL
 * 
 * Returns: 0 if the brackets are balanced. 1 otherwise. 
 */
int Parser::checkCurlBalance(int curl_balance) {
    // check for closing RCURL
    if (token != RCURL) {
        printErrorMessage("Missing '}'");
        return 1;
    }
    // check for trailing RCURL(s)
    while (token == RCURL) {
        curl_balance--;
        token = yylex();
    }
    if (curl_balance < 0) {
        printErrorMessage("Unexpected '}'");
        return 1;
    }
    return 0;
}

/* ===================================================================
 *           * * *  PARSING ROUTINES START HERE * * * 
 * ===================================================================
 * Note: For all parse routines return 1 if a parse error occurred and
 *       0 if the parse routine was successful.
 * -------------------------------------------------------------------
 * Parse Imports: parseImports() 
 *--------------------------------------------------------------------
 * Purpose: parses the first two lines in the source code checking for 
 *          the syntax:
 *
 *               1. SUBJECT COLON SUBJECT_TYPE NEWLINE 
 *               2. GRAPH COLON GRAPH VALUE NEWLINE 
 */
int Parser::parseImports() {
    // parse first line of source code: the subject type
    if (parseSubject() != 0) {
        return 1;
    }
    // make sure new line separates import statements
    if (getNextToken(NEWLINE, "new line") != 0) {
        return 1;
    }
    // parse next line of source code: the graph value 
    if (parseGraphValue() != 0) {
        return 1;
    }
    // if we're here, parse imports successful 
    return 0;
}

/* -------------------------------------------------------------------
 * Parse 'subject' line: parseSubject() 
 * -------------------------------------------------------------------
 * Purpose: 1. Parse SUBJECT COLON SUBJECT_TYPE 
 *          2. Save SUBJECT_TYPE value to subject data member
 */
int Parser::parseSubject() {
    if (getNextToken(SUBJECT, "Missing keyword: 'subject'") != 0) {
        return 1;
    }
    if (getNextToken(COLON, "Expected ':'") != 0) {
        return 1;
    }
    if (getNextToken(SUBJECT_TYPE, "Invalid subject type. Use: 'mechanics' or 'stats'") != 0) {
        return 1;
    }
    // if we're here: current token is SUBJECT_TYPE 
    setSubject(yylval.strval);
    return 0; 
}
// ===================================================================
/* -------------------------------------------------------------------
 * Parse 'graph' line: parseGraphImport()
 * -------------------------------------------------------------------
 * Purpose : 1. Parse GRAPH COLON GRAPH_VALUE
 *           2. Save GRAPH_VALUE value to graph data member
 */
int Parser::parseGraphValue() {
    if (getNextToken(GRAPH, "Missing keyword: 'graph'") != 0) {
        return 1;
    }
    if (getNextToken(COLON, "Expected ':'") != 0) {
        return 1;
    }
    if (getNextToken(GRAPH_VALUE, "Invalid graph value. Use: 'true' or 'false'.") != 0) {
        return 1;
    }
    // if we're here, current token is GRAPH_VALUE
    setGraphValue(yylval.strval);
    return 0;
}


/* -------------------------------------------------------------------
 * Parse 'equation' line: parseEquation() 
 * -------------------------------------------------------------------
 * Purpose:  1. Parse EQUATION COLON INT
 *           2. Save INT value to equation data member 
 */
int Parser::parseEquationID() {
    if (getNextToken(EQUATION, "Missing keyword: 'equation'") != 0) {
        return 1;
    }
    if (getNextToken(COLON, "Expected ':'") != 0) {
        return 1;
    }
    if (getNextToken(INT, "Expected equation number (integer value)") != 0) {
        return 1;
    }   
    // if we're here, current token is INT
    setEquationID(yylval.intval);
    return 0;
}

/* -------------------------------------------------------------------
 * Parse 'data' line: parseData() 
 * -------------------------------------------------------------------
 * Purpose: 1. Parse DATA COLON LCURL <data values> RCURL
 *                  
 *          2. The correct syntax for the data withing LCURL RCURL 
 *             depends on the subject, so check the subject type and 
 *             call the appropriate parse routine to get the data. 
 * 
 *          3. Balance check LCURL and RCURL
 * 
 * Returns: The data parsing routine will return 1 if a syntax error was 
 *          discovered, or an unreachable code section was entered due to
 *          an invalid value returned from a parser's data member. (This
 *          is not possible unless something went seriously wrong.)
 * 
 *          Otherwise, all data was successfully written into the data_values 
 *          data memeber, and 0 is returned.  
 */
int Parser::parseData() {
    if (getNextToken(DATA, "data") != 0) {
        return 1;
    }
    if (getNextToken(COLON, "Expected ':'") != 0) {
        return 1;
    }
    if (getNextToken(LCURL, "Expected '{'") != 0) {
        return 1;
    }
    int curl_balance = 1; // set curl balance flag

    // Get next token and call the parse routine for the current subject
    token = yylex();
    
    if (getSubject() == "mechanics") {
        if (parseMechanicsData() != 0) {
            return 1; 
        }
    }
    else if (getSubject() == "stats") {
        if (parseDataFile() != 0) {
            return 1;
        }  
    }
    else {
        // parseSubject() only accepts the values mechanics and stats, so this section should never be reached.
        printErrorMessage("parseData() in unreachable code section: getSubject() returned invalid subject type.");
        return 1;
    }
    // check that there is exactly 1 closing bracket 
    if (checkCurlBalance(curl_balance) != 0) {
        return 1;
    }
    return 0;
}

// -------------------------------------------------------------------
/* ===================================================================
 *                   Mechanics Domain
 * ===================================================================
 * -------------------------------------------------------------------
 * Parse int or float: parseIntOrFloat()                     INT|FLOAT
 * -------------------------------------------------------------------
 * Purpose:  1. Parse INT|FLOAT
 *           2. Add it to the end of the VARIABLE_NAME ASSIGN statement
 * 
 * Returns 1 if the assignment value is not an int or a float. 
 *   
 */
int Parser::parseIntOrFloat(std::string& assign_str) {
    token = yylex();

    // INT
    if (token == INT) {
        assign_str += std::to_string(yylval.intval);
        return 0;
    }
    // FLOAT
    else if (token == FLOAT) {
        assign_str += std::to_string(yylval.floatval);
        return 0;
    }
    // Invalid assignment value
    else {
        return 1;
    }
    printToken("I shouldn't be here (parseIntOrFloat() in unreachable section).");
    return 1;    
}

// -------------------------------------------------------------------
/* -------------------------------------------------------------------
 * Function: parseVariable()                    VARIABLE ASSIGNMENT
 * -------------------------------------------------------------------
 * Purpose: Parse variable assignments within the data section in the 
 *          case that the subject is mechanics and adds each assignment   
 *          to the data_values data member. Uses recursive calls to 
 *          parse each assignment. 
 *  
 *   Grammar: 
 *
 *  VARIABLE_NAME ASSIGN INT|FLOAT (COMMA VARIABLE_NAME ASSIGN INT|FLOAT)*
 *
 */
int Parser::parseVariable() {

    while (token != RCURL && token != EOF_TOKEN && token != NEWLINE) {

        // VARIABLE_NAME
        if (token != VARIABLE_NAME) {
            printErrorMessage("Expected variable name.");
            return 1;
        }
        std::string assignment_str = yylval.strval;  // assignment statement string
        std::string variable = assignment_str + "'"; // used for printing error 

        // ASSIGN
        token = yylex();
        if (token != ASSIGN) {
            printErrorMessage("Expected '='");
            return 1;
        }
        assignment_str += " = ";
   
        // INT|FLOAT
        if (parseIntOrFloat(assignment_str) != 0) {
            printErrorMessage("Missing/invalid assignment value for '" + variable);
            return 1;
        }
        addVariable(assignment_str); // add variable assignment to data_values
        token = yylex();             // get next token

        // (COMMA VARIABLE_NAME ASSIGN INT|FLOAT)* 
        if(token == COMMA) {
            token = yylex();
            if (token != VARIABLE_NAME) {
                printErrorMessage("Misplaced ',' or missing variable name in data.");
                return 1;
            }
            // get next variable with recursive call 
            if (parseVariable() != 0) {
                return 1;
            }
        }
    }
    return 0;
}

/* -----------------------------------------------------------------
 * Function: parseMechanicsData()                   MECHANICS DATA
 * -----------------------------------------------------------------
 * Purpose: Parse the variable assignment statements in data section.
 */
int Parser::parseMechanicsData() {

    switch(token) {

        case(RCURL): 
            return 0;

        case(VARIABLE_NAME):
            if (parseVariable() != 0) {
                return 1;
            }
            return 0;
    
        default: 
            printErrorMessage("Invalid data input for subject: mechanics. Expected variable assignment.");
            return 1;
    }
    printErrorMessage("I shouldn't be here (parseMechanicsData() in unreachable code section).");
    return 1;
}

// -------------------------------------------------------------------
/* ===================================================================
 *                   Stats Domain
 * ==================================================================*/
/* -------------------------------------------------------------------
 * Function: parseDataFile()                               DATA FILE 
 * ------------------------------------------------------------------
 * Purpose: Parse data file for stats equation. If we're here, current token
 *          should be DATA_FILE. If it is not, return an error. Otherwise, 
 *          save the filename to the data_values data member. 
 * 
 **/
int Parser::parseDataFile() {

    // stats MUST include data within the {}'s
    if (token == RCURL) {
        printErrorMessage("No file provided. Use: include '<filepath>.csv' in data section.");
        return 1;
    }
    if (token != DATA_FILE) {
        printErrorMessage("Invalid file type. Use: file with extension '.csv'.");
        return 1; 
    }
    data_values += yylval.strval; // add file name to data
    token = yylex();
    return 0;
}

