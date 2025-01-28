#ifndef PARSER_HPP
#define PARSER_HPP

#include <iostream>
#include <string>

class Parser {

public:
    Parser();
    int parse(FILE* input_file);
    
    std::string getSubject() const;
    std::string getGraphValue() const;
    std::string getData() const; 
    int getEquationID() const;
    

private:
    void setSubject(const std::string& subject);
    void setGraphValue(const std::string& graph);
    void setEquationID(int equation_id);

    int getNextToken(int expected_type, const std::string& type_error);
    int lineNumber();

    int parseImports();
    int parseSubject();     // subject
    int parseGraphValue();  // graph
    int parseEquationID();  // equation
    int parseData();        // data 

    // data methods
    int parseDataValues();
    int parseMechanicsData(); // mechanics parse methods --
    int parseVariable();
    int parseIntOrFloat(std::string& assign_str);
    void addVariable(const std::string& value); // --------
    int parseDataFile(); // parse stats data 

    int checkCurlBalance(int curl_balance);
    void printErrorMessage(const std::string& error_message);
    
    int token;
    int equation_id;
    std::string subject;
    std::string graph;
    std::string data_values;  
};


#endif // PARSER_HPP
