#ifndef WRITER_HPP
#define WRITER_HPP

#include <fstream>
#include "parser.hpp"

class Writer {

public: 
    Writer(std::ofstream& out);
    int write();    
    void setParser(const Parser& parser); 

private: 
    void writeToScript(const char* python_code);
    int writeGraphImport();
    void writeGraphGen(const char* plot_generator);
    int writeSubjectScript();
    void writeEquationNumber();

    // == mechanics domain ==
    void writeMechanicsData();
    int writeMechanicsScript();
    int writeMechDataAndEq(int expected_equation);
    int writeMechanicsEquation();
    void writeGravityAndDrag();
    int startMechanicsMain();
    
    int writeFreeFallEquationID();
    int writeFreeFallSolution();   // equations 1-4
    int writeProjectileSolution(); // equation 5

    // == stats domain ==
    int writeStatsScript();
    int writeStatsEquation();
    void writeFileName();
    int writeDescriptiveStats();  // equation 6

    Parser parser;
    std::ofstream& py_script;
};

#endif // WRITER_HPP
