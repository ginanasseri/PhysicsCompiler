/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * Physics Wizard: PhysLang compiler 
 *
 * A source-to-source compiler which executes PhysLang programs by translating
 * their contents into an equivalent Python program. 
 * 
 * Takes a .phys file as an argument and generates a Python script of the 
 * same filename with a .py extension, containing an equivalent tranlation
 * of the PhysLang program into Python. 
 * 
 * If the .phys file contains a valid PhysLang program, then the results 
 * are printed to the screen and the generated Python script will be in 
 * the directory. 
 * 
 * If the .phys file contains any syntax errors, then an error message with the
 * error type and line number is printed and no script is generated. 
 * 
 * The the .phys file contains no syntax errors but is missing required data,
 * a .py script is still created, but no solution is generated. Instead, An error 
 * message indicating which values are missing from the .phys file is printed.
 * Once corrected, running the program again will overwrite the .py script with 
 * the missing values and the solution will be generated. 
 *
 * Note on error handling: The compilation is performed in three phases:
 * 
 *  1. Lexical analysis
 *  2. Parsing
 *  3. Translating (any translation errors are due to Parsing errors)
 *  4. Execution 
 *
 * Errors are handled phase-to-phase. Phases 1-3 errors result in a value 
 * of 1 being returned to the main function. No script is generated. 
 * 
 * Phase 4 errors are due to missing data. If this occurs, the execution is 
 * halted and an error message indicating which data is missing is printed.
 * Note that the error is handled by the PhysWiz compiler (NOT Python).
 * 
 */

#include <cstdlib>
#include "writer.hpp"

/* -----------------------------------------------------------
 * Function: get_file_info()
 * -----------------------------------------------------------
 * Purpose: Extracts the name and extension type of a file.
 *
 * Args: 
 *       - const std::string& file: input file passed as argument
 *       - std::string& name: reference to a string to hold the name
 *       - std::string& type: reference to string to hold the extension
 *       
 */
void get_file_info(const std::string& file, std::string& name, std::string& type) {
    size_t dot_pos = file.find_last_of(".");

    if (dot_pos != std::string::npos) {
        name = file.substr(0,dot_pos);
        type = file.substr(dot_pos+1);
    }
    else {
        name = file;
        type = "";
    }
}

int main(int argc, char *argv[]) {

    // Check if a filename was provided
    if (argc < 2) {
        std::cerr << "Usage: ./main <filename>.phys" << std::endl;
        return 1;
    }
    // Verify the file has the extension .phys
    std::string file_name, file_type;
    get_file_info(argv[1], file_name, file_type);
    if (file_type != "phys") {
        std::cerr << "Invalid file type: '" << file_type << "' - Required: '.phys'" 
                  << std::endl;
        return 1;
    }
    // And exists
    FILE *input_file = fopen(argv[1], "r");
    if (!input_file) {
        std::cerr << "Failed to open input file." << std::endl;
        return 1;
    }

    // Create new parser instance to parse the .phys file
    Parser parser; 
    if (parser.parse(input_file) != 0) {
        return 1;
    }
    // If we're here, the parse was successful.


    // Create and open a file of the same name with a .py extension 
    file_name += ".py";
    std::ofstream out_file;
    out_file.open(file_name);


    // Make sure file is opened - this file will hold the PhysLang to Python translation
    if (!out_file.is_open()) {
        std::cerr << "unable to open file: " << file_name << std::endl;
        return 1;
    }
   
    // Create new writer instance to write the PhysLang to Python translation into the .py file
    Writer writer(out_file);
    writer.setParser(parser);
    if (writer.write() != 0) {
        out_file.close();
        return 1;
    }
    // If we're here, write was successful
    out_file.close();  

    // Execute the Python script
    std::string command = "python " + file_name;
    int result = system(command.c_str());

    if (result != 0) {
        std::cerr << "Error executing Python script: " << file_name << std::endl;
        return 1;
    } 
    return 0;
}
