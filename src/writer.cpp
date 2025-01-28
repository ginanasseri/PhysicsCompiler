#include "writer.hpp"
#include "outputs.hpp"

Writer::Writer(std::ofstream& out) : py_script(out) {}

/* =================== public section ========================= 
 * Function: write()
 * ------------------------------------------------------------
 * Purpose: Serves as the entry way into the writing process. 
 *          Handles the graph import section then passes control
 *          to the subject handler.  
 *
 * Returns: int - 0 if write was successful 
 *                1 otherwise. 
 */
int Writer::write() {
    if (writeGraphImport() != 0) {
        std::cerr << "Parser Error: getGraphValue() returned invalid graph value." << std::endl;
        return 1;
    }
    if (writeSubjectScript() != 0) {
        std::cerr << "Solution script generator failed." << std::endl;
        return 1;
    }
    return 0;
}

/* function: setParser() sets the parser data member */
void Writer::setParser(const Parser& parser) {
    this->parser = parser;
}

/* ================ private seciton ========================== 
 * -----------------------------------------------------------
 * Function: writeToScript()
 * -----------------------------------------------------------
 * Purpose: writes Python code into the writer instance's 
 *          std::ofstream py_script
 *
 * Args: 
 *       const char* python_code: string of Python code
 *       
 */
void Writer::writeToScript(const char* python_code) {
    py_script << python_code << std::endl;
}

/* ===================================================================
 *           * * *  WRITING ROUTINES START HERE * * * 
 * ===================================================================
 * 
 * Note on unreachable code sections:
 *       
 * If an unreachable code section is entered at any point, then and
 * invalid value was returned from the parser, meaning the parser
 * accepted an invalid token. This should never happen unless something
 * went seriously wrong. 
 *------------------------------------------------------------------- */

/* ------------------------------------------------------------------ 
 * Function: writeGraphImport()                         GRAPH
 * ------------------------------------------------------------------ 
 * Purpose: writes the matplotlib.pyplot import statement to 
 *          the script if the graph value is `true`
 * 
 *          If the graph value is neither true nor false, then 
 *          an invalid graph token was parsed.
 *
 * Returns: 1 if an invalid graph token was parsed. 
 *          0 otherwise.
 */
int Writer::writeGraphImport() {
    std::string graph_value = parser.getGraphValue();
    if (graph_value == "true") {
        writeToScript(graph_imports);
        return 0;
    }
    if (graph_value == "false") {
        return 0;
    }
    std::cerr << "Error: writeGraphImport() in unreachable code section:" << std::endl;
    return 1;
}
/* ------------------------------------------------------------------------- 
 * Function: writeSubjectScript()                     SUBJECT
 * ------------------------------------------------------------------------- 
 * Purpose: Validates the subject token value and dispatches to the 
 *          writer routine for the subject.
* 
 * Returns: 1 if an invalid subject token was parsed. Or
 *          if an error occurred in writing the script.
 *          0 otherwise.
 */
int Writer::writeSubjectScript() {
    std::string subject = parser.getSubject();

    if (subject == "mechanics") {
        if (writeMechanicsScript() != 0) {
            return 1;
        }
        return 0;
    }
  if (subject == "stats") {
        if (writeStatsScript() != 0) {
            return 1;
        }
        return 0;
    }
    std::cerr << "Error: writeSubjectScript() in unreachable code section: getSubject() returned invalid subject type." << std::endl;
    return 1;
}

/* ------------------------------------------------------------------ 
 * Function: writeGraphGen()                           GRAPHS 
 * ------------------------------------------------------------------ 
 * Purpose: checks if the graph setting is true and writes the python code 
 *          which produces the plots into the script.
 * 
 *    - const char* plot_generator: Python code segment of either a Python 
 *      function which generates a plot or the function calling statement. 
 */
void Writer::writeGraphGen(const char* plot_generator) {
    if (parser.getGraphValue() == "true") {
        writeToScript(plot_generator);
    }
}

/* ==================================================================
 *                   Mechanics Domain
 * ==================================================================
 * ------------------------------------------------------------------ 
 * Function: writeMechanicsScript()
 * ------------------------------------------------------------------ 
 * Purpose: Imports the Python libraries used in the mechanics 
 *          solutions, and adds functions used to find missing 
 *          required data and assign default values to optional
 *          args without provided parameters to the script. 
 * 
 *  If writeMehchanicsEquation() returns 1, then the equation
 *  number on line 3 of the .phys file was invalid, or an 
 *  error occurred during the solution write.   
 * 
 */
int Writer::writeMechanicsScript() {
    writeToScript(mechanics_imports);
    writeToScript(missing_data_check);

    // Mechanics Equation Dispatch table:
    if (writeMechanicsEquation() != 0) {
        return 1;
    }
    return 0;    
    
}

/* -------------------------------------------------------------------- 
 * Function: writeMechanicsData()                 MECHANICS DATA
 * -------------------------------------------------------------------- 
 * Purpose: Writes the variable assignment statements from the 
 *          data section in the .phys file into equivalent 
 *          assignment statements in Python. E.g.:
 * 
 *          data: {x = 12, y = 13}  is translated into
 * 
 *          # Your data:
 *          x = 12
 *          y = 13
 * 
 */
void Writer::writeMechanicsData() {
    std::string values = parser.getData();
    writeToScript("    # Your data:");
    py_script << values << std::endl;
}

/* -------------------------------------------------------------------- 
 * Function: writeEquationNumber()                WRITE EQUATION
 * -------------------------------------------------------------------- 
 *  Adds the equation number to the Python script following the 
 *  data write above. E.g., using the same example above, if the 
 *  equation number is 4, then the following code is generated:
 * 
 *          # Your data:
 *          x = 12
 *          y = 13
 *          equation = 4
 *
 *  Note that exactly one indentation is used as the data and 
 *  equation number are written into the Python script's main
 *  function. 
 *  
 */
void Writer::writeEquationNumber() {
    int id = parser.getEquationID();
    std::string equation = "    equation = ";
    equation += std::to_string(id); 
    py_script << equation << std::endl;
}


/* -------------------------------------------------------------------  
 * Function: writeMechDataAndEq()             MECH DATA AND EQ 
 * ------------------------------------------------------------------- 
 * Purpose: Triple checks that we're where we're supposed to be. 
 *
 *  First checks that the current equation ID is the equation ID
 *  we're expecting (otherwise, we should not be here), and then 
 *  writes the data and equation number into the script. 
 * 
 *  Returns: 0 if equation ID is valid 
 *           1 if equation ID is not valid. (This should not happen.)
 *
 * This fuction is only used for single value evaluation. The free
 * fall equations have a separate validation function as the range 
 * of the equations is 1-4 and not a single value. 
 *             
 */
int Writer::writeMechDataAndEq(int expected_equation) {
    int id = parser.getEquationID();
    if (id != expected_equation) {
        std::cerr << "Write Error: writeMechData: unexpected equation number: " << id
                  << ". Expected: " << expected_equation;
        return 1;
    }
    writeMechanicsData();
    writeEquationNumber();
    return 0;
}

/* ------------------------------------------------------------------- 
 * Function: writeFreeFallEquationID()    FREE FALL EQUATION ID
 * ------------------------------------------------------------------- 
 * Purpose: Triple checks that we're where we're supposed to be. 
 *
 * Ensures that writeFreeFallEquation hasn't been entered in error and that we 
 * indeed have a free fall equation ID number before writing it in to the script. 
 * 
 * If we're here, we've been called from a function that was entered only if 
 * the equation ID is 1,2,3,or 4. So the invalid equation ID error statement
 * should never be reached. 
 *
 *  Equation:
 *            1 - constant drag and gravitational forces
 *            2 - variable drag and constat gravitational forces 
 *            3 - constant drag and variable gravitational forces
 *            4 - variable drag and gravitational forces
 * 
 *  Returns: 0 if equation ID is valid 
 *           1 if equation ID is not valid. (This should not happen) 
 *
 */
int Writer::writeFreeFallEquationID() {
    int id = parser.getEquationID();
    if (id > 0 && id < 5) {
        writeEquationNumber();
        return 0;
    }
    std::cerr << "Error: writeFreeFallSolution() entered with invalid equation ID: " << id << std::endl;
    return 1;
}



/* ------------------------------------------------------------------- 
 * Function: writeMechanicsEquation()        MECHANICS EQUATION
 * ------------------------------------------------------------------- 
 * Purpose: Verifies that the equation number in the equation: 
 *          line (line 3) from the input file is valid, and 
 *          dispatches to the routine corresponding to the  
 *          equation ID. 
 *  
 * Returns: 1 if the equation ID is invalid or if a write error
 *          occured in the dispatch routine. 0 otherwise. 
 * 
 */
int Writer::writeMechanicsEquation() {

    // get the differential equation to solve
    int id = parser.getEquationID();

    // free fall differential equations
    if (id > 0 && id < 5) {
        if (writeFreeFallSolution() != 0) {
            return 1;
        }
        return 0; 
    }
    // projectile differential equation
    if (id == 5 || id == 6) {
        if (writeProjectileSolution() != 0 ) {
            return 1; 
        }
        return 0;
    }
    // pendulum 

    std::cerr << "Error: Invalid mechanics equation number: " << id  
              << "(line 3). See manual for table of equations." << std::endl;
    return 1;   
}
/* ------------------------------------------------------------------ 
 * Function: writeGravityAndDrag()           GRAVITY AND DRAG
 * ------------------------------------------------------------------ 
 *  Purpose: Routine that calculates the drag force constant,
 *           sets the parameters for the planet to their default 
 *           Earth values or values provided inthe file, and 
 *           calculates gravity using these parameters.
 *
 *           Also checks that all parameters needed to calculate
 *           the drag force constant were provided. 
 */
void Writer::writeGravityAndDrag() {
    writeToScript(drag_data_check); 
    writeToScript(set_gravitational_constant);
    writeToScript(calculate_drag_constant);
    writeToScript(set_planet_params);
}


/* ------------------------------------------------------------------ 
 * Function: startMechanicsMain()      MECHANICS MAIN FUCNTION
 * ------------------------------------------------------------------ 
 * Purpose: writes the def main(): into the script then initializes 
 *          parameters common to all mechanics problems (y0, v0, g, etc.)
 *          The parameters are all initialized with a value of None and 
 *          then written over by any values provided in the data: {}
 *          section in the .phys file. (Used to check which values were
 *          provided.) and double checks that the writer is where it 
 *          should be.
 * 
 */
int Writer::startMechanicsMain() {
    writeToScript(start_main_function);
    if (parser.getSubject() != "mechanics") {
        std::cerr << "Error: writeMechanicsMain() entered with invalid subject type." << std::endl;
        return 1;
    }
    writeToScript(initialize_mechanics_params);
    return 0;
}

/* ------------------------------------------------------------------ 
 * Function: writeFreeFallSolution()                      FREE FALL
 * ------------------------------------------------------------------ 
 *  Purpose: Writes python script which solves the free-fall problem 
 *  for a body falling through an  atmosphere under various
 *  conditions (constanct/non-constat drag and/or gravitational
 *  forces) using the initial conditions and paramters provided 
 *  in the source code.
 * 
 *  Returns: 1 if the routine was entered in error (invalid equation ID)
 *           0 otherwise. (Write was succesful)
 * 
 */
int Writer::writeFreeFallSolution() { 
    
    // -- functions for generating solution -- 
    writeGravityAndDrag();               // set planet parameters, calculate the drag and grav. constants
    writeToScript(free_fall_ODE_solver); // contains ODE and solver functions 
    writeGraphGen(free_fall_plots);      // writes graphing function if graph is true
    writeToScript(free_fall_printf);     // for printing results to table 
        
    // --- main function begins here ----
    if (startMechanicsMain() != 0) {     // start main function and initialize all variables as None
        return 1;
    }
    writeMechanicsData();                   // write the data (the variable assignment statements) from .phys file
    if (writeFreeFallEquationID() != 0) {   // and the free fall case number.
        return 1;
    }    
    writeToScript(get_free_fall_solution); // update variables with the user's data and call the ODE solver routine
    writeGraphGen(generate_ff_plots);      // call plot generation function if graph is true 
    writeToScript(main_entry_point);       // entry point to the main function (if __name__ == "__main__...)

    return 0; 
}


/* ------------------------------------------------------------ 
 * Function: writeProjectileSolution()           PROJECTILE
 * ------------------------------------------------------------
 * Purpose: Writes the Python script which solves the projectile motion
 *          with drag problem using the parameters provided in the 
 *          .phys file's data section. 
 * 
 *  Returns: 1 if the routine was entered in error (invalid equation ID)
 *           0 otherwise. (Write was succesful)
 */
int Writer::writeProjectileSolution() {
    // -- functions for generating solution -- 
    writeGravityAndDrag();                // set planet parameters, calculate the drag and grav. constants
    writeToScript(get_init_velocites);    // calculate the initial velocity components from the initial conditons
    writeToScript(projectile_ODE_solver); // contains ODE and solver functions 
    writeGraphGen(projectile_plot);       // writes graphing function if graph is true
    writeToScript(projectile_printf);     // for printing results to table 
  
    // -- main function begins here --
    if (startMechanicsMain() != 0) {     // start main function and is mechanics initialize all variables as None
        return 1;
    }
    writeToScript(init_proj_params);      // initialize projectile specific variables: `theta` and `v0`
    if (writeMechDataAndEq(parser.getEquationID()) != 0) { // write user data and equation number into script 
        return 1;                                          
    }
    writeToScript(get_projectile_solution); // update variables with the user's data and call the ODE solver routine
    writeGraphGen(generate_proj_plots);     // call plot generation function if graph is true
    writeToScript(main_entry_point);        // entry point to the main function (if __name__ == "__main__...)
    
    return 0; 
}

/* ==================================================================
 *                       Stats Domain
 * ==================================================================
 * Function: writeFileName()                          STATS DATA
 * ------------------------------------------------------------------
 * Purpose: Writes the file name provided in the data section in
 *          the .phys file into the Python script. E.g., 
 *          
 *    data: {my_file.csv} is translated into
 *       
 *    # Your data file:
 *    file = 'my_file.csv'
 */
void Writer::writeFileName() {
    std::string file_name = "    data_file = \'";
    file_name += parser.getData();
    file_name += "\'";
    writeToScript("    # Your data file:");
    py_script << file_name << std::endl;
}

/* ------------------------------------------------------------------
 * Function: writeStatsScript()                   
 * ------------------------------------------------------------------
 * Purpose: Import Python libraries needed to generate the stats 
 * solution (pandas, stats, etc) and ensure the equation number is valid.
 */
int Writer::writeStatsScript() {
    py_script << stats_imports << std::endl; // import libraries 
    if (writeStatsEquation() != 0) {
        return 1;
    }
    return 0;
}

/* ------------------------------------------------------------------
 * Function: writeStatsScript()                    DISPATCH TABLE
 * ------------------------------------------------------------------
 * Purpose: Call the Python script corresponding to the eqaution number
 *          (currently only one is implemented).
 */
int Writer::writeStatsEquation() {
    int id = parser.getEquationID();
    if (id != 1) { 
        std::cerr << "Error: Invalid stats equation number: " << id
        << " (line 3). See manual for table of equations." << std::endl;
        return 1;
    }
    if (writeDescriptiveStats() != 0) {
        return 1;
    }   
    return 0;
}

/* ------------------------------------------------------------------
 * Function: writeDescriptiveStats()              DESCRIPTIVE STATS
 * ------------------------------------------------------------------
 * Purpose: Writes script to print and, if requested, plot the desc-
 * criptive stats of a data file with columns x,y. The descriptive
 * stats are the standard deviation, max, and mean of each column along
 * with a histogram of the frequency with an overlaying probability 
 * distribution function and the correlation and covariance coeff-
 * cient between the two data sets and, if requested, a scatter plot
 * of the data sets with y(x).
 */    
int Writer::writeDescriptiveStats() {
    // --- functions for generating solution ---
    writeToScript(descriptive_stats);       // std. dev, mean, max, covariance, correlation
    writeGraphGen(descriptive_stats_plots); // histogram, PDF, scatter plot if graph requested
    
    // ----- main function -----
    writeToScript(start_main_function);
    writeFileName();                 // write data file to script
    writeToScript(get_desc_stats);   // call descriptive stats function with file 
    writeGraphGen(get_desc_plots);   // call plotting functions if graph is true
    writeToScript(main_entry_point); // entry point to the main function (if __name__ == "__main__...)
    return 0;
}
