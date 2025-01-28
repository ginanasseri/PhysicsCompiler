#ifndef OUTPUTS_HPP
#define OUTPUTS_HPP

#include "imports.hpp"
#include "verify_user_data.hpp"
#include "initialize_parameters.hpp"
#include "mechanics_solutions.hpp"
#include "stats_solutions.hpp"
#include "solution_plots.hpp"
#include "print_results.hpp"


//#include "parameters.hpp"
//#include "initialize_parameters.hpp"
//#include "print_results.hpp"
//#include "solution_scripts.hpp"
//#include "solution_plots.hpp"

const char* start_main_function = R"(
def main():)";

const char* main_entry_point = R"(if __name__ == "__main__":
    main())";

#endif // OUTPUTS_HPP
