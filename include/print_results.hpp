#ifndef PRINT_RESULTS_HPP
#define PRINT_RESULTS_HPP

const char* free_fall_printf = R"(
# ------------------------------------------------------------
#   SOLUTION SECTION ENDS HERE 
# ------------------------------------------------------------
#   PRINTING SECTION: PRINT DATA AND RESULTS 
# ------------------------------------------------------------
def print_parameters(local_vars):
    """
    Prints initial conditions and parameters used to calculate the solution in a table.
    Only relevant data is printed to the table: optional parameters are printed only if an argument was provided for them.
    """
    equation = local_vars['equation']
    forces_map = {
        1: "Drag: constant | Gravity: constant",
        2: "Drag: non-constant | Gravity: constant",
        3: "Drag: constant | Gravity: non-constant",
        4: "Drag: non-constant | Gravity: non-constant"
    }
    forces = forces_map.get(equation, "Unknown case")
    print(f"\nFree Fall Solution: Case {equation}.")
    print("-" * len(forces))
    print(forces)
    print("-" * len(forces))

    # Print data table
    print("Data:")
    keys = ['y0', 'm', 'g', 'c']
    units = ['m', 'kg', 'm/s^2', 'kg-^1']

    params = ['H', 'R', 'M', 'rho']
    defaults = [8000, 6.3781e6, 5.97219e24, 1.225]
    default_units = ['m', 'm', 'kg', 'kg/m^3']

#    if local_vars['using_k']:
#        keys.append('k')
#        units.append('kg/m')
#    else:
    if not local_vars['using_k']:
        keys += ['Cd', 'A']
        units += ['', 'm^2']

    for param, default, unit in zip(params, defaults, default_units):
        if local_vars[param] != default:
            keys.append(param)
            units.append(unit)

#    for key, unit in zip(keys, units):
#        value = local_vars[key]
#        print(f"{key:>9} = {value:<9.2f} {unit}")
#    print("")

    for key, unit in zip(keys, units):
        value = local_vars[key]
        if value > 99999 or value < 0.01:
            print(f"{key:>9} = {value:<9.2e} {unit}") # print large values in exponential form
        else:
            print(f"{key:>9} = {value:<9.2f} {unit}")
    print("")
        

def print_free_fall_results(data):
    """
    Prints the results of the solution: the time of impact and velocity data. The velocity data
    that is printed depends on if the terminal velocity was found. 
    
    If terminal velocity was found by the time of impact, then the maximum and terminal velocity 
    are printed. If the terminal velocity was not found, thenn the maximum and final velocity 
    are printed. 
    
    """
    # Table headers 
    headers = ["Event", "Value", "Time (s) "]

    # Column widths
    col_widths = [max(len(str(value)) for value in column) for column in zip(*([headers] + data))]

    # Print table: headers 
    print("-" * (sum(col_widths) + 6))
    print("   * * * * * * Results * * * * * *")
    print("-" * (sum(col_widths) + 6))
    header_str = " | ".join(f"{headers[i]:<{col_widths[i]}}" for i in range(len(headers)))
    print(header_str)
    print("-" * (sum(col_widths) + 6))

    # Print table: rows
    for row in data:
        row_str = " | ".join(f"{row[i]:<{col_widths[i]}}" for i in range(len(row)))
        print(row_str)
    print("-" * (sum(col_widths) + 6))
    
    
def print_velocity_message(c, converged):
    """
    Prints message on why terminal velocity was not found (if not found). 
    """
    print("Note  : -> Terminal velocity not found.")
    if c == 0:
        print("Reason: -> Not applicable in the absence of drag force.")
        
    elif not converged:
        print("Reason: -> No acceleration found within +/- 0.02 of zero.")
            
    else:
        print("Reason: -> Not reached by time of impact")


# ------------------------------------------------------------
#   PRINTING SECTION ENDS HERE
# ------------------------------------------------------------
#   MAIN FUNCTION STARTS HERE
# ------------------------------------------------------------)";

const char* print_free_fall_params = R"(
def print_parameters(local_vars):
    """
    Prints initial conditions and parameters used to calculate the solution in a table.
    Only relevant data is printed to the table: optional parameters are printed only if an argument was provided for them.
    """
    equation = local_vars['equation']
    forces_map = {
        1: "Drag: constant | Gravity: constant",
        2: "Drag: non-constant | Gravity: constant",
        3: "Drag: constant | Gravity: non-constant",
        4: "Drag: non-constant | Gravity: non-constant"
    }
    forces = forces_map.get(equation, "Unknown case")
    print(f"\nFree Fall Solution: Case {equation}.")
    print("-" * len(forces))
    print(forces)
    print("-" * len(forces))

    # Print data table
    print("Data:")
    keys = ['y0', 'm', 'g']
    units = ['m', 'kg', 'm/s^2']

    params = ['H', 'R', 'M', 'rho']
    defaults = [8000, 6.3781e6, 5.97219e24, 1.225]
    default_units = ['m', 'm', 'kg', 'kg/m^3']

    if local_vars['using_k']:
        keys.append('k')
        units.append('kg/m')
    else:
        keys += ['Cd', 'A']
        units += ['', 'm^2']

    for param, default, unit in zip(params, defaults, default_units):
        if local_vars[param] != default:
            keys.append(param)
            units.append(unit)

    for key, unit in zip(keys, units):
        value = local_vars[key]
        if value > 99999 or value < 0.01:
            print(f"{key:>9} = {value:<9.2e} {unit}") # print large values in exponential form
        else:
            print(f"{key:>9} = {value:<9.2f} {unit}")
    print("")



)";

const char* projectile_printf = R"(
# ------------------------------------------------------------
#   PRINTING SECTION: PRINT DATA AND RESULTS 
# ------------------------------------------------------------
def print_parameters(local_vars, using_c):
    title = "Projectile Motion Solution"
    print()
    print(title)
    print("-"*len(title))
    
    # Print data table
    print("Data:")
    keys = ['v0', 'theta', 'g', 'c']
    units = ['m/s','degrees', 'm/s^2', 'kg^-1']
    
#    if using_c:
#        keys.append('c')
#        units.append('kg^-1')

    # print m, Cd, A if used to calculate c
    if not using_c:
        keys += ['m', 'Cd', 'A']
        units += ['kg', '', 'm^2']
    

    params = ['x0', 'y0', 'R', 'M']
    defaults = [0, 0, 6.3781e6, 5.97219e24]
    default_units = ['m', 'm', 'm', 'kg']

    # print M, R, x0, y0 only if values were provided
    for param, default, unit in zip(params, defaults, default_units):
        if local_vars[param] != default:
            keys.append(param)
            units.append(unit)
            
    for key, unit in zip(keys, units):
        value = local_vars[key]
        if value > 99999 or value < 0.01:
            print(f"{key:>9} = {value:<9.2e} {unit}") # use exp. form for very large and very small values
        else:
            print(f"{key:>9} = {value:<9.2f} {unit}")
    print("")


def print_projectile_results(data):
    """
    Prints the results of the solution: the time of impact, max height time, max distance, and max height.
    """
    # Table headers
    headers = ["Event", "Value (m) ", "Time (s) "]

    # Convert the data dictionary to a list of lists
    data_list = [
        ["Time of Impact", f"N/A", f"{data['time_of_impact']:.2f}"],
        ["Max Height", f"{data['max_height']:.2f}", f"{data['max_h_t']:.2f}"],
        ["Max Distance", f"{data['max_distance']:.2f}", f"{data['time_of_impact']:.2f}"]
    ]

    # Column widths
    col_widths = [max(len(str(value)) for value in column) for column in zip(*([headers] + data_list))]

    # Print table: headers
    print("-" * (sum(col_widths) + 6))
    print("  * * * * * * Results * * * * * *")
    print("-" * (sum(col_widths) + 6))
    header_str = " | ".join(f"{headers[i]:<{col_widths[i]}}" for i in range(len(headers)))
    print(header_str)
    print("-" * (sum(col_widths) + 6))

    # Print table: rows
    for row in data_list:
        row_str = " | ".join(f"{row[i]:<{col_widths[i]}}" for i in range(len(row)))
        print(row_str)
    print("-" * (sum(col_widths) + 6))

# ------------------------------------------------------------
#   PRINTING SECTION ENDS HERE
# ------------------------------------------------------------
#   MAIN FUNCTION STARTS HERE
# ------------------------------------------------------------)"; 

#endif // PRINT_RESULTS_HPP
