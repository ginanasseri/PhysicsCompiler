#ifndef VERIFY_USER_DATA_HPP
#define VERIFY_USER_DATA_HPP

/* -------------
 * In this file: 
 *  
 *  missing_data_check : returns default values for keys which have a value of None
 *  find_missing_values: returns set difference between required and provided paramters
 *  check_initial_conditions: checks if any initial condition keys needed for the equation have a value of None
 *  drag_data_check: checks if sufficient data to calculate drag was provided, also includes entry into data checking 
 *                   routine and generates error messages for any values missing.
 */

// General functions to check if any data required for a mechanics Python solution scripts is missing
const char *missing_data_check= R"(
# --------------------------------------------------------------
#  DATA VERIFICATION SECTION: CHECKING FOR MISSING VALUES
# -------------------------------------------------------------- 
def get_value(variable_name, default, variables):
    """
    Returns the value of the key `variable_name` in the `variables` dictionary. If the value of the key 
    is None, then the value `default` is returned. 
    
    Args: 
        - `variable_name` (str)    : Name of key to check  
        - `default` (float or None): Value to return if the key `variable_name` has a value of None. 
        - `variables` (dict)       : Dictionary of the local variables in the main function and their values
    Returns:
        float : value of the key if the key has a non-None value OR 
        float or None: default value if the key has a value of None (keys may have a default value of None)
                       
    """
    value = variables.get(variable_name)
    return value if value is not None else default


# --------  FIND MISSING VALUES -------- 
def find_missing_values(required_variables, user_data):
    """
    Checks that all values needed to generate the solution were provided in the source code's data section.
    
    Creates a set `provided_variables`, containing the keys (variable names) in the `user_data` dict with  
    non-None values, and returns the set difference between the `required_variables` and `provided_variables` 
    sets as a list (i.e., returns a list of the required variables which were not provided). 
    
    If all required variables were provided, the list will be empty. 

    Args: 
        `required_variables` (set): contains the names of variables needed to generate the solution. 
        `user_data` (dict): dictionary of required variables and their provided values (None or float) 

    Returns:
        list: A list of the missing variable names. 
    """
    provided_variables = {parameter for parameter, value in user_data.items() if value is not None}
    return list(required_variables - provided_variables)


# -------- FIND MISSING INITIAL CONDITIONS -------- 
def check_initial_conditions(user_data, equation):
    """
    Checks if the initial conditions for the initial conditions needed for computing the integration
    constants in the differential equation solution were provided in the source code.

    If any initial conditions are missing then there is insufficient data to generate the solution. 
    
    Args: 
        `user_data` (dict): dictionary of required variables and their provided values (None or float) 
        
    Returns:
        list: A list of of variables with None values. An empty list is returned is neither are missing. 
    """
    
    if equation > 0 and equation < 5:
        missing_values = find_missing_values({'y0', 'm'}, user_data)

    elif equation == 5 or equation == 6:
        missing_values = find_missing_values({'v0', 'theta'}, user_data)
    
    # future implementations
    
    else:
        # This section should never be reached unless compiler error occurred (ignore this)
        print(f"Script generator error: Invalid equation number '{equation}' written to script.")
        sys.exit(1)
        
    if missing_values:
        return list(missing_values)
    
    # no missing values, return an empty list
    return []
    
)";


// functions which test for missing values unique to the free-fall problem 
const char* drag_data_check = R"(
# --------  FIND MISSING DRAG DATA -------- 
def check_drag_data(user_data, equation):
    """
    Checks if all variables needed to calculate the drag force component were provided the source code's 
    data section.
    
    Free-Fall:
    To calculate the drag force component, the drag constant `k` or both the drag coefficient `Cd` and 
    cross-sectional area `A` must be specified. If `k` has a value of None, and `Cd` and/or `A` have a 
    a None value, in the `user_data` dictionary, then there is insufficient data to calculate the drag 
    acceleration component, and the solution cannot be generated.

    Projectile:
    The same logic applies, but if 'c' was not specified, then 'Cd', 'A', and 'm' must be specified to 
    calculate the drag constant. 
    
    Args: 
        `user_data` (dict): dictionary of required variables and their provided values (None or float)
        
    Returns:
        list: empty or containing missing data.  
        
    Logic: Free-Fall
        1.  If `k` was provided, then an empty list is returned. Likewise, if `k` was not provided, 
            but both `Cd` and `A` were provided, an empty list is returned. 
        2.  If `k` was not provided and `Cd` and/or `A` were not provided, then a list containing 
            `k` along with the other missing parameter(s) (`Cd` and/or `A`) is returned. 
    
    Logic: Projectile - same logic applies; if 'c' was not provided, then 'Cd', 'A', 'm' must be provided.
    """
    # ----- projectile motion -----
    if equation == 5 or equation == 6:
        
        if user_data.get('c') is not None:
            return []
        
        missing_values = find_missing_values({'Cd', 'A', 'm'}, user_data) 
        if missing_values: 
            missing_values.append('c')
            return missing_values

    # -------- free-fall ---------
    if user_data.get('k') is not None:
        return []
    
    missing_values = find_missing_values({'Cd', 'A'}, user_data)
    if missing_values: 
        missing_values.append('k')
        return missing_values
       
    return []


# --------  ENTRY TO DATA CHECKING FUNCTIONS -------- 
def missing_data_check(user_data, equation):
    """
    Checks if any parameters required to generate the solution were not provided in the source code by
    looking up their keys in the `user_data` dictionary. If any initial conditions or drag data have a 
    value of None, an error message is generated indicating what values are missing, and the program 
    exits. 
    
    Note: This is a data entry error from the PhysLang program, not a Python error. It is handled by the 
    PhysWiz data verification section, which generates the error messages. Therefore, return 0, is used 
    to exit the program, rather than return 1, as no Python errors occurred
    
    Args: 
        `user_data` (dict): dictionary of required variables and their provided values (None or float)                            
        
    Returns:
        list: A list of the required variable names with missing values. If no missing values are found, 
              then an empty list is returned.
    """
    missing_drag = check_drag_data(user_data, equation)
    missing_inits = check_initial_conditions(user_data, equation)
    missing_params = missing_drag + missing_inits
    provided_params = list({parameter for parameter, value in user_data.items() if value is not None})
    
    if equation > 4:
        drag_constant = 'c'
        needed_values = "Cd, A, and m" # projectile 
    else:
        drag_constant = 'k'
        needed_values = "Cd and A"  # free fall 
    
    if missing_params:
        errors = []

        if missing_inits:
            errors.append(f"Missing initial condition(s) : {', '.join(missing_inits)}.",)

        if missing_drag:
            errors.append(f"Missing drag data: {', '.join(missing_drag)}\n      Note: if {drag_constant} is not known, then values for {needed_values} must be provided")

        error_message = "\nError: Solution generation failed:\n\n  Details:\n"

        for error in errors:
            error_message += f"    - {error}\n"
        error_message += f"\n  Summary:\n    - Missing values for required parameter(s): {', '.join(missing_params)}\n"
        error_message += f"    - Data you provided: {', '.join(provided_params)}\n"
        error_message += f"\n   Suggested: Add missing values and run the program again (see manual for additional help)."
        
        print(error_message)
        return missing_params

    return []

# -------------------------------------------------------------
#   DATA VERIFICATION SECTION ENDS HERE 
# -------------------------------------------------------------
#  PARAMETER INITIALIZATION SECTION: ASSIGNING PARAMETER VALUES
# -------------------------------------------------------------)";

#endif // GET_USER_DATA_HPP
