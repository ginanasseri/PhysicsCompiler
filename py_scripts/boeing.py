import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
import scipy.constants as const
import sys


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
# -------------------------------------------------------------
# --------  CALCULATE GRAVITY WITH GIVEN PARAMETERS -------- 
def set_gravitational_constant(M, R, g=None):
    """
    Returns the correct value for g based on the optional parameters provided in the source code. 
    
    The gravitational constant g, planet mass M, and planet radius R are optional parameters with 
    default values set to Earth's values. The default value for g can be overwritten EITHER by 
    providing a value for g, or by providing values for M and R.
    
    If a value for g is provided, that value is returned. Otherwise, if values for M and R (other than 
    Earth's values) were provided, then the gravitational constant for a planet with mass M and radius R 
    is automatically calulated and returned. Otherwise, the default value of Earth's g is returned
    
    Args:
        `M` (float): Mass of the planet (kg).
        `R` (float): Radius of the planet (m).
        `g` (float or None): Value for g provided in source code.
    
    Returns:
        float: The gravitational constant. The absolute value is returned to ensure the sign is 
               consistent with the equation of motion. 
        
    Logic:
        1. if `g` is None and M and R are equal to the default Earth values, return Earth's 
           gravitational constant. 
        2. if `g` is None and M and R do not equal the default Earth values, calculate and return
           the gravitational constant for a planet with mass M and radius R. 
        3. if `g` is not None, return that value. 
    """
    # Earth values:
    M_earth = 5.97219e24  # mass (kg) 
    R_earth = 6.3781e6    # radius (m) 
    g_earth = 9.81        # gravity (m/s^2)
    G = const.G           # universal gravitational constant (m^3 kg^-1 s^-2)
    
    if g is None:
        if M == M_earth and R == R_earth:
            return g_earth
        else:
            return G * M / R**2 # gravity of planet with mass M and radius R  
    else: 
        return abs(g)


# --------  CALCULATE DRAG CONSTANT FROM GIVEN PARAMETERS -------- 
def check_H_value(local_vars, g, H):
    """
    If the equation case specifies non-constant drag and (non-Earth) parameters were provided for M, R
    or g, but no parameter was provided for the atmospheric scale height H, then a message is printed 
    indicating that an Earth value was used in the solution and to provide the correct H value for the
    planet. (The solution will still be generated)
    """
    g_earth = 9.81
    H_earth = 8000  # Earth atmospheric scale height (m)
    
    equation = local_vars['equation']
    
    # constant drag, H doesn't matter 
    if equation < 3 or equation == 5:
        return
    
    # we're on earth
    if g == g_earth:
        return
    
    # we're not on earth, check if H for the planet was provided 
    if H == H_earth:
        print(f"\n - Note: No atmospheric scale height H was provided for planet with g = -{g:.2f}\n         Using atmospheric scale height of Earth: H = 8 km")
        
        
def check_rho_value(local_vars, g, rho):
    """
    If (non-Earth) parameters were provided for M, R, or g, and the drag constant was calculated using 
    Cd, A, m, and rho, but no parameter was provided for rho, then a message is printed indicating that
    indicating that an Earth value was used in the solution and to provide the correct rho value for the
    planet. (The solution will still be generated)
    """
    g_earth = 9.81
    rho_0 = 1.225   # Earth atmospheric density at sea level (kg/m^3)
    
    # we're on earth
    if g == g_earth:
        return
    
    # we're not on earth, check if rho for the planet was provided 
    if rho == rho_0:
        print(f"\n - Note: No air density rho for planet with g = -{g:.2f} was provided.\n         Using surface density of Earth: rho = 1.225 kg/m^3.")


def calculate_drag_constant(k, Cd, A, m, rho, H, g, local_vars):
    """
    Calculates the drag constant from the provided data. The drag force constant:
    
                                c = 1/2 * Cd * A * rho / m
            
    Where Cd and A are the drag force coefficient and cross-sectional area of the body, and rho is the 
    surrounding air density. The drag force constant is a measure of how aerodynamic a body is. The 
    higher the drag force constant, the less aerodynamic the body is. 
    
    Note: rho is an optional parameter with a default value of Earth's value (1.225 kg/m^3)
    
    Note: The constant is divided by the mass m, as it is used to to compute the drag acceleration 
          component in the solution (F/m = a). 
    
    Args:
        `k` (float or None) : Drag constant (kg/m)
        `Cd` (float or None): Drag coefficient (dimensionless)
        `A` (float or None) : Cross-sectional area of body in free-fall (m^2)
        `m` (float)         : Mass of body in free-fall (kg).
        `rho` (float)       : Air density. Default value is Earth unless value provided (kg/m^3)
    
    Returns:
        float: The drag force constant divided by the mass 
        boool: flag which indicates if k was provided or Cd and A were provided (used for printing results) 
    """

    using_k = True
    
    if get_value('k', None, local_vars) is not None:
        c = k/m
        
    else:
        using_k = False # for printing paramters
        c = (1/2 * Cd * rho * A) / m
        check_rho_value(local_vars, g, rho) 
        
    check_H_value(local_vars, g, H) 
    return c, using_k


# -------- ASSIGN PLANET PARAMETER VALUES WITH EARTH OR GIVEN VALUES -------- 
def get_planet_parameters(local_vars):
    """
    Sets optional parameters to their default value if no value was provided. 
    
    Args:
        - `local_vars` (dict): dictionary containing the local variables and their associated values. 
    
    Returns: 
        (float) values: H, R, M, rho, g as their default values or the values provided.
    """
    # Default Earth values for optional data
    H_earth = 8000        # Atmospheric scale height (m)
    R_earth = 6.3781e6    # Radius (m) 
    M_earth = 5.97219e24  # Mass (kg) 
    rho_0 = 1.225         # Earth atmospheric density at sea level (kg/m^3)
    G = const.G           # Gravitational constant 

    # Set optional parameters to default values if no value provided in source code
    H = get_value('H', H_earth, local_vars)
    R = get_value('R', R_earth, local_vars)
    M = get_value('M', M_earth, local_vars)
    rho = get_value('rho', rho_0, local_vars)
    
    # Get the gravitational constant from source code or calculate from provided parameters
    g_provided = get_value('g', None, local_vars)
    g = set_gravitational_constant(M, R, g_provided) # returns g consistent with source code parameters
    
    return H, R, M, rho, g


# ------------------------------------------------------------
#   PARAMETER INITIALIZATION SECTION ENDS HERE
# ------------------------------------------------------------
#   SOLUTION SECTION: SOLVING THE DIFFERENTIAL EQUATION
# ------------------------------------------------------------
# --------  EQUATION OF MOTION -------- 
def equation_of_motion(t, y, g, c, H, R, equation):
    """
    Returns the equation of motion for a body in free-fall with air resistance as a system of two 
    first order ODE's, where y1 = y(t) and y2 = v(t), and with initial conditions y(0) = y0 and 
    v(0) = v0. 
    
    The equation number is a setting which determines if the drag and/or gravitational forces are 
    constant or non-constant. (Settings listed below in list of Args, see manual for more details). 
    
    This function is defined to be used with solve_ivp() (solve initial value problem) which computes the 
    integral of the system over a defined time interval to obtain the position y(t) and velocity v(t) 
    solutions, using the initial conditions to solve for the integration constants. 
    
    See manual for more details on the differential equations. 
    
    Args:
        `t` (float or array-like): Single time value or time span to evaluate over. 
        `y` (list): Initial conditions vector containing [y0,v0].
        `g` (float): Gravitational constant. 
        `c` (float): Drag force constant divided by mass.
        `H` (float): Scale height of atmosphere (used in equation 2 and 4).
        `R` (float): Planet radius (used in equation 3 and 4).
        
        equation (int): Specifies the type of gravitational and drag forces;
            1: Constant gravity, constant drag.  
            2: Constant gravity, variable drag. 
            3: Variable gravity, constant drag.
            4: Variable gravity, variable drag.
            
    Returns:
        list: the system of equations as a list [velocity y'(t), acceleration y''(t)]      
    """
    y1, y2 = y
    
    if equation == 1:
        dydt = [y2, -g - c * y2 * abs(y2)]
    elif equation == 2:
        dydt = [y2, -g - c * np.exp(-y1/H) * y2 * abs(y2)]
    elif equation == 3:
        dydt = [y2, -g / (1 + (y1 / R))**2 - c * y2 * abs(y2)]
    elif equation == 4:
        dydt = [y2, -g / (1 + (y1 / R))**2 - c * np.exp(-y1/H) * y2 * abs(y2)]
    else:
        # This section should never be reached unless compiler error occurred (ignore this)
        print(f"Script generator error: Invalid equation number '{equation}' written to script.")
    return dydt


# --------  GET SOLUTION -------- 
def get_free_fall_solution(equation_of_motion, y_initial, t_span, t_eval, g, c, H, R, equation):
    """
    Generates the solutions y(t) and v(t) to the system of equations defined above in the 
    equation_of_motion() function over the time interval t_span with intital conditions [y0,v0].
    
    This function is defined to be used with solve_ivp() (solve initial value problem) which computes the 
    integral of the system over a defined time interval to obtain the position y(t) and velocity v(t) 
    solutions, using the initial conditions to solve for the integration constants. 

    The solve_ivp() (solve initial value problem) function generates the solution to the ODE defined in 
    the equation_of_motion() above by integrating the system at each point in `t_eval` over the time
    interval `t_span`, and using the intial conditions [y0, v0] to compute the integration constants.

    Args:
        `equation_of_motion` (function): system of ODE's describing a body in free fall 
        `y_initial` (list): Initial conditions vector [y0, v0]
        `t_span` (array-like): Time span to evaluate over.
        `t_eval` (array-like): Time points to evaluate 
        `g` (float): Gravitational constant
        `c` (float): Drag force constant
        `H` (float): Atmospheric scale height (used in equation 2 and 4)
        `R` (float): Radius of planet(used in equation 3 and 4)
        
        equation (int): specifies the type of gravitational and drag forces;
            1: Constant gravity, constant drag.
            2: Constant gravity, variable drag. 
            3: Variable gravity, constant drag.
            4: Variable gravity, variable drag.
        
    Returns:
        tuple: (time_of_impact, solution) where time_of_impact is the time at which the object impacts the ground,
               and solution is the full solution of the ODE system.
    """
    while True:
        # Solve the ODE system
        solution = solve_ivp(equation_of_motion, t_span, y_initial, t_eval=t_eval, args=(g, c, H, R, equation), method='RK45', rtol=1e-8, atol=1e-10)
        
        # Extract the solution
        t = solution.t     # time
        y = solution.y[0]  # altitude
        
        # Check if the time of impact is captured
        if np.any(y <= 0):
            impact_index = np.where(y <= 0)[0][0]
            time_of_impact = t[impact_index]
            return time_of_impact, impact_index, solution
        else:
            # Double the time span and try again if not captured
            t_span = (0, t_span[1] * 2)
            t_eval = np.linspace(t_span[0], t_span[1], 10000)
            
            
# --------  CALCULATE ACCELERATION -------- 
def calculate_acceleration(g, c, v, y, H, R, equation):
    """
    Computes the net acceleration as a function of time for a body in free fall subject to gravitational
    and drag forces specified by the equation number, using the y(t) and v(t) solutions obtained from 
    the get_free_fall_solution() defined above. 
    
    Args:
        `g` (float): Gravitational acceleration
        `c` (float): Drag force constant (divided by mass)
        `v` (list) : velocity solution v(t)  
        `y` (list) : position soltuition y(t)
        `H` (float): Atmospheric scale height
        `R` (float): planet radius
        `equation` (int): specifies if gravity and/or drag is non-constant (see above function)

    Returns:
        list: the net acceleration, a(t), of the body in free fall subject to the forces specified by 
              the equation number. 
    """
    if equation == 1:
        return -g - c * v * abs(v)
    elif equation == 2:
        return -g - c * np.exp(-y/H) * v * abs(v)
    elif equation == 3:
        return -g / (1 + (y / R))**2 - c * v * abs(v)
    elif equation == 4:
        return -g / (1 + (y / R))**2 - c * np.exp(-y/H) * v * abs(v)
    else:
        # This section should never be reached unless compiler error occurred (ignore this)
        print(f"Script generator error: Invalid equation number '{equation}' written to script.")
    return -g


# -------- TERMINAL VELOCITY --------- 
def find_terminal_velocity(v, a):
    """
    Searches the net acceleration, a(t), to find an index i where the acceleration is approximately 0, 
    within a defined threshold. 
    
    The threshold is intialized as 0.001, and the function loops, increasing the threshold until an index
    i is found where a(i) is within the max threshold 0.02 of 0. 

    If no index where a(t) is within |0.1| of 0, then None is returned, indicating that the net acceleration
    did not come within the threshold of 0 at any point in the free fall. 
    
    Note that the value where a is approximately 0 could either be the terminal velocity or the maximum 
    (most negative) value of the velocity function. Hence why "terminal" is in quotations below.  

    Args:
        v (array-like): the velocity function v(t)
        a (array-like): the acceleration function a(t) = v'(t)

    Returns:
        float, int: the terminal velocity and its corresponding index if found. None, None otherwise
    """
    deviation_from_zero = 1e-6
    terminal_velocity = None
    
    while terminal_velocity is None and deviation_from_zero <= 2e-1:
        terminal_velocity_i = np.where(np.abs(a) < deviation_from_zero)[0]
        
        if terminal_velocity_i.size > 0:
            return v[terminal_velocity_i[0]], terminal_velocity_i[0] 
        else:
            deviation_from_zero *= 10
            
    return None, None


# -------- CHECK IF TERMINAL VELOCITY WAS FOUND -------- 
def check_if_terminal(t, v, a, impact_index, returned_velocity, threshold=0.2):
    """
    Checks if the velocity returned from find_terminal_velocity is the terminal velocity or a critical point
    in the velocity solution.
    
    The velocity is accepted as the terminal velocity if it is within a threshold (default value = 0.2 m/s)
    distance from the final velocity. 
    
    Args:
        `t` (float or array-like): Solution time span
        `v (array-like): the velocity function v(t).
        `a` (array-like): the acceleration function a(t) = v'(t)
        `impact_index` (int): index i where y[i] = 0. 
        `returned_velocity` (float): the velocity returned from the find_terminal_velocity function
        `threshold` (float): accepted difference between returned and final velocity (optional)
        
    Returns:
       (bool): True if terminal velocity is within the threshold of the final velocity. False otherwise. 
    """
    difference = abs(returned_velocity - v[impact_index]) 
    if difference < threshold:
        return True
    return False


# ------- GET TERMINAL VELOCITY DATA: NON-ZERO DRAG ------ 
def get_terminal_data(t, v, a, impact_i, velocity_data, term_v, term_v_index, reached=False):
    """
    Adds either the terminal or final velocity to the data table (described in the get_velocity_data
    function below) depending on if the terminal velocity was reached by the time of impact. 
    
    Returns:
        (string, bool): The data table and reached = True if terminal velocity was reached. False otherwise. 
    """

    if check_if_terminal(t, v, a, impact_i, term_v):
        reached = True
        velocity_data.append(["Terminal Velocity", f"{term_v:.2f} m/s", f"{t[term_v_index]:.2f}"])
    else:
        reached = False
        velocity_data.append(["Final Velocity", f"{v[impact_i]:.2f} m/s", f"{t[impact_i]:.2f}"])
    
    return velocity_data, reached


# -------- VELOCITY DATA: MAX AND FINAL *OR* TERMINAL -------- 
def get_velocity_data(t, v, a, impact_index, c, reached=False, converged=False):
    """
    Creates data table of trajectory events: ground impact time, maximum impact, and either the 
    terminal or final velocity depending on if the terminal velocity was reached by the time of impact. 

    Final velocity is added to the data table if:
    
        - There is zero drag (terminal velocity N/A).        reached = False, converged = False
        - find_terminal_velocity failed to converge to zero. reached = False, converged = False
        - find_terminal_velocity returned a maximum value    reached = False, converged = True
        
    Terminal velocity is added to the data table if find_terminal_velocity returned the terminal
    velocity: reached = True, converged = True 
    
        Args:
        `t` (float or array-like): Solution time span
        `v (array-like): the velocity function v(t).
        `a` (array-like): the acceleration function a(t) = v'(t)
        `impact_index` (int): index of impact with ground
        `c` (float): drag force constant divided by mass
        `reached` (bool): True if terminal velocity was found, False otherwise (the max was found). 
        `converged` (bool): False if find_terminal_velocity returned value of None (no critical point found)
        
    Returns:
        (string, bool, bool): the velocity data, and reached, converged as desribed above.
    """
    velocity_data = [["Ground Impact", " N/A", f"{t[impact_index]:.2f}"]]
    velocity_data.append(["Maximum Velocity", f"{min(v):.2f} m/s", f"{t[np.argmin(v)]:.2f}"])
    returned_velocity = None

    # if no drag, terminal velocity not applicable 
    if c != 0:
        returned_velocity, returned_index = find_terminal_velocity(v,a)

    # find_terminal_velocity converged to zero, check if terminal velocity was found
    if returned_velocity is not None:
        converged = True
        velocity_data, reached = get_terminal_data(t, v, a, impact_index, velocity_data, returned_velocity, returned_index)
    
    # find_terminal_velocity did not converge to zero (or drag is zero). Add final velocity and return 
    else:
        velocity_data.append(["Final Velocity", f"{v[impact_index]:.2f} m/s", f"{t[impact_index]:.2f}"])
        
    return velocity_data, reached, converged

# -------------------------------------------------------------
#   PLOTTING SOLUTION 
# -------------------------------------------------------------
def generate_plots(t, y, v, a, y0):
    """
    Plots the trajectory y(t), velocity v(t), and total acceleration a(t) of a body in free-fall from 
    the start of the fall to the time of impact.

    Args:
        t (array-like): The time interval of the solution from 0 to the time of impact.
        y (list): The altitude as a function of time y(t).
        v (list): The velocity as a function of time y'(t).
        a (list): The net acceleration as a function of time y''(t).
        time_of_impact (float): The time t where y(t) = 0.
        y0 (float): The initial height of the body.

    Generates four images: a separate image for y(t), v(t), and a(t), along with a subplot combining 
    each case.
    
    Note: if using this script without PhysWiz, comment out the plt.savefig() lines and/or change the 
          path to your path name. 
    """
    # Altitude vs. Time plot
    plt.figure(figsize=(10, 5))
    plt.plot(t, y)
    plt.title("Altitude vs. Time")
    plt.xlabel('Time (s)')
    plt.ylabel('Height (m)')
    plt.xlim(0, t[-1])
    plt.ylim(0, y0 + (y0/10))
    plt.grid()
    plt.savefig('/home/vagrant/mirror/plots/free_fall_y_vs_t') # comment out or change to your own path name 

    # Velocity vs. Time plot
    plt.figure(figsize=(10, 5))
    plt.plot(t, v)
    plt.title("Velocity vs. Time")
    plt.xlabel('Time (s)')
    plt.ylabel('Velocity (m/s)')
    plt.xlim(0, t[-1])
    plt.grid()
    plt.savefig('/home/vagrant/mirror/plots/free_fall_v_vs_t') # comment out or change to your own path name

    # Acceleration vs. Time plot
    plt.figure(figsize=(10, 5))
    plt.plot(t, a)
    plt.title("Acceleration vs. Time")
    plt.xlabel('Time (s)')
    plt.ylabel('Acceleration (m/s²)')
    plt.xlim(0, t[-1])
    plt.grid()
    plt.savefig('/home/vagrant/mirror/plots/free_fall_a_vs_t') # comment out or change to your own path name

    # Subplot of each case 
    fig, axes = plt.subplots(3, 1, figsize=(9, 15))

    # Altitude vs. Time subplot
    axes[0].plot(t, y)
    axes[0].set_title("Altitude vs. Time")
    axes[0].set_xlabel('Time (s)')
    axes[0].set_ylabel('Height (m)')
    axes[0].set_xlim(0, t[-1])
    axes[0].set_ylim(0, y0 + (y0/10))
    axes[0].grid()

    # Velocity vs. Time subplot
    axes[1].plot(t, v) 
    axes[1].set_title("Velocity vs. Time")
    axes[1].set_xlabel('Time (s)')
    axes[1].set_ylabel('Velocity (m/s)')
    axes[1].set_xlim(0, t[-1])
    axes[1].grid()

    # Acceleration vs. Time subplot
    axes[2].plot(t, a)
    axes[2].set_title("Acceleration vs. Time")
    axes[2].set_xlabel('Time (s)')
    axes[2].set_ylabel('Acceleration (m/s²)')
    axes[2].set_xlim(0, t[-1])
    axes[2].grid()

    fig.tight_layout()
    plt.savefig('/home/vagrant/mirror/plots/free_fall_plots_combined') # comment out or change to your own path name


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
# ------------------------------------------------------------

def main():

    # -----------------------------------------------------
    #              INITIALIZE PARAMETERS
    # -----------------------------------------------------
    # Initialize all potential variables as None
    y0 = None           # initial height (m)
    m = None            # mass (kg)
    k = None            # drag constant (kg/m)
    Cd = None           # drag coefficient (dimensionless)
    A = None            # cross-sectional area (m^2)
    H = None            # scale height of the atmosphere (m)
    R = None            # planetary radius (m)
    M = None            # planetary mass (kg)
    rho = None          # air density (kg/m^3)
    g = None            # gravitational acceleration constant (m/s^2)
    v0 = None           # initial velocity (m/s)
    equation = 0        # equation ID (0 used as equation ID None value)
    # Your data:
    y0 = 4200
    m = 290
    A = 1.500000
    Cd = 1
    g = -3.710000
    rho = 0.020000
    H = 11100.000000
    equation = 2

    # Create a dictionary of the local variables and their values (float, int, or None)
    local_vars = {
        'y0': y0,
        'm': m,
        'k': k,
        'Cd': Cd,
        'A': A,
        'H': H,
        'R': R,
        'M': M,
        'rho': rho,
        'g': g,
        'v0': v0,
        'equation': equation
    }

    # Create dictionary of required variables and their values (None, int, or float) from Your data 
    data = {
        'y0': get_value('y0', None, local_vars), 
        'm': get_value('m', None, local_vars),   
        'k': get_value('k', None, local_vars),    # Note: must provide a value for k 
        'Cd': get_value('Cd', None, local_vars),  # OR both Cd and A
        'A': get_value('A', None, local_vars)
    }

    # Check if any required data is missing. If missing values are found, an error message is printed 
    # and program exits. (return 0 is used to exit as this is a PhysLang error, not a Python error)
    if missing_data_check(data, equation):
        return 0
        
    # Set any missing optional parameters to their default Earth values (see User Manual for details)
    H, R, M, rho, g = get_planet_parameters(local_vars) 
    v0 = get_value('v0', 0, local_vars)

    # calculate drag constant (using_k tells printer whether to print k or Cd and A) 
    c, using_k = calculate_drag_constant(k, Cd, A, m, rho, H, g, local_vars)


    # -----------------------------------------------------
    #               GENERATE SOLUTION                  
    # -----------------------------------------------------        
    # Define the initial time span for the integration 
    t_span = (0, 1000)
    t_eval = np.linspace(t_span[0], t_span[1], 10000)

    y_initial = [y0, v0] # initial conditions vector
    
    # Get the time of impact and solution for a body in free-fall with the defined parameters
    time_of_impact, impact_index, solution = get_free_fall_solution(equation_of_motion, y_initial, t_span, t_eval, g, c, H, R, equation)

    # Extract solution components up to the impact index
    t = solution.t[:impact_index+1]    # time span
    y = solution.y[0][:impact_index+1] # position y(t)
    v = solution.y[1][:impact_index+1] # velocity y'(t)

    # Net acceleration y''(t)
    accel = calculate_acceleration(g, c, v, y, H, R, equation)
    
    # Get the velocity data: reached = True if terminal velocity found, converged = False if no critial point
    # found in acceleration function.
    velocity_data, reached, converged = get_velocity_data(t, v, accel, impact_index, c)

    # --------------------------------------------------
    #                 PRINT RESULTS                    
    # --------------------------------------------------
    # update variables for printing values 
    local_vars = {
        'equation': equation,
        'y0': y0,
        'm': m,
        'k': k,
        'Cd': Cd,
        'A': A,
        'H': H,
        'R': R,
        'M': M,
        'rho': rho,
        'g': g,
        'v0': v0,
        'c': c,
        'using_k': using_k,
        'equation': equation
    }
    print_parameters(local_vars)
    
    # print velocities
    print_free_fall_results(velocity_data)
    
    # if terminal velocity was not found, print the reason why
    if not reached:
        print_velocity_message(c, converged)
    print("") 


    # --------------------------------------------------
    #                 PLOT RESULTS                    
    # --------------------------------------------------
    # Plot position, velocity, and acceleration solutions
    generate_plots(t, y, v, accel, y0)

if __name__ == "__main__":
    main()
