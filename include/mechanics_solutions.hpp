#ifndef MECHANICS_SOLUTIONS_HPP
#define MECHANICS_SOLUTIONS_HPP

/* -------------
 * In this file:
 *
 *    *_ODE_solver: Contain the differential equation of motion, integration functions, as well as additional 
 *                  functions for computing events relative to the problem being solved. E.g.,: 
 * 
 *                - free_fall_ODE_solver examines the velocity data: was the terminal velocity reached? (If so,
 *                  what time? If not, why not?), and builds a velocity data table with the maximum and terminal 
 *                  (if reached) or final (if not reached) velocities and the time they occurred, as well as the
 *  
 *                - projectile_ODE_solver finds the max height, time of max height, total distance and total time
 *                  of trajectory. 
 * 
 *                - pendulum_ODE_solver would examine the period and amplitude of the pendulum.
 *
 *                  
 *  get_*_soltuion: call the ODE solver functions once all the parameters have been intialized with the User's 
 *                  data or default parameters.
 *  
 */


const char* free_fall_ODE_solver = R"(# ------------------------------------------------------------
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
)";

const char* get_free_fall_solution = R"(
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
)";



const char* projectile_ODE_solver = R"(# ------------------------------------------------------------
#   PARAMETER INITIALIZATION SECTION ENDS HERE
# ------------------------------------------------------------
#   SOLUTION SECTION: SOLVING THE DIFFERENTIAL EQUATION
# ------------------------------------------------------------
# --------  EQUATION OF MOTION -------- 
def eq_of_motion(t, z, g, c, H, equation):
    """
    Returns the equation equation of motion for the x and y components of a projectile with 
    air resitance (constant on non-contsant) as a system of first order ODE's as a list. 
    Defined to be used over an integration time span t with initial conditions z to compute 
    the x(t), vx(t), y(t), and vy(t) solutions to the system. 

        Args:
        `t` (float or array-like): Single time value or time span to evaluate over. 
        `z` (list): Initial conditions vector containing [x0, vx0, y0. vy0]
        `g` (float): Gravitational constant. 
        `c` (float): Drag force constant.

    Returns:
        list: the system of equations as a list [x'(t), x''(t), y'(t), y''(t)]
    """
    x, vx, y, vy = z
    v = np.hypot(vx, vy)
    
    # non-constant drag
    if equation == 6:
        dvx_dt = -c * np.exp(-y/H) * v * vx
        dvy_dt = -g - c * np.exp(-y/H) * v * vy
    
    # constant drag 
    else:
        dvx_dt = -c * v * vx
        dvy_dt = -c * v * vy - g
    return vx, dvx_dt, vy, dvy_dt


# ------- INTEGRATION EVENTS -------- 
def max_distance(t, z, g, c, H, equation):
    """
    Returns y(t)

    Used as event in solve_ivp integration to find t when y(t) = 0 (time of max distance).
    """
    return z[2]

max_distance.terminal = True # stop integration once max_distance hits
max_distance.direction = -1  # ensures direction is decreasing to 0 (impact with ground)


def max_height(t, z, g, c, H, equation):
    """
    Returns vy(t).

    Used as event in solve_ivp to find t when vy(t) = 0 (time of max height)
    """
    return z[3]


# ------- GET SOLUTION --------
def get_projectile_solution(eq_of_motion, z0, g, c, H, equation, t_span=[0,1000]):
    """
    Generates the solutions x(t), vx(t), y(t) and vy(t) to the system of ODE's defined in the 
    proejectile motion eq_of_motion() function by integrating the system from 0 until the 
    max_distance event hits. The constants of integrations are computed from the initial 
    donditions z0 = [x0, vx0, y0, vy0]. 
    
    The while loop is used to ensure that the event y(t) = 0 is captured within the time span.  
    The maximum integration time is initialized as 1000 s. If the max_distance (the impact event)
    does not hit, then the time span is doubled until the event is captured. 
    """
    while True:
        solution = solve_ivp(eq_of_motion, t_span, z0, dense_output=True, events=(max_distance, max_height), args=(g, c, H, equation))
        
        if solution.t_events[0].size > 0:
            time_of_impact = solution.t_events[0][0]  # max_distance hit
            break
        else:
            t_span = (t_span[0], t_span[1] * 2)
    
    if solution.t_events[1].size > 0:
        time_of_max_height = solution.t_events[1][0]  # max_height hit 
    else:
        time_of_max_height = None
    
    return time_of_impact, time_of_max_height, solution


# --- SOLUTION EVALUATION GRANULARITY ---------
def get_time_evaluation_points(time_of_impact, t_eval_max=10000):
    """
    Determines the number of evaluation points to use based on the total trajectory time 
    and returns an array of time evaluation points. 
    
    If the trajectory is less than 100 seconds, then t_eval is set to 1000 points. Otherwise, 
    t_eval is set to 10 points per second, up to a maximum of 10000 points. 
    
    Args:
        time_of_impact (float): The total trajectory time.
        t_eval_max (int): Maximum number of evaluation points

    Returns:
        np.ndarray: Array of time evaluation points.
    """
    if time_of_impact < 100:  
        t_eval_points = 1000
    else:
        t_eval_points = int(min(10 * time_of_impact, t_eval_max))
    
    return np.linspace(0, time_of_impact, t_eval_points))";


const char* get_projectile_solution = R"(
    # Create a dictionary of the local variables and their values (float, int, or None)
    local_vars = {
        'm': m,     
        'Cd': Cd,   
        'A': A,     
        'rho': rho, 
        'k': None,  # ensures correct constant used in drag calculation
        'H': H, 
        'v0': v0,
        'theta': theta,
        'c': c,
        'y0': y0,
        'x0': x0,
        'M': M,
        'R': R,
        'g': g,
        'equation': equation
    }    
    # Create dictionary of required variables and their values (None, int, or float) from Your data 
    data = {
        'v0': get_value('v0', None, local_vars), 
        'theta': get_value('theta', None, local_vars),   
        'c': get_value('c', None, local_vars),
        'm': get_value('m', None, local_vars),     
        'Cd': get_value('Cd', None, local_vars),   
        'A': get_value('A', None, local_vars)      
    }
    
    # Check if any required data is missing. If missing values are found, an error message is printed 
    # and program exits. (return 0 is used to exit as this is a PhysLang error, not a Python error)
    if missing_data_check(data,equation):
        return 0
    

    # initial velocities:
    v0_x, v0_y = get_init_velocities(v0, theta)

    # set planet parameters 
    H, R, M, rho, g = get_planet_parameters(local_vars)
    
    # Set any missing optional parameters to their default values
    x0 = get_value('x0', 0, local_vars)
    y0 = get_value('y0', 0, local_vars)

    # calculate drag constant from given parameter values
    if get_value('c', None, local_vars) is not None:
        using_c = True # for printing parameters 
    else:
        c, using_c = calculate_drag_constant(k, Cd, A, m, rho, H, g, local_vars)
    

    # -----------------------------------------------------
    #               GENERATE SOLUTION                  
    # -----------------------------------------------------     
    z_initial = [x0, v0_x, y0, v0_y] # initial conditions vector
    
    t_span = [0,500] # initial time span estimate
    z0 = z_initial
    
    # Get the time of impact, time of max height, and solution for projectile with air resistance 
    time_of_impact, max_h_t, solution = get_projectile_solution(eq_of_motion, z0, g, c, H, equation)
    
    # Get number of evaluation points based on the trajectory time 
    t_eval = get_time_evaluation_points(time_of_impact)

    # Extract solution components: 
    z = solution.sol(t_eval)
    x = z[0]
    y = z[2]
    
    # Update local variables
    local_vars = {
        'm': m,     
        'Cd': Cd,   
        'A': A,     
        'rho': rho, 
        'H': H,     
        'v0': v0,
        'theta': theta,
        'c': c,
        'y0': y0,
        'x0': x0,
        'M': M,
        'R': R,
        'g': g
    }
    
    results = {
        'time_of_impact': time_of_impact,
        'max_h_t': max_h_t,
        'max_distance': x[-1],
        'max_height': max(y)
    }
    
    print_parameters(local_vars, using_c)
    print_projectile_results(results)
)";

#endif // MECHANICS_SOLUTIONS_HPP
