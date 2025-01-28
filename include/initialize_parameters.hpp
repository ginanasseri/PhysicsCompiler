#ifndef INITIALIZE_PARAMETERS_HPP
#define INITIALIZE_PARAMETERS_HPP

/* -------------
 * In this file: 
 * 
 *  set_gravitational_constant : return g based on source code optional args
 *  calculate_drag_constant    : return c or 1/2 * Cd * A * m depending what was provided in source code
 *  set_planet_params          : assign default values for omitted optional args
 *  initialize_mechanics_params: initialize all potential variables as None before data write
 *  get_init_velocities        : return x,y initial velocity components from v0, theta
 *  init_proj_params           : initialize extra parameters (theta,c) for projectile motion as None
 * 
 */

const char* set_gravitational_constant = R"(# --------  CALCULATE GRAVITY WITH GIVEN PARAMETERS -------- 
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
)";


const char* calculate_drag_constant = R"(
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
)";


const char* set_planet_params = R"(
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

)";


const char* initialize_mechanics_params = R"(
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
    equation = 0        # equation ID (0 used as equation ID None value))";



const char* get_init_velocites = R"(# ---- INITIAL VELOCITIES ----- 
def get_init_velocities(v0, theta):
    """
    Returns initial velocity x,y components from initial velocity v0 and initial angle theta.
    """
    angle = np.radians(theta)  # Convert angle to radians
    v0_x = v0 * np.cos(angle)
    v0_y = v0 * np.sin(angle)
    return v0_x, v0_y

)";


const char* init_proj_params = R"(    x0 = None           # x distance (m)
    theta = None        # initial angle (degrees)
    c = None        # drag force constant divided by mass
)";


#endif // INITIALIZE_PARAMETERS_HPP
