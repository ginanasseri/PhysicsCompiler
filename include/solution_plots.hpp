#ifndef SOLUTION_PLOTS_HPP
#define SOLUTIONS_PLOTS_HPP

const char* generate_ff_plots = R"(
    # --------------------------------------------------
    #                 PLOT RESULTS                    
    # --------------------------------------------------
    # Plot position, velocity, and acceleration solutions
    generate_plots(t, y, v, accel, y0)
)";

const char* generate_proj_plots = R"(
    # --------------------------------------------------
    #                 PLOT RESULTS                    
    # --------------------------------------------------
    plot_trajectory(x,y)
)";

const char* projectile_plot = R"(
# -------------------------------------------------------------
#   PLOTTING SOLUTION 
# -------------------------------------------------------------
def plot_trajectory(x,y):
    """
    Plot the trajectory of the projectile.
    """
    plt.figure(figsize=(10, 5))
    plt.plot(x, y)
    plt.title('Projectile Motion Trajectory')
    plt.xlabel('Distance (m)')
    plt.ylabel('Altitude (m)')
    plt.grid(True)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.savefig('/home/vagrant/mirror/projectile_plot//projectile_motion_trajectory') # comment out or change to your own path name   
)";

const char* free_fall_plots = R"(# -------------------------------------------------------------
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
)";

#endif // SOLUTION_PLOTS__HPP
