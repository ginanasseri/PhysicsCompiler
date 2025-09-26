import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd 
from scipy import stats

def get_cov(x,y):
    """
    Returns the covariance between the data sets x and y
    """
    cov_xy = np.cov(x, y)[0, 1]
    print(f'Covariance coefficient : {cov_xy:.3f}')
    return cov_xy
    
def get_corr(x,y):
    """
    Returns the correlation coefficient between the data sets x and y
    """
    corr_xy = np.corrcoef(x, y)[0, 1]
    print(f'Correlation coefficient: {corr_xy:.3f}\n')
    return corr_xy 


def print_descriptive_stats(x,y):
    # Calculate statistics for X
    mean_x = np.mean(x)
    var_x = np.var(x)
    stdev_x = np.std(x)

    # Calculate statistics for Y
    mean_y = np.mean(y)
    var_y = np.var(y)
    stdev_y = np.std(y)

    # Print results into table
    print("\nDescriptive Statistics Results:")
    title = "DATASET           X         Y   "
    print("-" * len(title))
    print(title)
    print("-" * len(title))
    table_format = "{:<16} {:<8.3f} {:<8.3f}"
    print(table_format.format("Mean", mean_x, mean_y))
    print(table_format.format("Variance", var_x, var_y))
    print(table_format.format("Std. Deviation", stdev_x, stdev_y))
    print("-" * len(title))
    
    # Print the covariance and correlation coeffificent of the data sets 
    get_cov(x,y)
    get_corr(x,y)

def normal_hist(data, label, filepath):
    """
    Creates a histogram and overlaying PDF of the data set `data`. The number of bins is 
    calculated using the square root choice. 
    
    The resolution for the PDF is set at bin_width/10. Future implementations could have 
    customized resolution. 
    """
    # create histogram
    min_x = min(data)
    max_x = max(data)
    mean_x = np.mean(data)
    stdev_x = np.std(data)
    N = int(np.ceil(np.sqrt(len(data)))) # np.ceil rounds up
    bin_width = (max_x - min_x)/N
    bin_edges = list(np.arange(min_x, max_x, bin_width))  
    plt.hist(data, bins=bin_edges, edgecolor='black', alpha=0.6, color='blue')
    
    # Scaled PDF
    n = len(data)
    fineness = 10 
    x_range = np.arange(min_x, max_x, bin_width/fineness) # finer range of x values
    scaled_pdfX = stats.norm.pdf(x_range, loc=mean_x, scale=stdev_x)*n*bin_width
    plt.plot(x_range, scaled_pdfX, 'r-', label='Scaled PDF')
    title = "Histogram and PDF for Dataset " + label
    plt.title(title)
    plt.xlabel(label)  
    plt.ylabel("Frequency or Probability Density")
    plt.legend()
    plt.savefig(filepath)
    plt.clf()
    #plt.show()      

 
def scatter_plot(x,y):    
    plt.scatter(x, y, s=1)
    plt.xlabel('X Values')
    plt.ylabel('Y Values')
    plt.title('Scatter Plot of Datasets X,Y')
    plt.savefig('/home/vagrant/mirror/plots/scatter')
    #plt.show()


def main():
    # Your data file:
    data_file = 'sample_data.csv'

    # Read the CSV file into the dataframe df
    df = pd.read_csv(data_file)

    # Check if file exists 
    if os.path.exists(data_file):
        # Read the file into the dataframe df
        df = pd.read_csv(data_file)
    else:
        print(f"\nError: Could not find data file: {file_path}. File path does not exist.\n")
        print("Suggested: Fix file path name and run program again (see manual for additional help).")

    # make sure file is properly formatted (column titles are x and y)
    if 'x' not in df.columns or 'y' not in df.columns:
        print("\nError reading file:",file)
        print("\nDetails: column headers in data file must be 'x' and 'y'")
        print("\nSuggested: update column titles and run program again (see manual for additonal help).")
        return 0

    # Extract the x and y data
    x = np.array(df['x'])
    y = np.array(df['y'])

    print_descriptive_stats(x,y)

    # plot results
    scatter_plot(x,y)
    normal_hist(x, "X", '/home/vagrant/mirror/plots/x_hist')
    normal_hist(y, "Y", '/home/vagrant/mirror/plots/y_hist')

if __name__ == "__main__":
    main()
