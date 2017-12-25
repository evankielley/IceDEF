"""
This module does analyis on the output of an Iceberg simulation.
"""

import matplotlib.pyplot as plt                                         
import scipy.io as sio                                                  

def main():

    diff_mat = calc_diff_mat()
    plot_model_diff(diff_mat[0][:], diff_mat[1][:], show=True)
                                                                    
def calc_diff_mat():                  

    root_dir = '/home/evankielley/IceDEF/TestCase/'
    output_dir = root_dir + 'Outputs/'
    base_file = output_dir + 'NewWagnerTestCaseOutput.mat'
    base_vec_list = ['XIL', 'YIL', 'mL', 'mW', 'mH']
    test_file = output_dir + 'NewPyWagnerTestCaseOutput.mat'
    test_vec_list = ['x_arr', 'y_arr', 'l_arr', 'w_arr', 'h_arr']
    diff_mat = []

    for i in range(len(base_vec_list)):
        base_vec = sio.loadmat(base_file)[base_vec_list[i]]
        test_vec = sio.loadmat(test_file)[test_vec_list[i]]
        diff_vec = base_vec - test_vec
        diff_mat.append(diff_vec)

    return diff_mat 


def plot_model_diff(xDiff, yDiff, show=False):

    f = plt.figure()
    plt.subplot(121)
    plt.plot(xDiff.transpose())
    plt.ylabel('Difference in Location')
    plt.xlabel('Timestep')
    plt.title('X')
    plt.subplot(122)
    plt.plot(yDiff.transpose())
    #plt.ylabel('Difference in Y Location')
    plt.xlabel('Timestep')
    plt.title('Y')

    if show:
        plt.show()

    else:
        return f        


def plot_berg_location(xil, yil):
    """
    This function plots the xy position of each iceberg over time.

    Args:
        xil (list): X location of the iceberg.
        yil (list): Y location of the iceberg.

    Generates:
        A plot showing the xy position of icebergs over time.
    """
    plt.plot(xil,yil)                                                       
    plt.xlabel('X Location')                                                   
    plt.ylabel('Y Location')                                                   
    plt.title('Iceberg Drift')                                          
    plt.show() 

if __name__ == "__main__":
	main()
