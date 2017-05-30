import matplotlib.pyplot as plt

def make_plots(x):
    plt.plot(x.transpose())
    #plt.plot(x)
    plt.xlabel('Timestep')
    plt.ylabel('XIL')
    #plt.ylim([-10,10])
    plt.title('Iceberg Drift')
    plt.show()

def make_plots2(x,y):
    plt.plot(x,y)
    plt.xlabel('XIL')
    plt.ylabel('YIL')
    #plt.ylim([-10,10])
    plt.title('Iceberg Drift')
    plt.show()
