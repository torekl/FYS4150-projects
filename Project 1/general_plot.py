import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0,1,10000)
u_exact = (np.exp(-10.)-1)*x+1-np.exp(-10*x)

for m in range(1,4):
    xm = []
    u_num = []
    err = []
    infile = open('general_%d' % m)
    infile.readline()
    for line in infile:
        data = line.split()
        xm.append(eval(data[0]))
        u_num.append(eval(data[1]))
        err.append(eval(data[3]))
    infile.close()
    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(x,u_exact,'C0',label='Exact solution')
    plt.plot(xm,u_num,'C1',label = 'Numerical solution (n = $10^{%d}$)' % m)
    plt.xlabel('x')
    plt.ylabel('u(x)')
    plt.subplot(2,1,2)
    plt.plot(xm,err,'C2',label = 'Error')
    plt.xlabel('x')
    plt.ylabel('Numerical-exact solution')
    plt.tight_layout()
    plt.figlegend()
    plt.savefig('general_plot_%d.pdf' % m)
