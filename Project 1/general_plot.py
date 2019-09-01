import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0,1,10000)
u_exact = (np.exp(-10.)-1)*x+1-np.exp(-10*x)

fig = plt.figure(figsize=(10,5))
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
    ax = plt.subplot(2,3,m)
    plt.title('n = $10^%d$' % m)
    if m == 1:
        plt.plot(x,u_exact,'C0',label='Exact solution')
        plt.plot(xm,u_num,'C1',label = 'Numerical solution')
    else:
        plt.plot(x,u_exact,'C0')
        plt.plot(xm,u_num,'C1')
    plt.xlabel('x')
    plt.ylabel('u(x)')
    box = ax.get_position()
    ax.set_position([box.x0,box.y0+0.2*box.height,box.width, box.height*1.2])
    ax = plt.subplot(2,3,m+3)
    if m == 1:
        plt.plot(xm,err,'C2',label = 'Error')
    else:
        plt.plot(xm,err,'C2')
    plt.xlabel('x')
    plt.ylabel('Numerical-exact solution')
    box = ax.get_position()
    ax.set_position([box.x0,box.y0,box.width, box.height*1.2])
plt.tight_layout()
plt.figlegend(loc='upper right')
plt.savefig('general_plot.pdf')
