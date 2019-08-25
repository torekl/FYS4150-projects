import numpy as np

b = np.array([2.]*10)
a = np.array([-1.]*10)
c = np.array([-1.]*10)
x = np.linspace(0,1,12)
g = (x[1]-x[0])**2*100*np.exp(-10*x[1:11])
n = len(g)
for i in range(1,n):
    a_b = a[i-1]/b[i-1]
    b[i] = b[i] - a_b*c[i-1]
    g[i] = g[i] - a_b*g[i-1]
g[-1] = g[-1]/b[-1]
for i in range(1,n):
    ind = n-1-i
    g[ind] = (g[ind]-c[ind]*g[ind+1])/b[ind]
print g
