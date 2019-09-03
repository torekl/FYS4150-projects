import numpy as np
import matplotlib.pyplot as plt

n = []
rel_err = []
infile = open('special_error','r')
infile.readline()
for line in infile:
    data = line.split()
    n.append(float(eval(data[0])))
    rel_err.append(eval(data[1]))
infile.close()

n = np.array(n)
rel_err = np.array(rel_err)
h = 1/(n+1)

plt.figure()
plt.plot(np.log10(h),np.log10(rel_err))
plt.xlabel('log10(h)')
plt.ylabel('log10(Relative error)')
plt.savefig('error_plot.pdf')

outfile = open('rel_error_log','w')
outfile.write('      log(h)   log(rel_error)\n')
for i in range(len(h)):
    outfile.write('%12.8f%17.8f\n' % (np.log10(h[i]), np.log10(rel_err[i])))
outfile.close()
