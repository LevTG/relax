#!/usr/bin/python3 -u

import numpy as np
import matplotlib.pyplot as plt


def individual_exps(t, p):
    timevalues = p[1:-1:2]
    arr = np.exp(np.divide(-t, np.vstack(timevalues)))
    return arr.transpose()

def sum_of_exps(t, p):
    """A powerlaw: a * e**(t*b)"""
    amplitudes = p[:-1:2]
    ie = individual_exps(t, p)
    arr = np.zeros(t.shape[0])

    #print("arr={}".format(arr))
    #print("arr.shape = ".format(arr.shape))

    for i,a in enumerate(amplitudes):
        arr+= np.multiply(ie[:,i], a)
    arr+=p[-1]

    #arr = np.sum(np.multiply(np.vstack(amplitudes),
    #                         np.exp(np.divide(-t, np.vstack(timevalues))))) + \
    #      p[-1]
    return arr


params = np.array([0.1, 1, 0.1, 10, 0.1, 100, 0.1, 200, 0.1, 500, 0.3])
time = np.linspace(0, 1000, 2000)

amplitudes = params[:-1:2]
timevalues = params[1:-1:2]
print("amplitudes = {}".format(amplitudes))
print("timevalues = {}".format(timevalues))

individual_exp_data = individual_exps(time, params)
plt.title("Individual exponents of shape {}".format(individual_exp_data.shape))
plt.plot(individual_exp_data)
plt.show()

plt.title("The first exponent indexed as [:,0]")
plt.plot(individual_exp_data[:,0])
plt.show()

data_exps = sum_of_exps(time, params)
plt.title("data of shape {} as returned from sum_of_exps(time, params)".format(data_exps.shape))
plt.plot(data_exps)
plt.show()

print(data_exps)

#all_data = np.zeros((2000, 2000))
all_data_noise = np.zeros((2000, 2000))
for i in range(2000):
    noise = np.random.normal(size=2000,scale=0.01)
    #all_data[i] = data_exps
    all_data_noise[i] = data_exps + noise

n_slice=10
plt.title("The generated data with noise, the {}-th slice indexed as [{},:]".format(n_slice,n_slice))
plt.plot(all_data_noise[n_slice,:])
plt.show()

plt.title("The generated data with noise, the {}-th slice indexed as [:,{}]".format(n_slice,n_slice))
plt.plot(all_data_noise[:,n_slice])
plt.show()

np.save('data_noise.npy', all_data_noise)
print("array of size {} saved to data_noise.npy".format(all_data_noise.shape))

