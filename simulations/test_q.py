import matplotlib.pyplot as plt
import numpy as np
import math
sig = 1
Deltas = np.linspace(0.1, 10, 100, endpoint = False)
bound = np.zeros(len(Deltas))
val = np.zeros(len(Deltas))
for i_delta in range (1, len(Deltas)) :
	delta = Deltas[i_delta]
	bound[i_delta] = math.pow(delta, 2) / 4 + math.pow(sig, 2) + (sig * delta)/math.sqrt(2 * math.pi) + 4 * math.sqrt(2/math.pi) * math.pow(sig, 3) / (3 * delta)
	summa = 0

	for z in range (1, pow(10,7)) :
		summa = summa + math.pow(delta,2) * math.erfc( (delta / (2 * math.sqrt(2) * sig) ) * z)  * math.pow((z + 1), 2) / 4	
	
	val[i_delta] = summa

one = plt.plot(Deltas,val,'rs',Deltas, bound,'bs')
plt.show()
plt.legend([one], ['simulated', 'bound'])

