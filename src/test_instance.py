from rar_class import Rar
import numpy as np


#parameters = np.array([56.0, 36.07038498980966, 63.4204354193215, 1.252654348e-5])      # GD-1
parameters = np.array([56.0, 37.765595, 66.34067, 1.1977342e-05])                        # Becerra-Vergara
halo = Rar(parameters, temperature_func=True, chemical_func=True, core_func=True)

print(halo.T(halo.r[-2]))