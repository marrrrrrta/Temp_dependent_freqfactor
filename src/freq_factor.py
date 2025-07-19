import numpy as np

def entropy(temperature, a, b, c):
    """
    Calculates the entropy of a system given its temperature and parameters a, b, c.
    
    Parameters:
    temperature (float or array-like): Temperature in Kelvin
    a (float or array-like): Parameter a
    b (float or array-like): Parameter b
    c (float or array-like): Parameter c
    
    Returns:
    float or array-like: Entropy in J/K
    """
    kB = 8.6173335e-5           # Boltzmann constant (eV/K)
    
    first = c * np.log(1 + c * np.exp(-a/(kB * temperature)))
    second = 1/(kB * temperature) * a * c**2 * (np.exp(-a/(kB * temperature)) / (1 + c * np.exp(-a/(kB * temperature))))
    third = -1/(kB * temperature) * b * temperature**(3/2)
    return kB * (first + second + third)

def freq_factor(temperature, value):
    """
    Computes the frequency factor that follows the expression:
                 S = nu * K * exp(DeltaS / kB)
    
    * nu: lattice phonon vibration frequency (s-1)  
    * K: transition probability constant
    * kB: Boltzmann constant (eV/K)

    !!!!!!!!
    
    Args:
        temperature: temperature range of the process
    """
    kB = 8.6173335e-5           # Boltzmann constant (eV/K)
    K = 1                       # Transition probability constant (dimensionless)
    nu_m = 90000                # from B.H. Bransden book (m-1)
    nu_s = 2.53e+16             # (s-1)
    
    # Add one more temperature value at the end
    temperature = np.append(temperature, temperature[-1]+1)
    
    # Calculate the entropy change
    
    S_values = entropy(temperature, value.a, value.b, value.c)
    deltaS = np.diff(S_values)
    deltaS = np.concatenate(([0], deltaS))
    
    
    freq_factor = nu_s * K * np.exp(deltaS / kB)
    return freq_factor