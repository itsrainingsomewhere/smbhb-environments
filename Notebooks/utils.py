import numpy as np
import astropy.units as u
from astropy.constants import G
from enterprise import constants as const

def logGWFrequencyDerivative(mass, gwFreq):
    """
    :param mass: chirp mass in solar masses, NOT logarithmic, not astropy Quantity
    :param gwFreq: gravitational wave frequency in Hz, NOT logarithmic, not astropy Quantity
    
    Returns the LOG of the evolution rate due to gravitational radiation in s^-2
    """
    orb = np.pi * gwFreq
    M = mass * const.Tsun # find mass expressed in seconds, rather than solar masses
    dOrb = (96/5) * (M**(5/3)) * (orb**(11/3))
    evolutionRate = dOrb / np.pi
    return np.log10(evolutionRate)

def gwFrequencyDerivative(mass, gwFreq):
    """
    :param mass: chirp mass in solar masses, not logarithmic, not astropy Quantity
    :param gwFreq: gravitational wave frequency in Hz, not logarithmic, not astropy Quantity
    
    Returns the rate of change of GW frequency due to gravitational radiation in s^-2
    """
    freq = gwFreq * np.pi #orbital frequency
    M = mass * 5 * (10 ** -6) # find mass expressed in seconds, rather than solar masses
    return (-96 / 5 / np.pi) * (M**(5/3)) * (freq**(11/3))

def gwFrequencyChangeApprox(mass, freq, yrs):
    """
    :param mass: chirp mass in solar masses, not logarithmic, not astropy Quantity
    :param freq: graviational wave frequency in Hz, not logarithmic, not astropy Quantity
    :param yrs: years over which we want the evolution
    
    Returns the total change in GW frequency due to gravitational radiation in s^-2 using an approximation that assumes constant derivative
    """
    seconds = yrs * const.yr # convert the years to seconds
    return -gwFrequencyDerivative(mass, freq) * seconds     # as if change rate per second * seconds = total change

def ssFrequencyDerivative(mass, sigma, rho, f, q=1):
    """
    :param mass: total mass; astropy Quantity object
    :param sigma: velocity dispersion; astropy Quantity object
    :param rho: mass density; astropy Quantity object
    :param f: GW frequency; astropy Quantity object
    
    Returns the rate of change of GW frequency due to stellar scattering in s^-2, as an astropy Quantity object
    """
    #TO CALCULATE H:
    #derive separation from frequency
    separation = (G * mass / (np.pi * f) ** 2) ** (1/3)
    separation = separation.to(u.m)
    
    #find ah
    M2 = q * mass / (1 + q)
    ah = G * M2 / 4 / (sigma ** 2)
    ah = ah.to(u.m)
    
    a = separation/ah   # define a to be orbital separation in units of ah
    A = 14.55
    a0 = 3.48
    gamma = -0.95

    H = A*( 1 + (a/a0) ) ** gamma

    #the actual calculation of df/dt
    dfdt = 3 * (2 * np.pi) ** (5/6) * G ** (4/3) * rho * H * mass ** (1/3) * f ** (1/3) / 2 / sigma
    return dfdt.to(u.s ** -2)

class Pulsar:
    def __init__(self, line):
        line = line.split()
        
        self.name = line[0]
        
        temp = line[1].split('+')
        self.meanDistance = float(temp[0])
        temporary = temp[1].split('-')
        self.plusError = float(temporary[0])
        self.minusError = float(temporary[1])
        
        self.rAsc = line[2]
        dividedUp = self.rAsc.split(':')
        hours = float(dividedUp[0])
        minutes = float(dividedUp[1])
        seconds = float(dividedUp[2])
        hours += (minutes / 60) + (seconds / 3600)
        self.phi = np.pi * hours / 12 #THIS IS IN RADIANS; DO WE WANT IT IN RADIANS?
        
        self.dec = line[3]
        divvyUp = self.dec.split(':')
        sign = divvyUp[0][0:1]
        divvyUp[0] = divvyUp[0][1:]
        degrees = float(divvyUp[0])
        minutes = float(divvyUp[1])
        seconds = float(divvyUp[2])
        degrees += (minutes / 60) + (seconds / 3600)
        if (sign.endswith('+')):
            self.theta = (np.pi/2) - (np.pi * degrees / 180)
        else:
            self.theta = (np.pi/2) + (np.pi * degrees / 180)
            
        self.pos = [np.sin(self.phi) * np.cos(self.theta), np.sin(self.phi) * np.sin(self.theta), np.cos(self.phi)]
    
    def __str__(self):
        return "{0}, located at {1}{2}: {3} +{4} -{5} kpc".format(self.name, self.rAsc, self.dec, self.meanDistance, self.plusError, self.minusError)
    
    