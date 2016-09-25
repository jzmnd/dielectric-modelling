#! /usr/bin/env python
"""
imp_model_lib.py
Library of functions for impedance modeling

Created by Jeremy Smith on 2015-09-24
University of California, Berkeley
j-smith@ecs.berkeley.edu
"""

import os
import sys
import myfunctions as mf
import numpy as np

__author__ = "Jeremy Smith"
__version__ = "1.2"


# Models for impedance
def zt2rx(zmag, theta):
	"""Converts Z-theta to Zreal-Zim"""
	zre = zmag*np.cos(theta)
	zim = zmag*np.sin(theta)
	return zre, zim


def rx2zt(zre, zim):
	"""Converts Zreal-Zim to Z-theta"""
	zmag = np.sqrt(zre**2 + zim**2)
	theta = np.arctan(zim/zre)
	return zmag, theta


def cpg2rx(cp, g, freq):
	"""Converts Cp-G to Zreal-Zim"""
	w = 2*np.pi*freq
	r = 1.0/g
	zre = r/(1 + (w*r*cp)**2)
	zim = -w*cp*(r**2)/(1 + (w*r*cp)**2)
	return zre, zim


def csg2rx(cs, g, freq):
	"""Converts Cs-G to Zreal-Zim"""
	w = 2*np.pi*freq
	r = 1.0/g
	zre = r*np.ones(len(freq))
	zim = -1.0/(w*cs)
	return zre, zim


# Materials models for impedance
def cpemodel2rx(r0, rinf, t0, alpha, freq):
	"""CPE model with series and parallel R"""
	w = 2*np.pi*freq
	b = 1 + (1j*w*t0)**(1 - alpha)
	z = rinf + (r0 - rinf)/b
	return z.real, z.imag


def cpemodelbeta2rx(r0, rinf, t0, alpha, beta, freq):
	w = 2*np.pi*freq
	b = (1 + (1j*w*t0)**(1 - alpha))**beta
	z = rinf + (r0 - rinf)/b
	return z.real, z.imag


def hnmodel2rx(r0, rinf, t0, alpha, beta, freq):
	w = 2*np.pi*freq
	b = (1 + (1j*w*t0)**(1 - alpha))**beta
	z = rinf + (r0 - rinf)/((r0/rinf)*(b - 1) + 1)
	return z.real, z.imag


def cpemodel2rx_withinterlayer(r0, rinf, cinter, t0, alpha, freq):
	w = 2*np.pi*freq
	b = 1 + (1j*w*t0)**(1 - alpha)
	z = rinf + (r0 - rinf)/b + 1.0/(1j*w*cinter)
	return z.real, z.imag

def oxidemodel_1(rs, q1, n1, r1, q2, n2, r2, freq):
	w = 2*np.pi*freq
	cperp1 = r1/(r1*(1j*w*q1)**n1 + 1)
	cperp2 = r2/(r2*(1j*w*q2)**n2 + 1)
	z = rs + cperp1 + cperp2
	return z.real, z.imag


# LCR models for extracting R and C
def seriesRC(zre, zim, freq):
	w = 2*np.pi*freq
	r = zre
	c = -1.0/(w*zim)
	return r, c


def parallelRC(zre, zim, freq):
	w = 2*np.pi*freq
	r = zre*(1 + (zim/zre)**2)
	c = -(zim/(zre*w*r))
	return r, c


default_starting_params = [1.0000000e+03, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00,
						   2.5000000e+04, 1.0000000e-09, 1.0000000e+00, 7.0000000e-01, 4.0000000e+00,
						   1.0000000e+09, 1.0000000e-09, 1.0000000e+00, 9.5000000e-01, 4.0000000e+00,
						   0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00,
						   0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00,
						   0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00,
						   0.0000000e+00, 1.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00,
						   0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00]

# File output functions
def levmwOutput(filename, path, zre, zim, freq, modelcircuit='A', parammask='1000011010010100000000000000000000000000', maxfev=1000, igacc='3', irch='2', paramlist=default_starting_params):
	if "results" not in os.listdir(path):
		os.mkdir(os.path.join(path, "results"))

	with open(os.path.join(path, "results", filename), 'w') as outfile:
		outfile.write("{:s}\n".format(filename))
		outfile.write("   9  ZZRRF {:s} .0000D+00C         0 .0000D+00     0     0\n".format(modelcircuit))
		outfile.write("   {:d}   40 {:d}    0    {:s}    0    0    1    {:s} .1000D+01   .5967D-01   .5973D-01\n".format(len(freq), maxfev, irch, igacc))

		for i in range(8):
			p = []
			for j in range(5):
				p.append("." + "{:.7E}".format(default_starting_params[j + i*5] * 10).replace("E", "D").replace(".", ""))
			outfile.write("  {:s}  {:s}  {:s}  {:s}  {:s}\n".format(p[0], p[1], p[2], p[3], p[4]))

		outfile.write("{:s}\n".format(parammask))

		for i, f in enumerate(freq):
			freqstr = "       ." + "{:.12E}".format(f * 10).replace("E", "D").replace(".", "")
			zrestr = "       ." + "{:.12E}".format(zre[i] * 10).replace("E", "D").replace(".", "")
			zimstr = "      -." + "{:.12E}".format(abs(zim[i]) * 10).replace("E", "D").replace(".", "")
			outfile.write("{:5d}{:s}{:s}{:s}\n".format(i + 1, freqstr, zrestr, zimstr))
	return
