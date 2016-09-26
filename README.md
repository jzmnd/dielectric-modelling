# dielectric-modelling

Script for processing impedance spectroscopy data and outputing LEVM file format.

Library of impedance spectroscopy functions:
zt2rx(zmag, theta)
rx2zt(zre, zim)
cpg2rx(cp, g, freq)
csg2rx(cs, g, freq)
cpemodel2rx(r0, rinf, t0, alpha, freq)
cpemodelbeta2rx(r0, rinf, t0, alpha, beta, freq)
hnmodel2rx(r0, rinf, t0, alpha, beta, freq)
cpemodel2rx_withinterlayer(r0, rinf, cinter, t0, alpha, freq)
oxidemodel_1(rs, q1, n1, r1, q2, n2, r2, freq)
seriesRC(zre, zim, freq)
parallelRC(zre, zim, freq)
levmwOutput(filename, path, zre, zim, freq, modelcircuit='A', parammask='1000011010010100000000000000000000000000', maxfev=1000, igacc='3', irch='2', paramlist=default_starting_params)
