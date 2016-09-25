#! /usr/bin/env python
"""
imp_converter_batch.py
Converts impedance measurement to Zre and Zim
Outputs standard text file and LEVM input file

Created by Jeremy Smith on 2016-07-19
University of California, Berkeley
j-smith@eecs.berkeley.edu
"""

import os
import sys
import xlrd
import myfunctions as mf
import imp_model_lib as imp
import numpy as np


__author__ = "Jeremy Smith"
__version__ = "1.2"

data_path = os.path.dirname(__file__)
files = os.listdir(data_path)   # All files in directory
data_extract = '1'              # Sheet number to extract data from
data_type = 'cpg'
skipini = 1                     # Data points to skip for LEMV output
skipend = 0


def main():
	print "\nBatch importing and processing C-V .xlsx files..."
	print "  PATH:  {:s}\n".format(data_path)
	for f in files:
		if ".xlsx" in f:
			workbook = xlrd.open_workbook(f, logfile=open(os.devnull, 'w'))
		else:
			continue
		if data_extract not in workbook.sheet_names():
			continue
		print "    FILE:  {:s}".format(f)
		datasheet = workbook.sheet_by_name(data_extract)
		data = {}
		if data_type == 'cpg':
			headers = ['freq', 'cp', 'g']
		elif data_type == 'ztd':
			headers = ['freq', 'zmag', 'theta']
		elif data_type == 'csg':
			headers = ['freq', 'cs', 'g']
		else:
			headers = ['freq', 'A', 'B']
		for h in headers:
			data[h] = []
		for row in range(datasheet.nrows):
			for col, h in enumerate(headers):
				data[h].append(float(datasheet.cell_value(row, col)))

		if data_type == 'ztd':
			freq = np.array(data['freq'])
			zmag = np.array(data['zmag'])
			tr = np.pi*np.array(data['theta'])/180
			zre, zim = imp.zt2rx(zmag, tr)

		if data_type == 'cpg':
			freq = np.array(data['freq'])
			cp = np.array(data['cp'])
			g = np.array(data['g'])
			zre, zim = imp.cpg2rx(cp, g, freq)

		if data_type == 'csg':
			freq = np.array(data['freq'])
			cp = np.array(data['cs'])
			g = np.array(data['g'])
			zre, zim = imp.csg2rx(cs, g, freq)


		mf.dataOutputHead("{:s}_{:s}.csv".format(f[:-5], data_type), data_path, [data[i] for i in headers], [[data_type.upper()], [','.join(headers)]],
			format_d="%.0f, %e, %e\n", format_h="%s\n %s")

		mf.dataOutput("{:s}_output.txt".format(f[:-5]), data_path, [freq,zre,zim],
			format="%.1e\t %e\t %e\n")
		if skipend == 0:
			imp.levmwOutput("{:s}_outputlevm".format(f[:-5]), data_path, zre[skipini:], zim[skipini:], freq[skipini:])
		else:
			imp.levmwOutput("{:s}_outputlevm".format(f[:-5]), data_path, zre[skipini:-skipend], zim[skipini:-skipend], freq[skipini:-skipend])

	return


if __name__ == "__main__":
    sys.exit(main())
