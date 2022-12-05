#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import argparse as ap
from prepin import dumpin1, dumpin2

def main_run(espdat,Istage,IIstage,\
				ptype,dtype,\
				nmol,charge,\
				QWTw,QWTs,PWTw,PWTs,\
				EXC12,EXC13,\
				DEPTH,\
				verbose):
	if Istage and IIstage:
		dumpin1(espdat,Istage,ptype,dtype,nmol,charge,QWTw,PWTw,EXC12,EXC13,DEPTH,verbose)
		dumpin2(espdat,IIstage,ptype,dtype,nmol,charge,QWTs,PWTs,EXC12,EXC13,DEPTH,0)
	elif len(Istage) !=0 and len(IIstage) == 0:
		dumpin1(espdat,Istage,ptype,dtype,nmol,charge,QWTw,PWTw,EXC12,EXC13,DEPTH,verbose)
	elif len(Istage) ==0 and len(IIstage) != 0:
		dumpin2(espdat,IIstage,ptype,dtype,nmol,charge,QWTs,PWTs,EXC12,EXC13,DEPTH,verbose)
	else:
		print('Please specify one or two stage input file name with option --Istage or --IIstage')


parser = ap.ArgumentParser(description='PyRESP Generation')
parser.add_argument('--espdat','-i', help='Input file of esp data',required=True)
parser.add_argument('--Istage','-f1',help='Output file for 1st stage',default='pyrespgen.1st')
parser.add_argument('--IIstage','-f2',help='Output file for 2nd stage',default='pyrespgen.2nd')
parser.add_argument('--ptype','-p',help='Polarization type: x chg| x ind | x perm | x perm-v',default='chg')
parser.add_argument('--dtype','-d',help='Damping Function type: x additive | x applequist | x tinker | x exp | x linear', default='additive')
parser.add_argument('--nmol','-n',type=int,help='Number of conformations',default=1)
parser.add_argument('--charge','-q',type=int,help='Total charge for this structure or conformer',default=0)
parser.add_argument('--QWTw','-qwtw',type=float,help='Charge Constraint 1st stage',default=0.0005)
parser.add_argument('--QWTs','-qwts',type=float,help='Charge Constraint 2nd stage',default=0.001)
parser.add_argument('--PWTw','-pwtw',type=float,help='Permanent dipoles Constraint',default=0.0005)
parser.add_argument('--PWTs','-pwts',type=float,help='Permanent dipoles Constraint',default=0.001)
parser.add_argument('--EXC12','-exc12',type=int,help='include (0) or exclude (1) 1-2 interaction',default=0)
parser.add_argument('--EXC13','-exc13',type=int,help='include (0) or exclude (1) 1-3 interaction',default=0)
parser.add_argument('--DEPTH','-depth',type=int,help='Maximum depth for searching equilvalance atoms',default=3)
parser.add_argument('--verbose','-v',type=int,help='Print verbose information',default=0)


args = parser.parse_args()

if __name__ == '__main__':
	try:
		main_run(args.espdat,args.Istage,args.IIstage,\
						args.ptype,args.dtype,\
						args.nmol,args.charge,\
						args.QWTw,args.QWTs,args.PWTw,args.PWTs,\
						args.EXC12,args.EXC13,\
						args.DEPTH,\
						args.verbose)
	except Exception as e:
		print(e)
		
