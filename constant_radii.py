#!/bin/python
import numpy as np
import pandas as pd


tbl = [['H',1,0.40],\
		['D',1,0.40],\
		['C',6,0.85],\
		['N',7,0.78],\
		['O',8,0.85],\
		['F',9,0.53],\
		['P',15,1.08],\
		['S',16,1.04],\
		['Cl',17,0.92],\
		['Br',35,1.09],\
		['I',53,1.50],\
		['Si',14,1.11],\
		['Be',4,0.56],\
		['B',5,0.81],\
		['He',2,0.32],\
		['Ne',10,0.69],\
		['Ar',18,0.97],\
		['Kr',36,1.10],\
		['Te',52,1.35],\
		['Xe',54,1.30],\
		['Na',11,0.97],\
		['Mg',12,0.66],\
		['K',19,1.33],\
		['Ca',20,0.99],\
		['Mn',25,0.80],\
		['Fe',26,0.74],\
		['Cu',29,0.96],\
		['Zn',30,0.74],\
		['Al',13,0.64],\
		['Li',3,0.85],\
		['Sc',21,0.84],\
		['Ti',22,0.77],\
		['V',23,0.74],\
		['Cr',24,0.72],\
		['Co',27,0.79],\
		['Ni',28,0.79],\
		['Ga',31,0.72],\
		['Ge',32,0.82],\
		['As',33,1.13],\
		['Se',34,1.12],\
		['Sr',38,1.25],\
		['Ba',56, 1.41],\
		['Ru',44,0.78],\
		['Rh',45,0.76],\
		['Pd',46,0.95],\
		['Ag',47,1.02],\
		['Cd',48,1.03],\
		['Pt',78,0.89],\
		['Au',79,1.42],\
		['Hg',80,1.26],\
		['Tl',81,1.55],\
		['Pb',82,1.26],\
		['Rb',37,2.11],\
		['Y',39,1.62],\
		['Zr',40,1.48],\
		['Nb',41,1.37],\
		['Mo',42,1.45],\
		['Tc',43,1.56],\
		['In',49,1.44],\
		['Sn',50,1.41],\
		['Sb',51,1.38],\
		['Cs',55,2.25],\
		['La',57,1.69],\
		['Lu',71,1.60],\
		['Hf',72,1.50],\
		['Ta',73,1.38],\
		['W',74,1.46],\
		['Re',75,1.59],\
		['Os',76,1.28],\
		['Ir',77,1.37],\
		['Bi',83,1.46],\
		['Rn',86,1.45]]

if __name__ =='__main__':
	tbl = pd.DataFrame(data=tbl,columns=['elem_syl','elem_num','elem_radii'])
	print(tbl)