#  PyRESP GEN PROGRAM version 1.0
---
This program is developed for automatically generating the input file for program PyRESP (a python version) or pol\_resp (fortran version)

Author info: Qiang Zhu  PostDoc _AT_ UC Irvine  [Qiang's Blog](https://csu1505110121.github.io/about/)

# Architecture
- constant.py        : definition of constant variables
- pyresp\_gen.py     : main program for pyresp gen
- prepin.py          : function of dumping 1st stage and 2nd stage input file
- zmatrix.py         : function for getting information from esp data, building distance matrix, finding equivalent atoms and bonds, and dump esp data to xyz for visualization.a
- pyresp_gen_env.yml : a environment library file
- examples           : a bundle of examples for test

# How to intall
## Dependencies
This program runs with python 3. To run this program, make sure the following dependencies are installed:
- [numpy](https://numpy.org/): A Python library supporting for large, multi-dimensional arrays and matrices.
- [scipy](https://scipy.org/): A Python library used for scientific computing and technical computing.

You can manually install the aformationed library one by one. However, for your convenience, a environment file named `pyresp_gen_env.yml` was provided. A corresponding environment could be created using following command:
```
conda env create -f pyresp_gen_env.yml
```

# How to use
## Step 1:
Downloading the source code and put them to your desired location `\your\path`;

## Step 2:
Installing the required libraries and activate them. If you are utilizing the [anaconda](https://www.anaconda.com), you can easily install the libraries with the environment file `pyresp_gen_env.yml` we supplied. Then, just activate the environment you created with the following command:

```
conda activate your_env_name
```

## Step 3:
Setting the environment path:
- For Mac User:
Editing the `~/.zshrc` and adding the following line:
```
export PATH=/your/path/to/pyresp_gen:$PATH
```

- For Linux User:
Editing the `~/.bashrc` file and adding the same line

## Executing the pyresp_gen command line
```
usage:
pyresp_gen.py [-h] --espdat ESPDAT [--Istage ISTAGE] [--IIstage IISTAGE] [--ptype PTYPE] [--dtype DTYPE] [--nmol NMOL] [--charge CHARGE] \
			  [--QWTw QWTW] [--QWTs QWTS] [--PWT PWT] [--EXC12 EXC12] [--EXC13 EXC13] [--DEPTH DEPTH] [--verbose VERBOSE]
```

For `help` just type `pyresp_gen.py -h` which will show the explanation of each arguments. Here, we explained them in detail.

- `-h`  or `--help`          : show help message and exit
- `-i`  or `--espdat`		 : input file of esp data
- `-f1` or `--Istage`        : name for the 1st stage input file, if you not specified, default name `pyrespgen.1st` will be applied
- `-f2` or `--IIstage`       : name for the 2nd stage input file, if you not specified, default name `pyrespgen.2nd` will be applied
- `-p`  or `--ptype`         : polarization type, in this version, only support 
								`chg`:point charge model                     <- default one
								`ind`: induced dipole model 
								`perm`: permanent dipole model
								`perm-v`: permanent dipole model virtual
- `-d`  or `--dtype`         : Damping function type, in this version supproted types are listed below:
								`additive`                                    <- default one
								`applequist`
								`tinker`
								`exp`
								`linear`

- `-nmol` or `--nmol`        : Number of conformations. In this version could only handle one conformation
- `-q`    or `--charge`      : Total charge for this structure or conformer
- `-qwtw` or `--QWTw`        : Weak charge constraint in the 1st stage: default is set to be 0.0005
- `-qwts` or `--QWTs`        : Strong charge constraint in the 2nd stage: defalt is set to be 0.001
- `-pwtw`  or `--PWTw`       : Weak permanent dipoles constraint
- `-pwts`  or `--PWTs`       : Strong permanent dipoles constraint
- `-exc12` or `--EXC12`      : including (0) or excluding (1) 1-2 interactions
- `-exc13` or `--EXC13`      : including (0) or excluding (1) 1-3 interactions
- `-depth` or `--DEPTH`      : Maximum depth for searching equivalance atoms
- `-v`     or `--verbose`    : Print verbose info.


