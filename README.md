hybrid-coal

hybrid-coal is used to compute gene tree probabilities given species network under coalescent process.


INSTALLATION

To install hybrid-coal, type the following commands:
./bootstrap
make;

Version             | Branch  | Travis CI Build Status                                                                                                                   | Circle CI Build Status
------------------- | ------- | ---------------------------------------------------------------------------------------------------------------------------------------- |--------------------------
Development Version | dev     | [![Build Status](https://travis-ci.com/hybridLambda/hybrid-coal.svg?branch=dev)](https://travis-ci.com/hybridLambda/hybrid-coal)         | [![Circle CI](https://circleci.com/gh/hybridLambda/hybrid-coal.svg?style=svg)](https://circleci.com/gh/hybridLambda/hybrid-coal)


##DOCUMENTATION
[Download](https://github.com/hybridLambda/hybrid-coal/raw/doc/doc/manual.pdf)

##INSTALLATION

### For developers
To install hybrid-coal, first install the following packages and libraries

on Debian/Ubuntu based systems:
```bash
apt-get install git-core build-essential autoconf autoconf-archive libcppunit-dev graphviz
```
on Mac OS:
```bash
port install git cppunit automake autoconf autoconf-archive graphviz
```

then type the following commands:
```bash
./bootstrap
make
```

##ASSUMPTION
Input network files are written in extended newick format.

##LICENCE
You can freely use all code in this project under the conditions of the GNU
GPL Version 3 or later.

##HOW IT WORKS

Program parameters and options:

Options              | Useage |
:-------------------:| ------------------------------- |
-h or -help          | Help. List the following content. |
             -gt STR | Input the gene tree string string through command line or a file.
             -sp STR | Input the species network/tree string through command line or a file.
              -gtopo | To generate the gene tree topologies of a given set of taxa.
 -plot/-dot [option] | Use LaTEX(-plot) or Dot (-dot) to draw the input (defined by -spcu) network(tree).
           [-branch] | Branch lengths will be labelled in the figure.
             -o STR  | Specify the file name prefix of the output.

##Examples:
```bash
hybrid-coal -sp '((((B:1,C:1)s1:1)h1#.5:1,A:3)s2:1,(h1#.5:1,D:3)s3:1)r;'
hybrid-coal -sp trees/4_tax_sp_nt1_para -gt '(((A,D),C),B);'
hybrid-coal -sp trees/4_tax_sp_nt1_para -gt trees/4_tax_gt4.tre -latex
hybrid-coal -sp trees/4_tax_sp_nt1_para -plot
hybrid-coal -sp trees/4_tax_sp_nt1_para -plot -branch
hybrid-coal -sp trees/4_tax_sp_nt1_para -plot -label
hybrid-coal -sp trees/4_tax_sp_nt1_para -dot
hybrid-coal -sp trees/4_tax_sp_nt1_para -gtopo
```
