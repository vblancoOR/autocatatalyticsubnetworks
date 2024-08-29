# The geometry and combinatorics of an autocatalytic ecology in chemical and cluster chemical reaction networks
## Python Codes for the Generation of Minimal Autocatalytic Subnetworks

Algorithm developed in:

Gagrani, P., Blanco, V., Smith, E., & Baum, D. (2024). Polyhedral geometry and combinatorics of an autocatalytic ecosystem. Journal of Mathematical Chemistry, 62(5), 1012-1078. https://link.springer.com/article/10.1007/s10910-024-01576-x


__________________________
File **autocatalytic_cores_lib** contains all the functions for the computation of minimal autocatalytic subnetworks from a given stoichiometric matrix.

Requirements:

- *pandas*: to export results to excel files and handle data frames
- *numpy*: to manipulate matrices
- *gurobipy*: for solving the optimization problems (it can be downloaded and installed from www.gurobi.com)

Functions:

- _OptModel_AutocatalyticCores(SM)_: Given the stoichiometric matrix SM computes the gurobi Integer Programming model to iterate in the construction of all the minimal autocatalytic subnetworks. The matrix SM is given such that its rows represent species and columns represent reactions.
- _ComputeAutocatalyticCores(SM, Excelfile: str, txtfile="", namesSp=[], namesRe=[], NumReact=0)_: Compute the entire list of minimal-support autocatalytic subnetworks. The stoichiometric matrix SM is given. The MAS are stored in an excel file "Excelfile" and optionally in a txt file. The names of the species and reactions can be optionally provided (otherwise the species are named C1, ...Cn, and the reactions R1, ... Rm). The function outputs a pandas dataframe with all the information about the cycles, namely:
  - AC (the label of the autocatalytic core, numbered)
numReact (number of reactions)
  - F1, ... Fn (taking value 1 if species n is food species for the core, and 0 otherwise).
  - M1, ... Mn (taking value 1 if species n is core species for the core, and 0 otherwise).
  - W1, ... Wn (taking value 1 if species n is waste species for the core, and 0 otherwise).
  - EM1, ... EMn (taking value 1 if species n is a member non core species for the core, and 0 otherwise).
  - X1, ..., Xm: a feasible flux for the core (SM*X_i > 0 for core species)
  - Zeros: Number of zeros in the core-restricted stoichiometric matrix (helpful to determine the type of the core).

Datasets:

In folder "data" we provide different files with stochiometric matrices for CCRN with 1-constituent (for L=4, ..., 7) called "ccrn4.txt",  "ccrn5.txt",  "ccrn6.txt", and  "ccrn7.txt".

**Usage** (example available in the Jupyter notebook _autocatalytic_subnetworks_example_ 

### Load the library
from autocatalytic_cores_lib import

### Load the dataset of 1-constituent CCRN with 4 species (files for more species are available in data folder)
SM=np.loadtxt("data/binarycrn_n4.txt", dtype='i', delimiter=' ')

### Compute minimal autocatalytic Subnetworks for that network
df = ComputeAutocatalyticCores(SM, "AC_n4.xlsx")



