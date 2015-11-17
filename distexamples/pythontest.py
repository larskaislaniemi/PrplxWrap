from ctypes import *
import sys
import time

p_pname_len = 15 # const, length of phase name, includes space for NUL
p_size_sysprops = 28 # const, num of system properties
p_size_phases = 25  # const, max num of phases
p_cname_len = 6   # const, includes space for NUL

components = []   # will hold the oxide names

lib_meemum = cdll.LoadLibrary('/data/home/lkaislan/software/perplex672/libmeemum.so') # this has to be compiled with gfortran, not gcc
libpx = cdll.LoadLibrary('/data/home/lkaislan/software/PrplxWrap/libperplex.so')

libpx.ini_phaseq.restype = c_int
libpx.ini_phaseq.argtypes = [c_char_p]
libpx.ini_phaseq(b'crustmod')

# libpx.print_comp_order.restype = None
# libpx.print_comp_order.argtypes = []
# libpx.print_comp_order()

libpx.number_of_components.restype = c_int
libpx.number_of_components.argtypes = []
n = libpx.number_of_components()

compordertxt = " " * (n * (p_cname_len-1)) 

libpx.get_comp_order_list.restype = c_int
libpx.get_comp_order_list.argtypes = [c_char_p]
libpx.get_comp_order_list(compordertxt)

print(n)

# print (compordertxt)

for i in range(0,n):
    this_component = compordertxt[(i*(p_cname_len-1)):((i+1)*(p_cname_len-1))]
    components.append(this_component.strip())

print ("Components are: ")

for i in range(0,n):
    print("[" + components[i] + "]")
    
comps = {
    'NA2O ' : 0.0327,
    'MGO  ' : 0.0248,
    'AL2O3' : 0.1540,
    'SIO2 ' : 0.6662,
    'K2O  ' : 0.0280,
    'CAO  ' : 0.0359,
    'TIO2 ' : 0.0064,
    'MNO  ' : 0.0010,
    'FEO  ' : 0.0504,
    'H2O  ' : 0.0000
}

P = 6000.0
T = 900.0
ncomp = n
dbgprint = 1
composition = [0.00,    2.48,   15.40,  66.62,  2.80,   3.59,   0.64,   0.10,   5.04,   3.27]
# [H2O]
# [MGO]
# [AL2O3]
# [SIO2]
# [K2O]
# [CAO]
# [TIO2]
# [MNO]
# [FEO]
# [NA2O]

# NA2O  1  3.27000      0.00000      0.00000     weight amount
# MGO   1  2.48000      0.00000      0.00000     weight amount
# AL2O3 1  15.4000      0.00000      0.00000     weight amount
# SIO2  1  66.6200      0.00000      0.00000     weight amount
# K2O   1  2.80000      0.00000      0.00000     weight amount
# CAO   1  3.59000      0.00000      0.00000     weight amount
# TIO2  1 0.640000      0.00000      0.00000     weight amount
# MNO   1 0.100000      0.00000      0.00000     weight amount
# FEO   1  5.04000      0.00000      0.00000     weight amount
# H2O   1  0.00000      0.00000      0.00000     weight amount

comp_type = c_double * ncomp
wtphases_type = c_double * p_size_phases
cphases_type = c_double * (p_size_phases * ncomp)
sysprop_type = c_double * p_size_sysprops
namephases_type = c_char * (p_size_phases * p_pname_len)
nphases_type = c_int * 1

wtphases = wtphases_type()
cphases = cphases_type()
sysprop = sysprop_type()
namephases = namephases_type()
nphases = nphases_type()

libpx.phaseq.restype = c_int
libpx.phaseq.argtypes = [c_double, c_double, c_int, POINTER(c_double), POINTER(c_int), 
                         POINTER(c_double), POINTER(c_double), POINTER(c_double),
                         POINTER(c_char), c_int]

print "Running..."
starttime = time.time()
for loops in range(0, 1000):
    retval = libpx.phaseq(P, T, ncomp, comp_type(*composition), nphases, wtphases,
                          cphases, sysprop, namephases, dbgprint)
endtime = time.time()
print "Done."

#int phaseq(double P, double T, int ncomp, double *comp, int *nphases, 
#    double *wtphases, double *cphases, double *sysprop, char *namephases, int dbgprint) {

nphases = nphases[0]
print (nphases)

print ("** Weights:")
for i in range(0, nphases):
    sys.stdout.write(namephases[(i*p_pname_len):((i+1)*p_pname_len)])
    sys.stdout.write("{:4.2f}".format(wtphases[i]))
    sys.stdout.write("\n")
    
print ("** Compositions:")
sys.stdout.write("Phase         \t")
for j in range(0, ncomp):
    sys.stdout.write(components[j])
    sys.stdout.write("\t");
sys.stdout.write("\n")

for i in range(0, nphases):
    sys.stdout.write(namephases[(i*p_pname_len):((i+1)*p_pname_len)])
    sys.stdout.write("\t")
    for j in range(0, ncomp):
        sys.stdout.write("{:4.2f}".format((cphases[i * ncomp + j])))
        sys.stdout.write("\t")
    sys.stdout.write("\n")
    
print("Values of sysprops:")
for i in range(0, p_size_sysprops):
    sys.stdout.write("{:4.2f}".format(sysprop[i]))
    sys.stdout.write("\n")

print (endtime-starttime)
