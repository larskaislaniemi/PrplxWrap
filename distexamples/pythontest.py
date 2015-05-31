from ctypes import *
import sys
import time

p_pname_len = 15 # const, length of phase name, includes space for NUL
p_size_sysprops = 28 # const, num of system properties
p_size_phases = 25  # const, max num of phases
p_cname_len = 6   # const, includes space for NUL

components = []   # will hold the oxide names

lib_meemum = cdll.LoadLibrary('/home/lars/Progs/Code/perplex668/libmeemum.so') # this has to be compiled with gfortran, not gcc
libpx = cdll.LoadLibrary('/home/lars/Progs/Code/PrplxWrapWT/libperplex.so')

libpx.ini_phaseq.restype = c_int
libpx.ini_phaseq.argtypes = [c_char_p]
libpx.ini_phaseq(b'/home/lars/Data/melto4/melto4')

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
    

P = 10000.0
T = 1623.0
ncomp = n
dbgprint = 1
composition = [45.0, 43.0, 8.0, 2.0, 0.0, 2.0]

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
