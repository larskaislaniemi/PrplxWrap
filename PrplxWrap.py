from ctypes import *
import sys

class perplex:
    """..."""

    p_pname_len = 15 # const, length of phase name, includes space for NUL
    p_size_sysprops = 28 # const, num of system properties
    p_size_phases = 25  # const, max num of phases
    p_cname_len = 6   # const, includes space for NUL

    components = []   # will hold the oxide names
    c_clr2prp = []    # mapping from caller component ordering to perplex component ordering
    c_prp2clr = []    # ... and vice versa
    ncomp = 0

    #lib_meemum = cdll.LoadLibrary('/Users/larskaislaniemi/Code/perplex668/libmeemum.so') # this has to be compiled with gfortran, not gcc
    #libpx = cdll.LoadLibrary('/Users/larskaislaniemi/Code/PrplxWrap/libperplex.so')
    libpx = cdll.LoadLibrary('/data/home/lkaislan/software/PrplxWrap/libperplex.so')

    syspropnum = {
         'V' : 0,
         'H' : 1,
         'Gruneisen_T' : 2,
         'Ks' : 3,
         'Mu' : 4,
         'V0' : 5,
         'Vp' : 6,
         'Vs' : 7,
         'Vp/Vs' : 8,
         'rho' : 9,
         'unknown10' : 10,
         'Cp' : 11,
         'alpha' : 12,
         'beta' : 13,
         'S' : 14,
         'unknown15' : 15,
         'N' : 16,
         'unknown17' : 17,
         'unknown18' : 18,
         'unknown19' : 19,
         'unknown20' : 20,
         'unknown21' : 21,
         'unknown22' : 22,
         'unknown23' : 23,
         'unknown24' : 24,
         'unknown25' : 25,
         'unknown26' : 26,
         'Cp/Cv' : 27,
         'unknown28'  : 28

    }
    syspropidx = [
         'V',
         'H',
         'Gruneisen_T',
         'Ks bar',
         'Mu bar',
         'V0 km/s',
         'Vp km/s',
         'Vs km/s',
         'Vp/Vs',
         'rho',
         'unknown10',
         'Cp',
         'alpha',
         'beta',
         'S',
         'unknown15',
         'N',
         'unknown17',
         'unknown18',
         'unknown19',
         'unknown20',
         'unknown21',
         'unknown22',
         'unknown23',
         'unknown24',
         'unknown25',
         'unknown26',
         'Cp/Cv',
         'unknown28'
    ]

    syspropidxunits = [
         'J/bar/mol',
         'J/mol',
         '',
         'bar',
         'bar',
         'km/s',
         'km/s',
         'km/s',
         '',
         'kg/m3',
         'unknown',
         'J/K/mol',
         '1/K',
         '1/bar',
         'J/K/mol',
         'unknown',
         'g/mol',
         'unknown',
         'unknown',
         'unknown',
         'unknown',
         'unknown',
         'unknown',
         'unknown',
         'unknown',
         'unknown',
         'unknown',
         '',
         'unknown'
    ]


    def phaseq(self, P, T, composition, debug = False):

        if debug:
            dbgprint = 1
        else:
            dbgprint = 0

        comp_type = c_double * self.ncomp
        wtphases_type = c_double * self.p_size_phases
        cphases_type = c_double * (self.p_size_phases * self.ncomp)
        phsysprop_type = c_double * (self.p_size_phases * self.p_size_sysprops)
        sysprop_type = c_double * self.p_size_sysprops
        namephases_type = c_char * (self.p_size_phases * self.p_pname_len)
        nphases_type = c_int * 1

        wtphases = wtphases_type()
        cphases = cphases_type()
        sysprop = sysprop_type()
        phsysprop = phsysprop_type()
        namephases = namephases_type()
        nphases = nphases_type()

        use_composition = [composition[i] for i in self.c_prp2clr]
        if dbgprint:
            print "Calling phaseq:", P, T, use_composition
        retval = self.libpx.phaseq(P, T, self.ncomp, comp_type(*use_composition), nphases, wtphases,
                                   cphases, sysprop, phsysprop, namephases, dbgprint)
        nphases = nphases[0]

        py_wtphases = []
        py_namephases = []
        for i in range(0, nphases):
            py_wtphases.append(wtphases[i])
            py_namephases.append((namephases[(i*self.p_pname_len):((i+1)*self.p_pname_len-1)]).strip())

        py_cphases = []
        for i in range(0, nphases):
            py_cphases.append([])
            for j in range(0, self.ncomp):
                py_cphases[i].append(cphases[i * self.ncomp + self.c_clr2prp[j]])

        py_phsysprop = []
        for i in range(0, nphases):
            py_phsysprop.append([])
            for j in range(0, self.p_size_sysprops):
                py_phsysprop[i].append(phsysprop[i * self.p_size_sysprops + j])

        py_sysprop = []
        for i in range(0, self.p_size_sysprops):
            py_sysprop.append(sysprop[i])

        ret = {
            'WTPHASES' : py_wtphases,
            'CPHASES' : py_cphases,
            'SYSPROP' : py_sysprop,
            'PHSYSPROP' : py_phsysprop,
            'NAMEPHASES' : py_namephases,
            'NPHASES' : nphases,
            'RETVAL' : retval
        }

        self.lastResult = ret

        return ret



    def __init__(self, inputfile, debug = False, callerCompnames = []):

        self.libpx.ini_phaseq.restype = c_int
        self.libpx.ini_phaseq.argtypes = [c_char_p]
        self.libpx.ini_phaseq(inputfile)   # b'in26klb_new'

        self.libpx.print_comp_order.restype = None
        self.libpx.print_comp_order.argtypes = []
        #self.libpx.print_comp_order()

        self.libpx.number_of_components.restype = c_int
        self.libpx.number_of_components.argtypes = []

        n = self.libpx.number_of_components()
        if n <= 0:
            raise Exception("Number of components from input file equals zero")
        self.ncomp = n
        compordertxt = " " * (n * (self.p_cname_len-1))

        self.libpx.get_comp_order_list.restype = c_int
        self.libpx.get_comp_order_list.argtypes = [c_char_p]
        self.libpx.get_comp_order_list(compordertxt)


        self.libpx.phaseq.restype = c_int
        self.libpx.phaseq.argtypes = [c_double, c_double, c_int, POINTER(c_double), POINTER(c_int),
                                      POINTER(c_double), POINTER(c_double), POINTER(c_double),
                                      POINTER(c_double), POINTER(c_char), c_int]


        for i in range(0,n):
            this_component = compordertxt[(i*(self.p_cname_len-1)):((i+1)*(self.p_cname_len-1))]
            self.components.append(this_component.strip())

        self.c_clr2prp = [-1] * self.ncomp
        self.c_prp2clr = [-1] * self.ncomp

        findCompFrom = self.components[:]
        if len(callerCompnames) == self.ncomp:
            self.callerComponents = callerCompnames[:]
            for i in range(0,self.ncomp):
                match = False
                for j in range(0,self.ncomp):
                    if callerCompnames[i].strip() == findCompFrom[j].strip():
                        findCompFrom[j] = " " * (self.p_cname_len-1)
                        self.c_clr2prp[i] = j
                        self.c_prp2clr[j] = i
                        match = True
                        break
                if not match:
                    raise Exception("Matching caller component names to PerpleX component names failed")
        else:
            for i in range(0,self.ncomp):
                self.c_clr2prp[i] = i
                self.c_prp2clr[i] = i
                self.callerComponents = self.components[:]


        if debug:
            print ("Components are: ")

            for i in range(0,n):
                print("[" + self.components[i] + "]")

        self.lastResult = None

