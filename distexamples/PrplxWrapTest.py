import PrplxWrap
import sys

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

perplex = PrplxWrap.perplex("crustmod", callerCompnames=comps.keys())
#print comps.keys()
#print comps.values()
#print perplex.components
res = perplex.phaseq(7000., 900., comps.values(), debug=False)

sys.stdout.write("phase\twt%")
for i in range(perplex.ncomp):
    sys.stdout.write("\t" + comps.keys()[i])
sys.stdout.write("\n")
for i in range(res['NPHASES']):
    sys.stdout.write(res['NAMEPHASES'][i] + "\t")
    sys.stdout.write('{:.4}'.format(res['WTPHASES'][i])+"\t")
    for j in range(perplex.ncomp):
        sys.stdout.write('{:.4}'.format(res['CPHASES'][i][j]) + "\t")
    sys.stdout.write("\n")
