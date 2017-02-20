import PrplxWrap
import sys

perplex = PrplxWrap.perplex("crustmod")
n = perplex.ncomp
comps = []
P = float(raw_input("P: "))
T = float(raw_input("T: "))
for i in range(0, n):
    x = raw_input("Enter " + perplex.components[i] + ":")
    comps.append(float(x))
res = perplex.phaseq(P, T, comps, debug=True)

sys.stdout.write("phase\twt%")
for i in range(perplex.ncomp):
    sys.stdout.write("\t" + perplex.components[i])
sys.stdout.write("\n")
for i in range(res['NPHASES']):
    sys.stdout.write(res['NAMEPHASES'][i] + "\t")
    sys.stdout.write('{:.4}'.format(res['WTPHASES'][i])+"\t")
    for j in range(perplex.ncomp):
        sys.stdout.write('{:.4}'.format(res['CPHASES'][i][j]) + "\t")
    sys.stdout.write("\n")
