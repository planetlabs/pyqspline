import os
import math
import unittest
import pyqspline

class PyQsplineTest(unittest.TestCase):

    def test(self):

        samplespath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'samples.dat')
        nominalpath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'nominal.dat')

        datafile = open(samplespath,'r')
        wi = [float(w) for w in datafile.readline().split(' ')]
        wf = [float(w) for w in datafile.readline().split(' ')]
        nout = int(datafile.readline())
        x = []
        y = []
        for row in datafile.readlines():
            t,qx,qy,qz,qs = [float(q) for q in row.strip().split(' ') if q!='']
            x.append(t)
            y.append([qx,qy,qz,qs])

        rt,rq,romega,ralpha = pyqspline.interpolate(len(x),0,0.5,10,0.000001,wi,wf,x,y)

        ndatafile = open(nominalpath,'r')
        rows =  ndatafile.readlines()
        nt = []
        nq = []
        nomega = []
        nalpha = []
        for i in range(len(rows)/3):
            t,qx,qy,qz,qs = [x for x in rows[i*3].strip().split(' ') if x!='']
            nt.append(float(t))
            nq.append([float(qx),float(qy),float(qz),float(qs)])
            nomega.append([float(x) for x in rows[i*3+1].strip().split(' ') if x!=''])
            nalpha.append([float(x) for x in rows[i*3+2].strip().split(' ') if x!=''])

        testcount = 0
        tol = 1e-5
        for tt in rt:
            if tt in nt:
                testcount += 1
                ni = nt.index(tt)
                nqstar = nq[ni]
                nomegastar = nomega[ni]
                nalphastar = nalpha[ni]
                ri = rt.index(tt)
                rqstar = rq[ri]
                romegastar = romega[ri]
                ralphastar = ralpha[ri]
                self.assertLess(math.sqrt(sum([(qa-qb)**2 for qa, qb in zip(rqstar,nqstar)])),tol)
                self.assertLess(math.sqrt(sum([(oa-ob)**2 for oa, ob in zip(romegastar,nomegastar)])),tol)
                self.assertLess(math.sqrt(sum([(aa-ab)**2 for aa, ab in zip(ralphastar,nalphastar)])),tol)
        self.assertGreater(testcount,0,'Could not find the initial sample points in the resulting interpolation. No testing possible.')

if __name__=='__main__':
    unittest.main()