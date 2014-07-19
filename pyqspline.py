import qspline
import numpy
import matplotlib.pyplot as plt

def flatten(nestedlist):
    return [item for sub in nestedlist for item in sub] 

def pyqspline(n,ns,ds,maxit,tol,wi,wf,x,y):
    qspline.qspline(n,ns,ds,maxit,tol,wi,wf,x,flatten(y))
    exit(0)

def test():
    datafile = open('in.dat','r')
    wi = [float(w) for w in datafile.readline().split(' ')]
    wf = [float(w) for w in datafile.readline().split(' ')]
    nout = int(datafile.readline())
    x = []
    y = []
    for row in datafile.readlines():
        t,qx,qy,qz,qs = [float(q) for q in row.strip().split(' ') if q!='']
        x.append(t)
        y.append([qx,qy,qz,qs])

    rt,rq,romega,ralpha = pyqspline(len(x),nout,-1,10,0.000001,wi,wf,x,y)

    ndatafile = open('out.dat','r')
    rows =  ndatafile.readlines()
    nt = []
    nq = []
    nomega = []
    nalpha = []
    for i in range(len(rows)/3):
        t,qx,qy,qz,qs = [x for x in rows[i*3].strip().split(' ') if x!='']
        nt.append(t)
        nq.append([qx,qy,qz,qs])
        nomega.append([x for x in rows[i*3+1].strip().split(' ') if x!=''])
        nalpha.append([x for x in rows[i*3+2].strip().split(' ') if x!=''])

    plt.plot(nt, numpy.array(nq)[:,3], 'go')
    plt.plot(rt, numpy.array(rq)[:,3], 'ro')
    plt.show()
    
if __name__=='__main__':
    test()
