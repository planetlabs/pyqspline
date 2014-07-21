import numpy
import matplotlib.pyplot as plt
from qspline import qspline


def pyqspline(n,ns,ds,maxit,tol,wi,wf,x,y):
    t, q, omega, alpha = qspline(n,ns,ds,maxit,tol,wi,wf,x,flatten(y))
    q_out = unflatten(q,ns,4)
    omega_out = unflatten(omega,ns,3)
    alpha_out = unflatten(alpha,ns,3)
    return [t,q_out,omega_out,alpha_out]

def flatten(nestedlist):
    return [item for sub in nestedlist for item in sub] 

def unflatten(flatlist,m,n):
    out = []
    for r in range(m):
        row = []
        for c in range(n):
            row.append(flatlist[n*r+c])
        out.append(row)
    return out

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

    plt.plot(nt, numpy.array(nq), 'gD')
    plt.plot(rt, numpy.array(rq), 'rx')
    plt.show()
    
if __name__=='__main__':
    test()
