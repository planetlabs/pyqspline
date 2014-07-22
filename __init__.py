from qspline import qspline

def pyqspline(n,ns,ds,maxit,tol,wi,wf,x,y):

    if n<4:
        raise ValueError('n must be greater or equal to 4')
    if ns<2:
        raise ValueError('ns must be greater or equal to 2')
    if len(wi)!=3:
        raise ValueError('wi must be size 3')
    if len(wf)!=3:
        raise ValueError('wf must be size 3')
    if len(x)!=n:
        raise ValueError('x must have the same length as n')
    if len(y)!=n:
        raise ValueError('y must have the same length as y')
    for i in range(n):
        if len(y[i])!=4:
            raise ValueError('y must have 4 columns')

    t, q, omega, alpha = qspline(n,ns,ds,maxit,tol,wi,wf,x,flatten(y))
    
    return [t,unflatten(q,ns,4),unflatten(omega,ns,3),unflatten(alpha,ns,3)]

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
    import os
    import numpy
    import matplotlib.pyplot as plt

    indatapath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'in.dat')
    outdatapath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'out.dat')

    datafile = open(indatapath,'r')
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

    ndatafile = open(outdatapath,'r')
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
