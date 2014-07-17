import qspline
import numpy

def pyqspline(n, ns, ds, maxit, tol, wi, wf, x, y):

    # Should probably perform some checks on the inputs
    # note that wi,wf,x and y should be numpy array
    
    print str(y)
    output = qspline.qspline(n,ns,ds,maxit,tol,wi,wf,x,y)
    print "OUTPUT: %s" % str(output)

    t = None
    q = None
    omega = None
    alpha = None

    return [t,q,omega,alpha]

def test():
    datafile = open('in.dat','r')

    wi = numpy.array([float(w) for w in datafile.readline().split(' ')],dtype="double")
    wf = numpy.array([float(w) for w in datafile.readline().split(' ')],dtype="double")

    nout = int(datafile.readline())

    data = []
    for row in datafile.readlines():
        data.append([q for q in row.strip().split(' ') if q!=''])
    data = numpy.array(data,dtype="double")

    [t,q,omega,alpha] = pyqspline(data.shape[0],nout,-1,10,1.0e-6,wi,wf,data[:,0],data[:,1:])

if __name__=='__main__':
    test()
