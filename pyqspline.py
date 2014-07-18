import qspline
import numpy

def pyqspline(n, ns, ds, maxit, tol, wi, wf, x, y):

    # Should probably perform some checks on the inputs
    # note that wi,wf,x and y should be numpy array
    
    [t,q,omega,alpha] = qspline.qspline(n,ns,ds,maxit,tol,wi,wf,x,y)
    result = {'t':t,
              'q':q,
              'omega': omega,
              'alpha': alpha,
              }
    return result

def test():
    datafile = open('in.dat','r')

    wi = numpy.array([float(w) for w in datafile.readline().split(' ')],dtype="double")
    wf = numpy.array([float(w) for w in datafile.readline().split(' ')],dtype="double")

    nout = int(datafile.readline())

    data = []
    for row in datafile.readlines():
        data.append([q for q in row.strip().split(' ') if q!=''])
    data = numpy.array(data,dtype="double")

    result = pyqspline(data.shape[0],nout,-1,10,1.0e-6,wi,wf,data[:,0],data[:,1:])

    print result.get('q')

if __name__=='__main__':
    test()
