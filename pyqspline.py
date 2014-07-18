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

    nomfile = open('out.dat','r')
    nominal = { 't':[]),
                'q':[]),
                'omega':[]),
                'alpha':[]),
                }
    rows =  nomfile.readlines()
    m = len(rows)/3
    for i in range(m):
        t,qx,qy,qz,qs = [x for x in rows[i*3].strip().split(' ') if x!='']
        omega = [x for x in rows[i*3+1].strip().split(' ') if x!='']
        alpha = [x for x in rows[i*3+2].strip().split(' ') if x!='']
        nominal['t'].append(t)
        nomdial['q'].append([qx,qy,qz,qs])
        nominal['omega'].append(omega)
        nominal['alpha'].append(alpha)
    nominal['t'] = numpy.array(nominal.get('t'))
    nominal['q'] = numpy.array(nominal.get('q'))
    nominal['omega'] = numpy.array(nominal.get('omega'))
    nominal['alpha'] = numpy.aaray(nominal.get('alpha'))
    import pdb;psb.set_trace()

if __name__=='__main__':
    test()
