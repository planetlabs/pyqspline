import qspline

def interpolate(n,ns,ds,maxit,tol,wi,wf,x,y):

    if ds>0:
        ns = int((x[-1]-x[0])/ds)+1

    # tuples throw off qspline, check sizes and also
    # make sure they are list objects
    if n<4:
        raise ValueError('n must be greater or equal to 4')
    if ns<2:
        raise ValueError('ns must be greater or equal to 2')
    wi = list(wi)
    if len(wi)!=3:
        raise ValueError('wi must be size 3')
    wf = list(wf)
    if len(wf)!=3:
        raise ValueError('wf must be size 3')
    x = list(x)
    if len(x)!=n:
        raise ValueError('x must have the same length as n')
    y = list(y)
    if len(y)!=n:
        raise ValueError('y must have the same length as y')
    for i in range(n):
        y[i] = list(y[i])
        if len(y[i])!=4:
            raise ValueError('y must have 4 columns')

    t,q,omega,alpha = qspline.pytocqspline(n,ns,ds,maxit,tol,wi,wf,x,flatten(y))
    
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