Pyqspline
=========
A python interface for the [qspline project](http://qspline.sourceforge.net/) that provides quaternion interpolation. This is a project originally authored by James McEnnan. Planet Labs only added the python interface, and is not actively supporting it in any way.

Installation
------------
    sudo python setup.py install

Usage
-----
    import pyqspline
    
    x = range(100) # a list of times
    y = [[1,0,0,0] for t in x] # a list of quaternions for these times (qx,qy,qz,qs)
    wi = [0,0,0] # initial angular velocity
    wf = [0,0,0] # final angular velocity
    
    times,quats,omega,alpha = pyqspline.interpolate(len(x),0,0.5,10,0.000001,wi,wf,x,y)
    
    time_series = zip(times,quats)
    
Contact
-------
* James McEnnan, jmcennan@mailaps.org
* Planet Labs, pyqspline@planet.com
