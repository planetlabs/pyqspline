#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define DZERO (double **)0
#define ZERO (double *)0

static PyObject* pyqspline(PyObject *self, PyObject *args)
{

  int n, ns, maxit;
  double ds, tol;
  PyObject *wi, *wf, *x, *y;

  if (!PyArg_ParseTuple(args, "iididOOOO", &n, &ns, &ds, &maxit, &tol, &wi, &wf, &x, &y))
    return NULL;

  wi = PyArray_FROM_OTF(wi, NPY_DOUBLE, NPY_IN_ARRAY);
  wf = PyArray_FROM_OTF(wf, NPY_DOUBLE, NPY_IN_ARRAY);
  x = PyArray_FROM_OTF(x, NPY_DOUBLE, NPY_IN_ARRAY);
  if (wi==NULL || wf==NULL || x==NULL){
    Py_XDECREF(wi);
    Py_XDECREF(wf);
    Py_XDECREF(x);
    return NULL;
  }
  double* wiptr = (double *) PyArray_DATA(wi);
  double* wfptr = (double *) PyArray_DATA(wf);
  double* xptr = (double *) PyArray_DATA(x);

  double **yptr;
  npy_intp ydims[2];
  PyArray_Descr *descr = PyArray_DescrFromType(NPY_DOUBLE);
  if (PyArray_AsCArray(&y,(void *)&yptr,(npy_intp *)&ydims,2,descr) < 0) {
    Py_XDECREF(wi);
    Py_XDECREF(wf);
    Py_XDECREF(x);
    return NULL;
  }

  /*
  double* t = (double *) malloc(ns*sizeof(double));
  double** q = (double **) malloc(ns * sizeof(double*));
  double** omega = (double **) malloc(ns * sizeof(double*));
  double** alpha = (double **) malloc(ns * sizeof(double*));
  int row;
  for (row=0; row<ns; row++) {
    q[row] = (double *) malloc(4 * sizeof(double));
  }
  for (row=0; row<ns; row++) {
    omega[row] = (double *) malloc(3 * sizeof(double));
  }
  for (row=0; row<ns; row++) {
    alpha[row] = (double *) malloc(3 * sizeof(double));
  }
  */
  double t[10000];
  double q[10000][4];
  double omega[10000][3];
  double alpha[10000][3];
  
  qspline(n,ns,ds,maxit,tol,wiptr,wfptr,xptr,yptr,t,q,omega,alpha);

  Py_DECREF(wi);
  Py_DECREF(wf);
  Py_DECREF(x);
  Py_DECREF(y);
  
  /*
  free(t);
  for (row=0; row<ns; row++) {
    free(q[row]);
  }
  for (row=0; row<ns; row++) {
    free(omega[row]);
  }
  for (row=0; row<ns; row++) {
    free(alpha[row]);
  }
  free(q);
  free(omega);
  free(alpha);
  */

  //return Py_BuildValue("OOOO", t, q, omega, alpha);
  return Py_None;

}

static PyMethodDef QsplineMethods[] = {
  {"qspline", pyqspline, METH_VARARGS,
   "Produces a quaternion spline interpolation of sparse data."},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initqspline(void)
{
  (void) Py_InitModule("qspline", QsplineMethods);
  import_array();
}

qspline(n,ns,ds,maxit,tol,wi,wf,x,y,t,q,omega,alpha)
/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

purpose

Subroutine qspline produces a quaternion spline interpolation of sparse data.
The method is based on a third-order polynomial expansion of the
rotation angle vector.

calling sequence

variable     i/o     description
--------     ---     -----------

n             i      number of input points (n >= 4).

ns            i      number of output points (ns >= 2).

ds            i      interval between output points (ds > 0.0).

maxit         i      maximum number of iterations.

tol           i      convergence tolerance (rad/sec) for iteration termination.

wi            i      initial angular rate vector.

wf            i      final angular rate vector.

x             i      pointer to input vector of time values.

y             i      pointer to input vector of quaternion values.

t             o      pointer to output vector of time values.

q             o      pointer to output array of interpolated quaternion values.

omega         o      pointer to output array of interpolated
                     angular rate values (rad/sec).

alpha         o      pointer to output array of interpolated
                     angular acceleration values (rad/sec^2).

return value

 ns >= 2 -> normal return, no error
 -1      -> insufficient input/output data (n < 4 or ns < 2)
 -2      -> independent variable is not monotonic increasing
 -3      -> memory allocation failure

external references

freeall
getang
rates
slew3_init
slew3

programming

J. J. McEnnan, April, 2003.

COPYRIGHT (C) 2003 by James McEnnan

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
int n, ns, maxit;
double ds, tol, wi[], wf[], x[], y[][4], t[], q[][4], omega[][3], alpha[][3];
{
  int i, j;
  double dx, xi, *h, *a, *b, *c, *dtheta, **e, **w, **wprev;
  double dum1[3], dum2[3], getang();

  /* error checking. */

  if(n < 4)
  {
    fprintf(stderr,"qspline: insufficient input data.\n");

    return -1;
  }

  if(ds > 0.0)
  {
    dx = ds;

    ns = (int) ((x[n - 1] - x[0])/ds) + 1;
  }
  else
    dx = (x[n - 1] - x[0])/(ns - 1);

  if(ns < 2)
  {
    fprintf(stderr,"qspline: too few output points.\n");

    return -1;
  }

  i = n*sizeof(double);

  if((w = (double **) malloc((unsigned) i)) == DZERO)
  {
    fprintf(stderr,"qspline: memory allocation failure.\n");

    return -3;
  }

  if((wprev = (double **) malloc((unsigned) i)) == DZERO)
  {
    fprintf(stderr,"qspline: memory allocation failure.\n");

    free(w);

    return -3;
  }

  i = (n - 1)*sizeof(double);

  if((a = (double *) malloc((unsigned) i)) == ZERO)
  {
    fprintf(stderr,"qspline: memory allocation failure.\n");

    free(w);
    free(wprev);

    return -3;
  }

  if((b = (double *) malloc((unsigned) i)) == ZERO)
  {
    fprintf(stderr,"qspline: memory allocation failure.\n");

    free(w);
    free(wprev);
    free(a);

    return -3;
  }

  if((c = (double *) malloc((unsigned) i)) == ZERO)
  {
    fprintf(stderr,"qspline: memory allocation failure.\n");

    free(w);
    free(wprev);
    free(a);
    free(b);

    return -3;
  }

  if((h = (double *) malloc((unsigned) i)) == ZERO)
  {
    fprintf(stderr,"qspline: memory allocation failure.\n");

    free(w);
    free(wprev);
    free(c);
    free(b);
    free(a);

    return -3;
  }

  if((dtheta = (double *) malloc((unsigned) i)) == ZERO)
  {
    fprintf(stderr,"qspline: memory allocation failure.\n");

    free(w);
    free(wprev);
    free(c);
    free(b);
    free(a);
    free(dtheta);

    return -3;
  }

  if((e = (double **) malloc((unsigned) i)) == DZERO)
  {
    fprintf(stderr,"qspline: memory allocation failure.\n");

    free(w);
    free(wprev);
    free(e);
    free(dtheta);
    free(c);
    free(b);
    free(a);

    return -3;
  }

  for(i = 0;i < n;i++)
    w[i] = ZERO;

  for(i = 0;i < n;i++)
    wprev[i] = ZERO;

  for(i = 0;i < n - 1;i++)
    e[i] = ZERO;

  j = 3*sizeof(double);

  for(i = 0;i < n;i++)
    if((w[i] = (double *) malloc((unsigned) j)) == ZERO)
    {
      fprintf(stderr,"qspline: memory allocation failure.\n");

      freeall(n,h,a,b,c,dtheta,e,w,wprev);

      return -3;
    }

  for(i = 0;i < n;i++)
    if((wprev[i] = (double *) malloc((unsigned) j)) == ZERO)
    {
      fprintf(stderr,"qspline: memory allocation failure.\n");

      freeall(n,h,a,b,c,dtheta,e,w,wprev);

      return -3;
    }

  for(i = 0;i < n - 1;i++)
    if((e[i] = (double *) malloc((unsigned) j)) == ZERO)
    {
      fprintf(stderr,"qspline: memory allocation failure.\n");

      freeall(n,h,a,b,c,dtheta,e,w,wprev);

      return -3;
    }

  for(i = 0;i < n;i++)
    for(j = 0;j < 3;j++)
      w[i][j] = 0.0;

  for(i = 0;i < n - 1;i++)
  {
    h[i] = x[i + 1] - x[i];

    if(h[i] <= 0.0)
    {
      fprintf(stderr,"qspline: x is not monotonic.\n");

      freeall(n,h,a,b,c,dtheta,e,w,wprev);

      return -2;
    }
  }

  /* compute spline coefficients. */

  for(i = 0;i < n - 1;i++)
    dtheta[i] = getang(y[i],y[i + 1],e[i]);

  rates(n,maxit,tol,wi,wf,h,a,b,c,dtheta,e,w,wprev);

  /* interpolate and output results. */

  i = 0;
  xi = x[0];

  slew3_init(h[0],dtheta[0],e[0],w[0],dum1,w[1],dum2);

  for(j = 0;j < ns;j++)
  {
    while(xi >= x[i + 1] && i < n - 2)
    {
      i++;

      slew3_init(h[i],dtheta[i],e[i],w[i],dum1,w[i + 1],dum2);
    }

    t[j] = xi;

    slew3(xi - x[i],h[i],y[i],q[j],omega[j],alpha[j],dum1);

    xi += dx;
  }

  freeall(n,h,a,b,c,dtheta,e,w,wprev);

  return ns;
}

freeall(n,h,a,b,c,dtheta,e,w,wprev)
/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
purpose

subroutine freeall frees previously allocated memory.

calling sequence

variable     i/o     description
--------     ---     -----------

n             i      number of input data points.

h             i      pointer to vector of x-interval values.

a             i      pointer to intermediate work space.

b             i      pointer to intermediate work space.

c             i      pointer to intermediate work space.

dtheta        i      pointer to vector of rotation angles.

e             i      pointer to array of rotation axis vectors.

w             i      pointer to output intermediate angular rates.

wprev         i      pointer to previous intermediate angular rates.

return value

none

external references

none

programming

J. J. McEnnan, April, 2003.

COPYRIGHT (C) 2003 by James McEnnan

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
int n;
double *h, *a, *b, *c, *dtheta, **e, **w, **wprev;
{
  int i;

  for(i = 0;i < n;i++)
    if(w[i] != ZERO)
      free(w[i]);

  for(i = 0;i < n;i++)
    if(wprev[i] != ZERO)
      free(wprev[i]);

  for(i = 0;i < n - 1;i++)
    if(e[i] != ZERO)
      free(e[i]);

  free(w);
  free(wprev);
  free(e);
  free(dtheta);
  free(c);
  free(b);
  free(a);
  free(h);
}

double getang(qi,qf,e)
/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

purpose

Subroutine getang computes the slew angle and axis between the input initial and
final states.

calling sequence

variable     i/o     description
--------     ---     -----------

qi            i      initial attitude quaternion.

qf            i      final attitude quaternion.

e             o      unit vector along slew eigen-axis.

return value

slew angle in radians

external references

unvec

programming

J. J. McEnnan, May, 2000.

COPYRIGHT (C) 2003 by James McEnnan

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
double qi[], qf[], e[];
{
  double dtheta, sa, ca, temp[3], unvec();

  temp[0] = qi[3]*qf[0] - qi[0]*qf[3] - qi[1]*qf[2] + qi[2]*qf[1];
  temp[1] = qi[3]*qf[1] - qi[1]*qf[3] - qi[2]*qf[0] + qi[0]*qf[2];
  temp[2] = qi[3]*qf[2] - qi[2]*qf[3] - qi[0]*qf[1] + qi[1]*qf[0];

  ca =  qi[0]*qf[0] + qi[1]*qf[1] + qi[2]*qf[2] + qi[3]*qf[3];

  sa = unvec(temp,e);

  dtheta = 2.0*atan2(sa,ca);

  return dtheta;
}

rates(n,maxit,tol,wi,wf,h,a,b,c,dtheta,e,w,wprev)
/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

purpose

subroutine rates computes intermediate angular rates for interpolation.

calling sequence

variable     i/o     description
--------     ---     -----------

n             i      number of input data points.

maxit         i      maximum number of iterations.

tol           i      convergence tolerance (rad/sec) for iteration termination.

wi            i      initial angular rate vector.

wf            i      final angular rate vector.

h             i      pointer to vector of time interval values.

a             i      pointer to intermediate work space.

b             i      pointer to intermediate work space.

c             i      pointer to intermediate work space.

dtheta        i      pointer to vector of rotation angles.

e             i      pointer to array of rotation axis vectors.

w             o      pointer to output intermediate angular rate values.

wprev         o      pointer to previous intermediate angular rate values.

return value

none

external references

bd
rf

programming

J. J. McEnnan, April, 2003.

COPYRIGHT (C) 2003 by James McEnnan

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
int n, maxit;
double tol, wi[], wf[], h[], *a, *b, *c, *dtheta, **e, **w, **wprev;
{
  int i, j, iter;
  double dw, temp1[3], temp2[3];

  iter = 0;

  do                                                 /* start iteration loop. */
  {
    for(i = 1;i < n - 1;i++)
      for(j = 0;j < 3;j++)
        wprev[i][j] = w[i][j];

    /* set up the tridiagonal matrix. d initially holds the RHS vector array;
       it is then overlaid with the calculated angular rate vector array. */

    for(i = 1;i < n - 1;i++)
    {
      a[i] = 2.0/h[i - 1];
      b[i] = 4.0/h[i - 1] + 4.0/h[i];
      c[i] = 2.0/h[i];
  
      rf(e[i - 1],dtheta[i - 1],wprev[i],temp1);

      for(j = 0;j < 3;j++)
        w[i][j] = 6.0*(dtheta[i - 1]*e[i - 1][j]/(h[i - 1]*h[i - 1]) +
                       dtheta[i    ]*e[i    ][j]/(h[i    ]*h[i    ])) -
                  temp1[j];
    }
  
    bd(e[0    ],dtheta[0    ],1,wi,temp1);
    bd(e[n - 2],dtheta[n - 2],0,wf,temp2);
  
    for(j = 0;j < 3;j++)
    {
      w[1    ][j] -= a[1    ]*temp1[j];
      w[n - 2][j] -= c[n - 2]*temp2[j];
    }
  
    /* reduce the matrix to upper triangular form. */
  
    for(i = 1;i < n - 2;i++)
    {
      b[i + 1] -= c[i]*a[i + 1]/b[i];
  
      for(j = 0;j < 3;j++)
      {
        bd(e[i],dtheta[i],1,w[i],temp1);
        
        w[i + 1][j] -= temp1[j]*a[i + 1]/b[i];
      }
    }
  
    /* solve using back substitution. */
  
    for(j = 0;j < 3;j++)
      w[n - 2][j] /= b[n - 2];
  
    for(i = n - 3;i > 0;i--)
    {
      bd(e[i],dtheta[i],0,w[i + 1],temp1);
  
      for(j = 0;j < 3;j++)
        w[i][j] = (w[i][j] - c[i]*temp1[j])/b[i];
    }
  
    dw = 0.0;

    for(i = 1;i < n - 1;i++)
      for(j =  0;j < 3;j++)
        dw += (w[i][j] - wprev[i][j])*(w[i][j] - wprev[i][j]);

    dw = sqrt(dw);
  }
  while(iter++ < maxit && dw > tol);

  /* solve for end conditions. */
  
  for(j = 0;j < 3;j++)
  {
    w[0    ][j] = wi[j];
    w[n - 1][j] = wf[j];
  }
}

#define EPS 1.0e-6
bd(e,dtheta,flag,xin,xout)
/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

purpose

Subroutine bd performs the transformation between the coefficient vector and
the angular rate vector.

calling sequence

variable     i/o     description
--------     ---     -----------

e             i      unit vector along slew eigen-axis.

dtheta        i      slew angle (rad).

flag          i      flag determining direction of transformation.
                      = 0 -> compute coefficient vector from
                      angular rate vector
                      = 1 -> compute angular rate vector from
                      coefficient vector

xin           i      input vector.

xout          o      output vector.

return value

 0 -> no error
-1 -> transformation direction incorrectly specified.

external references

crossp

programming

J. J. McEnnan, April, 2003.

COPYRIGHT (C) 2003 by James McEnnan

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
int flag;
double dtheta, e[], xin[], xout[];
{
  int i;
  double sa, ca, b0, b1, b2, temp1[3], temp2[3];

  if(dtheta > EPS)
  {
    ca = cos(dtheta);
    sa = sin(dtheta);

    if(flag == 0)
    {
      b1 = 0.5*dtheta*sa/(1.0 - ca);
      b2 = 0.5*dtheta;
    }
    else if(flag == 1)
    {
      b1 = sa/dtheta;
      b2 = (ca - 1.0)/dtheta;
    }
    else
      return -1;

    b0 = xin[0]*e[0] + xin[1]*e[1] + xin[2]*e[2];

    crossp(e,xin,temp2);

    crossp(temp2,e,temp1);

    for(i = 0;i < 3;i++)
      xout[i] = b0*e[i] + b1*temp1[i] + b2*temp2[i];
  }
  else
  {
    for(i = 0;i < 3;i++)
      xout[i] = xin[i];
  }

  return 0;
}

rf(e,dtheta,win,rhs)
/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

purpose

Subroutine rf computes the non-linear rate contributions to the final
angular acceleration.

calling sequence

variable     i/o     description
--------     ---     -----------

e             i      unit vector along slew eigen-axis.

dtheta        i      slew angle (rad).

win           i      input final angular rate vector.

rhs           o      output vector containing non-linear rate contributions
                     to the final acceleration.

return value

none

external references

crossp

programming

J. J. McEnnan, May, 2003.

COPYRIGHT (C) 2003 by James McEnnan

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
double dtheta, e[], win[], rhs[];
{
  int i;
  double sa, ca, dot, mag, c1, r0, r1, temp1[3], temp2[3];

  if(dtheta > EPS)
  {
    ca = cos(dtheta);
    sa = sin(dtheta);

    crossp(e,win,temp2);

    crossp(temp2,e,temp1);

    dot = win[0]*e[0] + win[1]*e[1] + win[2]*e[2];

    mag = win[0]*win[0] + win[1]*win[1] + win[2]*win[2];

    c1 = (1.0 - ca);

    r0 = 0.5*(mag - dot*dot)*(dtheta - sa)/c1;

    r1 = dot*(dtheta*sa - 2.0*c1)/(dtheta*c1);

    for(i = 0;i < 3;i++)
      rhs[i] = r0*e[i] + r1*temp1[i];
  }
  else
  {
    for(i = 0;i < 3;i++)
      rhs[i] = 0.0;
  }
}

static double a[3][3], b[3][3], c[2][3], d[3];
slew3_init(dt,dtheta,e,wi,ai,wf,af)
/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

purpose

Subroutine slew3_init computes the coefficients for a third-order polynomial
interpolation function describing a slew between the input initial and
final states.

calling sequence

variable     i/o     description
--------     ---     -----------

dt            i      slew time (sec).

dtheta        i      slew angle (rad).

e             i      unit vector along slew eigen-axis.

wi            i      initial body angular rate (rad/sec).

ai            i      initial body angular acceleration (rad/sec^2)
                     (included for compatibility only).

wf            i      final body angular rate (rad/sec).

af            i      final body angular acceleration (rad/sec^2)
                     (included for compatibility only).

return value

none

external references

none

programming

J. J. McEnnan, March, 2003.

COPYRIGHT (C) 2003 by James McEnnan

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
double dt, dtheta, e[], wi[], ai[], wf[], af[];
{
  int i;
  double sa, ca, c1, c2;
  double b0, bvec1[3], bvec2[3], bvec[3];

  if(dt <= 0.0)
    return;

  sa = sin(dtheta);
  ca = cos(dtheta);

  /* final angular rate terms. */

  if(dtheta > EPS)
  {
    c1 = 0.5*sa*dtheta/(1.0 - ca);

    c2 = 0.5*dtheta;

    b0 = e[0]*wf[0] + e[1]*wf[1] + e[2]*wf[2];

    crossp(e,wf,bvec2);

    crossp(bvec2,e,bvec1);

    for(i = 0;i < 3;i++)
      bvec[i] = b0*e[i] + c1*bvec1[i] + c2*bvec2[i];
  }
  else
  {
    for(i = 0;i < 3;i++)
      bvec[i] = wf[i];
  }

  /* compute coefficients. */

  for(i = 0;i < 3;i++)
  {
    b[0][i] = wi[i];
    a[2][i] = e[i]*dtheta;
    b[2][i] = bvec[i];

    a[0][i] =  b[0][i]*dt;
    a[1][i] = (b[2][i]*dt - 3.0*a[2][i]);

    b[1][i] = (2.0*a[0][i] + 2.0*a[1][i])/dt;
    c[0][i] = (2.0*b[0][i] +     b[1][i])/dt;
    c[1][i] = (    b[1][i] + 2.0*b[2][i])/dt;

       d[i] = (    c[0][i] +     c[1][i])/dt;
  }
}

slew3(t,dt,qi,q,omega,alpha,jerk)
/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

purpose

Subroutine slew3 computes the quaternion, body angular rate, acceleration and
jerk as a function of time corresponding to a third-order polynomial
interpolation function describing a slew between initial and final states.

calling sequence

variable     i/o     description
--------     ---     -----------

t             i      current time (seconds from start).

dt            i      slew time (sec).

qi            i      initial attitude quaternion.

q             o      current attitude quaternion.

omega         o      current body angular rate (rad/sec).

alpha         o      current body angular acceleration (rad/sec^2).

jerk          o      current body angular jerk (rad/sec^3).

return value

none

external references

unvec

programming

J. J. McEnnan, March, 2003.

COPYRIGHT (C) 2003 by James McEnnan

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
double t, dt, qi[], q[], omega[], alpha[], jerk[];
{
  int i;
  double x, ang, sa, ca, u[3], x1[2], unvec();
  double th0[3], th1[3], th2[3], th3[3], temp0[3], temp1[3], temp2[3];
  double thd1, thd2, thd3, w2, td2, ut2, wwd;
  double w[3], udot[3], wd1[3], wd1xu[3], wd2[3], wd2xu[3];

  if(dt <= 0.0)
    return;

  x = t/dt;

  x1[0] = x - 1.0;

  for(i = 1;i < 2;i++)
    x1[i] = x1[i - 1]*x1[0];

  for(i = 0;i < 3;i++)
  {
    th0[i] = ((x*a[2][i] + x1[0]*a[1][i])*x + x1[1]*a[0][i])*x;

    th1[i] = (x*b[2][i] + x1[0]*b[1][i])*x + x1[1]*b[0][i];

    th2[i] = x*c[1][i] + x1[0]*c[0][i];

    th3[i] = d[i];
  }

  ang = unvec(th0,u);

  ca = cos(0.5*ang);
  sa = sin(0.5*ang);

  q[0] = ca*qi[0] + sa*( u[2]*qi[1] - u[1]*qi[2] + u[0]*qi[3]);
  q[1] = ca*qi[1] + sa*(-u[2]*qi[0] + u[0]*qi[2] + u[1]*qi[3]);
  q[2] = ca*qi[2] + sa*( u[1]*qi[0] - u[0]*qi[1] + u[2]*qi[3]);
  q[3] = ca*qi[3] + sa*(-u[0]*qi[0] - u[1]*qi[1] - u[2]*qi[2]);

  ca = cos(ang);
  sa = sin(ang);

  if(ang > EPS)
  {
    /* compute angular rate vector. */

    crossp(u,th1,temp1);

    for(i = 0;i < 3;i++)
      w[i] = temp1[i]/ang;

    crossp(w,u,udot);

    thd1 = u[0]*th1[0] + u[1]*th1[1] + u[2]*th1[2];

    for(i = 0;i < 3;i++)
      omega[i] = thd1*u[i] + sa*udot[i] - (1.0 - ca)*w[i];

    /* compute angular acceleration vector. */

    thd2 = udot[0]*th1[0] + udot[1]*th1[1] + udot[2]*th1[2] +
              u[0]*th2[0] +    u[1]*th2[1] +    u[2]*th2[2];

    crossp(u,th2,temp1);

    for(i = 0;i < 3;i++)
      wd1[i] = (temp1[i] - 2.0*thd1*w[i])/ang;

    crossp(wd1,u,wd1xu);

    for(i = 0;i < 3;i++)
      temp0[i] = thd1*u[i] - w[i];

    crossp(omega,temp0,temp1);

    for(i = 0;i < 3;i++)
      alpha[i] = thd2*u[i] + sa*wd1xu[i] - (1.0 - ca)*wd1[i] +
      thd1*udot[i] + temp1[i];

    /* compute angular jerk vector. */

    w2 = w[0]*w[0] + w[1]*w[1] + w[2]*w[2];

    thd3 = wd1xu[0]*th1[0] + wd1xu[1]*th1[1] + wd1xu[2]*th1[2] -
           w2*(u[0]*th1[0] + u[1]*th1[1] + u[2]*th1[2]) +
           2.0*(udot[0]*th2[0] + udot[1]*th2[1] + udot[2]*th2[2]) +
           u[0]*th3[0] + u[1]*th3[1] + u[2]*th3[2];

    crossp(th1,th2,temp1);

    for(i = 0;i < 3;i++)
      temp1[i] /= ang;

    crossp(u,th3,temp2);

    td2 = (th1[0]*th1[0] + th1[1]*th1[1] + th1[2]*th1[2])/ang;

    ut2 = u[0]*th2[0] + u[1]*th2[1] + u[2]*th2[2];

    wwd = w[0]*wd1[0] + w[1]*wd1[1] + w[2]*wd1[2];

    for(i = 0;i < 3;i++)
      wd2[i] = (temp1[i] + temp2[i] - 2.0*(td2 + ut2)*w[i] -
      4.0*thd1*wd1[i])/ang;

    crossp(wd2,u,wd2xu);

    for(i = 0;i < 3;i++)
      temp2[i] = thd2*u[i] + thd1*udot[i] - wd1[i];

    crossp(omega,temp2,temp1);

    crossp(alpha,temp0,temp2);

    for(i = 0;i < 3;i++)
      jerk[i] = thd3*u[i] + sa*wd2xu[i] - (1.0 - ca)*wd2[i] +
      2.0*thd2*udot[i] + thd1*((1.0 + ca)*wd1xu[i] - w2*u[i] - sa*wd1[i]) -
      wwd*sa*u[i] + temp1[i] + temp2[i];
  }
  else
  {
    crossp(th1,th2,temp1);

    for(i = 0;i < 3;i++)
    {
      omega[i] = th1[i];
      alpha[i] = th2[i];
       jerk[i] = th3[i] - 0.5*temp1[i];
    }
  }
}

double unvec(a,au)
/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

purpose

subroutine unvec unitizes a vector and computes its magnitude.

calling sequence

variable     i/o     description
--------     ---     -----------

a             i      input vector.

au            o      output unit vector.

return value

magnitude of vector a.

external references

none

programming

J. J. McEnnan, December, 1987.

COPYRIGHT (C) 2003 by James McEnnan

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
double a[], au[];
{
  double amag;

  amag = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);

  if(amag > 0.0)
  {
    au[0] = a[0]/amag;
    au[1] = a[1]/amag;
    au[2] = a[2]/amag;
  }
  else
  {
    au[0] = 0.0;
    au[1] = 0.0;
    au[2] = 0.0;
  }

  return amag;
}

crossp(b,c,a)
/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

purpose

subroutine crossp computes the vector cross product b x c.

calling sequence

variable     i/o     description
--------     ---     -----------

b             i      input vector.

c             i      input vector.

a             o      output vector = b x c.

return value

none

external references

none

programming

J. J. McEnnan, February, 1988.

COPYRIGHT (C) 2003 by James McEnnan

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
double a[], b[], c[];
{
  a[0] = b[1]*c[2] - b[2]*c[1];
  a[1] = b[2]*c[0] - b[0]*c[2];
  a[2] = b[0]*c[1] - b[1]*c[0];
}
