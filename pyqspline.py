"""
COPYRIGHT (C) 2014 by Planet Labs

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
    02110-1301 USA.
"""

import qspline


def interpolate(n, ns, ds, maxit, tol, wi, wf, x, y):
    """ Interpolate quaternions

    Args:
        n: number of input points (n >= 4).
        ns: number of output points (ns >= 2).
        ds: interval between output points
            (if ds>0.0, it is used over ns).
        maxit: maximum number of iterations.
        tol: convergence tolerance (rad/sec) for iteration termination.
        wi: initial angular rate vector [wx,wy,wz].
        wf: final angular rate vector [wx,wy,wz].
        x: input list of time values.
        y: input list of quaternion lists [qx,qy,qz,qs].
    Returns:
        t: output list of time values.
        q: output list of interpolated quaternion values.
        omega: output list of interpolated angular rate [wx,wy,wz]
            (rad/sec).
        alpha: output list of interpolated angular acceleration
            [ax,ay,az] (rad/sec^2).

    """

    if ds > 0:
        ns = int((x[-1] - x[0]) / ds) + 1

    # tuples throw off qspline, check sizes and also
    # make sure they are list objects
    if n < 4:
        raise ValueError('n must be greater or equal to 4')
    if ns < 2:
        raise ValueError('ns must be greater or equal to 2')
    wi = list(wi)
    if len(wi) != 3:
        raise ValueError('wi must be size 3')
    wf = list(wf)
    if len(wf) != 3:
        raise ValueError('wf must be size 3')
    x = list(x)
    if len(x) != n:
        raise ValueError('x must have the same length as n')
    y = list(y)
    if len(y) != n:
        raise ValueError('y must have the same length as y')
    for i in range(n):
        y[i] = list(y[i])
        if len(y[i]) != 4:
            raise ValueError('y must have 4 columns')

    t, q, omega, alpha = qspline.pytocqspline(
        n, ns, ds, maxit, tol, wi, wf, x, flatten(y))

    return [
        t,
        unflatten(q, ns, 4),
        unflatten(omega, ns, 3),
        unflatten(alpha, ns, 3),
    ]


def flatten(nestedlist):
    """ Flattens a list

    Args:
        nestedlist: a list of lists
    Returns:
        a list containing all the elements of nestedlist

    """

    return [item for sub in nestedlist for item in sub]


def unflatten(flatlist, m, n):
    """ Construct a nested list from a flat list

    Args:
        flatlist: a list of elements
        m: the number of lists inside the final list
        n: the number of element in each sublist
    Returns:
        a nested list of size m, each with n elements

    """

    out = []
    for r in range(m):
        row = []
        for c in range(n):
            row.append(flatlist[n * r + c])
        out.append(row)

    return out
