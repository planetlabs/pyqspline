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

from __future__ import division

import os
import math
import unittest
import pyqspline


class PyQsplineTest(unittest.TestCase):

    def test(self):

        samplespath = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), 'samples.dat')
        nominalpath = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), 'nominal.dat')

        with open(samplespath) as f:
            wi = [float(w) for w in f.readline().split(' ')]
            wf = [float(w) for w in f.readline().split(' ')]

            x = []
            y = []
            for row in f.readlines():
                t, qx, qy, qz, qs = [
                    float(q) for q in row.strip().split(' ') if q != '']
                x.append(t)
                y.append([qx, qy, qz, qs])

        rt, rq, romega, ralpha = pyqspline.interpolate(
            len(x), 0, 0.5, 10, 0.000001, wi, wf, x, y)

        with open(nominalpath) as f:
            rows = f.readlines()
        nt = []
        nq = []
        nomega = []
        nalpha = []

        rng_stop = len(rows) / 3
        assert rng_stop.is_integer()
        rng_stop = int(rng_stop)
        for i in range(rng_stop):
            t, qx, qy, qz, qs = [
                x for x in rows[i * 3].strip().split(' ') if x != '']
            nt.append(float(t))
            nq.append([float(qx), float(qy), float(qz), float(qs)])
            nomega.append([
                float(x)
                for x in rows[i * 3 + 1].strip().split(' ') if x != ''])
            nalpha.append([
                float(x)
                for x in rows[i * 3 + 2].strip().split(' ') if x != ''])

        testcount = 0
        tol = 1e-5
        for tt in rt:
            if tt in nt:
                testcount += 1
                ni = nt.index(tt)
                nqstar = nq[ni]
                nomegastar = nomega[ni]
                nalphastar = nalpha[ni]
                ri = rt.index(tt)
                rqstar = rq[ri]
                romegastar = romega[ri]
                ralphastar = ralpha[ri]
                self.assertLess(
                    math.sqrt(sum(
                        [(qa - qb) ** 2 for qa, qb in zip(rqstar, nqstar)])),
                    tol)
                self.assertLess(
                    math.sqrt(sum([
                        (oa - ob) ** 2
                        for oa, ob in zip(romegastar, nomegastar)])),
                    tol)
                self.assertLess(
                    math.sqrt(sum([
                        (aa - ab) ** 2
                        for aa, ab in zip(ralphastar, nalphastar)])),
                    tol)
        self.assertGreater(
            testcount, 0,
            ('Could not find the initial sample points in the resulting '
             'interpolation. No testing possible.'))

if __name__ == '__main__':
    unittest.main()
