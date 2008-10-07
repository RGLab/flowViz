/*

  Code based, very heavily in some parts on the cpp Binner class from
  the scagnostics package.

which includes this statement:
****************
 * Binner
 *
 * Leland Wilkinson (SPSS, Inc.)
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, THE AUTHORS MAKE NO
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
****************
*/

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h>

#include <stdlib.h>

void binHex(int *n, double *x, double *y, int *xBins, int *nBin,
    int *count, double *xbin, double *ybin) {

 // scaling constants

  double con1 = .25;
  double con2 = 1. / 3.;
  double c1 = (double) (*xBins - 1);
  double c2 = c1 / sqrt(3.);
  int jinc = *xBins;
  int iinc = 2 * (*xBins);
//  int nBin = (nBins + 20) * (nBins + 20);

//  int *count = int[nBin];
//  double *xbin = double[nBin];
//  double *ybin = double[nBin];

  // fill bins
  for(int i = 0; i < *nBin; ++i) {
    count[i] = 0;
    xbin[i] = 0;
    ybin[i] = 0;
  }

  for (int i = 0; i < *n; ++i) {
    double sx = c1 * x[i];
    double sy = c2 * y[i];
    int i1 = (int) (sy + .5);
    int j1 = (int) (sx + .5);
    double dy = sy - ((double) i1);
    double dx = sx - ((double) j1);
    double dist1 = dx * dx + 3. * dy * dy;
    int m = 0;
    if (dist1 < con1) {
      m = i1 * iinc + j1;
    } else if (dist1 > con2) {
      m = ((int) sy) * iinc + ((int) sx) + jinc;
    } else {
      int i2 = (int) sy;
      int j2 = (int) sx;
      dy = sy - ((double) i2) - .5;
      dx = sx - ((double) j2) - .5;
      double dist2 = dx * dx + 3. * dy * dy;
      if (dist1 <= dist2) {
        m = i1 * iinc + j1;
      } else {
        m = i2 * iinc + j2 + jinc;
      }
    }
    ++count[m];
    xbin[m] += (x[i] - xbin[m]) / count[m];
    ybin[m] += (y[i] - ybin[m]) / count[m];
  }
}

