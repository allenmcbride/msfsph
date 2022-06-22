// fparams.h
// Allen McBride
// April 20, 2022
//
// The Fparams class maintains size information to go along with a
// Fields object. It also encapsulated interconversions between physical
// coordinates and grid coordinates.


#ifndef FPARAMS_H
#define FPARAMS_H

#include "point.h"

#include <cmath>
#include <iostream>

struct Fparams {
   int ypixels; // range of y in pixels
   double yphys; // range of y in physical units is 0.0 to this number
   double hy; // height of a grid cell
   int xpixels; // range of x in pixels
   double xphys; // range of x in physical units is 0.0 to this number
   double hx; // width of a grid cell

   // Construct Fparams object based on requested physical dimensions
   // and grid cell size. Enforce square grid cells for now, adjusting
   // physical grid size in the x direction as needed. to achieve square
   // cells. Solute diffusion would need to be rewritten to support
   // non-square cells.
   Fparams(const double xphysRequested, const double yphys, double h) : 
      ypixels{static_cast<int>(std::ceil(yphys / h))},
      yphys{yphys},
      hy{yphys / ypixels},
      xpixels{static_cast<int>(std::ceil(xphysRequested / hy))},
      xphys{hy * xpixels},
      hx{hy}
   {
      if (yphys / h != ypixels)
         std::cerr << "For whole number of grid cells, h reduced from " << h << " to " << hx << std::endl;
      if (xphysRequested / hy != xpixels)
         std::cerr << "For square grid cells, xphys modified from " << xphysRequested << " to " << xphys << std::endl;
   }

   // Unfinished deserializing constructor.
   Fparams(std::istream& is) {
      is.read(reinterpret_cast<char *>(this), sizeof(*this));
   }

   // Return row or column numbers corresponding to given physical
   // coordinates.
   int yToRow(double y) const { return static_cast<int>(std::floor(y / hy)); }
   int xToCol(double x) const { return static_cast<int>(std::floor(x / hx)); }

   // These two return indices to grid cell boundaries rather than cells,
   // so their range is one more than for yToRow and xToCol.
   int yToRowBound(double y) const { return static_cast<int>(std::round(y / hy)); }
   int xToColBound(double x) const { return static_cast<int>(std::round(x / hx)); }

   // Return point representing lower-indexed corner of given cell
   Point boundToPt(int r, int c) const { return { c * hx, r * hy }; }

   // Return point representing center of given cell
   Point coordsToPt(int r, int c) const { return { (c + 0.5) * hx, (r + 0.5) * hy }; }

   // Return point representing center of given edge (indexed as in
   // EdgeGrid class).
   Point horizEdgeToPt(int r, int c) const { return { (c + 0.5) * hx, (r + 1) * hy }; }
   Point vertEdgeToPt(int r, int c) const { return { (c + 1) * hx, (r + 0.5) * hy }; }

   // Return area and diagonal length of a grid cell.
   double cellArea() const { return hx * hy; }
   double cellDiag() const { return std::hypot(hx, hy); }

   // Minimum grid cell dimension, in case rectangular pixels are allowed
   // in future.
   double minH() const { return std::min(hx, hy); }

   // Unfinished serializer
   void serialize(std::ostream& os) const {
      os.write(reinterpret_cast<const char *>(this), sizeof(*this));
   }
};

#endif
