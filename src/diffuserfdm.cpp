// diffuserfdm.cpp
// Allen McBride
// April 20, 2022
//
// Details in diffuserfdm.h

#include "diffuserfdm.h"
#include "fields.h"
#include "substances.h"

#include <algorithm>
#include <cmath>
#include <iostream> // DEBUG
#include <limits>
#include <stdexcept>
#include <utility>


double DiffuserFDM::perCellTransfer(int rOut, int cOut, const Grid<double>& grid, double giveEachWay) const {
   if (!fields.inBounds(rOut, cOut) || fields.fullyObstructed(rOut, cOut)) 
      return 0.0;
   return grid.at(rOut, cOut) - giveEachWay;
}


double DiffuserFDM::laplacian(int r, int c, const Grid<double>& grid) const {
   const double giveEachWay = grid.at(r, c);
   double laplacianNumerator = 
      sideCoeff * 
      (perCellTransfer(r    , c + 1, grid, giveEachWay) +
       perCellTransfer(r - 1, c    , grid, giveEachWay) +
       perCellTransfer(r    , c - 1, grid, giveEachWay) +
       perCellTransfer(r + 1, c    , grid, giveEachWay)) +
      cornerCoeff *
      (perCellTransfer(r - 1, c + 1, grid, giveEachWay) +
       perCellTransfer(r - 1, c - 1, grid, giveEachWay) +
       perCellTransfer(r + 1, c - 1, grid, giveEachWay) +
       perCellTransfer(r + 1, c + 1, grid, giveEachWay));
   return laplacianNumerator / fields.fp.cellArea();
}


DiffuserFDM::DiffuserFDM(
      Fields& fields, 
      Substances subs,
      double targetLaplacian,
      double radius
      ) :
   fields{fields}, 
   buffer{fields.fp, 0.0}, 
   subs{subs},
   targetLaplacian{targetLaplacian},
   radius{radius}
{
   if (fields.fp.hx != fields.fp.hy)
      throw std::invalid_argument("hx and hy must be same for this diffuser class");
}


void DiffuserFDM::diffuse(double dt) { 
   if (dt > maxDt())
      throw std::runtime_error("dt too large for DiffuserFDM::diffuse");

   Grid<double>& grid = fields.grids.at(DiffStates::one);

#pragma omp parallel for 
   for (int r = 0; r < fields.fp.ypixels; ++r) {
      for (int c = 0; c < fields.fp.xpixels; ++c) {
         const Point center{fields.fp.xphys / 2.0, fields.fp.yphys / 2.0};
         const auto point = fields.fp.coordsToPt(r, c);
         if (radius >= 0.0 && (point - center).mag() > radius) {
            buffer.at(r, c) = 0.0;
         } else {
            const auto l = laplacian(r, c, grid);
            //if (l > 0.000001) std::cerr << l << std::endl;
            buffer.at(r, c) = grid.at(r, c) 
               + dt * subs.at(DiffStates::one).difRate * (laplacian(r, c, grid) - targetLaplacian);
            //std::cerr << l << " -- " << dt * subs.at(DiffStates::one).difRate << " -- " << targetLaplacian << " -- " << buffer.at(r, c) - grid.at(r, c) << std::endl;
         }

         if (buffer.at(r, c) < 0.0) {
            if (buffer.at(r, c) < -std::numeric_limits<double>::epsilon() * grid.at(r, c)) {
               std::cerr << "negative value: " << buffer.at(r, c) << std::endl;
               std::cerr << grid.at(r - 1, c - 1) << "      " << grid.at(r - 1, c) << "      " << grid.at(r - 1, c + 1) << std::endl;
               std::cerr << grid.at(r, c - 1) << "      " << grid.at(r, c) << "      " << grid.at(r, c + 1) << std::endl;
               std::cerr << grid.at(r + 1, c - 1) << "      " << grid.at(r + 1, c) << "      " << grid.at(r + 1, c + 1) << std::endl;
               std::cerr << subs.at(DiffStates::one).difRate << std::endl;
               throw std::logic_error( // DEBUG
                     "Negative substance value in DiffuserFDM::diffuse. isub: " + std::to_string(DiffStates::one)
                     + "; r,c: " + std::to_string(r) + "," + std::to_string(c)
                     + "; val: " + std::to_string(buffer.at(r, c)));
            } else {
               buffer.at(r, c) = 0.0;
            }
         }
      }
   }

   grid.swap(buffer);
}


double DiffuserFDM::maxDt(const Substances& subs, double h) {
   double maxDiffusionRate = 0.0;
   for (const auto& s : subs) maxDiffusionRate = std::max(maxDiffusionRate, s.difRate);
   return h * h / (centralCoeff * maxDiffusionRate);
}
