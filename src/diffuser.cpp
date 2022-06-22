// diffuser.cpp
// Allen McBride
// March 10, 2022
//
// Details in diffuser.h

#include "diffuser.h"
#include "fields.h"
#include "fluidvel.h"
#include "substances.h"

#include <algorithm>
#include <cmath>
#include <iostream> // DEBUG
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>


inline double decayFactor(double decayRate, double dt);


// When out-cell is partly or fully out of bounds, enforce homogeneous Neumann boundary conditions.
double Diffuser::perCellTransferNopass(int rOut, int cOut, const Grid<double>& grid, double giveEachWay, double fpIn) const {
   if (!fields.inBounds(rOut, cOut) || fields.fullyObstructed(rOut, cOut)) 
      return 0.0;
   const double fpOut = 1.0 - fields.propObstructed(rOut, cOut);
   return (grid.at(rOut, cOut) - giveEachWay) * std::min(fpIn, fpOut);
}


double Diffuser::perCellTransferPassObs(int rOut, int cOut, const Grid<double>& grid, double giveEachWay) const {
   if (!fields.inBounds(rOut, cOut) || fields.fullyObstructed(rOut, cOut)) 
      return 0.0;
   return grid.at(rOut, cOut) - giveEachWay;
}


double Diffuser::perCellTransferPassNoObs(int rOut, int cOut, const Grid<double>& grid, double giveEachWay) const {
   if (!fields.inBounds(rOut, cOut)) 
      return 0.0;
   return grid.at(rOut, cOut) - giveEachWay;
}


void Diffuser::diffuseOneSubNopass(int isub, double dt) {

   Grid<double>& grid = fields.grids.at(isub);
   const double gridDifRate = subs.at(isub).difRate / fields.fp.cellArea();

#pragma omp parallel for 
   for (int r = 0; r < fields.fp.ypixels; ++r) {
      for (int c = 0; c < fields.fp.xpixels; ++c) {
         if (!fields.fullyObstructed(r, c)) {
            const double giveEachWay = grid.at(r, c);
            const double fpIn = 1.0 - fields.propObstructed(r, c);
            buffer.at(r, c) = grid.at(r, c) + dt * gridDifRate *
               (sideCoeff * 
                (perCellTransferNopass(r    , c + 1, grid, giveEachWay, fpIn) +
                 perCellTransferNopass(r - 1, c    , grid, giveEachWay, fpIn) +
                 perCellTransferNopass(r    , c - 1, grid, giveEachWay, fpIn) +
                 perCellTransferNopass(r + 1, c    , grid, giveEachWay, fpIn)) +
                cornerCoeff *
                (perCellTransferNopass(r - 1, c + 1, grid, giveEachWay, fpIn) +
                 perCellTransferNopass(r - 1, c - 1, grid, giveEachWay, fpIn) +
                 perCellTransferNopass(r + 1, c - 1, grid, giveEachWay, fpIn) +
                 perCellTransferNopass(r + 1, c + 1, grid, giveEachWay, fpIn)));
         } else {
            buffer.at(r, c) = grid.at(r, c);
         }
         buffer.at(r, c) *= decayFactor(subs.at(isub).decayRate, dt);
      }
   }

   grid.swap(buffer);
}


void Diffuser::diffuseOneSubPassObs(int isub, double dt) {

   Grid<double>& grid = fields.grids.at(isub);
   const double gridDifRate = subs.at(isub).difRate / fields.fp.cellArea();

#pragma omp parallel for 
   for (int r = 0; r < fields.fp.ypixels; ++r) {
      for (int c = 0; c < fields.fp.xpixels; ++c) {
         if (!fields.fullyObstructed(r, c)) {
            const double giveEachWay = grid.at(r, c);
            buffer.at(r, c) = grid.at(r, c) + dt * gridDifRate *
               (sideCoeff * 
                (perCellTransferPassObs(r    , c + 1, grid, giveEachWay) +
                 perCellTransferPassObs(r - 1, c    , grid, giveEachWay) +
                 perCellTransferPassObs(r    , c - 1, grid, giveEachWay) +
                 perCellTransferPassObs(r + 1, c    , grid, giveEachWay)) +
                cornerCoeff *
                (perCellTransferPassObs(r - 1, c + 1, grid, giveEachWay) +
                 perCellTransferPassObs(r - 1, c - 1, grid, giveEachWay) +
                 perCellTransferPassObs(r + 1, c - 1, grid, giveEachWay) +
                 perCellTransferPassObs(r + 1, c + 1, grid, giveEachWay)));
         } else {
            buffer.at(r, c) = grid.at(r, c);
         }
         buffer.at(r, c) *= decayFactor(subs.at(isub).decayRate, dt);
      }
   }

   grid.swap(buffer);
}


void Diffuser::diffuseOneSubPassNoObs(int isub, double dt) {

   Grid<double>& grid = fields.grids.at(isub);
   const double gridDifRate = subs.at(isub).difRate / fields.fp.cellArea();

#pragma omp parallel for 
   for (int r = 0; r < fields.fp.ypixels; ++r) {
      for (int c = 0; c < fields.fp.xpixels; ++c) {
         const double giveEachWay = grid.at(r, c);
         buffer.at(r, c) = grid.at(r, c) + dt * gridDifRate *
            (sideCoeff * 
             (perCellTransferPassNoObs(r    , c + 1, grid, giveEachWay) +
              perCellTransferPassNoObs(r - 1, c    , grid, giveEachWay) +
              perCellTransferPassNoObs(r    , c - 1, grid, giveEachWay) +
              perCellTransferPassNoObs(r + 1, c    , grid, giveEachWay)) +
             cornerCoeff *
             (perCellTransferPassNoObs(r - 1, c + 1, grid, giveEachWay) +
              perCellTransferPassNoObs(r - 1, c - 1, grid, giveEachWay) +
              perCellTransferPassNoObs(r + 1, c - 1, grid, giveEachWay) +
              perCellTransferPassNoObs(r + 1, c + 1, grid, giveEachWay)));
         //if (buffer.at(r, c) < 0.0) {
         //   std::cerr.precision(std::numeric_limits<double>::max_digits10);
         //   for (int i = -2; i <= 2; ++i) {
         //      for (int j = -2; j <= 2; ++j) {
         //         std::cerr << grid.at(r, c) << " , ";
         //      }
         //      std::cerr << std::endl;
         //   }
         //   std::cerr << std::endl;
         //   std::cerr << buffer.at(r, c) << std::endl;
         //   exit(1);
         //}
         buffer.at(r, c) *= decayFactor(subs.at(isub).decayRate, dt);
      }
   }

   grid.swap(buffer);
}


void Diffuser::diffuseOneSubcycle(int isub, double dt) {

   if (!passthrough)
      diffuseOneSubNopass(isub, dt);
   else if (fields.obstacles)
      diffuseOneSubPassObs(isub, dt);
   else
      diffuseOneSubPassNoObs(isub, dt);
}


void Diffuser::transferToHigh(const Grid<double>& grid, int rowLow, int colLow, int rowHigh, int colHigh, double component, double h, double dt) {
   const bool flowTowardHigh = component >= 0.0;
   const int rowUpwind = flowTowardHigh ? rowLow : rowHigh;
   const int colUpwind = flowTowardHigh ? colLow : colHigh;
   const double transfer = dt * component * grid.at(rowUpwind, colUpwind) / h;
   buffer.at(rowLow, colLow) -= transfer;
   buffer.at(rowHigh, colHigh) += transfer;
}


Diffuser::Diffuser(
      Fields& fields, 
      Substances subs,
      bool passthrough,
      double viscosity,
      double viscosityTolerance,
      double pressureTolerance
      ) :
   fields{fields}, 
   buffer{fields.fp, 0.0}, 
   subs{subs}, 
   passthrough{passthrough}, 
   fluidVel{passthrough ? std::nullopt : std::make_optional<FluidVel>(fields, viscosity, viscosityTolerance, pressureTolerance)}
{
   if (fields.fp.hx != fields.fp.hy)
      throw std::invalid_argument("hx and hy must be same for this diffuser class");
}


Diffuser::Diffuser(Fields& fields, std::istream& is) :
   fields{fields},
   buffer{fields.fp, 0.0}, 
   subs{readObjectFromStream<Substances>(is)},
   passthrough{static_cast<bool>(is.get())},
   fluidVel{std::nullopt} // TODO finish fluidVel serialization
{
}


void Diffuser::diffuse(double dt) { 
   for (int isub = 0; isub < nSubs; ++isub) 
      if (subs.at(isub).include)
         diffuseOne(isub, dt); 
}


void Diffuser::diffuseOne(int isub, double dt) { 
   if (subs.at(isub).difRate > 0.0) {
      const auto maxdt = fields.fp.cellArea() / (centralCoeff * subs.at(isub).difRate);
      const auto subcycles = static_cast<int>(std::ceil(dt / maxdt));
      const auto subdt = dt / subcycles;
      if (subcycles > 1) std::cerr << "  dt: " << dt << "; maxdt: " << maxdt << "; diffuseOne subcycles: " << subcycles << std::endl; // DEBUG
      for (int icycle = 0; icycle < subcycles; ++icycle) 
         diffuseOneSubcycle(isub, subdt);
   }
}


void Diffuser::advect(double dt) {
   //for (int row = 0; row < fields.fp.ypixels - 1; ++row)
   //   for (int col = 0; col < fields.fp.xpixels - 1; ++col)
   //      if (fields.grids.at(1).at(row, col) > 0.0) { std::cerr << "adv: r: " << row << "; c: " << col << "; val: " << fields.grids.at(1).at(row, col) << std::endl; exit(1); } // DEBUG

   if (!passthrough) {

      updateFlow(dt);
      
      for (int isub = 0; isub < nSubs; ++isub) {

         if (subs.at(isub).include) {
            // The grid must be copied because we don't know in advance which direction the upwinding happens across the various boundaries
            buffer.copy(fields.grids.at(isub));
#pragma omp parallel for 
            for (int row = 0; row < fields.fp.ypixels - 1; ++row) {
               for (int col = 0; col < fields.fp.xpixels - 1; ++col) {
                  transferToHigh(fields.grids.at(isub), row, col, row, col + 1, 
                        fluidVel->highX(row, col), fields.fp.hx, dt);
                  transferToHigh(fields.grids.at(isub), row, col, row + 1, col, 
                        fluidVel->highY(row, col), fields.fp.hy, dt);
               }
            }

            for (int row = 0; row < fields.fp.ypixels; ++row)
               for (int col = 0; col < fields.fp.xpixels; ++col)
                  if (buffer.at(row, col) < 0.0)
                     buffer.at(row, col) = 0.0;

            fields.grids.at(isub).swap(buffer);
         }
      }
   }
}


double Diffuser::maxDt(const Substances& subs, double h) {
   double maxDiffusionRate = 0.0;
   for (const auto& s : subs) maxDiffusionRate = std::max(maxDiffusionRate, s.difRate);
   const auto ideal = h * h / (centralCoeff * maxDiffusionRate);
   return ideal * (1.0 - 8.0 * std::numeric_limits<double>::epsilon());
}


inline double decayFactor(const double decayRate, const double dt) {
   //return std::exp(-dt * decayRate);
   return 1.0 - dt * decayRate;
}

