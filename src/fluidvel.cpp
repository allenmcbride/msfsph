// fluidvel.cpp
// Allen McBride
// April 19, 2022
//
// Descriptions in fluidvel.h

#include "fluidvel.h"
#include "fields.h"
#include "fparams.h"
#include "grid.h"

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <vector>


// Build matrix for backward Euler momentum diffusion
// Use ghost points to deal with staggered grid.
void FluidVel::Component::setViscosityCoeff(double dt) {
   std::vector<Eigen::Triplet<double> > coeffs;
   // Column-major order is used in terms of how grid indices are converted to matrix indices,
   // which seems standard in CFD.
   // But this is independent of the storage format of the Laplacian matrix.
   for (int col = 0; col < xsize; ++col) {
      for (int row = 0; row < ysize; ++row) {
         const int coeffRow = index(row, col);
         const double valx = viscosity * dt / (fp.hx * fp.hx);
         const double valy = viscosity * dt / (fp.hy * fp.hy);
         coeffs.emplace_back(coeffRow, coeffRow, 1.0 + 2.0 * (valx + valy));

         if (col == xsize - 1) {
            if (dim == FluidVel::Dim::y) 
               coeffs.emplace_back(coeffRow, coeffRow, valx);
         } else {
            coeffs.emplace_back(coeffRow, index(row, col + 1), -valx);
         }
         if (col == 0) {
            if (dim == FluidVel::Dim::y) 
               coeffs.emplace_back(coeffRow, coeffRow, valx);
         } else {
            coeffs.emplace_back(coeffRow, index(row, col - 1), -valx);
         }
         if (row == ysize - 1) {
            if (dim == FluidVel::Dim::x) 
               coeffs.emplace_back(coeffRow, coeffRow, valy);
         } else {
            coeffs.emplace_back(coeffRow, index(row + 1, col), -valy);
         }
         if (row == 0) {
            if (dim == FluidVel::Dim::x) 
               coeffs.emplace_back(coeffRow, coeffRow, valy);
         } else {
            coeffs.emplace_back(coeffRow, index(row - 1, col), -valy);
         }
      }
   }
   viscMatrix.setFromTriplets(coeffs.begin(), coeffs.end());
}


FluidVel::Component::Component(
      Fparams fp, 
      const EdgeGrid<EdgeType>& edgeFreepart, 
      const EdgeGrid<EdgeType>& edgeBodyVelProp, 
      double viscosity, 
      FluidVel::Dim dim, 
      double viscosityTolerance
      ) : 
   fp{fp}, 
   viscosity{viscosity},
   dim{dim},
   xsize{dim == FluidVel::Dim::x ? fp.xpixels - 1 : fp.xpixels},
   ysize{dim == FluidVel::Dim::y ? fp.ypixels - 1 : fp.ypixels},
   data{EVec::Zero(xsize * ysize)},
   edgeFreepart{edgeFreepart},
   edgeBodyVelProp{edgeBodyVelProp},
   viscMatrix{xsize * ysize, xsize * ysize}
{
   viscSolver.setTolerance(viscosityTolerance); // A loose tolerance increases speed dramatically but only with large grid sizes
   setViscosityCoeff(1.0);
   viscSolver.analyzePattern(viscMatrix);
   if (viscSolver.info() != Eigen::Success) throw std::runtime_error("viscSolver initialization failed");
}


void FluidVel::Component::diffuseMomentum(double dt) {
   if (dt != prevDt) {
      setViscosityCoeff(dt);
      viscSolver.factorize(viscMatrix);
      if (viscSolver.info() != Eigen::Success) throw std::runtime_error("viscSolver factorize failed");
      prevDt = dt;
   }
   data = viscSolver.solve(data).eval();
   if (viscSolver.info() != Eigen::Success) throw std::runtime_error("viscSolver solve failed");
}


void FluidVel::Component::applyPressure(const EVec& pressureVec) {
#pragma omp parallel for
   for (int col = 0; col < xsize; ++col) {
      for (int row = 0; row < ysize; ++row) {
         const double h = dim == FluidVel::Dim::x ? fp.hx : fp.hy;
         const double pressureLow = pressureVec(FluidVel::index(lowCellRow(row), lowCellCol(col), fp));
         const double pressureHigh = pressureVec(FluidVel::index(highCellRow(row), highCellCol(col), fp));
         const double pressureGradComponent = (pressureHigh - pressureLow) / h;

         const int i = index(row, col);
         data(i) -= pressureGradComponent;
      }
   }
}


void FluidVel::Component::reimposeBodyVel() { 
#pragma omp parallel for
   for (int col = 0; col < xsize; ++col)
      for (int row = 0; row < ysize; ++row)
         data(index(row, col)) = data(index(row, col)) * edgeFreepart.at(row, col) + edgeBodyVelProp.at(row, col);
}


void FluidVel::Component::printRaw(std::ostream& s) const {
   for (int row = 0; row < ysize; ++row) {
      for (int col = 0; col < xsize; ++col) {
         s << data(index(row, col)) << std::endl;
      }
   }
}


// Build RHS of Poisson equation for pressure based on divergence of
// intermediate velocity.  As with Laplacian matrix, column-major order
// is used to convert grid indices to vector indices.
EVec FluidVel::negativeDivergence() const {
   EVec negDivergence(fields.fp.xpixels * fields.fp.ypixels);
#pragma omp parallel for
   for (int col = 0; col < fields.fp.xpixels; ++col) {
      for (int row = 0; row < fields.fp.ypixels; ++row) {
         const double divergenceX = (highX(row, col) - lowX(row, col)) / fields.fp.hx;
         const double divergenceY = (highY(row, col) - lowY(row, col)) / fields.fp.hy;
         const double divergence = divergenceX + divergenceY;
         negDivergence(index(row, col)) = -divergence;
      }
   }
   return negDivergence;
}


FluidVel::FluidVel(
      const Fields& fields, 
      double viscosity, 
      double viscosityTolerance, 
      double pressureTolerance
      ) : 
   fields{fields},
   xcomponents{fields.fp, fields.edgeFreepartVertical.value(), fields.edgeBodyVelPropVertical.value(), viscosity, FluidVel::Dim::x, viscosityTolerance},
   ycomponents{fields.fp, fields.edgeFreepartHorizontal.value(), fields.edgeBodyVelPropHorizontal.value(), viscosity, FluidVel::Dim::y, viscosityTolerance},
   laplMatrix{fields.fp.xpixels * fields.fp.ypixels, fields.fp.xpixels * fields.fp.ypixels}
{
   std::vector<Eigen::Triplet<double> > laplCoeffs;
   // Column-major order is used in terms of how grid indices are converted to matrix indices,
   // which seems standard in CFD.
   // But this is independent of the storage format of the Laplacian matrix.
   for (int col = 0; col < fields.fp.xpixels; ++col) {
      for (int row = 0; row < fields.fp.ypixels; ++row) {
         const int coeffRow = index(row, col);
         const double valx = -1.0 / (fields.fp.hx * fields.fp.hx);
         const double valy = -1.0 / (fields.fp.hy * fields.fp.hy);

         laplCoeffs.emplace_back(coeffRow, coeffRow, -2.0 * (valx + valy));
         if (col == fields.fp.xpixels - 1)
            laplCoeffs.emplace_back(coeffRow, coeffRow, valx);
         else
            laplCoeffs.emplace_back(coeffRow, index(row, col + 1), valx);
         if (col == 0)
            laplCoeffs.emplace_back(coeffRow, coeffRow, valx);
         else
            laplCoeffs.emplace_back(coeffRow, index(row, col - 1), valx);
         if (row == fields.fp.ypixels - 1)
            laplCoeffs.emplace_back(coeffRow, coeffRow, valy);
         else
            laplCoeffs.emplace_back(coeffRow, index(row + 1, col), valy);
         if (row == 0)
            laplCoeffs.emplace_back(coeffRow, coeffRow, valy);
         else
            laplCoeffs.emplace_back(coeffRow, index(row - 1, col), valy);
      }
   }
   laplMatrix.setFromTriplets(laplCoeffs.begin(), laplCoeffs.end());
   pressureSolver.compute(laplMatrix);
   pressureSolver.setTolerance(pressureTolerance); // A loose tolerance increases speed dramatically but only with large grid sizes
   if (pressureSolver.info() != Eigen::Success) throw std::runtime_error("pressureSolver compute failed");
}


// Based on Kajishima et al., 2001
//void FluidVel::update(const Grid<Point>& bodyVelProp, const Grid<double>& freepart, double dt) {
void FluidVel::update(double dt) {
   reimposeBodyVel();

   // Adjust velocity based on viscosity, ignoring advection of momentum as consistent with Stokes flow
   xcomponents.diffuseMomentum(dt);
   ycomponents.diffuseMomentum(dt);
   
   // reimpose agent bodies again for greater accuracy
   reimposeBodyVel();
   
   // If we needed physical pressure we would multiply by density/dt in the RHS.
   const EVec& pressureVec = pressureSolver.solve(negativeDivergence());
   if (pressureSolver.info() != Eigen::Success) throw std::runtime_error("pressureSolver solve failed");

   // Update velocity based on gradient of newly computed pressure, and interpolate between
   // this assuming-free-fluid velocity and body velocity based on cell occupancies.
   // If we were using physical pressure we would divide by density/dt in these updates.
   xcomponents.applyPressure(pressureVec);
   ycomponents.applyPressure(pressureVec);

   //std::cerr << xcomponents << std::endl;
   //std::cerr << ycomponents << std::endl;
}
