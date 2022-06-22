// diffuserfdm.h
// Allen McBride
// April 20, 2022
//
// The DiffuserFDM class is mostly a copy of the diffusing part of
// the Diffuser class, modified to solve a Poisson equation instead of
// Laplace's equation, for use in computing ideal fields for comparison
// with MSF-SPH simulations.

#ifndef DIFFUSERFDM_H
#define DIFFUSERFDM_H

#include "fields.h"
#include "substances.h"

#include <ostream>

class DiffuserFDM {

   private:
      Fields& fields;
      Grid<double> buffer;
      const Substances subs;
      const double targetLaplacian;
      const double radius;

      double perCellTransfer(int rOut, int cOut, const Grid<double>& grid, double giveEachWay) const;
      double laplacian(int r, int c, const Grid<double>& grid) const;
   public:
      static constexpr int difWays = 8;
      static constexpr double cornerCoeff = 1.0 / 6.0;
      static constexpr double sideCoeff = 2.0 / 3.0;
      static constexpr double centralCoeff = difWays == 4 ? 4.0 : cornerCoeff * 4 + sideCoeff * 4;
      static double maxDt(const Substances& subs, double h);

      DiffuserFDM(
            Fields& fields, 
            Substances subs,
            double targetLaplacian,
            double radius
            ); 
      void diffuse(double dt);
      double maxDt() const { return maxDt(subs, fields.fp.minH()); }
};

#endif
