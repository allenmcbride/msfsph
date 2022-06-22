// diffuser.h
// Allen McBride
// March 10, 2022
//
// The Diffuser class handles diffusion and degradation of morphogens.
// When the "swimming" model is used, a FluidVel object, contained by
// Diffuser, computes a fluid velocity field.  This field is then used
// by Diffuser to compute advection of morphogen.
//
// A forward Euler method with a nine-point stencil is used to compute
// diffusion.  For details, see Appendix A.2, "Modeling Diffusion".


#ifndef DIFFUSER_H
#define DIFFUSER_H

#include "fields.h"
#include "fluidvel.h"
#include "substances.h"

#include <optional>
#include <ostream>


class Diffuser {

   private:
      Fields& fields; // Reference to the Fields class passed to the constructor
      Grid<double> buffer; // Buffer copy of grid, to make updates synchronous [TODO test if slowdown if this is put on stack]
      const Substances subs; // Copy of Substances object passed to constructor
      const bool passthrough; // Copy of flag passed to constructor; true for crawling, false for swimming
      std::optional<FluidVel> fluidVel; // Optional object to compute fluid velocity for swimming model

      // Compute amount to transfer to a neighboring cell. Different
      // versions are for efficiency; these are called in an inner loop
      // so every conditional matters. Nopass is for swimming; PassObs
      // is for crawling when obstacles are present; PassNoObs is for
      // crawling when no obstacles are present
      double perCellTransferNopass(int rOut, int cOut, const Grid<double>& grid, double giveEachWay, double fpIn) const;
      double perCellTransferPassObs(int rOut, int cOut, const Grid<double>& grid, double giveEachWay) const;
      double perCellTransferPassNoObs(int rOut, int cOut, const Grid<double>& grid, double giveEachWay) const;

      // Loop through cells, computing diffusion and degradation for a
      // single morphogen. Three versions as above for perCellTransfer.
      void diffuseOneSubNopass(int isub, double dt);
      void diffuseOneSubPassObs(int isub, double dt);
      void diffuseOneSubPassNoObs(int isub, double dt);

      // Call the appropriate diffuseOneSub version above.
      void diffuseOneSubcycle(int isub, double dt);

      // For a given pair of cells, advect morphogen between them in an
      // upwind fashion by applying the appropriate velocity component
      // at their border.
      void transferToHigh(const Grid<double>& grid, int rowLow, int colLow, int rowHigh, int colHigh, double component, double h, double dt);

   public:
      // Corner, side, and central values for the nine-point stencil used.
      static constexpr double cornerCoeff = 1.0 / 6.0;
      static constexpr double sideCoeff = 2.0 / 3.0;
      static constexpr double centralCoeff = cornerCoeff * 4 + sideCoeff * 4;

      // This class seems like a logical place to define this quantity; 6
      // pixels across seems like a reasonable lower limit for resolution
      // relative to agent bodies for swimming. Crawling simulations do
      // not need as high resolution because there we are not trying to
      // resolve the effect of agent bodies on fluid velocity.
      static constexpr double minimumRecommendedGridCellWidthsPerMinorAxisSwimming = 6.0; 

      // Compute maximum timestep allowable, across all substances,
      // based on diffusion stability criterion.
      static double maxDt(const Substances& subs, double h);

      // Normal constructor.
      Diffuser(
            Fields& fields, 
            Substances subs, 
            bool passthrough, 
            double viscosity, 
            double viscosityTolerance, 
            double pressureTolerance
            ); 

      // Deserializing constructor (unfinished).
      Diffuser(
            Fields& fields, 
            std::istream& is
            );

      // Perform diffusion for all substances.
      void diffuse(double dt);

      // Perform diffusion for one substance. If needed for diffusion
      // stability criterion, call appropriate diffuseOneSub method
      // several times for subdivided timesteps.
      void diffuseOne(int isub, double dt);

      // Get updated fluid velocity field from FluidVel object and apply
      // it to advect morphogens.
      void advect(double dt);

      // Public version of maxDt; no parameters required.
      double maxDt() const { return maxDt(subs, fields.fp.minH()); }

      // Ask FluidVel to update fluid velocity field.
      void updateFlow(double dt) { if (!passthrough) fluidVel->update(dt); }

      // Print fluid velocity field to given stream, for debugging.
      void printVel(std::ostream& s) const { if (!passthrough) fluidVel->printVel(s); }
};

#endif
