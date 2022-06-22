// fluidvel.h
// Allen McBride
// April 14, 2022
//
// The FluidVel class encapsulates code for computing the velocity of the
// fluid medium due to pressure changes imposed by moving agents. The
// algorithm solves for unsteady Stokes flow with no-slip boundary
// conditions at walls and agent bodies. See Appendix A of McBride 2022
// for details. The Eigen library is used to solve matrix equations.

#ifndef FLUIDVEL_H
#define FLUIDVEL_H

#include "fields.h"
#include "fparams.h"
#include "grid.h"

#include <Eigen/Sparse>
#include <iomanip>
#include <ostream>
#include <vector>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::VectorXd EVec;

class FluidVel {
   private:

      // The Dim enum class is for indicating whether a given function
      // should operate on the x- or y-component of velocity.
      enum class Dim { x, y };

      // The Component subclass stores a matrix representing the x-
      // or y-component of the fluid velocity field.
      class Component {

         private:

            // Data members
            const Fparams fp;
            const double viscosity;
            const Dim dim;
            const int xsize;
            const int ysize;
            EVec data;
            const EdgeGrid<EdgeType>& edgeFreepart;
            const EdgeGrid<EdgeType>& edgeBodyVelProp;
            SpMat viscMatrix;
            Eigen::BiCGSTAB<SpMat> viscSolver;
            double prevDt = -1.0;

            // In-line utility methods to hide details of indexing
            int index(int row, int col) const { return col * ysize + row; }
            int lowCellRow(int row) const { return row; }
            int lowCellCol(int col) const { return col; }
            int highCellRow(int row) const { return dim == Dim::x ? row : row + 1; }
            int highCellCol(int col) const { return dim == Dim::y ? col : col + 1; }
            double pointComponent(const Point& p) const { return dim == Dim::x ? p.x : p.y; }

            // Set up matrix to solve for viscosity
            void setViscosityCoeff(double dt);


         public:

            // Provide size of grid
            int provideXsize() const { return xsize; };
            int provideYsize() const { return ysize; };
            
            // Provide component of fluid velocity at a given grid cell
            double at(int row, int col) const { return data(index(row, col)); }

            // Construct component grid, copying in key data (using
            // references for externally-owned grid-based data)
            Component(Fparams fp, const EdgeGrid<EdgeType>& edgeFreepart, const EdgeGrid<EdgeType>& edgeBodyVelProp, double viscosity, Dim dim, double viscosityTolerance);

            // Perform one timestep of momentum diffusion (that is,
            // solve for effect of viscosity)
            void diffuseMomentum(double dt);

            // After pressure is computed, update velocity accordingly
            void applyPressure(const EVec& pressureVec);
            
            // At agent bodies, force velocity field to equal agent body
            // velocity, interpolating for grid edges partially inside
            // agent body
            void reimposeBodyVel();

            // Print data to given stream, unformatted
            void printRaw(std::ostream& s) const;

            // Print data to given stream in grid format
            friend std::ostream& operator<<(std::ostream& s, const Component& a) {
               s << "ysize: " << a.ysize << "; xsize: " << a.xsize << std::endl;
               for (int r = 0; r < a.ysize; ++r) {
                  for (int c = 0; c < a.xsize; ++c)
                     s << std::setprecision(3) << a.at(r, c) << " ";
                  s << std::endl;
               }
               return s;
            }
      };

      // Data members
      const Fields& fields;
      Component xcomponents;
      Component ycomponents;
      // pressureSolver keeps a reference to laplMatrix
      SpMat laplMatrix;
      Eigen::BiCGSTAB<SpMat> pressureSolver;

      // Encapsulate indexing
      static int index(int row, int col, const Fparams& fparams) { return col * fparams.ypixels + row; }
      int index(int row, int col) const { return index(row, col, fields.fp); }

      // Call reimposeBodyVel for each component
      void reimposeBodyVel() { xcomponents.reimposeBodyVel(); ycomponents.reimposeBodyVel(); }

      // Return matrix of negative divergence of velocity, for use in
      // Poisson equation for pressure
      EVec negativeDivergence() const;

   public:

      // Constructor sets up matrix to use in solving Poisson equation
      // for pressure
      FluidVel(const Fields& fields, double viscosity, double viscosityTolerance, double pressureTolerance);

      // Advance fluid velocity field one timestep. Agent body information
      // comes from stored references in the Component subclass.
      void update(double dt);

      // Inline methods to encapsulate indexing logic
      double highX(int row, int col) const { return col == fields.fp.xpixels - 1 ? 0.0 : xcomponents.at(row, col); }
      double lowX(int row, int col) const { return col == 0 ? 0.0 : xcomponents.at(row, col - 1); }
      double highY(int row, int col) const { return row == fields.fp.ypixels - 1 ? 0.0 : ycomponents.at(row, col); }
      double lowY(int row, int col) const { return row == 0 ? 0.0 : ycomponents.at(row - 1, col); }

      // Print velocity field, one component at a time, unformatted
      void printVel(std::ostream& s) const { xcomponents.printRaw(s); ycomponents.printRaw(s); }
};

#endif
