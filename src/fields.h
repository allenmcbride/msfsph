// fields.h
// Allen McBride
// March 23, 2022
//
// The Fields class holds most grid-based information for a simulation.
// This includes morphogen fields, agent bodies, agent body velocies,
// and obstacles, among others.  The exception is fluid velocity, which
// is stored in the FluidVel class.
//
// Fields also holds the direction of gravity, because it can be thought
// of as a constant field.
//
// Fields also encapsulates some key grid operations, such as checking
// whether a location is in bounds or a grid cell is obstructed, or
// injecting a given amount of morphogen by mass.
//
// This class started as a struct and accreted; at this point, it ought
// to be split into header and implementation files.

#ifndef FIELDS_H
#define FIELDS_H

#include "fparams.h"
#include "grid.h"
#include "point.h"
#include "serialization.h"
#include "substances.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>
#include <optional>
#include <string>
#include <vector>


typedef double AgentBodyType;
typedef int InstrumentType;
typedef double FrontType;
typedef char ObstacleType;
typedef double EdgeType;
typedef double ConcType;


// Exception to throw if fields would, if created, take up more RAM than
// we want to allow.
class FieldsTooBig : public std::runtime_error {
   public: using runtime_error::runtime_error;
};


class Fields {
   private:
      static void checkSize(Fparams fp, bool passthrough, bool obstaclesFlag, bool cuesFlag, bool visual) {
         const long long bytesPerPix
            = sizeof(AgentBodyType)
            + nSubs * sizeof(ConcType) // for grids
            + (obstaclesFlag ? sizeof(ObstacleType) : 0)
            + (cuesFlag ? sizeof(CueType) : 0)
            + (visual ? sizeof(FrontType) : 0)
            + (passthrough ? 0 : 4 * sizeof(EdgeType)); // for the four edge grids
         const long long npix = fp.xpixels * fp.ypixels;
         if (npix * bytesPerPix > 7'000'000'000ll) {
            std::cerr << "Fields constructor called with " << fp.xpixels << " by " << fp.ypixels << " pixels, and roughly " << bytesPerPix << " bytes per pixel." << std::endl;
            throw FieldsTooBig("Fields would have been too big.");
         }
      }

   public:
      const Fparams fp;
      const bool passthrough;

      Point gravityDir; // This must be a unit vector
      Grid<AgentBodyType> agentBodies;
      std::vector<Grid<ConcType> > grids;
      
      // obstacles should be Grid<bool>, but we need to avoid the
      // std::vector<bool> specialization because it doesn't allow
      // returning references using [].
      std::optional<Grid<ObstacleType> > obstacles; 
      std::optional<Grid<CueType> > cues;
      std::optional<Grid<FrontType> > front;
      std::optional<Grid<InstrumentType> > instruments;
      std::optional<EdgeGrid<EdgeType> > edgeFreepartHorizontal;
      std::optional<EdgeGrid<EdgeType> > edgeFreepartVertical;
      std::optional<EdgeGrid<EdgeType> > edgeBodyVelPropHorizontal;
      std::optional<EdgeGrid<EdgeType> > edgeBodyVelPropVertical;

      void assertObstacles() {
         if (!passthrough) {
            for (int r = 0; r < fp.ypixels; ++r) {
               for (int c = 0; c < fp.xpixels; ++c) {
                  if (obstacles && obstacles->at(r, c)) {
                     if (r > 0) edgeFreepartHorizontal->at(r - 1, c) = 0.0;
                     if (r < fp.ypixels - 1) edgeFreepartHorizontal->at(r, c) = 0.0;
                     if (c > 0) edgeFreepartVertical->at(r, c - 1) = 0.0;
                     if (c < fp.xpixels - 1) edgeFreepartVertical->at(r, c) = 0.0;
                  }
               }
            }
         }
      }

      Fields(Fparams fp, bool passthrough, bool obstaclesFlag = false, bool cuesFlag = false, bool visual = false) : 
         fp{(checkSize(fp, passthrough, obstaclesFlag, cuesFlag, visual), fp)}, // Use comma operator to check potential size before initializing fields
         passthrough{passthrough},
         gravityDir{Point{1.0, 0.0}.rotated(M_PI * 1.5)},
         agentBodies{fp, 0.0},
         grids{nSubs, Grid<ConcType>(fp, 0.0)},
         obstacles{obstaclesFlag ? std::make_optional<Grid<ObstacleType> >(fp, false) : std::nullopt},
         cues{cuesFlag ? std::make_optional<Grid<CueType> >(fp, CueType{}) : std::nullopt},
         front{visual ? std::make_optional<Grid<FrontType> >(fp, 0.0) : std::nullopt},
         instruments{passthrough ? std::nullopt : std::make_optional<Grid<InstrumentType> >(fp, 0)},
         edgeFreepartHorizontal{passthrough ? std::nullopt : std::make_optional<EdgeGrid<EdgeType> >(fp, 1.0, true)},
         edgeFreepartVertical{passthrough ? std::nullopt : std::make_optional<EdgeGrid<EdgeType> >(fp, 1.0, false)},
         edgeBodyVelPropHorizontal{passthrough ? std::nullopt : std::make_optional<EdgeGrid<EdgeType> >(fp, 0.0, true)},
         edgeBodyVelPropVertical{passthrough ? std::nullopt : std::make_optional<EdgeGrid<EdgeType> >(fp, 0.0, false)}
      { 
         assertObstacles();
      }

      Fields(std::istream& is) :
         fp{is},
         passthrough{static_cast<bool>(is.get())},
         gravityDir{readObjectFromStream<Point>(is)},
         agentBodies{is},
         obstacles{static_cast<bool>(is.get()) ? std::make_optional<Grid<ObstacleType> >(is) : std::nullopt},
         cues{static_cast<bool>(is.get()) ? std::make_optional<Grid<CueType> >(is) : std::nullopt},
         front{static_cast<bool>(is.get()) ? std::make_optional<Grid<FrontType> >(is) : std::nullopt},
         instruments{passthrough ? std::nullopt : std::make_optional<Grid<InstrumentType> >(is)},
         edgeFreepartHorizontal{passthrough ? std::nullopt : std::make_optional<EdgeGrid<EdgeType> >(is)},
         edgeFreepartVertical{passthrough ? std::nullopt : std::make_optional<EdgeGrid<EdgeType> >(is)},
         edgeBodyVelPropHorizontal{passthrough ? std::nullopt : std::make_optional<EdgeGrid<EdgeType> >(is)},
         edgeBodyVelPropVertical{passthrough ? std::nullopt : std::make_optional<EdgeGrid<EdgeType> >(is)}
      {}

      void clearAgentInfo() {
         agentBodies.fill(0.0);
         if (front) front->fill(0.0);

         if (!passthrough) {
            instruments->fill(0);
            edgeFreepartHorizontal->fill(1.0);
            edgeFreepartVertical->fill(1.0);
            assertObstacles();
            edgeBodyVelPropHorizontal->fill(0.0);
            edgeBodyVelPropVertical->fill(0.0);
         }
      }

#ifdef XWRAP
      bool inBoundsXLow(double x) const { return true; }
      bool inBoundsXHigh(double x) const { return true; }
      bool inBounds(int row, int col) const {
         return row >= 0
            && row < fp.ypixels;
      }
#else
      bool inBoundsXLow(double x) const { return x >= 0.0; }
      bool inBoundsXHigh(double x) const { return x < fp.xphys; }
      bool inBounds(int row, int col) const {
         return row >= 0
            && col >= 0
            && row < fp.ypixels
            && col < fp.xpixels;
      }
#endif

      bool inBoundsYLow(double y) const { return y >= 0.0; }
      bool inBoundsYHigh(double y) const { return y < fp.yphys; }

      bool inBounds(const Point& p) const {
         return inBoundsXLow(p.x)
            && inBoundsXHigh(p.x)
            && inBoundsYLow(p.y)
            && inBoundsYHigh(p.y);
      }

      double propOccupied(int row, int col) const {
         double propOccupied = agentBodies.at(row, col);
         if (obstacles) propOccupied += obstacles->at(row, col);
         return propOccupied;
      }
      double propOccupied(const Point& p) const {
         return propOccupied(fp.yToRow(p.y), fp.xToCol(p.x));
      }

      double propObstructed(int row, int col) const { 
         if (obstacles && obstacles->at(row, col)) 
            return 1.0;
         if (passthrough) 
            return 0.0;
         return agentBodies.at(row, col);
         // Commented-out version below is for if I make obstacle field non-boolean
         //double prop = 0.0;
         //if (obstacles) prop += obstacles->at(row, col);
         //if (!passthrough) prop += agentBodies.at(row, col);
         //return prop;
      }
      double propObstructed(const Point& p) const { 
         return propObstructed(fp.yToRow(p.y), fp.xToCol(p.x)); 
      }

      bool fullyObstructed(int row, int col, double additional = 0.0) const {
         return propObstructed(row, col) + additional >= 1.0;
      }
      bool fullyObstructed(const Point& p, double additional = 0.0) const {
         return fullyObstructed(fp.yToRow(p.y), fp.xToCol(p.x), additional);
      }

      ConcType concentration(int isub, int row, int col) const {
         return grids.at(isub).at(row, col);
      }

      ConcType concentration(int isub, const Point& p) const {
         return concentration(isub, fp.yToRow(p.y), fp.xToCol(p.x));
      }

      void injectMass(int isub, const Point& p, double mass) {
         grids.at(isub).at(p) += mass / fp.cellArea();
      }

      // TODO this should just take a stream; I should have a helper function in main() that does the file stuff
      void outputField(std::string basefn, int isub) const {
         std::string subfilename{basefn + "-fields-" + std::to_string(isub)};
         std::ofstream subfile(subfilename);
         subfile << fp.xpixels << std::endl;
         subfile << fp.ypixels << std::endl;
         subfile << fp.xphys << std::endl;
         subfile << fp.yphys << std::endl;
         for (int r = 0; r < fp.ypixels; ++r) {
            for (int c = 0; c < fp.xpixels; ++c) {
               if (!fullyObstructed(r, c)) {
                  subfile << grids.at(isub).at(r, c) << std::endl;
               } else {
                  subfile << -1.0 << std::endl;
               }
            }
         }
         subfile << std::endl;
         subfile.close();
         std::system(("bzip2 -f " + subfilename).c_str());
      }

      void serialize(std::ostream& os) const {
         os.write(reinterpret_cast<const char *>(&fp), sizeof(fp));
         os.put(static_cast<char>(passthrough));
         os.write(reinterpret_cast<const char *>(&gravityDir), sizeof(gravityDir));
         agentBodies.serialize(os);
         for (const auto& g : grids) g.serialize(os);
         os.put(static_cast<char>(obstacles.has_value()));
         if (obstacles)
            obstacles->serialize(os);
         os.put(static_cast<char>(cues.has_value()));
         if (cues)
            cues->serialize(os);
         os.put(static_cast<char>(front.has_value()));
         if (front)
            front->serialize(os);
         if (!passthrough) {
            instruments->serialize(os);
            edgeFreepartHorizontal->serialize(os);
            edgeFreepartVertical->serialize(os);
            edgeBodyVelPropHorizontal->serialize(os);
            edgeBodyVelPropVertical->serialize(os);
         }
      }
};


inline bool whetherToUpdate(double seconds, double intervalInSeconds, double dt) {
   if (intervalInSeconds < 0.0) return false;
   if (intervalInSeconds == 0.0) return true;

   return std::fmod(seconds, intervalInSeconds) < std::nextafter(dt, 0.0);
}


#endif
