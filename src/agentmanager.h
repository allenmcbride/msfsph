// agentmanager.h
// Allen McBride
// March 21, 2022
//
// The AgentManager class encapsulates the code that relates to
// a particular agent but that needs to know more than the agents
// themselves should. Most notably this knowledge includes a reference
// to the Fields object that holds all actual grids, and it includes an
// agent's location and heading. The most important behaviors defined
// in this class are agent initialization, motion, secretion, and sensing.

#ifndef AGENTMANAGER_H
#define AGENTMANAGER_H

#include "agent.h"
#include "fields.h"
#include "managerinterface.h"
#include "point.h"
#include "substances.h"

#include <array>
#include <exception>
#include <random>
#include <string>
#include <variant>
#include <vector>


// Type for number of particles bounds to a sensor.
typedef long long NBoundType;


// Just a convenient way to refer to which of the two corners of a grid
// cell, along a given side, are inside the body of an agent. "high"
// (or "low") means the corner with the higher-valued (or lower-valued)
// coordinate inside the agent, and the other is not. "span" means
// the body does cut through this side of the grid cell, but both
// intersections are between the two corners, so both lie outside.
enum class CornersInside { both, neither, high, low, span };


// Values to return to the World class to indicate whether, for a given
// update, this agent has decided to reproduce or apoptose.
enum class LifeEvent { none, reproduce, die };


// For a given grid line, horizontal or vertical, flag whether the agent
// intersects with it (called "real" based on the relevant roots being
// real numbers), and, if so, the lower and higher valued coordinate
// (in the relevant dimension) of these intersections.
struct Intersections {
   bool real;
   double low;
   double high;
   friend std::ostream& operator<<(std::ostream& s, const Intersections& a) { s << "<" << a.low << "," << a.high << ">" << (a.real ? "" : "none"); return s; }
};


// Map to hold the information in Intersections for all grid lines of
// a given dimension (horizontal or vertical)
typedef std::map<double, Intersections> IntersectionsByLine;


// Agent could not be placed at requested location because it would at
// least partly intersect with an obstruction.
class ObstructionException : public std::runtime_error { 
   public: using runtime_error::runtime_error;
};


// New agent would be out of bounds at requested location.
class OOBException : public std::runtime_error { 
   public: using runtime_error::runtime_error;
};


// Struct to gather parameters relating to the overall simulations,
// not just a single agent.
struct SimulationParams {
   bool passthrough;
   double leakRate;
   double maxSpeed;
   double maxAngularSpeed;
   double speedupFactor;
};


// For all noise parameters, zero should turn off that type of noise.
struct AgentManagerParams {
   double semimajorAxis;
   double semiminorAxis;
   double stalkLen;

   // Noise-type parameters:
   double sensorBiasLogSD;
   // A non-positive value for sensorEventRatePerConcentration implies
   // a continuous model, or infinite resolution for a sensor.
   double sensorEventRatePerConcentration; 
   double sensorMeanBindingTime;
   double oriBiasSD;
   double speedBiasLogSD;
   double motionNoiseRelSD;
   double sensorMTTF;
   double propulsionMTTF;
   double productionMTTF;

   // Assume density of agent is 1.0, as would be roughly true for
   // a biological cell, so that mass density has the same value as
   // area/volume.
   double mass() const {
      return M_PI * semimajorAxis * semiminorAxis;
   }

   // Return a copy of self but with all noise-type parameters set
   // to zero.
   AgentManagerParams noNoise() const {
      AgentManagerParams newParams{};
      newParams.semimajorAxis = this->semimajorAxis;
      newParams.semiminorAxis = this->semiminorAxis;
      newParams.stalkLen = this->stalkLen;
      newParams.sensorMeanBindingTime = this->sensorMeanBindingTime;
      return newParams;
   }
};


// Agent's current location and heading
struct AMExternalState {
   Point pos;
   double ori;
};


class AgentManager {

   // ManagerInterface mediates the Agent class's access to AgentManager
   // so that the former can know only what it logically should. Most
   // notably, the Agent should not have access to the grid except
   // through its instruments, and it should not know the agent's location
   // or heading.
   friend class ManagerInterface;


   // The next two lines, marking the Nucleus or Agent classes as
   // friends, are useful for debugging but should be commented out
   // to demonstrate that agents do not need localization or global
   // information in MSF-SPH.
   //friend class Nucleus; // DEBUG
   //friend class Agent; // DEBUG


   private:
      static constexpr double propOfOriginalToAllow = 0.01;
      static int nextid; // DEBUG


      // Sensor defines each of the four sensors around an agent's body.
      struct Sensor {
         Point location{};
         std::array <NBoundType, nSubs> nBound{};
         std::array <double, nSubs> amountBound{};
         std::array <double, nSubs> reading{};
         double bias = 0.0; // TODO Make bias a per-substance trait
         double lifespan = 0.0;
      };


      // The Record and Memory classes are only for debugging. Records
      // of various information are written into a circular buffer so
      // that the history of an agent can be recalled when some future
      // condition is met.
      struct Record {
         int id;
         Point pos;
         double orientation;
         Point absVel;
         Point densGradRel;
         double speed;
         Point relVel;
         //Point acc;
         NuclearOutput nucOut;
         SensorVals freeparts;
      };


      class Memory {
         private:
            std::array <Record, 200> records;
            size_t head = 0;
            size_t nUsed = 0;
         public:
            void remember(const Record& r);
            void report(int n) const;
      };


      // ---------- Data members ----------

      // Persistent data: constant or updated incrementally (at least conceptually)
      const AgentManagerParams amParams;
      AMExternalState extState;
      Fields& fields;
      std::default_random_engine& gen;
      const double maxSpeed;
      const bool passthrough;
      const double maxSettleTime;
      Agent agent;

      // Persistent data set or modified by init()
      int id = 0; // Global ID for debugging; note Agent itself does not have access to this
      std::array <Sensor, 4> sensors{};
      std::array <Point, 4> prodLocations{};
      double oriBias = 0.0;
      double speedBias = 0.0;
      double productionLifespan = 0.0;
      double propulsionLifespan = 0.0;

      // Other persistent data (TODO: sensorLifetime, readinessCount,
      // and sensorReady might be more appropriate inside Sensor class.)
      double lifetime = 0.0;
      double sensorLifetime = 0.0;
      Memory memory{};
      bool sensorsReady = false;
      NBoundType readinessCount = 0;
      double postCollisionTimer = 0.0;

      // Data needed temporarily: within update or from one update to next
      Point vel{};
      double angularVel = 0.0;
      int cornerRow = 0;
      int cornerCol = 0;
      std::vector<std::vector<AgentBodyType> > proportions;
      std::vector<std::vector<EdgeType> > edgesOccludedVertical;
      std::vector<std::vector<EdgeType> > edgesOccludedHorizontal;
      

      // ---------- Methods accessible through ManagerInterface ----------
      
      // Report directional cue relative to heading. p should be a
      // unit vector.
      Point senseDir(const Point& p) const { return p.rotated(-extState.ori); }
                                    
      // Report sensor readings for a given morphogen
      SensorVals sense(int isub) const;

      // Secrete a given morphogen for one timestep
      void produce(int isub, double rate, double dt) const;

      // Attempt to move as requested for one timestep
      void move(double desiredSpeed, double desiredAngularVel, double collResSpeed, double dt);


      // ---------- Internal motion-related helper functions ----------

      // Calculate body velocity at given location. Argument should be
      // a point inside agent's body.
      Point bodyVelAtPt(const Point& p) const { return vel + angularVel * (p - extState.pos).rotated(M_PI / 2.0); }

      // Based on body position, update locations of sensors and secretion
      // effectors relative to grids.
      void updateInstrumentLoc();

      // Calculate intersection of agent's surface with ray from center at
      // given angle relative to heading. extraRadius adds more distance
      // from the center.
      Point surfacePoint(double angleParam, double extraRadius) const;

      // Compute "proportions" array, which tracks how much of each grid
      // cell near an agent's body is occupied by that agent. Compute
      // similar information with respect to grid edges for purposes of
      // "swimming" fluid dynamics model. Update cornerRow and cornerCol
      // as needed, which record position of proportions array in overall
      // grid. Call updateInstrumentLoc(). Check for collisions and
      // return the number of grid cells in which collisions are found.
      int occlusion(double leeway = 0.0, bool debug = false);

      // Compute maps of intersections of body surface with grid edges,
      // one for map for horizontal and one for vertical grid edges.
      // Given the x- or y-coordinate for a given grid edge, this map
      // contains an Intersection object recording the y- or x-coordinate,
      // respectively, of intersections (if any).
      std::pair<IntersectionsByLine, IntersectionsByLine> findBoundIntersections(int rowmin, int rowmax, int colmin, int colmax) const;

      // Based on cell and agent body intersect with a given grid edge,
      // determine which corners of grid cell (along given edge) are
      // inside agent body. Possibilities are: low or high (corner with
      // lower or higher coordinate value is inside agent), neither,
      // both, or "span" (agent intersects this grid cell between the
      // two corners but neither is inside agent)
      CornersInside findCornersInside(double cornerLow, double cornerHigh, const Intersections& intersections) const;

      // Similar to findCornersInside, but calculate length of grid cell
      // edge inside agent.
      double findEdgeOccluded(double cornerLow, double cornerHigh, const Intersections& intersections) const;

      // Calculate proportion of grid cell (identified by col and row)
      // inside agent body, approximating edge of body within grid cell
      // as a straight line segment.
      double proportionOccluded(int col, int row, const IntersectionsByLine& colIntersections, const IntersectionsByLine& rowIntersections, bool debug = false);

      // Modify shared grids as appropriate to reflect presence of agent
      // (assertOccupation) or remove it (eraseOccupation), such as is
      // needed to determine potential collisions.
      void assertOccupation(bool debug = false);
      void eraseOccupation(bool debug = false);

      // Like assertOccupation() but for information not needed to
      // determine collisions. Can be run just once per update cycle
      // rather than each time an agent is moved while looking for a
      // collision-free location.
      void assertTransientInfo();

      // Return a vector from the center of the agent's position toward
      // the high-valued corner of a grid-aligned box circumscribing
      // the agent's body.
      Point boundingBoxRelative() const;

      // If agent is not in bounds, move it to the nearest in-bounds
      // position.
      void moveInBounds(const Point& targetLocation);

      // If occlusion() returns a positive number for an agent's desired
      // new location, attempt to find a nearby position such that
      // occlusion() returns zero. Eventually resort to original position.
      bool resolveCollisions(const AMExternalState& origExtState, const Point& displacement, double oriChange, double collResDisp);


      // ---------- Internal sensing-related helper function ----------
      void sensorUpdate(double dt, bool trueValue = false);




      // ---------- Factored constructor code ----------
      void init();


   public:

      // ---------- Constructors, destructors, and related ----------
      
      // Constructor for new agents that are not offspring of another
      AgentManager(
            const AgentManagerParams& amParams,
            const AMExternalState& extState,
            Fields& fields, 
            const Substances& subs, 
            const SimulationParams& simParams,
            std::default_random_engine& gen,
            const NuclearOutput nucOut = NuclearOutput{1.0} // Initialize density to 1, all else to 0.
            );

      // Constructor for reproduction/recruitment. Copies some state,
      // with or without modification, from state of parent agent.
      AgentManager(const AgentManager& am, const Point& pos);

      // Unfinished deserializing constructor
      AgentManager(Fields& fields, std::default_random_engine& gen, std::istream& is);

      // Destructor removes agent from shared grids
      ~AgentManager();

      // Unfinished serializer
      void serialize(std::ostream& os) const;


      // ---------- Methods for use throughout simulation ----------

      // Main update function returns enum signalling if agent should
      // reproduce or apoptose or neither
      LifeEvent update(double dt, bool debug = false);

      // Check if agent is in bounds
      bool inBounds() const;
      

      // ---------- Inline methods ----------
      
      // Call methods needed to assert agent info on shared grids
      void assertAgentInfo() { assertOccupation(); assertTransientInfo(); }

      // Not for ordinary simulation, but to help determine maximum
      // angular speed
      void rotateOnly(double angle) { eraseOccupation(); extState.ori += angle; occlusion(); assertOccupation(); }

      // Run calibration procedure for one timestep
      CalibBySub calibrate(Point force, double multiplier, double dt, bool debug = false) { sensorUpdate(dt, true); return agent.calibrate(force, multiplier, dt, debug); }

      // Report position and heading
      AMExternalState provideState() const { return extState; }

      // Report agent size dimensions
      std::pair<double, double> provideAxes() const { return {amParams.semimajorAxis, amParams.semiminorAxis}; }


      // ---------- Setup and final reporting ----------
      
      // Report true mass density morphogen concentration
      double provideDensOneConc() const;

      // Perform one timestep of burn-in procedure (that is, secrete
      // morphogen based on field values, but don't react)
      void burnin(int isub, double dt) { sensorUpdate(dt); agent.burnin(isub, dt); };

      // To support burn-in procedure. Based on value for a given field,
      // estimate total steady state morphogen mass contributed by
      // this agent.
      double integral(int isub) { return agent.integral(isub); }

      // Report mass
      double provideMass() const { return agent.provideMass(); }


      // ---------- Debugging ----------
      
      void report(int n) const { memory.report(n); agent.printInfo(); } // DEBUG
      auto provideNucOut() const { return agent.provideNucOut(); } // DEBUG
      DiffStateVals provideDens() const { return agent.provideDens(); } // DEBUG
      Point provideDensGrad() const { return agent.provideDensGrad(); } // DEBUG
      Point provideColorGrad() const { return agent.provideColorGrad(); } // DEBUG
      Point provideVel() const { return vel; } // DEBUG
      int provideID() const { return id; } // DEBUG
};

#endif
