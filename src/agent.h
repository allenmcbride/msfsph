// agent.h
// Allen McBride
// April 14, 2022
//
// The Agent class encapsulates the code that is within an agent's
// control, with the exception of integrating the ODEs directly translated
// from Morphgen PDEs, which is encapsulated by Nucleus. In other words,
// Agent handles, in an interpretive manner, the logic described in Stages
// II and III of compilation as described in McBride 2022. Other than
// for debugging, key information is hidden from an Agent, most notably
// its location and any information about its neighbors not mediated by
// morphogen fields.

#ifndef AGENT_H
#define AGENT_H

#include "managerinterface.h"
#include "point.h"
#include "nucleus.h"
#include "substances.h"

#include <array>
#include <iostream>
#include <map>
#include <numeric>


class MaxSpeedExceeded : public std::runtime_error { 
   public: using runtime_error::runtime_error;
};


// Parameters to be passed from higher-level classes to new Agent
// instances.
struct AgentParams {
   double sensordistMajor;
   double sensordistMinor;
   double mass;
   double leakRate;
   double maxSpeed; // Values less than zero mean no limitation
   double maxAngularSpeed;
   double maxSettleTime;
   double sensorLagTime;
   double speedupFactor; // DEBUG
};


// Information and operations related to outcome (or current state)
// of a single calibration run.
struct CalibRecord {
   double val;
   Point grad;
   friend CalibRecord operator+(const CalibRecord& lhs, const CalibRecord& rhs) { 
      return {lhs.val + rhs.val, lhs.grad + rhs.grad}; 
   }
   friend CalibRecord operator-(const CalibRecord& lhs, const CalibRecord& rhs) { 
      return {lhs.val - rhs.val, lhs.grad - rhs.grad}; 
   }
   friend CalibRecord operator*(const double k, const CalibRecord& a) { return {k * a.val, k * a.grad}; }
   friend CalibRecord operator/(const CalibRecord& a, const double k) { return {a.val / k, a.grad / k}; }
   friend std::ostream& operator<<(std::ostream& s, const CalibRecord& a) { s << a.val << "," << a.grad; return s; }
};


typedef std::array<CalibRecord, nSPHSubs> CalibBySub;
typedef std::map<double, CalibBySub> CalibData;


// Look up secretion rate and speed in calibration table and return
// interpolated value.
CalibRecord self(int isub, double prod, double speed, bool debug = false);


class Agent {

   // Other than for debugging, I hide information in Agent from the
   // Nucleus class except as mediated by AgentInterface. Constraining
   // information shared between Agent and Nucleus helps clarify the
   // distinction between Stage I of compilation, the result of which
   // is reflected in Nucleus, and Stages II and III, which are handled
   // in interpretive form by Agent.
   friend class AgentInterface;
   //friend class Nucleus; // DEBUG

   private:

      static constexpr double propOfOriginalToAllow = 0.01;


      // ---------- Data members ----------
      
      // Constants
      const AgentParams aParams;
      const Substances subs;
      const bool isOffspring;
      const DiffStateVals parentalTargetDensities{}; // TODO I have four of these arrays. Maybe I should make an array of structs instead.
      const DiffStateVals parentalRecentDensities{};

      // ManagerInterface and Nucleus I think of metaphorically as
      // organelles.
      ManagerInterface mi;
      Nucleus nucleus;

      // Sensor values
      std::array<SensorVals, nSubs> sensors{};
      Point gravityDir{};
      
      // Epigenetic memory (to be inherited by offspring)
      FieldVals fieldValsTimesState{};
      FieldVals laplCorrection{};
      DiffStateVals leakyDensityIntegrals{}; // for findLambdaProportional
      DiffStateVals lambdas{}; // for findLambdaProportional
      
      // Miscellaneous memory, not inherited by offspring
      double lifetime = 0.0;
      SPHVals prod {}; // raw production levels. Zero-init because it's used in self() before it's assigned.
      double speed = 0.0; // most recent speed requested. Old speed used before assigning new, hence init here
      double angularVel = 0.0; // Not using the memory of this yet, but it could be useful to be remembered in future
      std::array<bool, nDiffStates> prevGrowthSwitches{}; // most recent growth switches, to detect change from false to true
      DiffStateVals targetDensities{};
      //Point accMem; // DEBUG


      // ---------- Miscellaneous methods ----------
      
      // If a sensor's raw reading is valid, update sensor values,
      // otherwise leave them unchanged.
      void interpretSensors(int isub);

      // Field values should often be multiplied by the value of the
      // corresponding differentiation state, so these products are stored
      // as data members and updated once per update by this function.
      void calcFieldValsTimesState();

      // In SPH, the value of a field divided by density is often
      // needed. Special cases are useful When the field is density
      // itself, or when differentiation states are involved, so I use
      // a function for this operation.
      double valOverRho(int isub) const;

      // Compute lambda, which determines probability of
      // reproduction/recruitment or apoptosis/inactivation, based
      // on proportional control method and integral control method,
      // respectively.
      double findLambdaProportional(double dt);
      double findLambdaIntegral(double dt);

      // Adjust velocity based on Morphgen program's specification for
      // how aggressively to counteract unsmoothness (that is, variations
      // with high spatial frequency) in the density and velocity fields.
      Point smoothnessAdjustments(const Point& desiredVel) const;

      // Adjust velocity based on linear and angular speed limits.
      void setSpeedsFromLimits(const Point& desiredVel, double dt);


      // ---------- Field value methods, scalar and vector ---------- 
      
      // Estimate value of a field without applying any corrections.
      double rawVal(int isub) const;

      // Compute correction term for field value based on speed and
      // calibration data.
      double speedValueCorrection(int isub) const;

      // Compute correction term for field value based on a new
      // agent's not having had a chance for its own secreted morphogen
      // contribution to have reached steady state.
      double offspringValueCorrection(int isub) const;

      // Compute corrected field value estimates, scalar and vector.
      double findVal(int isub) const;
      Point findVal(int isubx, int isuby) const;


      // ---------- Field gradient methods ---------- 

      // Estimate gradient of a field without applying any corrections.
      Point rawGrad(SensorVals sv, bool debug = false) const;

      // Compute corrected field gradient estimate.
      Point findGrad(int isub) const;

      // Estimate velocity . gradient(field) with upstream correction
      // term.
      double velDotGradUpstream(int isub, const Point& vel) const;


      // ---------- Laplacian methods, scalar and vector ---------- 
      double findLaplacian(int isub, bool debug = false) const;
      Point findLaplacian(int isubx, int isuby) const;


      // ---------- Other derivatives ----------
      double findDiv(int isubx, int isuby, Point refDirX) const; // Untested


      // ---------- Edge detection, density and color field versions ----------
      double edgeDetect(int isub) const;
      double colorEdge(int isub) const;


   public:

      //ManagerInterface mi; // DEBUG (should be private)
      static CalibData cdat;
      static void printCalibDat();


      // ---------- Methods for normal operation ----------

      // Construct new agents at beginning of simulation.
      Agent(AgentParams aParams, Substances subs, ManagerInterface mi, NuclearOutput nucOut);

      // Construct new agent based on a parent agent.
      Agent(const Agent& a, ManagerInterface mi);

      // Unfinished deserializing constructor.
      Agent(ManagerInterface mi, std::istream& is);

      // Update agent state one timestep, for ordinary simulation.
      double update(double dt, bool debug = false);


      // ---------- Utility methods for setup ----------
      
      // Execture one timestep of burn-in procedure, which performs
      // secretion but does not move or update agent state.
      void burnin(int isub, double dt);

      // To support burn-in procedure. Based on value for a given field,
      // estimate total steady state morphogen mass contributed by
      // this agent.
      double integral(int isub);

      // Run calibration procedure for one timestep.
      CalibBySub calibrate(Point force, double multiplier, double dt, bool debug = false);

      // Report all SPH density estimates.
      DiffStateVals provideDens() const;

      // Report values for all fields that the Nucleus class tracks.
      auto provideNucOut() const { return nucleus.provideVals(); } // DEBUG


      // ---------- Other utility ----------

      // Report agent mass.
      double provideMass() const { return aParams.mass; }

      // Unfinished serialization.
      void serialize(std::ostream& os) const;


      // ---------- Debugging methods ---------- 
      void printInfo() const {} // DEBUG
      Point provideDensGrad() const { return findGrad(DiffStates::one); } // DEBUG
#if PROGRAM == EDGE
      Point provideColorGrad() const { return findGrad(FieldVars::color); } // DEBUG
#else
      Point provideColorGrad() const { return {}; } // DEBUG
#endif
      double provideSpeed() const { return speed; } // DEBUG
      Point provideVel() const { return Point::polarPt(speed, angularVel); } // DEBUG
      //Point provideAcc() const { return accMem; } // DEBUG
      //double provideLifetime() const { return lifetime; } // DEBUG
};

#endif
