// world.h
// Allen McBride
// April 19, 2022
//
// The World class encapsulates information and methods relating to
// the simulation globally. It sets up various fields defining the
// domain, makes the initial placement of agents, and deals with agent
// reproduction/recruitment and apoptosis. It also manages the output
// of information at the end of a simulation.

#ifndef WORLD_H
#define WORLD_H

#include "agentmanager.h"
#include "diffuser.h"
#include "fields.h"
#include "point.h"
#include "substances.h"

#include <functional>
#include <iostream>
#include <list>
#include <random>
#include <set>
#include <string>


// AgentPlacement encapsulates information about how agents should be
// placed initially, which is passed in from main().
struct AgentPlacement {
   size_t nAgents;
   double initDens;
   std::function<bool(const Point&)> inShape;
};


class World {
   private:

      // Constants defining precision of fluid velocity calculations.
      static constexpr double viscosityTolerance = 0.05;
      static constexpr double pressureTolerance = 0.01;

      // Data members
      std::function<bool(const Point&)> inShape;
      Fields fields;
      Substances subs;
      Diffuser diffuser;
      std::list <AgentManager> agentManagers;
      std::default_random_engine& gen;
      //int nextid = 0; // DEBUG

      // Set any environmental cues in the domain.
      void setCues();
      
      // Place any obstacles in the domain.
      void placeObstacles();

      // Provide initial field values at a given location.
      NuclearOutput initVals(const Point& point) const;

      // Place a new agent near a given agent, if possible.
      void placeNearby(const AgentManager& am);

      // Use Monte Carlo integration to estimate area of shape in which
      // agents are to be placed.
      double areaInShape();


   public:

      // Constructor for setting up new World
      World(
            const AgentPlacement& agentPlacement,
            const Substances& subs,
            double visc,
            const AgentManagerParams& amParams,
            const SimulationParams& simParams,
            const Fparams& fparams,
            bool visual,
            std::default_random_engine& gen
           );

      // Unfinished deserializing constructor
      World(
            std::default_random_engine&,
            std::istream& is
            );

      // Unfinished serializer
      void serialize(std::ostream& os);

      // Advance World instance one timestep (agents and fluid medium)
      void update(double dt, bool debug = false);

      // Allow agents to secrete based on their initial values but not
      // update in any other way, until only propOfOriginalToAllow
      // proportion of original state remains. Allows for starting
      // simulation without an initial period of confusion
      void burnin(double propOfOriginalToAllow);

      // Destroys agents in a band at a particular time. For testing
      // healing with GROWTH program.
      void cut();

      // Getter for Fields object
      const Fields& provideFields() const { return fields; }

      // Provide number of agents remaining in domain
      int nAgentsRemaining() const { return agentManagers.size(); }

      // Print some statistics to the console
      void textupdate(double t);

      // Print velocity field to a file of given name
      void outputVel(std::string basefn) const;

      // Make Gaussian-based kernel density estimate, print to
      // file. O(n^2); could be made O(n ln n) with an octree, but
      // currently this is only called once at the end of a simulation
      // so it's not worth the trouble.
      void gaussianDens(std::string basefn) const;

      // Print to file list of agents and their values for field
      // corresponding to substance isub.
      void agentFieldVals(std::string basefn, int isub) const;

      // Print to file value of edgeDetect and colorEdge for each agent,
      // for testing these mechanisms.
      void edgeInfo(std::string basefn) const;

      // Print statistics about agents to console
      void printAgentStats() const;

      // Fill region of initial shape with given concentration of density
      // morphogen. For comparing MSF-SPH results with ideal results.
      void fillShape(double concentration);

      // Fill domain with Gaussian distribution of density morphogen,
      // for comparing MSF-SPH ADVECTION program with ideal results.
      void fillGaussian();

      // Find average orientation; for debugging
      void avgOri() const; // DEBUG
      
      // Fill domain with square dots of density morphogen, for
      // visualization in testing fluid velocity.
      void dots(); // DEBUG
      
      // Fill domain with density morphogen of constant concentration,
      // for visualization in testing fluid velocity.
      void fill(); // DEBUG
      
      //Print information about a particular agent. For debugging.
      //void reportOneAgent(int id, int nRec) const { agentManagers.at(id).report(nRec); }
};

#endif
