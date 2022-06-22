// managerinterface.h
// Allen McBride
// April 20, 2022
//
// The ManagerInterface class controls access by an Agent object to the
// AgentManager object it is owned by. Each Agent object owns one of
// these ManagerInterface objects, constructed and given by the owning
// AgentManager. ManagerInterface is a friend of AgentManager. The
// functions here are trivial, yet must be put in a separate .cpp file
// to avoid circular dependencies.

#ifndef MANAGERINTERFACE_H
#define MANAGERINTERFACE_H

#include "point.h"

#include <string>
#include <array>

class AgentManager;

typedef std::array <double, 4> SensorVals;

class ManagerInterface {

   private:

      // For debugging, it is helpful to give an Agent direct access
      // to its AgentManager, which can be accomplished by moving this
      // declaration to public.
      AgentManager& am;


   public:
      
      // Construct new ManagerInterface with reference to given
      // AgentManager.
      ManagerInterface(AgentManager& am) : am(am) {}

      // Provide direction of gravity as sensed by AgentManager.
      Point senseGravityDir() const;

      // Provide whether AgentManager is in contact with a given
      // environmental cue.
      bool senseCue(int icue) const;

      // Call AgentManager's sense() method.
      SensorVals sense(int isub) const;

      // Provide whether AgentManager's sensors have seen enough binding
      // events to render an acceptable concentration estimate.
      bool sensorsReady() const;

      // Call AgentManager's produce() method (e.g., secretion).
      void produce(int isub, double quantity, double dt) const;

      // Call AgentManager's move() method.
      void move(double desiredSpeed, double desiredAngularVel, double collResSpeed, double dt);

      //void report(int n = -1) const; // DEBUG
};

#endif
