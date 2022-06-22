// agentinterface.h
// Allen McBride
// April 20, 2022
//
// The AgentInterface class controls access by a Nucleus object to
// the Agent object it is owned by. Each Nucleus object owns one of
// these AgentInterface objects, constructed and given by the owning
// Agent. AgentInterface is a friend of Agent. The functions here are
// trivial, yet must be put in a separate .cpp file to avoid circular
// dependencies.

#ifndef AGENTINTERFACE_H
#define AGENTINTERFACE_H

#include "point.h"

class Agent;

class AgentInterface {

   private:

      // For debugging, it is helpful to give a Nucleus direct access to
      // its Agent, which can be accomplished by moving this declaration
      // to public.
      Agent &a;


   public:

      // Construct new AgentInterface with reference to given Agent
      AgentInterface(Agent &a) : a{a} {}

      // Get field value estimates.
      double findVal(int isub) const;
      Point findVal(int isubx, int isuby) const;

      // Get field gradient-related estimates.
      Point findGrad(int isub) const;
      double velDotGradUpstream(int isub, Point vel) const;

      // Get field Laplacian estimates.
      double findLaplacian(int isub, bool debug = false) const;
      Point findLaplacian(int isubx, int isuby) const;

      // Get divergence estimates (untested).
      double findDiv(int isubx, int isuby, Point refDirX) const;
      
      // Get edge detection estimates.
      double edgeDetect(int isub) const;
      double colorEdge(int isub) const;

      // Get direction of gravity.
      Point findGravityDir() const;

      // Check if Agent is in contact with given environmental cue.
      bool findCue(int icue) const;

      // Get most recent Agent speed.
      double speed() const;
};

#endif
