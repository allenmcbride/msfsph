// nucleus.h
// Allen McBride
// June 21, 2022
//
// The Nucleus class encapsulates the code that varies from one Morphgen
// program to the next. That is, it handles integration of the ODEs
// translated from Morphgen PDEs. This file also defines structures for
// the state that Nucleus contains and passes to its owning Agent. There
// is also a DualRail subclass, intended for representing signed
// quantities using two substances, that needs more testing.

#ifndef SPHEQS_H
#define SPHEQS_H

#include "agentinterface.h"
#include "point.h"
#include "serialization.h"
#include "substances.h"

#include <array>
#include <iostream>

typedef std::array<double, nFieldVars> FieldVals;
typedef std::array<double, nDiffStates> DiffStateVals;
typedef std::array<double, nSPHSubs> SPHVals;

struct NuclearOutput {
   std::array<double, nDiffStates> diffStateVals = {1.0};
   std::array<double, nFieldVars> fieldVals = {};
   std::array<bool, nDiffStates> growthSwitches = {};
   std::array<double, nDiffStates> growthRates = {};
   std::array<double, nDiffStates> smoothDensityGuards = {};
   std::array<double, nDiffStates> smoothVelocityGuards = {};
   // Velocity here is relative to the agent's orientation, which means that
   // heading changes are integrated. So a constant 'velocity' would generally
   // result in acceleration.
   Point velocity = {};
   void print() {
      for (auto& e : diffStateVals) std::cout << e << ", ";
      std::cout << ";     ";
      for (auto& e : fieldVals) std::cout << e << ", ";
      std::cout << ";     ";
      for (auto& e : growthSwitches) std::cout << (e ? "true" : "false") << ", ";
      std::cout << ";     ";
      for (auto& e : growthRates) std::cout << e << ", ";
      std::cout << std::endl;
   }
};


struct InternalVals {
#if PROGRAM == GROWTH || PROGRAM == MIXTURE || PROGRAM == GROWX
      double timer = 0.0;
      double lock = 0.0;
#elif PROGRAM == PATH
      double lock = 0.0;
#endif
};


class Nucleus {
   private:
      static int nextid; // DEBUG

      class DualRail {
         private:
            double &pos;
            double &neg;
         public:
            DualRail(double& pos, double& neg) : pos(pos), neg(neg) {}
            void operator=(double val) { if (val > 0.0) { pos = val; neg = 0.0; } else { pos = 0.0; neg = -val; } }
            explicit operator double() { return pos - neg; }
      };

      int id; // DEBUG
      AgentInterface ai;
      NuclearOutput outVals;
      InternalVals internal;

      double& fieldRef(int isub);
      void fieldSet(int isub, double val) { fieldRef(isub) = val; }
      void fieldChange(int isub, double rate, double dt) { fieldRef(isub) += dt * rate; }
   public:

      Point nucAdj; //DEBUG
      bool hit; //DEBUG

      Nucleus(AgentInterface ai, NuclearOutput outVals) : 
         ai{ai}, 
         outVals{outVals} 
      {
         //DualRail actDR(sVars.at(SubsID::activatorPos), sVars.at(SubsID::activatorNeg));
         id = nextid++; // DEBUG
         //std::cerr << id << std::endl; // DEBUG
      }
      Nucleus(AgentInterface ai, const Nucleus& parent) : 
         ai{ai},
         outVals{parent.outVals},
         internal{parent.internal}
      {
         id = nextid++; // DEBUG
         //std::cerr << id << std::endl; // DEBUG
      }
      Nucleus(AgentInterface ai, std::istream& is) :
         id{readObjectFromStream<int>(is)},
         ai{ai},
         outVals{readObjectFromStream<NuclearOutput>(is)},
         internal{readObjectFromStream<InternalVals>(is)}
      {
         if (id >= nextid) nextid = id + 1;
      }

      void update(double dt, bool debug = false);
      NuclearOutput provideVals() const { return outVals; }

      void serialize(std::ostream& os) const {
         os.write(reinterpret_cast<const char *>(&id), sizeof(id));
         os.write(reinterpret_cast<const char *>(&outVals), sizeof(outVals));
         os.write(reinterpret_cast<const char *>(&internal), sizeof(internal));
      }
};

#endif
