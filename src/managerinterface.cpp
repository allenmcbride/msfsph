// managerinterface.cpp
// Allen McBride
// April 20, 2022
//
// Descriptions in managerinterface.h

#include "managerinterface.h"
#include "agentmanager.h"

Point ManagerInterface::senseGravityDir() const {
   return am.senseDir(am.fields.gravityDir);
}

bool ManagerInterface::senseCue(int icue) const {
   return am.fields.cues->at(am.extState.pos).test(icue);
}

SensorVals ManagerInterface::sense(int isub) const { 
   return am.sense(isub); 
}

bool ManagerInterface::sensorsReady() const {
   return am.sensorsReady;
}

void ManagerInterface::produce(int isub, double quantity, double dt) const { 
   am.produce(isub, quantity, dt); 
}

void ManagerInterface::move(double desiredSpeed, double desiredAngularVel, double collResSpeed, double dt) {
   am.move(desiredSpeed, desiredAngularVel, collResSpeed, dt); 
}

// For debugging only. An agent should not logically have access to
// this information.
//void ManagerInterface::report(int n) const { 
//   am.report(n); 
//}
