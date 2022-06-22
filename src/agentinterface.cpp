// agentinterface.cpp
// Allen McBride
// April 20, 2022
//
// Descriptions in agentinterface.h

#include "agentinterface.h"
#include "agent.h"

double AgentInterface::findVal(int isub) const { return a.findVal(isub); }
Point AgentInterface::findVal(int isubx, int isuby) const { return a.findVal(isubx, isuby); }

Point AgentInterface::findGrad(int isub) const { return a.findGrad(isub); }
double AgentInterface::velDotGradUpstream(int isub, Point vel) const { return a.velDotGradUpstream(isub, vel); }

double AgentInterface::findLaplacian(int isub, bool debug) const { return a.findLaplacian(isub, debug); }
Point AgentInterface::findLaplacian(int isubx, int isuby) const { return a.findLaplacian(isubx, isuby); }

double AgentInterface::findDiv(int isubx, int isuby, Point refDirX) const { return a.findDiv(isubx, isuby, refDirX); }

double AgentInterface::edgeDetect(int isub) const { return a.edgeDetect(isub); }
double AgentInterface::colorEdge(int isub) const { return a.colorEdge(isub); }
Point AgentInterface::findGravityDir() const { return a.gravityDir; }
bool AgentInterface::findCue(int icue) const { return a.mi.senseCue(icue); }
double AgentInterface::speed() const { return a.speed; }
