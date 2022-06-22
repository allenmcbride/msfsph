// agent.cpp
// Allen McBride
// April 14, 2022
//
// Descriptions in agent.h

#include "agent.h"
#include "agentinterface.h"
#include "managerinterface.h"
#include "nucleus.h"
#include "substances.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>


CalibData Agent::cdat;
//DiffStateVals findEdgeConsts(Substances subs);
CalibData::iterator findAbove(double speed);
double speedProp(CalibData::iterator above, CalibData::iterator below, double speed, bool debug);
double limitedSpeed(double desired, double max);


// Interpret raw sensor data. Currently uses previous value for a given sensor if it is in error.
void Agent::interpretSensors(int isub) {
   SensorVals readings = mi.sense(isub);
   for (size_t i = 0; i < sensors.at(isub).size(); ++i)
      if (readings.at(i) >= 0.0) sensors.at(isub).at(i) = readings.at(i) / subs.at(isub).gain;
}


void Agent::calcFieldValsTimesState() {
   for (int iFieldVar = 0; iFieldVar < nFieldVars; ++iFieldVar) {
      fieldValsTimesState.at(iFieldVar) = nucleus.provideVals().fieldVals.at(iFieldVar) 
         * nucleus.provideVals().diffStateVals.at(diffVarForFieldVar(iFieldVar + nDiffStates));
   }
}


double Agent::valOverRho(int isub) const {
   if (isub < nDiffStates) {
      // For differentiation state variables, density is in both numer and denom, so it's divided out, leaving just differentiation state
      return nucleus.provideVals().diffStateVals.at(isub);
   } else {
      const auto numer = fieldValsTimesState.at(isub - nDiffStates);
      const auto denom = findVal(diffVarForFieldVar(isub));
      if (denom == 0.0) {
         //if (numer > 0.0) throw std::runtime_error("in valOverRho(), denom == 0 but numer > 0 for isub = " + std::to_string(isub) + ": " + std::to_string(numer) + " / " + std::to_string(denom));
         return 0.0;
      }
      return numer / denom;
   }
}


double Agent::findLambdaProportional(double dt) {
   for (int iDiff = 0; iDiff < nDiffStates; ++iDiff) {
      if (nucleus.provideVals().growthSwitches.at(iDiff)) {
         if (prevGrowthSwitches.at(iDiff) == false && leakyDensityIntegrals.at(iDiff) == 0.0) {
            // TODO the leakyDensityIntegrals needs to be done in the main update, not just if we're looking for lambda
            leakyDensityIntegrals.at(iDiff) = findVal(iDiff) / aParams.leakRate; 
         } else {
            leakyDensityIntegrals.at(iDiff) += findVal(iDiff) * dt;
            leakyDensityIntegrals.at(iDiff) *= std::exp(-dt * aParams.leakRate);
            if (findVal(iDiff) > 0.0) {
               const auto leakyIntegralEstimate = aParams.leakRate * findVal(iDiff) - leakyDensityIntegrals.at(iDiff) * std::pow(aParams.leakRate, 2.0);
               const auto densityRateDeficit = nucleus.provideVals().growthRates.at(iDiff) - leakyIntegralEstimate;
               const auto lambdaDeficit = densityRateDeficit * nucleus.provideVals().diffStateVals.at(iDiff) / findVal(iDiff);
               lambdas.at(iDiff) += dt * lambdaDeficit / subs.at(iDiff).timeToSteady(propOfOriginalToAllow);
               //std::cerr << "est: " << leakyIntegralEstimate << "; measured: " << findVal(iDiff) << std::endl;
               //std::cerr << "est: " << leakyDensityIntegrals.at(iDiff) * aParams.leakRate << "; measured: " << findVal(iDiff) << std::endl;
            }
         }
         prevGrowthSwitches.at(iDiff) = true;
      }
   }
   return std::accumulate(lambdas.begin(), lambdas.end(), 0.0);
}


// Members changed here: prevGrowthSwitches, targetDensities; all are only used here (19mar21)
double Agent::findLambdaIntegral(double dt) {
   DiffStateVals lambdas{};
   for (int iDiff = 0; iDiff < nDiffStates; ++iDiff) {
      // TODO: write code for if a growth switch is turned back off
      if (nucleus.provideVals().growthSwitches.at(iDiff)) {
         if (prevGrowthSwitches.at(iDiff) == false) {

            if (mi.sensorsReady()) {
               targetDensities.at(iDiff) = parentalTargetDensities.at(iDiff) 
#ifndef NOTARGCORR
                  + rawVal(iDiff) - parentalRecentDensities.at(iDiff)
#endif
                  ;

               lambdas.at(iDiff) = 0.0; // TODO: I think this is redundant

               prevGrowthSwitches.at(iDiff) = true;
            }

         } else {

            targetDensities.at(iDiff) += dt * nucleus.provideVals().growthRates.at(iDiff);
            const auto deltaRho = targetDensities.at(iDiff) - findVal(iDiff);

            if (findVal(iDiff) > 0.0) {
               const auto time = subs.at(iDiff).timeToSteady(propOfOriginalToAllow);
               lambdas.at(iDiff) = deltaRho * nucleus.provideVals().diffStateVals.at(iDiff) / (time * findVal(iDiff));

            } else {
               std::cerr << "in Agent::findLambdaDisc()... findVal(iDiff) nonpositive: " << findVal(iDiff) << std::endl;
            }
            if (lambdas.at(iDiff) < -100.0) {
               std::cerr << " l: " << lambdas.at(iDiff) << "; rawVal: " << rawVal(iDiff) << "; findVal: " << findVal(iDiff) 
                  << "; targ: " << targetDensities.at(iDiff) << "; propMiss: " << std::exp(-subs.at(iDiff).decayRate * lifetime) << "; self.val: "
                  << self(iDiff, nucleus.provideVals().diffStateVals.at(iDiff) * subs.at(iDiff).gain * subs.at(iDiff).decayRate * aParams.mass, 0.0).val
                  << "; speed: " << speed << std::endl;
            }
         }
      }
   }
   return std::accumulate(lambdas.begin(), lambdas.end(), 0.0);
}


Point Agent::smoothnessAdjustments(const Point& desiredVel) const {
   auto minMaxSpeed = std::numeric_limits<double>::max();
   auto maxVelGradDens = Point{};
   for (int iDiff = 0; iDiff < nDiffStates; ++iDiff) {
      if (nucleus.provideVals().smoothVelocityGuards.at(iDiff) > 0) {
         const auto logVelGuard = std::log(nucleus.provideVals().smoothVelocityGuards.at(iDiff));
         minMaxSpeed = std::min(minMaxSpeed, -nucleus.provideVals().diffStateVals.at(iDiff) * logVelGuard * std::sqrt(aParams.mass / findVal(iDiff)) / aParams.sensorLagTime);
      }
      if (nucleus.provideVals().smoothDensityGuards.at(iDiff) > 0) {
         const auto logDensGuard = std::log(nucleus.provideVals().smoothDensityGuards.at(iDiff));
         maxVelGradDens = std::max(
               maxVelGradDens, 
               desiredVel.mag() * nucleus.provideVals().diffStateVals.at(iDiff) * std::sqrt(aParams.mass / findVal(iDiff)) * findGrad(iDiff) / logDensGuard,
               [](const Point& p1, const Point& p2) { return p1.mag() < p2.mag(); }
               );
      }
   }
   const auto adjustedForGradDens = desiredVel + maxVelGradDens;
   //const auto mag = adjustedForGradDens.mag();
   //const auto newv = mag == 0.0 || minMaxSpeed == 0.0
   //   ? Point{}
   //   : minMaxSpeed * std::tanh(mag / minMaxSpeed) * adjustedForGradDens / mag;

   const auto newv = minMaxSpeed < adjustedForGradDens.mag() ? minMaxSpeed * adjustedForGradDens.normed() : adjustedForGradDens;
   //if (std::rand() % 100000 == 0 && adjustedForGradDens != Point{})
   //   std::cerr << newv.mag() / adjustedForGradDens.mag() << std::endl;

   return newv;
}


void Agent::setSpeedsFromLimits(const Point& desiredVel, double dt) {
   const auto turningDuration = std::max(aParams.maxSettleTime, aParams.sensorLagTime);
   auto angularSpeedLimit = std::numeric_limits<double>::max();
   if (turningDuration > 0.0) {
      const auto desiredHeadingChange = std::abs(desiredVel.angle());
      const auto angularSpeedLimitDueToSensorLag = desiredHeadingChange / turningDuration;
      //std::cerr << "desired change: " << desiredHeadingChange << " ; turndur: " << turningDuration << " ; senslag: " << aParams.sensorLagTime << " ; maxSettle: " << aParams.maxSettleTime << std::endl;
      angularSpeedLimit = std::min(angularSpeedLimit, angularSpeedLimitDueToSensorLag);
   }
   if (aParams.maxAngularSpeed > 0.0) {
      if (aParams.maxAngularSpeed < angularSpeedLimit) {
         angularSpeedLimit = aParams.maxAngularSpeed;
         std::cerr << "Angular velocitiy limited beyond sensor-time consideration" << std::endl; // DEBUG
      }
   }
   const auto desiredAngularVel = desiredVel.angle() / dt;
   angularVel = std::clamp(desiredAngularVel, -angularSpeedLimit, angularSpeedLimit);
   //std::cerr << "angvel: " << angularVel << " ; anglim: " << angularSpeedLimit << std::endl;

   if (std::abs(desiredVel.angle()) < M_PI / 2.0) {
      const auto expectedAngleAfterDt = angularVel * dt;
      const auto expectedAngleFromDesired = desiredVel.angle() - expectedAngleAfterDt;
      const auto speedConsideringAngle = std::cos(expectedAngleFromDesired) * desiredVel.mag();
      speed = limitedSpeed(speedConsideringAngle, aParams.maxSpeed);
   } else {
      speed = 0.0;
   }
   //speed = limitedSpeed(desiredVel.mag(), aParams.maxSpeed);
}


double Agent::rawVal(int isub) const { 
   return std::accumulate(sensors.at(isub).begin(), sensors.at(isub).end(), 0.0) / sensors.at(isub).size(); 
}


double Agent::speedValueCorrection(int isub) const {
   return self(isub, prod.at(isub), 0.0).val - self(isub, prod.at(isub), speed).val;
}


double Agent::offspringValueCorrection(int isub) const {
   return isOffspring
      ? std::exp(-subs.at(isub).decayRate * lifetime) * self(isub, prod.at(isub), speed).val
      : 0.0;
}


double Agent::findVal(int isub) const { 
   return rawVal(isub)
#ifndef NOCALIBVAL
      + speedValueCorrection(isub) 
#endif
#ifndef NOLIFECORR
      + offspringValueCorrection(isub)
#endif
      ; 
}


Point Agent::findVal(int isubx, int isuby) const { 
   return {findVal(isubx), findVal(isuby)}; 
}


Point Agent::rawGrad(SensorVals sv, bool debug) const { 
   return {(sv.at(0) - sv.at(2)) / aParams.sensordistMajor, (sv.at(1) - sv.at(3)) / aParams.sensordistMinor}; 
}


Point Agent::findGrad(int isub) const { 
#ifdef NOCALIBGRAD
   return rawGrad(sensors.at(isub)); 
#else
   return rawGrad(sensors.at(isub)) - self(isub, prod.at(isub), speed).grad; 
#endif
}


double Agent::velDotGradUpstream(int isub, const Point& vel) const {
   if (isub < nDiffStates || isub >= nDiffStates + nFieldVars) throw std::invalid_argument("Agent::velDotGradUpstream() only makes sense for field variables.");
   const auto dens = findVal(diffVarForFieldVar(isub));
   if (dens <= 0.0) {
      std::cerr << "Warning: in velDotGradUpstream, 0.0 measured for differentation state field " << diffVarForFieldVar(isub) << std::endl;
      return vel.dot(findGrad(isub)); // return uncorrected if we can't get the differentiation state reading needed for upstream
   } else {
      const auto numberDensity = dens / aParams.mass;
      const auto interparticleDist = 1.0 / std::sqrt(numberDensity);
      const auto secondDirectional = 0.5 * findLaplacian(isub);
      return vel.dot(findGrad(isub)) + 0.5 * vel.mag() * secondDirectional * interparticleDist;
   }
}


double Agent::findLaplacian(int isub, bool debug) const {
   const int iFieldVal = isub - nDiffStates;
   if (iFieldVal < 0 || iFieldVal >= nFieldVars) throw std::invalid_argument("Agent::findLaplacian() only makes sense for field variables.");

   const auto rawSubtraction = findVal(isub) - fieldValsTimesState.at(iFieldVal);

#ifdef NOLAPLCORR
   const auto finalSubtraction = rawSubtraction;
#else
   const auto potentialCorrectedSubtraction = rawSubtraction + laplCorrection.at(iFieldVal);
   const auto wouldCrossZero = (potentialCorrectedSubtraction > 0) != (rawSubtraction > 0);
   const auto finalSubtraction = wouldCrossZero ? 0.0 : potentialCorrectedSubtraction;
   if (debug)
      std::cout << " ... raw: " << findVal(isub) << " - " << fieldValsTimesState.at(iFieldVal) << ", corr: " << laplCorrection.at(iFieldVal) << ", sub: " << potentialCorrectedSubtraction << ", final: " << finalSubtraction << std::endl;
#endif

   return finalSubtraction * subs.at(isub).decayRate / subs.at(isub).difRate;
}


Point Agent::findLaplacian(int isubx, int isuby) const { 
   return {findLaplacian(isubx), findLaplacian(isuby)}; 
}


double Agent::findDiv(int isubx, int isuby, Point refDirX) const {
   return findGrad(isubx).scalarProjection(refDirX) 
      + findGrad(isuby).scalarProjection(refDirX.rotated(M_PI / 2));
}


double Agent::edgeDetect(int isub) const {
   if (isub >= nDiffStates) throw std::logic_error("edgeDetect() only makes sense for substances representing differentiation states.");
   if (subs.at(isub).difRate == 0.0) throw std::runtime_error("edgeDetect() called on substance with zero diffusion rate.");
   if (findVal(isub) == 0.0) return 0.0; // Assume it's not an edge if we can't detect any morphogen.
   const auto edgeToTest = findGrad(isub).mag() / findVal(isub);
   const auto edgeExpected = std::sqrt(subs.at(isub).decayRate / subs.at(isub).difRate);
   return edgeToTest / edgeExpected;
}


double Agent::colorEdge(int isub) const {
   if (isub < nDiffStates || isub >= nDiffStates + nFieldVars) throw std::logic_error("colorEdge() needs to work on a color field, which must be a regular abstract field.");
   if (subs.at(isub).difRate == 0.0) throw std::runtime_error("colorEdge() called on substance with zero diffusion rate.");
   if (findVal(isub) == 0.0) return 0.0; // Assume it's not an edge if we can't detect any morphogen.
   const auto edgeToTest = findGrad(isub).mag();
   const auto edgeExpected = std::sqrt(subs.at(isub).decayRate / subs.at(isub).difRate);
   return edgeToTest / edgeExpected;
}


Agent::Agent(
      AgentParams aParams,
      Substances subs, 
      ManagerInterface mi,
      NuclearOutput nucOut
      ) : 
   aParams{aParams},
   subs{subs},
   isOffspring{false},
   mi{mi},
   nucleus{AgentInterface{*this}, nucOut}
{ 
   calcFieldValsTimesState();
}


Agent::Agent(
      const Agent& a,
      ManagerInterface mi
      ) : 
   aParams{a.aParams},
   subs{a.subs},
   isOffspring{true},
   parentalTargetDensities{a.targetDensities},
   parentalRecentDensities{a.provideDens()},
   mi{mi},
   nucleus{AgentInterface{*this}, a.nucleus},
   fieldValsTimesState{a.fieldValsTimesState},
   laplCorrection{a.laplCorrection},
   leakyDensityIntegrals{a.leakyDensityIntegrals},
   lambdas{a.lambdas}
{ }


Agent::Agent(
      ManagerInterface mi,
      std::istream& is
      ) :
   aParams{readObjectFromStream<AgentParams>(is)},
   subs{readObjectFromStream<Substances>(is)},
   isOffspring{static_cast<bool>(is.get())},
   parentalTargetDensities{readObjectFromStream<DiffStateVals>(is)},
   parentalRecentDensities{readObjectFromStream<DiffStateVals>(is)},
   mi{mi},
   nucleus{AgentInterface{*this}, is},
   sensors{readObjectFromStream<std::array<SensorVals, nSubs> >(is)},
   gravityDir{readObjectFromStream<Point>(is)},
   fieldValsTimesState{readObjectFromStream<FieldVals>(is)},
   laplCorrection{readObjectFromStream<FieldVals>(is)},
   leakyDensityIntegrals{readObjectFromStream<DiffStateVals>(is)},
   lambdas{readObjectFromStream<DiffStateVals>(is)},
   lifetime{readObjectFromStream<double>(is)},
   prod{readObjectFromStream<SPHVals>(is)},
   speed{readObjectFromStream<double>(is)},
   angularVel{readObjectFromStream<double>(is)},
   prevGrowthSwitches{readObjectFromStream<std::array<bool, nDiffStates> >(is)},
   targetDensities{readObjectFromStream<DiffStateVals>(is)}
   //accMem{readObjectFromStream<Point>(is)}
{ }


// At the moment all these use the same dt, but the agent could decide to sense, update eqs, produce, and move
// all on different time intervals if it wanted
double Agent::update(double dt, bool debug) {

   lifetime += dt;
   
   // TODO maybe split these into separate functions (sense, think, move, produce)?

   // *** SENSE ***
   for (int isub = 0; isub < nSubs; ++isub) interpretSensors(isub);
   gravityDir = mi.senseGravityDir();

   // *** THINK ***
   FieldVals oldFieldValsTimesState{fieldValsTimesState};
   nucleus.update(dt, debug);
   const auto desiredVel = debugflow
      ? aParams.maxSpeed * gravityDir.rotated(M_PI / 2.0)
      : aParams.speedupFactor * smoothnessAdjustments(nucleus.provideVals().velocity);
   calcFieldValsTimesState();
   for (int iFieldVar = 0; iFieldVar < nFieldVars; ++iFieldVar) {
      laplCorrection.at(iFieldVar) += fieldValsTimesState.at(iFieldVar) - oldFieldValsTimesState.at(iFieldVar) 
         - subs.at(iFieldVar + nDiffStates).decayRate * laplCorrection.at(iFieldVar) * dt;
   }
   //accMem = nucleus.provideAcc(); //DEBUG
   //Point newVel = Point{speed, 0.0} + dt * nucleus.provideAcc(); // THIS ONE WORKS (note speed is still old speed at this point)

   // Figure actual speed and heading according to max speed and max angular speed
   if (!debugspin) setSpeedsFromLimits(desiredVel, dt);

   if (debugspin) {
      static auto flag = true;
      setSpeedsFromLimits(Point{1.0, 0.0}.rotated(M_PI * 1.0 / 2.0), dt);
      if (flag) {
         std::cerr << "angular velocity: " << angularVel << " ; speed: " << speed << std::endl;
         flag = false;
      }
   }

   // *** MOVE ***
   mi.move(speed, angularVel, std::min(desiredVel.mag(), aParams.maxSpeed), dt);

   // *** PRODUCE ***
   for (int isub = 0; isub < nSPHSubs; ++isub) {
      prod.at(isub) = valOverRho(isub) * subs.at(isub).gain * subs.at(isub).decayRate * aParams.mass;
      if (!debugflow) mi.produce(isub, prod.at(isub), dt);
   }
   
#ifdef DERIVCTRL
   return findLambdaProportional(dt);
#else
   return findLambdaIntegral(dt);
#endif
}


void Agent::printCalibDat() {
   for (const auto& rec : cdat) {
      for (size_t i = 0; i < rec.second.size(); ++i) {
         std::cout << "Speed: " << rec.first << ";   Substance: " << i << ";   Value: " << rec.second.at(i).val << ";   Gradient: " << rec.second.at(i).grad.x << ", " << rec.second.at(i).grad.y << std::endl;
      }
   }
}


void Agent::burnin(int isub, double dt) {
   if (isub >= nDiffStates) interpretSensors(diffVarForFieldVar(isub));
   mi.produce(isub, valOverRho(isub) * subs.at(isub).gain * subs.at(isub).decayRate * aParams.mass, dt);
}


// This allows the World object to compute a target field sum for burn-in.
double Agent::integral(int isub) {
   if (isub >= nSPHSubs) throw std::invalid_argument("Agent::integral() only makes sense for SPH variables.");
   if (isub >= nDiffStates) interpretSensors(diffVarForFieldVar(isub));
   const auto integral = valOverRho(isub) * subs.at(isub).gain * aParams.mass;

   // DEBUG
   if (!std::isfinite(integral)) {
      std::cerr << "Warning: in Agent::integral, integral not finite." << std::endl;
      std::cerr << "Sensors: ";
      for (auto e : sensors.at(isub)) std::cerr << e << ", ";
      std::cerr << std::endl;
   }
   if (std::isnan(integral) == true) std::cerr << "Integral tests as isnan" << std::endl;

   return integral;
}


// Order should be same as for update()
CalibBySub Agent::calibrate(Point force, double multiplier, double dt, bool debug) {
   for (int isub = 0; isub < nSPHSubs; ++isub)
      if (subs.at(isub).include) interpretSensors(isub);
   mi.move(force.mag(), force.angle(), 0.0, dt);
   CalibBySub calibBySub{};
   for (int isub = 0; isub < nSPHSubs; ++isub) {
      if (subs.at(isub).include) {
         mi.produce(isub, 1.0 * multiplier, dt);
         calibBySub.at(isub) = CalibRecord{rawVal(isub), rawGrad(sensors.at(isub), debug)};
      }
   }
   return calibBySub;
}


DiffStateVals Agent::provideDens() const {
   DiffStateVals dens;
   for (int iDiff = 0; iDiff < nDiffStates; ++iDiff) {
      dens.at(iDiff) = findVal(iDiff);
   }
   return dens;
}


void Agent::serialize(std::ostream& os) const {
   os.write(reinterpret_cast<const char *>(&aParams), sizeof(aParams));
   os.write(reinterpret_cast<const char *>(&subs), sizeof(subs));
   os.put(static_cast<char>(isOffspring));
   os.write(reinterpret_cast<const char *>(&parentalTargetDensities), sizeof(parentalTargetDensities));
   os.write(reinterpret_cast<const char *>(&parentalRecentDensities), sizeof(parentalRecentDensities));
   nucleus.serialize(os);
   os.write(reinterpret_cast<const char *>(&sensors), sizeof(sensors));
   os.write(reinterpret_cast<const char *>(&gravityDir), sizeof(gravityDir));
   os.write(reinterpret_cast<const char *>(&fieldValsTimesState), sizeof(fieldValsTimesState));
   os.write(reinterpret_cast<const char *>(&laplCorrection), sizeof(laplCorrection));
   os.write(reinterpret_cast<const char *>(&leakyDensityIntegrals), sizeof(leakyDensityIntegrals));
   os.write(reinterpret_cast<const char *>(&lambdas), sizeof(lambdas));
   os.write(reinterpret_cast<const char *>(&lifetime), sizeof(lifetime));
   os.write(reinterpret_cast<const char *>(&prod), sizeof(prod));
   os.write(reinterpret_cast<const char *>(&speed), sizeof(speed));
   os.write(reinterpret_cast<const char *>(&angularVel), sizeof(angularVel));
   os.write(reinterpret_cast<const char *>(&prevGrowthSwitches), sizeof(prevGrowthSwitches));
   os.write(reinterpret_cast<const char *>(&targetDensities), sizeof(targetDensities));
   //os.write(reinterpret_cast<const char *>(&accMem), sizeof(accMem));
}


CalibData::iterator findAbove(double speed) {
   auto above = Agent::cdat.upper_bound(speed);

   // Special cases to adjust "above" when speed is not between two calibration speeds,
   // so we can do a linear extrapolation beyond calibration data 
   if (above == Agent::cdat.end()) {
      above = std::prev(above); // --foo.end() would have been unsafe
      std::cerr << "Warning: speed (" << speed << ") above max used in calibration (" << above->first << ").\n";
   } else if (above == Agent::cdat.begin()) {
      std::cerr << "Warning: speed (" << speed << ") below min used in calibration (" << above->first << ").\n";
      above = std::next(above);
   }

   return above;
}


double speedProp(CalibData::iterator above, CalibData::iterator below, double speed, bool debug) {
   if (debug) std::cerr << speed << " " << Agent::cdat.size() << "\n"; // DEBUG
   if (debug) std::cerr << below->first << " " << above->first << "\n"; // DEBUG
   if (debug) std::cerr << (speed - below->first) / (above->first - below->first) << "\n"; // DEBUG

   return (speed - below->first) / (above->first - below->first);
}


// correct for speed by linear interpolation from table
CalibRecord self(const int isub, const double prod, const double speed, const bool debug) {
   // Special cases for small amounts of calibration data or if speed is exact match
   if (Agent::cdat.empty()) return {0.0, {0.0, 0.0}};
   if (Agent::cdat.size() == 1) return prod * Agent::cdat.begin()->second.at(isub);
   const auto match = Agent::cdat.find(speed);
   if (match != Agent::cdat.end()) return prod * match->second.at(isub);

   const auto above = findAbove(speed);
   const auto below = std::prev(above);

   // Linear interpolation (or extrapolation)
   if (debug) {
      std::cerr << "grad: " << below->second.at(isub).grad << ";  " << above->second.at(isub).grad << ";  " << prod << "\n"; // DEBUG
      std::cerr << "val: " << below->second.at(isub).val << ";  " << above->second.at(isub).val << ";  " << prod << "\n"; // DEBUG
   }
   return prod * (speedProp(above, below, speed, debug) * (above->second.at(isub) - below->second.at(isub)) + below->second.at(isub));
   //auto rv = prod * (speedProp(above, below, speed, debug) * (above->second.at(isub) - below->second.at(isub)) + below->second.at(isub));
   //rv.grad.y = 0.0;
   //return rv;
}


double limitedSpeed(double desired, double max) {
   if (max > 0.0 && desired - max > std::numeric_limits<double>::epsilon()) {
      std::cerr << "Speed requested: " << desired << "; max speed: " << max << ". Difference of " << desired - max << std::endl;
      throw MaxSpeedExceeded("Maximum speed exceeded");
   }
   return desired;

   //const auto limited = desired > max + std::numeric_limits<double>::epsilon() && max > 0.0 ? max : desired;
   //if (limited < desired) {
   //   std::cerr << "Limited speed from " << desired << " to " << max << ". Difference of " << desired - max << std::endl;
   //   throw MaxSpeedExceeded("Maximum speed exceeded");
   //}
   //mi.report(1); // TODO if needed again, just pass in mi as a reference
   //return limited;
}
