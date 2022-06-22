// calibrator.cpp
// Allen McBride
// June 21, 2022
//
// Descriptions in calibrator.h

#include "calibrator.h"

#include "agentmanager.h"
#include "diffuser.h"
#include "fields.h"
#include "picture.h"
#include "screen.h"
#include "substances.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <sys/file.h>
#include <unistd.h>
#include <valarray>
#include <vector>


std::map<double, CalibRecord> readFile(std::string filename);
void replaceSpeedIfMoreSamples(std::string filename, double speed, int nsamplesNew, std::string newLine);
void removeSpeed(std::string filename, double speedToRemove);
void appendStringToFile(std::string toAppend, std::string filename);
std::array<int, nSPHSubs> findDuplicates(Substances& subs);
CalibBySub calibBySubInvalid();
Substances subsNeededOnly(Substances subs, const CalibData& cdat, double speed);
double findDistToWall(const Substances& subs, double semimajorAxis);
double guessMaxDuration(const Substances& subs, double speed, double propOfOriginalToAllow, bool passthrough);
SampleStats calibStats(const std::vector<CalibSample>& samples);
double timeToDiffuseBeyond(double radius, double viscosity, double proportion);


Calibrator::Integrator::Integrator(const Substances& subs) : 
   subs{subs}
{
   int count = 0;
   for (size_t isub = 0; isub < nSPHSubs; ++isub)
      if (subs.at(isub).include) ++count;
   lastIntegral.resize(count, 0.0);
   subsDecayRates.resize(count, 0.0);
   candidateCutoffTimes.resize(count, 0.0);
   int i = 0;
   for (const auto& s : subs) 
      if (s.include) subsDecayRates[i++] = s.decayRate;
}


bool Calibrator::Integrator::addAndCheck(CalibBySub calibBySub, double t, double dt, bool debug) {
   if (Calibrator::multOfSteadyTimeToLook < 1.0) throw std::invalid_argument("multOfSteadyTimeToLook must be at least 1.0");
   CalibValBySub current(0.0, subsDecayRates.size());
   int i = 0;
   for (size_t isub = 0; isub < nSPHSubs; ++isub) 
      if (subs.at(isub).include)
         current[i++] = calibBySub.at(isub).val;
   lastIntegral += current * dt;
   calibIntegrals[t] = lastIntegral;
   const double candidateCutoffTime = t / Calibrator::multOfSteadyTimeToLook;
   const bool pulseAlreadyRegistered = current.min() > 0.0;
   bool allRatiosSmallEnough = false;
   if (pulseAlreadyRegistered) {
      const auto candidateRec = calibIntegrals.lower_bound(candidateCutoffTime);
      const CalibValBySub integralCurrentToInfApprox = current / subsDecayRates;
      const CalibValBySub integralZeroToInfApprox = lastIntegral + integralCurrentToInfApprox;
      const CalibValBySub integralCandidateToInfApprox = integralZeroToInfApprox - candidateRec->second;
      const CalibValBySub ratio = integralCandidateToInfApprox / integralZeroToInfApprox;
      candidateCutoffTimes[ratio >= Calibrator::propOfOriginalToAllow] = candidateCutoffTime;
      allRatiosSmallEnough = ratio.max() < Calibrator::propOfOriginalToAllow;
      if (debug) std::cerr << current[0] <<std::endl;
      if (debug) std::cerr << integralCandidateToInfApprox[0] << "; " << integralZeroToInfApprox[0] << std::endl;
      if (debug) std::cerr << "max ratio: " << ratio.max() << "; waiting for: " << Calibrator::propOfOriginalToAllow << std::endl;
   }
   return pulseAlreadyRegistered && allRatiosSmallEnough;
}


std::array<double, nSPHSubs> Calibrator::Integrator::provideCutoffTimes() const {
   std::array<double, nSPHSubs> cutoffTimesBySub{};
   int i = 0;
   for (size_t isub = 0; isub < nSPHSubs; ++isub)
      if (subs.at(isub).include)
         cutoffTimesBySub.at(isub) = candidateCutoffTimes[i++];
   return cutoffTimesBySub;
}


double Calibrator::calcYphys(const Substances& subsForSpeed, double speed, double duration) const {
   const double distToWall = findDistToWall(subsForSpeed, amParams.semimajorAxis);
   std::cerr << "distToWall: " << distToWall << "; 3*sqrt(semimajor): " << 3.0 * std::sqrt(amParams.semimajorAxis) << std::endl;
   return 2.0 * distToWall + speed * duration + 2.0 * (amParams.semimajorAxis + amParams.stalkLen);
}


AMExternalState Calibrator::calcStartState(double speed, Fparams fparams, double duration) const {
   AMExternalState startState;
   std::uniform_real_distribution<double> randUnit{0.0, 1.0};
   const double distanceOverDuration = duration * speed;
   const double macroOffsetDist = distanceOverDuration / 2.0;
   const Point center{fparams.xphys / 2.0, fparams.yphys / 2.0};
   const Point microOffset{(randUnit(gen) - 0.5) * fparams.hx, (randUnit(gen) - 0.5) * fparams.hy};
   startState.ori = 2 * M_PI * randUnit(gen);
   const auto macroOffset = Point{macroOffsetDist, 0.0}.rotated(startState.ori + M_PI);
   startState.pos = center + microOffset + macroOffset;
   return startState;
};


std::string Calibrator::calcFilename(int isub, double speed) const {
   std::ostringstream fnss;
   fnss.precision(4);
   fnss << "calibdat/calibdat";
   if (speed >= 0.0) 
      fnss << "samp-sp" << speed;
   if (isub >= 0)
      fnss << "-dif" << subs.at(isub).difRate
         << "-dec" << subs.at(isub).decayRate;
   fnss << "-rad" << amParams.semimajorAxis;
   if (amParams.semimajorAxis != amParams.semiminorAxis)
      fnss << "-min" << amParams.semiminorAxis;
   fnss << "-stk" << amParams.stalkLen 
      << "-visc" << visc
      //<< "-h" << gridh // TODO: Remember to uncomment before more test calibrations
      //<< "-tpc" << maxTimestepsPerCell // TODO: Remember to uncomment before more test calibrations
      //<< "-dt" << requestedDt
      << (passthrough ? "-pass" : "-nopass")
      << ".txt";
   return fnss.str();
}


void Calibrator::deleteOldData(const std::set<double>& speeds, bool onecalib) const {
   for (int isub = 0; isub < nSPHSubs; ++isub) {
      if (subs.at(isub).include) {
         const std::string filename{calcFilename(isub)};
         if (onecalib) {
            if (speeds.size() != 1) throw std::logic_error("onecalib set, but size of speeds set not 1.");
            removeSpeed(filename, *speeds.cbegin());
         } else {
            std::cerr << "Recalibrating. Old file " << filename;
            if (std::filesystem::remove(filename))
               std::cerr << " found and removed." << std::endl;
            else
               std::cerr << " not found." << std::endl;
         }
         for (double speed : speeds) {
            const std::string speedfile{calcFilename(isub, speed)};
            std::cerr << "File for calibration samples " << speedfile;
            if (std::filesystem::remove(speedfile))
               std::cerr << " found and removed." << std::endl;
            else
               std::cerr << " not found." << std::endl;
         }
      }
   }
}


// Read in any existing calibration data. Does not check whether speed read is in speeds list
std::map<double, CalibRecord> Calibrator::readFile(std::string filename) const {
   std::map<double, CalibRecord> cdatSub;

   std::ifstream prevDat{filename};
   if (!prevDat) {
      std::cerr << filename << " cannot be opened to read saved calibration data. Assuming it does not exist." << std::endl;
   } else {
      // Each line of file should be a calibration record, space delimited, in form "speed grad.x grad.y <optional comments...>"
      FILE* toLock = std::fopen(filename.c_str(), "r");
      flock(fileno(toLock), LOCK_EX);
      std::cerr << "Found " << filename << ". Attempting to read data." << std::endl;
      std::string line;
      while (std::getline(prevDat, line)) {
         
         //std::cerr << "Examining line; " << line << std::endl; // DEBUG

         // Prepare to read a line. Assume it's good until found otherwise.
         bool goodLine = true;

         std::istringstream liness(line);

         // Attempt to read speed
         double speed;
         if (!(liness >> speed)) {
            std::cerr << "Failed to read speed in line: " << line << std::endl;
            goodLine = false;
         }

         CalibRecord cr;
         if (!(liness >> cr.val)) {
            std::cerr << "Failed to read val in line: " << line << std::endl;
            goodLine = false;
         }
         if (!(liness >> cr.grad.x)) {
            std::cerr << "Failed to read grad.x in line: " << line << std::endl;
            goodLine = false;
         }
         if (!(liness >> cr.grad.y)) {
            std::cerr << "Failed to read grad.y in line: " << line << std::endl;
            goodLine = false;
         }

         size_t nsamples;
         if (!(liness >> nsamples)) {
            std::cerr << "Failed to read nsamples in line: " << line << std::endl;
            goodLine = false;
         }
         if (nsamples < nsamplesNeeded(speed)) {
            std::cerr << nsamples << " is less than the " << nsamplesNeeded(speed) << " required in line: " << line << std::endl;
            goodLine = false;
         }

         // If successful, add to map
         if (goodLine) cdatSub.emplace(speed, cr);
      }
      fclose(toLock);
   }

   return cdatSub;
}


void Calibrator::processSamples(double speed) const {
   for (int isub = 0; isub < nSPHSubs; ++isub) {
      if (subs.at(isub).include) {
         std::ifstream speedFile{calcFilename(isub, speed)};
         if (!speedFile) {
            std::cerr << calcFilename(isub, speed) << " cannot be opened to read saved calibration sample data. Assuming it does not exist." << std::endl;
         } else {
            std::vector<CalibSample> calSamples;
            std::vector<double> durationsOneSpeed;
            std::cerr << "Found " << calcFilename(isub, speed) << ". Attempting to read calibration sample data." << std::endl;
            std::string line;
            while (std::getline(speedFile, line)) {
               if (line.compare(0, 6, "pulse ") != 0) {
                  std::istringstream sampless{line};
                  CalibSample sample{};

                  bool goodLine = true;
                  if (!(sampless >> sample.record.val)) {
                     std::cerr << "Failed to read calibration value in line: " << line << std::endl;
                     goodLine = false;
                  }
                  if (!(sampless >> sample.record.grad.x)) {
                     std::cerr << "Failed to read x-component of calibration gradient in line: " << line << std::endl;
                     goodLine = false;
                  }
                  if (!(sampless >> sample.record.grad.y)) {
                     std::cerr << "Failed to read y-component of calibration gradient in line: " << line << std::endl;
                     goodLine = false;
                  }
                  if (!(sampless >> sample.duration)) {
                     std::cerr << "Failed to read calibration duration in line: " << line << std::endl;
                     goodLine = false;
                  }
                  if (!(sampless >> sample.yphys)) {
                     std::cerr << "Failed to read domain physical size in line: " << line << std::endl;
                     goodLine = false;
                  }
                  if (!(sampless >> sample.dt)) {
                     std::cerr << "Failed to read timestep in line: " << line << std::endl;
                     goodLine = false;
                  }

                  if (goodLine) calSamples.push_back(sample);
               }
            }

            if (calSamples.size() >= nsamplesNeeded(speed)) {
               const auto sampleStats = calibStats(calSamples);
               std::ostringstream newDatStr;
               // Write calibration data and (as comment) parameters that don't invalidate comparisons between files
               newDatStr << speed << " " << sampleStats.mean.val << " " << sampleStats.mean.grad.x << " " << sampleStats.mean.grad.y << " " << calSamples.size()
                  << " // dt: " << sampleStats.meanDt << "; duration: " << sampleStats.meanDuration << "; h: " << gridh << "; width: " << sampleStats.meanYphys 
                  << "; stdDev: " << sampleStats.stdDev << "; stdErr: " << sampleStats.stdErr << std::endl;
               replaceSpeedIfMoreSamples(calcFilename(isub), speed, calSamples.size(), newDatStr.str());
            }
         }
      }
   }
}


double Calibrator::findDuration(const Substances& subsForSpeed, const double speed) const {
   std::array<double, nSPHSubs> durationBySub{};
   auto subsNeedingDuration = subsForSpeed;
   bool someNeeded = true;
   int count = 0;
   while (someNeeded) {
      someNeeded = false;
      for (int isub = 0; isub < nSPHSubs; ++isub) {
         if (subsForSpeed.at(isub).include) {
            std::ifstream speedFile{calcFilename(isub, speed)};
            if (!speedFile) {
               std::cerr << calcFilename(isub, speed) << " cannot be opened to read saved pulse duration data. Assuming it does not exist." << std::endl;
            } else {
               std::vector<double> durationsOneSpeed;
               std::cerr << "Found " << calcFilename(isub, speed) << ". Attempting to read pulse duration data." << std::endl;
               std::string line;
               while (std::getline(speedFile, line)) {
                  if (line.compare(0, 6, "pulse ") == 0) {
                     std::istringstream durationss{line.substr(6)};
                     double duration = 0;
                     if (!(durationss >> duration)) 
                        std::cerr << "Failed to read duration in line: " << line << std::endl;
                     durationsOneSpeed.push_back(duration);
                  }
               }
               if (durationsOneSpeed.size() == 0.0) {
                  std::cerr << "No durations found in file" << calcFilename(isub, speed) << "." << std::endl;
               } else {
                  durationBySub.at(isub) = std::accumulate(durationsOneSpeed.begin(), durationsOneSpeed.end(), 0.0) / durationsOneSpeed.size();
                  std::cerr << "Average pulse duration for sub " << isub << ", speed " << speed 
                     << ", based on " << durationsOneSpeed.size() << " samples: " << durationBySub.at(isub) << std::endl;
                  subsNeedingDuration.at(isub).include = false;
               }
            }
         }
         if (subsNeedingDuration.at(isub).include) {
            someNeeded = true;
            std::cerr << "Still need duration for substance " << isub << ", speed " << speed << std::endl;
         }
      }

      if (someNeeded) {
         std::cerr << "Using pulse to find duration for speed " << speed << std::endl;
         calcDurationWithPulse(subsNeedingDuration, speed);
         ++count;
      }
   }
   if (count > 1) std::cerr << "calcDurationWithPulse should not have needed to be called more than once; was called " << count << " times." << std::endl;
   return *std::max_element(durationBySub.begin(), durationBySub.end());
}


// This looks like it does nothing, but that's because calibOneSpeed() will write the relevant info to file
void Calibrator::calcDurationWithPulse(const Substances& subsForSpeed, const double speed) const {
   const double arbitraryMultiplier = passthrough ? 1.0 : 0.7; // guessMaxDuration is just a rough heuristic anyway
   const auto maxDuration = guessMaxDuration(subsForSpeed, speed, propOfOriginalToAllow, passthrough);
   auto pulseDuration = arbitraryMultiplier * maxDuration * multOfSteadyTimeToLook;
   bool steadyNotFound = true;
   //constexpr int maxTries = 10;
   int nTries = 0;
   bool fieldsTooBig = false;
   do {
      const double yphysPulse = calcYphys(subsForSpeed, speed, pulseDuration);
      const Fparams fparamsPulse{yphysPulse, yphysPulse, gridh};
      std::cerr << "Looking for steady state, trying duration: " << pulseDuration << std::endl;
      const double duration = calibOneSpeed(speed, subsForSpeed, fparamsPulse, pulseDuration, true); 
      steadyNotFound = duration < 0.0;
      if (steadyNotFound) {
         pulseDuration *= 1.5;
         ++nTries;
         std::cerr << "Time to steady state not found. Trying again with duration = " << pulseDuration << std::endl;
      } else {
         std::cerr << "Time to steady state: " << duration << std::endl;
      }
   } while (steadyNotFound && !fieldsTooBig);
   //} while (steadyNotFound && nTries < maxTries && !fieldsTooBig);
   if (steadyNotFound) {
      std::cerr << "Failed to find steady state after " << nTries << " tries for speed " << speed << ". Giving up." << std::endl;
      throw std::runtime_error("Failure to find steady state");
   }
}


// TODO If I want, I think I can just do what I do in calibrate() in the middle of the main loop here, maybe every 500 iterations or whatever:
//processSamples(speed);
//const auto subsOneSpeed = subsNeededOnly(subs, readCalibDataFromFile(speeds), speed);
//enoughSamples = true;
//for (const auto& s : subsOneSpeed) if (s.include) enoughSamples = false;
double Calibrator::calibOneSpeed(double speed, const Substances& subsForSpeed, Fparams fparams, double duration, bool pulse) const {
   constexpr bool obstacles = false;
   constexpr bool cues = false;
   const bool visual = liveSecondsPerUpdate >= 0.0;
   Fields cfields{fparams, passthrough, obstacles, cues, visual};
   Diffuser cdiffuser{cfields, subsForSpeed, passthrough, visc, viscosityTolerance, pressureTolerance};
   auto startState = calcStartState(speed, fparams, duration);
   AgentManager am{amParams, startState, cfields, subsForSpeed, SimulationParams{passthrough, -1.0, -1.0}, gen};
   std::optional<Screen> screen = visual ? std::make_optional<Screen>(cfields, SubstanceVars{5000.0}, pmag, calcFilename(-1, speed)) : std::nullopt; // -1 to suppress particular substance info
   //std::optional<Picture> pic = visual && !pulse ? std::make_optional<Picture>(cfields, SubstanceVars{5000.0}, calcFilename(0, speed).substr(0, calcFilename(0, speed).length() - 4)) : std::nullopt;
   std::optional<Picture> pic = std::nullopt;
   Integrator integrator{subsForSpeed};

   const double dt = speed == 0.0 ? cdiffuser.maxDt() : std::min(cdiffuser.maxDt(), gridh / (maxTimestepsPerCell * speed)); 
   const long long iterationsNeeded = fparams.xpixels * fparams.ypixels * duration / dt;
   if (iterationsNeeded > 3'000'000'000'000ll) {
      std::cerr << "calibOneSpeed would take too long. Speed: " << speed 
         << "; Iterations needed: " << iterationsNeeded
         << "; Pixels: " << fparams.xpixels << " by " << fparams.ypixels 
         << ". dt: " << dt << "; duration: " << duration << std::endl;
      //throw TooManyIterationsNeeded("Too many iterations would be needed by calibOneSpeed");
   }

   std::cerr << "calibOneSpeed. speed " << speed << ", yphys " << fparams.yphys << ", gridh " << gridh << ", hx " << fparams.hx << ", dt " << dt << ", duration " << duration <<std::endl;
   if (speed > 0.0) std::cerr << "timesteps per traversal: " << fparams.hx / (speed * dt) << std::endl;

   CalibBySub calibBySub{};
   double t = 0.0;
   bool pulseDone = false;
   auto startTime = std::chrono::steady_clock::now();
   while (!pulseDone && t < duration) {
      calibBySub = am.calibrate(Point{speed, 0.0}, (pulse && t > 0.0) ? 0.0 : 1.0, dt, !pulse);
      cfields.clearAgentInfo();
      am.assertAgentInfo();
      if (speed > 0.0) cdiffuser.advect(dt);
      cdiffuser.diffuse(dt);
      if (whetherToUpdate(t, liveSecondsPerUpdate, dt) && screen) {
         screen->update();
         if (pic) pic->update();
      }
      //cfields.clearTransientInfo();
      if (pulse) pulseDone = integrator.addAndCheck(calibBySub, t, dt, false);
      t += dt;
   }
   std::chrono::duration<double> calibDuration{std::chrono::steady_clock::now() - startTime};
   std::cerr << "Calibration run took " << calibDuration.count() << " seconds for " 
      << iterationsNeeded << " pixel iterations." << std::endl;
   
   // -1 will signify that duration wasn't long enough for pulse to be undetectable
   const auto cutoffTimes = integrator.provideCutoffTimes();
   const auto pulseDuration = pulseDone ? *std::max_element(cutoffTimes.begin(), cutoffTimes.end()) : -1.0;

   for (int isub = 0; isub < nSPHSubs; ++isub) {
      if (subsForSpeed.at(isub).include) {
         if (pulse) {
            if (pulseDone) appendStringToFile("pulse " + std::to_string(cutoffTimes.at(isub)) + "\n", calcFilename(isub, speed));
         } else {
            std::cerr << "Calibration sample, substance " << isub << ": " << calibBySub.at(isub) << std::endl;
            std::ostringstream sampStr;
            sampStr << calibBySub.at(isub).val << " " << calibBySub.at(isub).grad.x << " " << calibBySub.at(isub).grad.y 
               << " " << duration << " " << fparams.yphys << " " << dt << std::endl;
            appendStringToFile(sampStr.str(), calcFilename(isub, speed));
         }
      }
   }

   return pulseDuration;
}


Calibrator::Calibrator(
      Substances subs, 
      double visc,
      AgentManagerParams amParams,
      double nCellWidthsPerDiameter,
      double maxTimestepsPerCell,
      double liveSecondsPerUpdate,
      int pmag,
      bool passthrough,
      std::default_random_engine& gen
      ) :
   subs{subs}, 
   visc{visc}, 
   amParams{amParams},
   maxTimestepsPerCell{maxTimestepsPerCell},
   liveSecondsPerUpdate{liveSecondsPerUpdate},
   pmag{pmag},
   passthrough{passthrough},
   gen{gen},
   gridh{amParams.semiminorAxis * 2.0 / nCellWidthsPerDiameter},
   dupLocations{findDuplicates(this->subs)}
{
   if (nCellWidthsPerDiameter < 3.0)
      std::cerr << "*** Warning: " << nCellWidthsPerDiameter << " grid cell widths per agent diameter; minimum recommended is 3.0. *** " << std::endl;
}


// Generate calibration data for any speeds not found in files
CalibData Calibrator::calibrate(std::set<double> speeds, bool recalib, bool onecalib) const {

   if (std::filesystem::create_directory("calibdat", "."))
      std::cerr << "Directory calibdat did not exist; created." << std::endl;

   if (recalib) deleteOldData(speeds, onecalib);

   std::vector<double> speedsShuffled{speeds.begin(), speeds.end()};
   //std::random_device rd;
   //std::shuffle(speedsShuffled.begin(), speedsShuffled.end(), std::default_random_engine{rd()});
   std::shuffle(speedsShuffled.begin(), speedsShuffled.end(), gen);
   for (auto speed : speedsShuffled) {
      //std::cerr << std::endl << "Calibrating for speed " << speed << std::endl; // DEBUG

      bool enoughSamples = false;
      while (!enoughSamples) {
         processSamples(speed);
         const auto subsOneSpeed = subsNeededOnly(subs, readCalibDataFromFile(speeds), speed);
         enoughSamples = true;
         for (const auto& s : subsOneSpeed) if (s.include) enoughSamples = false;
         if (!enoughSamples) {
            const auto duration = findDuration(subsOneSpeed, speed);
            const auto yphys = calcYphys(subsOneSpeed, speed, duration);
            const Fparams fparams{yphys, yphys, gridh};

            std::cerr << "  semimajor: " << amParams.semimajorAxis << "; semiminor: " << amParams.semiminorAxis << "; duration: " << duration << "; trav: " << speed * duration << "; phys: " << fparams.yphys << "; px: " << fparams.ypixels << "; h: " << fparams.hx << "; hreq: " << gridh << std::endl;

            calibOneSpeed(speed, subsOneSpeed, fparams, duration, false);
         }
      }
   }

   return readCalibDataFromFile(speeds);
}


// The only significance of the speeds parameter is that if data is not found for some somestances for one of these speeds,
// the returned map will have invalid entries for missing data for that speed, so that calibrate() knows to calibrate for that
// data. For any other speeds, the only entries in the returned map will be complete for all substances.
CalibData Calibrator::readCalibDataFromFile(const std::set<double>& speeds) const {
   CalibData cdat;
   for (auto speed : speeds)
      cdat[speed] = calibBySubInvalid();

   for (int isub = 0; isub < nSPHSubs; ++isub) {
      if (subs.at(isub).include) {
         const std::string filename{calcFilename(isub)};
         for (const auto& crecForSpeed : readFile(filename)) {
            const auto speed = crecForSpeed.first;
            const auto entryForSpeed = cdat.find(speed);
            if (entryForSpeed == cdat.end()) { // If we haven't seen this speed before, 
               CalibBySub calibBySub{};
               calibBySub.fill({-1.0, {0.0, 0.0}}); // -1.0 to indicate invalid entry
               cdat[speed] = calibBySub;
               cdat[speed] = calibBySubInvalid();
            }
            cdat[speed].at(isub) = crecForSpeed.second;
            //std::cerr << "Found saved data for speed " << crecForSpeed.first << ", isub " << isub << ": " << crecForSpeed.second << std::endl;
         }
      } else {
         //std::cerr << "Copying saved data for substance " << dupLocations.at(isub) << " to substance " << isub << std::endl;
         if (dupLocations.at(isub) >= 0) { // Check if this substance is marked as a duplicate of another
            for (auto& cdatEntry : cdat) {
               if (cdatEntry.second.at(dupLocations.at(isub)).val >= 0.0) {
                  cdatEntry.second.at(isub) = cdatEntry.second.at(dupLocations.at(isub));
                  //std::cerr << "Copied saved data for speed " << cdatEntry.first << ", isub " << dupLocations.at(isub) << " to isub " << isub << std::endl;
               }
            }
            //for (auto speed : speeds) {
            //   if (cdat.at(speed).at(dupLocations.at(isub)).val >= 0.0) {
            //      cdat.at(speed).at(isub) = cdat.at(speed).at(dupLocations.at(isub));
            //      //std::cerr << "Copied saved data for speed " << speed << ", isub " << dupLocations.at(isub) << " to isub " << isub << std::endl;
            //   }
            //}
         }
      }
   }

   // Only include data found in files for speeds not specified in speeds set if it exists for all substances for a given speed.
   std::set<double> speedsToErase;
   for (const auto& cdatEntry : cdat) {
      if (speeds.find(cdatEntry.first) == speeds.end()) {
         bool validForAllSubs = true;
         for (const auto& crec : cdatEntry.second)
            if (crec.val < 0.0) validForAllSubs = false;
         if (!validForAllSubs) speedsToErase.insert(cdatEntry.first);
      }
   }
   for (auto speedToErase : speedsToErase) cdat.erase(speedToErase);

   return cdat;
}


void replaceSpeedIfMoreSamples(std::string filename, double speedToReplace, int nsamplesNew, std::string newLine) {
   FILE* datafile = std::fopen(filename.c_str(), "a+");
   if (!datafile) {
      std::cerr << "Error: " << filename << " cannot be opened to replace speed " << speedToReplace << " with line: " << newLine << std::endl;
   } else {
      flock(fileno(datafile), LOCK_EX);

      std::vector<std::string> lines;
      rewind(datafile);
      char * line = nullptr;
      size_t len = 0;
      int matchCount = 0;
      bool appendNewLine = true;
      while (getline(&line, &len, datafile) != -1) {
         std::istringstream liness{line};
         double speed = -1.0;
         bool replaceThisLine = false;
         if (liness >> speed && speed == speedToReplace) {
            ++matchCount;
            CalibRecord cr;
            if (liness >> cr.val && liness >> cr.grad.x && liness >> cr.grad.y) {
               int nsamples;
               if (liness >> nsamples) {
                  if (nsamplesNew > nsamples) {
                     std::cerr << "For speed " << speed << ", replacing average of " << nsamples << " with average of " << nsamplesNew << " samples." << std::endl;
                     replaceThisLine = true;
                  } else {
                     //std::cerr << "For speed " << speed << ", average of " << nsamplesNew << " not greater than existing average of " << nsamples << "." << std::endl;
                     appendNewLine = false;
                  }
               }
            }
         }
         if (!replaceThisLine) lines.emplace_back(line);
      }
      if (matchCount > 1) std::cerr << "Error: More than one matching speed record found." << std::endl;
      free(line);

      if (appendNewLine) {
         rewind(datafile);
         ftruncate(fileno(datafile), 0);
         lines.emplace_back(newLine);
         for (auto& e : lines) fputs(e.c_str(), datafile);
      }
      fclose(datafile);
   }
}


void removeSpeed(std::string filename, double speedToRemove) {
   FILE* datafile = std::fopen(filename.c_str(), "r+");
   if (!datafile) {
      std::cerr << filename << " cannot be opened to delete previous calibration entry. Assuming it does not exist." << std::endl;
   } else {
      flock(fileno(datafile), LOCK_EX);
      std::vector<std::string> lines;

      char * line = nullptr;
      size_t len = 0;
      while (getline(&line, &len, datafile) != -1) {
         std::istringstream liness{line};
         double speed;
         if (liness >> speed && speed == speedToRemove) {
            std::cerr << "Removing previous calibration entry for speed " << speed << " from " << filename << std::endl;
         } else {
            lines.emplace_back(line);
         }
      }
      free(line);

      rewind(datafile);
      ftruncate(fileno(datafile), 0);
      for (auto& e : lines) fputs(e.c_str(), datafile);
      fclose(datafile);
   }
}


void appendStringToFile(std::string toAppend, std::string filename) {
   FILE* file = std::fopen(filename.c_str(), "a");
   if (!file) {
      std::cerr << "Error: Failed to open file " << filename << ". Cannot append string:\n" << toAppend << std::endl;
   } else {
      flock(fileno(file), LOCK_EX);
      fputs(toAppend.c_str(), file);
      fclose(file);
   }
}


// findDuplicates will map duplicate locations but also set duplicates' "include" member to false
// higher-indexed substances will always point to themselves or to lower-indexed substances, which
// makes the array more convenient to use
std::array<int, nSPHSubs> findDuplicates(Substances& subs) {
   typedef std::pair<double, double> SubPhysParams;

   std::map<SubPhysParams, int> locationsByParams;
   std::array<int, nSPHSubs> dupLocations{};
   dupLocations.fill(-1); // indicate invalid entry, in case a Substance was already marked not-include
   for (int isub = 0; isub < nSPHSubs; ++isub) {
      if (subs.at(isub).include) {
         const SubPhysParams params(subs.at(isub).difRate, subs.at(isub).decayRate);
         const auto locationOfParams = locationsByParams.find(params);
         if (locationOfParams == locationsByParams.end()) {
            locationsByParams[params] = isub;
            dupLocations.at(isub) = isub;
            std::cerr << "No duplicate found for isub " << isub << std::endl;
         } else {
            subs.at(isub).include = false;
            dupLocations.at(isub) = locationOfParams->second;
            std::cerr << "Duplicate found. New isub " << isub << " duplicates old isub " << locationOfParams->second << std::endl;
         }
      }
   }
   return dupLocations;
}


CalibBySub calibBySubInvalid() {
   CalibBySub calibBySub{};
   calibBySub.fill({-1.0, {0.0, 0.0}}); // -1.0 to indicate invalid entry
   return calibBySub;
}


Substances subsNeededOnly(Substances subs, const CalibData& cdat, double speed) {
   for (int isub = 0; isub < nSPHSubs; ++isub) {
      if (cdat.at(speed).at(isub).val >= 0.0) {
         subs.at(isub).include = false;
         //std::cerr << "Don't need new calibration data for substance " << isub << std::endl; // DEBUG
      } else {
         if (subs.at(isub).include)
            std::cerr << "Need new calibration data for substance " << isub << std::endl;
         else
            std::cerr << "Substance " << isub << " either a duplicate or otherwise to be ignored." << std::endl;
      }
   }
   return subs;
}


// At what minimum distance beyond an agent must walls be such that the
// resulting boundary conditions do not have more than 1% effect on the
// readings at the agent's sensors?
// With no walls, the steady-state smoothing function is:
//    BesselK[0,r/h]/(2 Pi h^2).
// With walls at radius a, the steady-state smoothing function is:
//    i(BesselJ[1,ia/h]BesselY[0,-ir/h]+BesselJ[0,ir/h]BesselY[1,-ia/h])
//       / (4 h^2 BesselI[1,a/h]).
// The numbers below are the result of taking the ratio of the latter to
// the former, setting it equal to 1.01, and solving for a/h for various
// values of r/h, where r is the distance from producing to sensing
// devices. (I treat agents as circles here with their semimajor axis
// as their radius.  This results in a conservative overestimate of r
// in the case of ellipses.)
double findDistToWall(const Substances& subs, double semimajorAxis) {
   const auto circumference = 2.0 * M_PI * semimajorAxis;
   const auto distToSensor = circumference / 8.0;
   auto maxH = double{};
   for (int isub = 0; isub < nSPHSubs; ++isub) {
      if (subs.at(isub).include) {
         const auto smoothingH = std::sqrt(subs.at(isub).difRate / subs.at(isub).decayRate);
         maxH = std::max(maxH, smoothingH);
      }
   }
   const auto distOverH = distToSensor / maxH;
   auto aOverH = double{};
   if (distOverH < 1.0 / 128.0)
      aOverH = 2.26;
   else if (distOverH < 1.0 / 64.0)
      aOverH = 2.33;
   else if (distOverH < 1.0 / 32.0)
      aOverH = 2.41;
   else if (distOverH < 1.0 / 16.0)
      aOverH = 2.51;
   else if (distOverH < 1.0 / 8.0)
      aOverH = 2.63;
   else if (distOverH < 1.0 / 4.0)
      aOverH = 2.81;
   else if (distOverH < 1.0 / 2.0)
      aOverH = 3.07;
   else
      aOverH = 3.54;

   return aOverH * maxH;
}


// Just a rough heuristic
double guessMaxDuration(const Substances& subs, double speed, double propOfOriginalToAllow, bool passthrough) {
   double maxDuration = 0.0;
   for (int isub = 0; isub < nSPHSubs; ++isub) {
      const auto& s = subs.at(isub);
      if (s.include) {
         const auto durationGlobalSettle = s.timeToSteady(propOfOriginalToAllow);
         const auto durationOutrunGuess = speed == 0.0 ? std::numeric_limits<double>::max() : -4.0 * subs.at(isub).difRate * std::log(0.1 * propOfOriginalToAllow) / std::pow(speed, 2.0);
         const auto maxDurationForSub = std::min(durationGlobalSettle, durationOutrunGuess);
         maxDuration = std::max(maxDuration, maxDurationForSub);

         //std::cerr << "Distance to global steady state: " << durationGlobalSettle * speed << std::endl;
         //std::cerr << "Duration to global steady state: " << durationGlobalSettle << std::endl;
         //std::cerr << "Distance to idealized outrunning: " << (speed == 0.0 ? -1 : -4.0 * subs.at(isub).difRate * std::log(propOfOriginalToAllow) / std::pow(speed, 1.0)) << std::endl;
         //std::cerr << "Duration to idealized outrunning: " << (speed == 0.0 ? -1 : -4.0 * subs.at(isub).difRate * std::log(propOfOriginalToAllow) / std::pow(speed, 2.0)) << std::endl;
         //std::cerr << "Distance for speed at which two should be same: " << -2.0 * smoothingH * std::log(propOfOriginalToAllow) << std::endl;
         //std::cerr << "Duration for this speed to go above distance: " << (speed == 0.0 ? -1 : durationOutrunGuess) << std::endl;
      }
   }
   return maxDuration;
}


SampleStats calibStats(const std::vector<CalibSample>& samples) {
   if (samples.size() < 1) return {}; 

   const CalibSample avg = std::accumulate(samples.begin(), samples.end(), CalibSample{}) / samples.size();
   if (samples.size() < 2) return {avg.record, {}, {}, avg.duration, avg.yphys, avg.dt};

   CalibRecord sumSqErr{};
   for (const auto& samp : samples) {
      const CalibRecord err = samp.record - avg.record;
      sumSqErr.grad.x += err.grad.x * err.grad.x;
      sumSqErr.grad.y += err.grad.y * err.grad.y;
      sumSqErr.val += err.val * err.val;
   }
   const auto variance = sumSqErr / (samples.size() - 1);
   const CalibRecord stdDev = { std::sqrt(variance.val), { std::sqrt(variance.grad.x), std::sqrt(variance.grad.y) }};
   const auto stdErr = stdDev / std::sqrt(samples.size());
   //std::cerr << "stdDev: " << stdDev << ";  stdErr: " << stdErr << ";  N: " << samples.size() << ";  avg: " << avg.record << std::endl; // DEBUG
   return {avg.record, stdDev, stdErr, avg.duration, avg.yphys, avg.dt};
}


double timeToDiffuseBeyond(double radius, double viscosity, double proportion) {
   return -std::pow(radius, 2.0) / (4.0 * viscosity * std::log(proportion));
}
