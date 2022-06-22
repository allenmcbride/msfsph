// world.cpp
// Allen McBride
// April 19, 2022
//
// Descriptions in world.h

#include "world.h"
#include "agent.h"
#include "agentmanager.h"
#include "fields.h"
#include "substances.h"
#include "point.h"

#include <algorithm>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <list>
#include <random>
#include <set>
#include <string>


template<int n> void printArrayStats(const std::vector<std::array<double, n> >& a);
double advectionGaussian(const Point& point, const Fparams& fp);


#if PROGRAM == GROWTH
static constexpr bool obstaclesFlag = false;
static constexpr bool cuesFlag = false;
void World::setCues() {}
void World::placeObstacles() {}
NuclearOutput World::initVals(const Point& point) const {
   NuclearOutput initV{};
   initV.diffStateVals = {1.0};
   initV.fieldVals.at(FieldVars::color - nDiffStates) = 1.0;
   return initV;
}


#elif PROGRAM == GROWX
static constexpr bool obstaclesFlag = false;
static constexpr bool cuesFlag = false;
void World::setCues() {}
void World::placeObstacles() {}
NuclearOutput World::initVals(const Point& point) const {
   NuclearOutput initV{};
   initV.diffStateVals = {1.0, 10.0 - point.x > point.y, point.y > 10.0 - point.x};
   initV.fieldVals.at(FieldVars::colora - nDiffStates) = 1.0;
   initV.fieldVals.at(FieldVars::colorb - nDiffStates) = 1.0;
   return initV;
}


#elif PROGRAM == DIFFUSION || PROGRAM == DIFFUSION_TRANSLATION
static constexpr bool obstaclesFlag = false;
static constexpr bool cuesFlag = false;
void World::setCues() {}
void World::placeObstacles() {}
NuclearOutput World::initVals(const Point& point) const {
   return {};
}


#elif PROGRAM == PATH
static constexpr bool obstaclesFlag = true;
static constexpr bool cuesFlag = true;
double scale(double x) { return 1.0 * (x - 5.0) + 5.0; };
void World::setCues() {
   constexpr auto squareSize = 1.0;
   constexpr auto sourceLowCorner = Point{1.5, 4.5};
   constexpr auto targetLowCorner = Point{7.5, 4.5};
   for (int r = fields.fp.yToRow(scale(sourceLowCorner.y)); r < fields.fp.yToRow(scale(sourceLowCorner.y + squareSize)); ++r)
      for (int c = fields.fp.xToCol(scale(sourceLowCorner.x)); c < fields.fp.xToCol(scale(sourceLowCorner.x + squareSize)); ++c)
         fields.cues->at(r, c).set(Cues::source);
   for (int r = fields.fp.yToRow(scale(targetLowCorner.y)); r < fields.fp.yToRow(scale(targetLowCorner.y + squareSize)); ++r)
      for (int c = fields.fp.xToCol(scale(targetLowCorner.x)); c < fields.fp.xToCol(scale(targetLowCorner.x + squareSize)); ++c)
         fields.cues->at(r, c).set(Cues::target);
}
void World::placeObstacles() {
   constexpr auto barrierDimensions = Point{0.5, 6.0};
   constexpr auto barrierLowCorner = Point{4.75, 2.0};
   for (int r = fields.fp.yToRow(scale(barrierLowCorner.y)); r < fields.fp.yToRow(scale((barrierLowCorner + barrierDimensions).y)); ++r)
      for (int c = fields.fp.xToCol(scale(barrierLowCorner.x)); c < fields.fp.xToCol(scale((barrierLowCorner + barrierDimensions).x)); ++c)
         fields.obstacles->at(r, c) = true;
   constexpr auto barrierBDimensions = Point{6.0, 0.5};
   constexpr auto barrierBLowCorner = Point{barrierLowCorner.x - 2.75, barrierLowCorner.y + barrierDimensions.y - barrierBDimensions.y};
   for (int r = fields.fp.yToRow(scale(barrierBLowCorner.y)); r < fields.fp.yToRow(scale((barrierBLowCorner + barrierBDimensions).y)); ++r)
      for (int c = fields.fp.xToCol(scale(barrierBLowCorner.x)); c < fields.fp.xToCol(scale((barrierBLowCorner + barrierBDimensions).x)); ++c)
         fields.obstacles->at(r, c) = true;
}
NuclearOutput World::initVals(const Point& point) const {
   NuclearOutput initV{};
   initV.diffStateVals = {1.0};
   if (fields.cues->at(point).test(Cues::source)) 
      initV.fieldVals.at(FieldVars::seeker - nDiffStates) = 1.0;
   return initV;
}


#elif PROGRAM == ADVECTION
static constexpr bool obstaclesFlag = false;
static constexpr bool cuesFlag = false;
void World::setCues() {}
void World::placeObstacles() {}
NuclearOutput World::initVals(const Point& point) const {
   NuclearOutput initV{};
   initV.fieldVals.fill(advectionGaussian(point, fields.fp));
   return initV;
}


#elif PROGRAM == MIXTURE
static constexpr bool obstaclesFlag = false;
static constexpr bool cuesFlag = true;
void World::setCues() {
   for (int r = fields.fp.yToRow(0.0); r < fields.fp.yToRow(3.0); ++r)
      for (int c = fields.fp.xToCol(0.0); c < fields.fp.xToCol(3.0); ++c)
         fields.cues->at(r, c).set(Cues::left);
   for (int r = fields.fp.yToRow(0.0); r < fields.fp.yToRow(3.0); ++r)
      for (int c = fields.fp.xToCol(6.0); c < fields.fp.xToCol(9.0); ++c)
         fields.cues->at(r, c).set(Cues::right);
}
void World::placeObstacles() {}
NuclearOutput World::initVals(const Point& point) const {
   NuclearOutput initV{};
   initV.diffStateVals = {
      1.0, 
      fields.cues->at(point).test(Cues::left) ? 1.0 : 0.0, 
      fields.cues->at(point).test(Cues::right) ? 1.0 : 0.0
   };
   return initV;
}


#elif PROGRAM == PARABOLOID
static constexpr bool obstaclesFlag = false;
static constexpr bool cuesFlag = true;
void World::setCues() {
   for (int r = 0; r < fields.fp.ypixels; ++r)
      for (int c = 0; c < fields.fp.xpixels; ++c)
         if ((fields.fp.coordsToPt(r, c) - Point{5.0, 5.0}).mag() < 3.0)
            fields.cues->at(r, c).set(Cues::circle);
}
void World::placeObstacles() {}
NuclearOutput World::initVals(const Point& point) const { return {}; }


#elif PROGRAM == EDGE
static constexpr bool obstaclesFlag = false;
static constexpr bool cuesFlag = false;
void World::setCues() {}
void World::placeObstacles() {}
NuclearOutput World::initVals(const Point& point) const { 
   NuclearOutput initV{};
   initV.diffStateVals = {1.0};
   initV.fieldVals = {1.0};
   return initV;
}


#else
static constexpr bool obstaclesFlag = false;
static constexpr bool cuesFlag = false;
void World::setCues() {}
void World::placeObstacles() {}
NuclearOutput World::initVals(const Point& point) const { return {}; }


#endif


void World::placeNearby(const AgentManager& am) {
   std::uniform_real_distribution<double> randAngle(0.0, 2.0 * M_PI);
   double r = 2.0 * am.provideAxes().second;
   const double maxR = 20.0 * am.provideAxes().first;
   int nObstructed = 0;
   int nOutOfBounds = 0;
   bool foundPlace = false;
   while (r < maxR && !foundPlace) {
      const double circumference = 2.0 * M_PI * r;
      constexpr double roughTriesPerCell = 2.0;
      const int nTriesAround = static_cast<int>(roughTriesPerCell * circumference / fields.fp.cellDiag());
      for (int i = 0; i < nTriesAround && !foundPlace; ++i) {
         const double angle = randAngle(gen);
         const double x = am.provideState().pos.x + r * std::cos(angle);
         const double y = am.provideState().pos.y + r * std::sin(angle);
         try {
            agentManagers.emplace_back(am, Point{x, y});
            //agentManagers.emplace_back(am, AMExternalState{Point{x, y}, randAngle(gen)});
            //agentManagers.emplace_back(am, AMExternalState{Point{x, y}, M_PI / 2.0});
            //std::cerr << "am #" << am.provideID() << " created am #" << agentManagers.back().provideID() << std::endl; // DEBUG
            foundPlace = true;
         } catch (ObstructionException& e) { 
            ++nObstructed; 
         } catch (OOBException& e) { 
            ++nOutOfBounds; 
         }
      }
      r += fields.fp.cellDiag() / roughTriesPerCell;
   }
   if (!foundPlace)
      std::cerr << "No place for new agent found. nObstructed: " << nObstructed << "; nOutOfBounds: " << nOutOfBounds << std::endl;
}


double World::areaInShape() {
   std::uniform_real_distribution<double> randX{0.0, fields.fp.xphys};
   std::uniform_real_distribution<double> randY{0.0, fields.fp.yphys};
   int nInShape = 0;
   //constexpr int trials = 200000000;
   constexpr int trials = 100000;
   for (int i = 0; i < trials; ++i) {
      const Point randomPoint{randX(gen), randY(gen)};
      if (inShape(randomPoint) && !(fields.obstacles && fields.obstacles->at(randomPoint))) 
         ++nInShape;
   }
   return fields.fp.xphys * fields.fp.yphys * static_cast<double>(nInShape) / trials;
}


// parameter nAgents ignored if agentPlacement.initDens >= 0.0
World::World(
      const AgentPlacement& agentPlacement,
      const Substances& subs,
      double visc,
      const AgentManagerParams& amParams,
      const SimulationParams& simParams,
      const Fparams& fparams,
      bool visual,
      std::default_random_engine& gen
      ) : 
   inShape{agentPlacement.inShape},
   fields{fparams, simParams.passthrough, obstaclesFlag, cuesFlag, visual}, 
   subs{subs}, 
   diffuser{fields, subs, simParams.passthrough, visc, viscosityTolerance, pressureTolerance},
   gen{gen}
{

   // FOR PLACING A SINGLE AGENT
   if (debugflow) {
      const Point location{(debugspin ? 0.5 : 0.1) * fparams.xphys, 0.5 * fparams.yphys};
      agentManagers.emplace_back(
            amParams, 
            AMExternalState{location, 0.0}, 
            fields, 
            subs, 
            simParams,
            gen,
            initVals(location)
            );
      return;
   }

   placeObstacles();
   setCues();

   std::uniform_real_distribution<double> randX{0.0, fields.fp.xphys};
   std::uniform_real_distribution<double> randY{0.0, fields.fp.yphys};
   std::uniform_real_distribution<double> randAngle{0.0, 2.0 * M_PI};
   bool startover = true;
   int nOutOfBounds = 0;
   const auto shapeArea = areaInShape();
   if (shapeArea <= 0.0) 
      throw std::runtime_error("areaInShape returned non-positive");

   auto lastReportedTime = std::chrono::steady_clock::now();
   while (startover) {
      agentManagers.clear();
      fields.clearAgentInfo();
      startover = false;
      int nCollision = 0;

      double curDensityRatio = 0.0;
      double massSum = 0.0;
      while (!startover 
            && (agentPlacement.initDens >= 0.0 
               ? agentPlacement.initDens > 0.0 && curDensityRatio < 1.0
               : agentManagers.size() < agentPlacement.nAgents)
            ) {
         const Point location{randX(gen), randY(gen)};
         if (inShape(location)) {

            const auto currentTime = std::chrono::steady_clock::now();
            const std::chrono::duration<double> sinceLastReported = currentTime - lastReportedTime;
            if (sinceLastReported.count() > 30.0) {
               lastReportedTime = currentTime;
               if (agentPlacement.initDens >= 0.0)
                  std::cerr << "Current proportion of target initial density achieved: " << curDensityRatio << std::endl;
               else
                  std::cerr << agentManagers.size() << " agents created successfully out of " << agentPlacement.nAgents << " desired" << std::endl;
            }

            try {
               agentManagers.emplace_back(
                     amParams, 
#if PROGRAM == ADVECTION || PROGRAM == DIFFUSION_TRANSLATION
                     AMExternalState{location, 0.0},
#else
                     AMExternalState{location, randAngle(gen)}, 
#endif
                     fields, 
                     subs, 
                     simParams,
                     gen,
                     initVals(location)
                     );
               //std::cerr << "number: " << agentManagers.size() << " ;  massSum: " << massSum << " ;  shapeArea: " << shapeArea << std::endl;
               massSum += agentManagers.back().provideMass();
               curDensityRatio = massSum / (shapeArea * agentPlacement.initDens);
               nCollision = 0;
            } catch (ObstructionException& e) { 
               ++nCollision; 
            } catch (OOBException& e) { 
               ++nOutOfBounds; 
            }
         }

         constexpr int maxFailures = 1e6;
         if (nCollision + nOutOfBounds > maxFailures) {
            startover = true; 
            std::cerr << "Too many failed placement attempts; starting over." << std::endl; 
         }
      }
   }
   std::cerr << "Number of agents not created because they would be out of bounds: " << nOutOfBounds << std::endl; // DEBUG
   std::cerr << "Final number of agents at start: " << agentManagers.size() << std::endl;

   /* FOR PLACING AGENTS IN A RECTANGULAR GRID
   std::uniform_real_distribution<double> randUnit(0.0, 1.0);
   const double dist = std::sqrt(M_PI * amParams.radius * amParams.radius / agentPlacement.initDens);
   std::cout << dist << std::endl;
   for (double x = dist; x < fields.fp.xphys - dist; x += dist) 
      for (double y = dist; y < fields.fp.xphys - dist; y += dist) 
         agentManagers.emplace_back(amParams, AMExternalState{Point{x, y}, randUnit(gen) * M_PI * 2}, fields, subs, gen,
               passthrough,
               fparams.minH() / dt, // maxSpeed
               maxAngularSpeed,
               //initV
               initVals(x, y)
               );
   */

   fields.clearAgentInfo();
   for (auto& am : agentManagers)
      am.assertAgentInfo();
}


World::World(
      std::default_random_engine& gen,
      std::istream& is
      ) :
   fields{is},
   subs{readObjectFromStream<Substances>(is)},
   diffuser{fields, is},
   gen{gen}
{
   const auto nAgentManagers = readObjectFromStream<std::list<AgentManager>::size_type>(is);
   for(std::list<AgentManager>::size_type iAM = 0; iAM < nAgentManagers; ++iAM)
      agentManagers.emplace_back(fields, gen, is);
}


void World::serialize(std::ostream& os) {
   fields.serialize(os);
   os.write(reinterpret_cast<const char *>(&subs), sizeof(subs));
   //diffuser.serialize(os);
   const auto nAgentManagers = agentManagers.size();
   os.write(reinterpret_cast<const char *>(&nAgentManagers), sizeof(nAgentManagers));
   for (const auto& am : agentManagers)
      am.serialize(os);
}


// As we move the agents, erase AgentManager's as they go out of bounds or die.
// This would need to be a separate loop if we want to parallelize the ami->update().
void World::update(double dt, bool debug) {
   //std::cerr << ".";
   fields.clearAgentInfo();
   for (auto& am : agentManagers)
      am.assertAgentInfo();
   diffuser.advect(dt);
   if (!debugflow) diffuser.diffuse(dt);

   // Using while loop instead of ranged for loop simplifies deleting agents
   auto am = agentManagers.begin();
   while (am != agentManagers.end()) {
      LifeEvent event = am->update(dt, debug);
      if (!am->inBounds() || event == LifeEvent::die) {
         //std::cerr << "am #" << am->provideID() << " destroyed" << std::endl;
         am = agentManagers.erase(am);
         //std::cerr << "Agents remaining: " << agentManagers.size() << "\n";
      } else {
         if (event == LifeEvent::reproduce) placeNearby(*am);
         ++am;
      }
   }
}


void World::burnin(double propOfOriginalToAllow) {
   for (int isub = 0; isub < nSPHSubs; ++isub) {

      // Estimate target sum
      double target = 0.0;
      for (auto& am : agentManagers) {
         target += am.integral(isub);
         if (!std::isfinite(am.integral(isub))) std::cerr << "Position: " << am.provideState().pos << std::endl; // DEBUG
      }
      target /= fields.fp.cellArea(); //integral() gives mass, which we convert to concentration here.
      std::cerr << "Burnin for substance " << isub << ", target: " << target << std::endl;

      // Set up loop variables
      const auto dt = diffuser.maxDt();
      const auto durationGlobalSettle = -std::log(propOfOriginalToAllow) / subs.at(isub).decayRate;
      std::cerr << "Global settling duration: " << durationGlobalSettle << "; dt: " << dt << std::endl;
      const auto startTime = std::chrono::steady_clock::now();
      auto lastReportedTime = startTime;
      double t = 0.0;

      if (target > 0.0) {
         while (t < durationGlobalSettle) {

            // Estimate time remaining for burnin
            const auto currentTime = std::chrono::steady_clock::now();
            const std::chrono::duration<double> sinceLastReported = currentTime - lastReportedTime;
            if (sinceLastReported.count() > 30.0) {
               lastReportedTime = currentTime;
               const std::chrono::duration<double> elapsedSoFar = currentTime - startTime;
               std::cerr << "Burn-in, substance " << isub << ". Minutes elapsed so far: " << elapsedSoFar.count() / 60.0 << "; estimated remaining: ";
               const auto proportionOfLimit = t / durationGlobalSettle;
               const auto estimatedRemaining = elapsedSoFar.count() * (1.0 / proportionOfLimit - 1.0);
               std::cerr << estimatedRemaining / 60.0 << std::endl;
            }

            // Actual burnin iteration
            for (auto& am : agentManagers) am.burnin(isub, dt);
            diffuser.diffuseOne(isub, dt);
            t += dt;
         }

         // Print info to confirm target was reached
         const auto sum = fields.grids.at(isub).sum();
         std::cerr << "Burnin complete; t: " << t << "; sum/target = " << sum / target 
            << "; proportion should be: " << 1.0 - propOfOriginalToAllow << std::endl;
      }
   }
}


void World::cut() {
   auto am = agentManagers.begin();
   while (am != agentManagers.end()) {
      const auto posy = am->provideState().pos.y;
      if (posy >= 2.5 && posy < 4.0)
         am = agentManagers.erase(am);
      else ++am;
   }
}


void World::textupdate(double t) {
   std::cerr << "Time: " << t << "; nAgents: " << agentManagers.size() << "; max agentBodies: " << fields.agentBodies.max() << "; min agentBodies: " << fields.agentBodies.min() << std::endl;
   for (int isub = 0; isub < nSubs; ++isub) {
      double maxc = -DBL_MAX;
      double minc = DBL_MAX;
      for (int r = 0; r < fields.fp.ypixels; ++r) {
         for (int c = 0; c < fields.fp.xpixels; ++c) {
            if (!fields.fullyObstructed(r, c) && fields.concentration(isub, r, c) > maxc) maxc = fields.concentration(isub, r, c);
            if (!fields.fullyObstructed(r, c) && fields.concentration(isub, r, c) < minc) minc = fields.concentration(isub, r, c);
         }
      }
      std::cerr << "  Substance " << isub 
         << ";  Max conc: " << maxc
         << ";  Min conc: " << minc
         << std::endl;
   }
}


void World::outputVel(std::string basefn) const {
   std::string velfilename{basefn + "-vel"};
   std::ofstream velfile(velfilename);
   diffuser.printVel(velfile);
   velfile.close();
   std::system(("bzip2 -f " + velfilename).c_str());
}


// First, calculate average number density of agents and use this to calculate a good smoothing
// length for a Gaussian smoothing function.
void World::gaussianDens(std::string basefn) const {

   // TODO comment based on notes from Aug 1, 2021
   double sumH = 0.0;
   std::vector<double> gaussianHByAgent;
   for (const auto& am : agentManagers) {
      constexpr int desiredNAgents = 21; // From Liu and Liu 2003, pg. 130
      std::vector<double> distances;
      for (const auto& neighbor : agentManagers)
         distances.push_back((am.provideState().pos - neighbor.provideState().pos).mag());
      std::sort(distances.begin(), distances.end());
      if (distances.size() >= desiredNAgents)
         gaussianHByAgent.push_back(distances.at(desiredNAgents - 1) / 2.0);
      else
         gaussianHByAgent.push_back(distances.back() * std::sqrt(static_cast<double>(desiredNAgents) / agentManagers.size()) / 2.0);
      sumH += gaussianHByAgent.back();
   }
   if (nAgentsRemaining() > 0) 
      std::cerr << "Average smoothing length: " << sumH / agentManagers.size() << std::endl;
   
   std::string densfilename{basefn + "-gaussianDens"};
   std::ofstream densfile(densfilename);
   for (int r = 0; r < fields.fp.ypixels; ++r) {
      for (int c = 0; c < fields.fp.xpixels; ++c) {
         double sumSPH = 0.0;
         const auto gridPos = fields.fp.coordsToPt(r, c);
         auto hIterator = gaussianHByAgent.begin();
         for (const auto& am : agentManagers) {
            const auto mass = am.provideMass();
            const auto agentPos = am.provideState().pos;
            const auto distance = (agentPos - gridPos).mag();
            const auto h = *hIterator;
            sumSPH += mass * std::exp(-std::pow(distance / h, 2.0)) / (M_PI * std::pow(h, 2.0));
            ++hIterator;
         }
         densfile << sumSPH << std::endl;
      }
   }
   densfile.close();
   std::system(("bzip2 -f " + densfilename).c_str());
}


void World::agentFieldVals(std::string basefn, int isub) const {
   const int iFieldVal = isub - nDiffStates;
   if (iFieldVal < 0 || iFieldVal >= nFieldVars) throw std::invalid_argument("Agent::agentFieldVals() only makes sense for field variables.");
   std::string affilename{basefn + "-agentField-" + std::to_string(iFieldVal)};
   std::ofstream affile(affilename);
   for (const auto& am : agentManagers) {
      const auto location = am.provideState().pos;
      affile << location.x << " " 
         << location.y << " "
         << am.provideNucOut().fieldVals.at(iFieldVal) << std::endl;
   }
   affile.close();
   std::system(("bzip2 -f " + affilename).c_str());
}


void World::edgeInfo(std::string basefn) const {
   std::string eifilename{basefn + "-edgeInfo"};
   std::ofstream eifile(eifilename);
   eifile << std::sqrt(subs.at(DiffStates::one).decayRate / subs.at(DiffStates::one).difRate) << std::endl;
   for (const auto& am : agentManagers) {
      const auto location = am.provideState().pos;
      eifile << location.x << " " << am.provideDens().at(DiffStates::one) << " " << am.provideDensGrad().mag() << " " << am.provideColorGrad().mag() << std::endl;
   }
   eifile.close();
   std::system(("bzip2 -f " + eifilename).c_str());
}


void World::printAgentStats() const {
   if (nAgentsRemaining() == 0) return;

   std::vector<std::array<double, nDiffStates> > densitiesAll;
   std::vector<std::array<double, 1> > densgradmagAll;
   std::vector<std::array<double, 1> > densgradmagHighmk;
   std::vector<std::array<double, nDiffStates> > diffStateValsAll;
   std::vector<std::array<double, nFieldVars> > fieldValsAll;
   for (const auto& am : agentManagers) {
      densitiesAll.push_back(am.provideDens());
      if (am.provideDens().at(0) > 0.0) densgradmagAll.push_back({am.provideDensGrad().mag() / am.provideDens().at(0)});
      //if (am.provideNucOut().fieldVals.at(2) > 0.95) densgradmagHighmk.push_back({am.provideDensGrad().mag() / am.provideDens().at(0)});
      diffStateValsAll.push_back(am.provideNucOut().diffStateVals);
      fieldValsAll.push_back(am.provideNucOut().fieldVals);
   }

   std::cout << "Agent-sensed densities: " << std::endl;
   printArrayStats<nDiffStates>(densitiesAll);
   std::cout << "Agent-sensed density gradient magnitude over density: " << std::endl;
   printArrayStats<1>(densgradmagAll);
   //std::cout << "Agent-sensed density gradient magnitude over density when marker over 0.95: " << std::endl;
   //printArrayStats<1>(densgradmagHighmk);
   std::cout << "Agent differentiation states: " << std::endl;
   printArrayStats<nDiffStates>(diffStateValsAll);
   if (nFieldVars > 0) {
      std::cout << "Agent field values: " << std::endl;
      printArrayStats<nFieldVars>(fieldValsAll);
   }
}


void World::fillShape(double concentration) {
   for (int r = 0; r < fields.fp.ypixels; ++r) {
      for (int c = 0; c < fields.fp.xpixels; ++c) {
         const Point truePt = fields.fp.coordsToPt(r, c);
         if (inShape(truePt)) {
            fields.grids.at(DiffStates::one).at(r, c) = concentration;
         }
      }
   }
}


void World::fillGaussian() {
   for (int r = 0; r < fields.fp.ypixels; ++r) {
      for (int c = 0; c < fields.fp.xpixels; ++c) {
         const Point truePt = fields.fp.coordsToPt(r, c);
         fields.grids.at(DiffStates::one).at(r, c) = advectionGaussian(truePt, fields.fp);
      }
   }
}


// find average orientation, originally for debugging but could be useful in future
void World::avgOri() const {
   double avgSin = 0.0;
   double avgCos = 0.0;
   for (const auto& am : agentManagers) {
      avgSin += sin(am.provideState().ori);
      avgCos += cos(am.provideState().ori);
   }
   avgSin /= agentManagers.size();
   avgCos /= agentManagers.size();
   std::cerr << "avgOri: " << atan2(avgSin, avgCos) << ";         avgRad: " << sqrt(avgSin * avgSin + avgCos * avgCos) << std::endl;
}


void World::dots() {
   static const int a = std::max(fields.fp.ypixels / debugndots, 1);
   for (int r = 0; r < fields.fp.ypixels; ++r) {
      for (int c = 0; c < fields.fp.xpixels; ++c) {
         if (((r / a) % 2 == 0 && (c / a) % 2 == 0)) {
            fields.grids.at(DiffStates::one).at(r, c) = 1.0;
         }
      }
   }
}


void World::fill() {
   for (int r = 0; r < fields.fp.ypixels; ++r)
      for (int c = 0; c < fields.fp.xpixels; ++c)
         fields.grids.at(DiffStates::one).at(r, c) = 1.0;
}


template<int n>
void printArrayStats(const std::vector<std::array<double, n> >& a) {
   if (a.size() == 0) {
      std::cout << "No agents present." << std::endl;
   } else {
      for (int iVal = 0; iVal < n; ++iVal) {
         double avg = 0.0;
         double max = 0.0;
         double min = std::numeric_limits<double>::max();
         for (size_t iAgent = 0; iAgent < a.size(); ++iAgent) {
            const double v = a.at(iAgent).at(iVal);
            avg += v;
            if (v > max) max = v;
            if (v < min) min = v;
         }
         avg /= a.size();
         std::cout << "  " << iVal << "; Avg: " << avg << "; Min: " << min << "; Max: " << max << std::endl;
      }
   }
}


double advectionGaussian(const Point& point, const Fparams& fp) {
   const Point center{Point{fp.xphys, fp.yphys} / 2.0};
   const Point relativeToCenter = point - center;
   return std::exp(-std::pow(relativeToCenter.mag() / 3.0, 2.0));
}
