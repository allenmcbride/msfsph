// agents.cpp
// Allen McBride
// March 10, 2022
//
// This is the main module for an MSF-SPH simulator. This program
// simulates an agent-based implementation of the Morphgen language for
// artificial morphogenesis. Morphgen describes morphogenetic processes
// in continuous terms using PDEs. This implementation uses MSF-SPH,
// a natural computating variant of smoothed particle hydrodynamics. In
// MSF-SPH, agents do not communicate directly, but rather through
// shared fields of dilute substances that diffuse and are degraded
// in an aqueous medium. MSF-SPH is, conceptually, a global-to-local
// compilation technique for Morphgen; however, this simulator takes a
// largely interpretive approach. Several Morphgen programs are hard-coded
// in nucleus.cpp.
//
// McBride 2022 presents in detail MSF-SPH as an approach to Morphgen
// compilation. The purpose, meaning, and overall structure of this
// simulator are detailed in the chapter entitled "Model." The reader
// should refer to that document as a starting point; comments in this
// code focus on implementation details.


#include "cxxopts.hpp"
#include "agent.h"
#include "agentmanager.h"
#include "calibrator.h"
#include "diffuser.h"
#include "diffuserfdm.h"
#include "fields.h"
#include "picture.h"
#include "screen.h"
#include "substances.h"
#include "world.h"

#include <link.h>

#include <cfenv>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>


class InsufficientCalibration : public std::runtime_error {
   public: using runtime_error::runtime_error;
};


// Number of seconds between each type of update: live (graphical, on
// screen), pic (graphical, to file), text (various statistics to console)
struct SecondsPerUpdate {
   double live;
   double pic;
   double text;
};


// Runtime parameters relating to calibration. Meanings are documented
// for the corresponding command-line options.
struct CalibInfo {
   bool nocalib;
   bool nonewcalib;
   bool recalib;
   double onecalib;
};


// Turn command line arguments into a string. This was more useful when
// an earlier version of the cxxopts parser modified its arguments.
std::string computeCmdLine(int argc, char** argv);


// This block of functions compute various internal simulation parameters
// from command-line parameters. Even simple logic for this task is
// wrapped in functions for consistency.
double computeStalkMult(double stalkMultPrelim, bool passthrough);
AgentManagerParams computeAgentManagerParams(const AgentManagerParams& amParamsPrelim, bool noNoise, double diam, double mass, double eccentricity, double stalkMult);
std::tuple<double, double> computeSpatialResolution(double nCellWidthsPerDiameterPrelim, double maxCellWidth, double cellWidthPrelim, double semiminorAxis, bool passthrough);
std::vector<double> computeDifRates(const std::vector<double>& difRatesPrelim, const std::vector<double>& decayRates, const std::vector<double>& smoothingH);
std::vector<double> computeEnvDifRates(const std::vector<double>& envDifRatesPrelim, const std::vector<double>& EnvDecayRates);
std::vector<double> computeDecayRates(const std::vector<double>& decayRatesPrelim, const std::vector<double>& halfLives, bool idealResults);
Substances computeSubstances(const std::vector<double>& difRates, const std::vector<double>& decayRates, const std::vector<double>& gains, const std::vector<double>& envDifRates, const std::vector<double>& envDecayRates);
std::tuple<double, double, double> computeTimeParams(double dt, double minTimestepsPerCell, double maxspeedPrelim, double tlimitPrelim, double speedupFactor, const Substances& subs, double h, bool idealResults);
double findMaxAngularSpeed(double gridh, const AgentManagerParams& amParams, double dt, double proportionAllowed, bool passthrough);
unsigned int computeSeed(std::string seedString);


// Print miscellaneous info about substances for user feedback.
void printSubstanceInfo(const Substances& subs, double initDens, AgentManagerParams amParams);


// Write to file info that would be needed to reproduce a given run:
// Code, command line, libraries, seed.
void writeInfoToReproduce(const std::string& picdir, const std::string& cmdLine, unsigned int seed);


// Adjust update timing based on optional artificial speedup parameter.
SecondsPerUpdate adjustSecondsPerUpdate(const SecondsPerUpdate& secondsPerUpdatePrelim, double speedupFactor);


// Set up Calibrator object and run calibration.
void calibrate(
      const CalibInfo& calibInfo, 
      const Substances& subs, 
      double visc, 
      const AgentManagerParams& amParams, 
      double nCellWidthsPerDiameter,
      double minTimestepsPerCell,
      double liveSecondsPerUpdate,
      int pmag,
      const SimulationParams& simParams,
      std::default_random_engine& gen
      );


// For those Morphgen programs where it makes sense to use an FDM
// simulation to find results, perform this simulation. At least that's
// the intention; in reality this is only used for the paraboloid-forming
// program. For the heat equation programs, it's simpler just to hack
// the main runSimulation() function below.
//
void runSimulationFDM(
      const Fparams& fparams,
      const Substances& subs,
      double targetLaplacian,
      double radius,
      double dt, 
      double tlimit, 
      const std::string& picdir, 
      const SecondsPerUpdate& secondsPerUpdate, 
      const SubstanceVars& fixedMaxVals, 
      int pmag
      );


// Contains main loop for simulation and deals with graphical output.
void runSimulation(World& world, double dt, double tlimit, const std::string& picdir, const SecondsPerUpdate& secondsPerUpdate, const SubstanceVars& fixedMaxVals, int pmag);


// Print version of shared libraries.
int dlPrint(dl_phdr_info* info, size_t, void* outfile);


// Extend a vector by repeating the last element until some tarted size
// reached. Used to avoid specifying repeated substance parameters.
template <typename T> T repeatBack(const std::vector<T>& vec, size_t i);


// For forward Euler method, compute number of subdivisions of a timestep
// needed to avoid violating diffusion stability criterion.
int subcyclesNeeded(double dt, double maxdt);


// Convert a half-life, in minutes, to the corresponding decay rate,
// in inverse seconds.
double halfMinToDecaySec(double half);


// Figure out how much an ellipse can rotate in one timestep without
// covering or uncovering an entire grid cell at once.  Negative return
// value indicates no limit with regards to covering/uncovering grid cells.
double findMaxRotationWithoutCoveringCellProportion(double gridh, const AgentManagerParams& amParams, double proportionAllowed);


// Helper function for above; return maximum proportion of cell newly
// covered by a given rotation.
bool rotationCoversCellProportion(double gridh, const AgentManagerParams& amParams, double angle, double proportion);


// Define shapes in which agents should be initially placed, for different
// programs, by returning true when an agent is inside the shapes. One
// of these will be partially evaluated in this module and then passed
// to World as a lambda taking only a location as parameter.
bool growthRectangle(const Point& point, const Fparams& fp, double height);
bool growxShape(const Point& point, const Fparams& fp);
bool circle(const Point& point, const Fparams& fp, double radius);
bool twoCircles(const Point& point, const Fparams& fp, double radius);


int main(int argc, char** argv) try {

   feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO);

   std::vector<std::string> helpGroups;
   cxxopts::Options options(argv[0], "SPH simulation using natural smoothing functions");

   helpGroups.emplace_back("Numerical simulation parameters");
   options.add_options(helpGroups.back())
      ("dt", "Time interval for updating (s)", cxxopts::value<double>()->default_value("-1.0")) // Invalid default to detect explicit setting
      ("tpc", "Minimum timesteps per grid cell traversal (inverse of Courant number)", cxxopts::value<double>()->default_value("1.0"))
      ("tlim", "Time to end simulation (s)", cxxopts::value<double>()->default_value("0.0")) // Set negative to continue until no agents remain
      ("nocalib", "Neither generate nor use calibration data. Overrides all other calibration options.")
      ("nonewcalib", "Use existing calibration data, but do not perform any new calibration. Overrides --recalib")
      ("recalib", "Force recalibration, overwriting any existing file.")
      ("onecalib", "Calibrate for just one speed, so a script can parallelize.", cxxopts::value<double>()->default_value("-1.0")) // invalid default to detect explicit setting
      ("swimming", "Substances flow around agent bodies.")
      ("hperdiam", "Minimum grid cell widths per agent diameter for calibration", cxxopts::value<double>()->default_value("3.0"))
      ("minxpix", "Minimum number of grid cells in x dimension.", cxxopts::value<int>()->default_value("0")) // invalid default to detect explicit setting 
      ("gridh", "Grid cell width (microns). Overrides --hperdiam and --minxpix.", cxxopts::value<double>()->default_value("0.0"))
      ("xphys", "Physical size of x dimension (mm)", cxxopts::value<double>()->default_value("10.0"))
      ("yphys", "Physical size of y dimension (mm)", cxxopts::value<double>()->default_value("10.0"))
      ("burnin", "Proportion of steady state amount to wait for during burnin", cxxopts::value<double>()->default_value("0.99"))
      ("seed", "Seed for pseudorandom number generator", cxxopts::value<std::string>()->default_value(""))
      ;

   helpGroups.emplace_back("Agent parameters");
   options.add_options(helpGroups.back())
      ("diam", "Diameter or major axis of agents (microns)", cxxopts::value<double>()->default_value("40"))
      ("mass", "Mass of agents (grams, assumes density is 1.0). Overrides --diam.", cxxopts::value<double>()->default_value("0.0"))
      ("ecc", "Eccentricity of agents (default zero)", cxxopts::value<double>()->default_value("0.0"))
      ("stk", "Stalk length, as multiple of maximum recommended grid cell diagonal", cxxopts::value<double>()->default_value("-1.0")) // invalid default to detect explicit setting
      ("sensbias", "Standard deviation of natural logarithm of log-normal distribution for sensor bias multiplier (unitless)", cxxopts::value<double>()->default_value("0.0"))
      ("bindrate", "Frequency of binding events for sensor at unit concentration (events/(micromoles*seconds)); set to zero for a continuous model.", cxxopts::value<double>()->default_value("0.0"))
      ("bindhalf", "Half life for molecules remaining bound to a sensor (minutes); set to zero for instantaneous and continuous measurements", cxxopts::value<double>()->default_value("16.0"))
      ("oribias", "Standard deviation of a wrapped normal distribution, a per-agent sample from which is added to each agent rotation (radians)", cxxopts::value<double>()->default_value("0.0"))
      ("speedbias", "Standard deviation of natural logarithm of log-normal distribution, a per-agent sample from which is multiplied by each speed (unitless)", cxxopts::value<double>()->default_value("0.0"))
      ("motionnoise", "Standard deviation of 2D normal distribution, a sample from which is scaled by agent speed and added to motion (unitless)", cxxopts::value<double>()->default_value("0.0"))
      ("sensormttf", "Mean time to failure of sensor system (seconds)", cxxopts::value<double>()->default_value("0.0"))
      ("propulsionmttf", "Mean time to failure of propulsion system (seconds)", cxxopts::value<double>()->default_value("0.0"))
      ("productionmttf", "Mean time to failure of production system (seconds)", cxxopts::value<double>()->default_value("0.0"))
      ("nonoise", "Turn off all noise in simulation, regardless of other parameters.")
      ("maxspeed", "Maximum agent speed (mm/s). If specified alone with --dt, --tpc is ignored.", cxxopts::value<double>()->default_value("-1.0")) // invalid default to detect explicit setting
      ("leakhalf", "Half life for leaky integral estimating density change rate (minutes)", cxxopts::value<double>()->default_value("3840.0"))
      ;

   helpGroups.emplace_back("Other physical simulation parameters");
   options.add_options(helpGroups.back())
      ("nagents", "Number of agents", cxxopts::value<size_t>()->default_value("0"))
      //("isf", "Use instantaneous smoothing functions instead of morphogen smoothing functions")
      ("visc", "Kinematic viscosity of substrate, (mm^2/s, or centiStokes)", cxxopts::value<double>()->default_value("1.0"))
      ("decay", "Degradation rates (proportion/s). Overrides --half.", cxxopts::value<std::vector<double> >()->default_value("-1.0")) // invalid default to detect explicit setting
      ("half", "Half lives (minutes)", cxxopts::value<std::vector<double> >()->default_value("16"))
      ("diff", "Diffusion rates (microns^2/s)", cxxopts::value<std::vector<double> >()->default_value("100"))
      ("smoothingh", "Smoothing length (h, in microns)", cxxopts::value<std::vector<double> >()->default_value("-1.0")) // invalid default to detect explicit setting
      ("gain", "Gain for each substance (dimensionless)", cxxopts::value<std::vector<double> >()->default_value("1.0"))
      ("envdiff", "Environmetnal diffusion rates (microns^2/s)", cxxopts::value<std::vector<double> >()->default_value("0.0"))
      ("envdecay", "Environmental decay rates (proportion/s)", cxxopts::value<std::vector<double> >()->default_value("0.0"))
      ("initdens", "Initial density (proportion of space). Overrides --nagents.", cxxopts::value<double>()->default_value("-1.0"))
      ("speedup", "Factor to multiply all Morphgen velocities by, for speed-related testing.", cxxopts::value<double>()->default_value("1.0"))
      ;

   helpGroups.emplace_back("Output parameters");
   options.add_options(helpGroups.back())
      ("live", "Seconds of simulation time between each live update (Zero for every iteration)", cxxopts::value<double>()->default_value("-1"))
      ("pic", "Seconds of simulation time between each picture file update (Zero for every iteration)", cxxopts::value<double>()->default_value("-1"))
      ("text", "Seconds of simulation time between each text update (Zero for every iteration)", cxxopts::value<double>()->default_value("-1"))
      ("pmag", "Magnification factor for live updates (integer)", cxxopts::value<int>()->default_value("1"))
      ("picdir", "Name associated with run, for various records.", cxxopts::value<std::string>()->default_value("agent_pics"))
      ("subdir", "Optional subdirectory of fields/ to put final output for agent simulations.", cxxopts::value<std::string>()->default_value("."))
      ("ideal", "Whether to output fields representing ideal results, if defined for a given Morphgen program.")
      ("h,help", "Print usage")
      ;

   const auto cmdLine = computeCmdLine(argc, argv);
   const auto result = options.parse(argc, argv);
   if (result.count("help")) {
      std::cout << options.help(helpGroups) << std::endl;
      return 0;
   }

   // Numerical simulation parameters
   const auto dtPrelim = result["dt"].as<double>();
   const auto minTimestepsPerCell = result["tpc"].as<double>();
   const auto tlimitPrelim = result["tlim"].as<double>();
   const CalibInfo calibInfo{
      result["nocalib"].as<bool>(),
      result["nonewcalib"].as<bool>(),
      result["recalib"].as<bool>(),
      result["onecalib"].as<double>(),
   };
   const auto nCellWidthsPerDiameterPrelim = result["hperdiam"].as<double>();
   const auto minimumXPixels = result["minxpix"].as<int>();
   const auto cellWidthPrelim = result["gridh"].as<double>();
   const auto xphys = result["xphys"].as<double>();
   const auto yphys = result["yphys"].as<double>();
   const auto passthrough = !result["swimming"].as<bool>();
   const auto burninProp = result["burnin"].as<double>();
   const auto seedString = result["seed"].as<std::string>(); // This is a string because any negative sentinel value necessitates a signed integral type, which might exclude some seeds

   // Agent parameters
   const auto diam = result["diam"].as<double>();
   const auto mass = result["mass"].as<double>();
   const auto eccentricity = result["ecc"].as<double>();
   const auto stalkMultPrelim = result["stk"].as<double>();
   const AgentManagerParams amParamsPrelim{
      0.0, 0.0, 0.0,
      result["sensbias"].as<double>(),
      result["bindrate"].as<double>(),
      60.0 * result["bindhalf"].as<double>() / std::log(2.0),
      result["oribias"].as<double>(),
      result["speedbias"].as<double>(),
      result["motionnoise"].as<double>(),
      result["sensormttf"].as<double>(),
      result["propulsionmttf"].as<double>(),
      result["productionmttf"].as<double>()
   };
   const auto noNoise = result["nonoise"].as<bool>();
   const auto maxspeedPrelim = result["maxspeed"].as<double>();
   const auto leakhalf = result["leakhalf"].as<double>();

   // Other physical simulation parameters
   const auto nagents = result["nagents"].as<size_t>();
   const auto visc = result["visc"].as<double>();
   const auto decayRatesPrelim = result["decay"].as<std::vector<double> >();
   const auto halfLives = result["half"].as<std::vector<double> >();
   const auto difRatesPrelim = result["diff"].as<std::vector<double> >();
   const auto smoothingH = result["smoothingh"].as<std::vector<double> >();
   const auto gains = result["gain"].as<std::vector<double> >();
   const auto envDifRatesPrelim = result["envdiff"].as<std::vector<double> >();
   const auto envDecayRates = result["envdecay"].as<std::vector<double> >();
   const auto initDens = result["initdens"].as<double>();
   const auto speedupFactor = result["speedup"].as<double>();

   // Output parameters
   const SecondsPerUpdate secondsPerUpdatePrelim{
      result["live"].as<double>(),
      result["pic"].as<double>(),
      result["text"].as<double>()
   };
   const auto pmag = result["pmag"].as<int>();
   const auto picdir = result["picdir"].as<std::string>();
   const auto subdir = result["subdir"].as<std::string>();
   const auto idealResults = result["ideal"].as<bool>();

   // Call helper functions to compute internal parameters from those given on the command line
   const auto stalkMult = computeStalkMult(stalkMultPrelim, passthrough);
   const auto amParams = computeAgentManagerParams(amParamsPrelim, noNoise, diam, mass, eccentricity, stalkMult);
   const auto maxCellWidth = minimumXPixels > 0 ? xphys / minimumXPixels : std::numeric_limits<double>::max();
   const auto [prelimGridH, nCellWidthsPerDiameter] = computeSpatialResolution(nCellWidthsPerDiameterPrelim, maxCellWidth, cellWidthPrelim, amParams.semiminorAxis, passthrough);
   const Fparams fparams{xphys, yphys, prelimGridH};
   const auto decayRates = computeDecayRates(decayRatesPrelim, halfLives, idealResults);
   const auto difRates = computeDifRates(difRatesPrelim, decayRates, smoothingH);
   const auto envDifRates = computeEnvDifRates(envDifRatesPrelim, envDecayRates);
   const auto subs = computeSubstances(difRates, decayRates, gains, envDifRates, envDecayRates);
   printSubstanceInfo(subs, initDens, amParams);
   const auto [dt, maxspeed, tlimit] = computeTimeParams(dtPrelim, minTimestepsPerCell, maxspeedPrelim, tlimitPrelim, speedupFactor, subs, fparams.minH(), idealResults);
   const auto maxAngularSpeed = findMaxAngularSpeed(fparams.minH(), amParams, dt, 0.01, passthrough);
   const SimulationParams simParams{passthrough, halfMinToDecaySec(leakhalf), maxspeed, maxAngularSpeed, speedupFactor};
   const auto seed = computeSeed(seedString);
   std::default_random_engine gen{seed};
   writeInfoToReproduce(picdir, cmdLine, seed);
   const auto secondsPerUpdate = adjustSecondsPerUpdate(secondsPerUpdatePrelim, speedupFactor);
   const bool visual = secondsPerUpdate.live >= 0.0 || secondsPerUpdate.pic >= 0.0;
   std::filesystem::create_directory("fields");
   std::filesystem::create_directory("fields/" + subdir);

   // Program-specific hard-coded parameters
#if PROGRAM == GROWTH
   constexpr double initHeight = 3.0;
#endif

   // If requested, compute ideal results for a given program, for comparison with MSF-SPH simulations
   if (idealResults) {
#if PROGRAM == GROWTH
      const double finalHeight = initHeight + 0.0001 * (tlimit - 10800.0);
      World ideal{{0, 0.0, [fparams, finalHeight](const Point& p){ return growthRectangle(p, fparams, finalHeight); }}, subs, 0.0, {}, {}, fparams, visual, gen};
      ideal.fillShape(initDens);
      ideal.provideFields().outputField("fields/" + picdir + "-idealFinal", DiffStates::one);
#elif PROGRAM == DIFFUSION || PROGRAM == DIFFUSION_TRANSLATION
      World ideal{{0, 0.0, [fparams](const Point& p){ return circle(p, fparams, 1.0); }}, subs, 0.0, {}, {true, 0.0, 0.0}, fparams, visual, gen};
      ideal.fillShape(initDens);
      const SubstanceVars fixedMaxVals{0.0};
      runSimulation(ideal, dt, tlimit, picdir, secondsPerUpdate, fixedMaxVals, pmag);
      ideal.provideFields().outputField("fields/" + picdir + "-idealFinal", DiffStates::one);
#elif PROGRAM == ADVECTION
      World ideal{{0, 0.0, [](const Point& p){ return true; }}, subs, 0.0, {}, {}, fparams, visual, gen};
      ideal.fillGaussian();
      ideal.provideFields().outputField("fields/" + picdir + "-idealFinal", DiffStates::one);
#elif PROGRAM == PARABOLOID
      //const SubstanceVars fixedMaxVals{2.5};
      const SubstanceVars fixedMaxVals{0.0};
      runSimulationFDM(fparams, subs, -1.0, 3.0, dt, tlimit, picdir + "-idealFinal", secondsPerUpdate, fixedMaxVals, pmag);
#else
      std::cerr << "No ideal results defined for this Morphgen program." << std::endl;
#endif

   } else {

      // Set up and run MSF-SPH simulation
#if PROGRAM == GROWTH
      const auto inShape = [fparams, initHeight](const Point& p) { return growthRectangle(p, fparams, initHeight); };
      const SubstanceVars fixedMaxVals{0.2, 1.0};
#elif PROGRAM == GROWX
      const auto inShape = [fparams](const Point& p) { return growxShape(p, fparams); };
      const SubstanceVars fixedMaxVals{0.0, 0.2, 0.2, 1.0, 1.0};
#elif PROGRAM == DIFFUSION || PROGRAM == DIFFUSION_TRANSLATION
      const auto inShape = [fparams](const Point& p) { return circle(p, fparams, 1.0); };
      const SubstanceVars fixedMaxVals{0.25};
#elif PROGRAM == PATH
      const auto inShape = [](const Point& p) { return true; };
      const SubstanceVars fixedMaxVals{1.0, 1.0, 1.0, 1.0};
#elif PROGRAM == ADVECTION
      const auto inShape = [](const Point& p) { return true; };
      const SubstanceVars fixedMaxVals{1.0, 1.0};
#elif PROGRAM == MIXTURE
      const auto inShape = [fparams](const Point& p) { return twoCircles(p, fparams, 1.5); };
      const SubstanceVars fixedMaxVals{0.0, 0.2, 0.2, 0.1, 0.1};
#elif PROGRAM == PARABOLOID
      const auto inShape = [](const Point& p) { return true; };
      const SubstanceVars fixedMaxVals{1.0, 2.5};
#elif PROGRAM == EDGE
      const auto inShape = [fparams](const Point& p) { return p.x > fparams.xphys / 2.0; };
      const SubstanceVars fixedMaxVals{0.0};
#else // e.g., CALIBRATE, where it doesn't matter because actual simulation shouldn't be run
      const auto inShape = [](const Point& p) { return true; };
      const SubstanceVars fixedMaxVals{};
#endif

      calibrate(calibInfo, subs, visc, amParams, nCellWidthsPerDiameter, minTimestepsPerCell, secondsPerUpdate.live, pmag, simParams, gen);
      const AgentPlacement agentPlacement{nagents, initDens, inShape};

      if (tlimit != 0.0) {
         std::cout << "Grid is " << fparams.ypixels << " rows by " << fparams.xpixels << " columns." << std::endl;
         World world{agentPlacement, subs, visc, amParams, simParams, fparams, visual, gen};
         if (debugflow) world.dots();
         else world.burnin(1.0 - burninProp);
         runSimulation(world, dt, tlimit, picdir, secondsPerUpdate, fixedMaxVals, pmag);
         std::cerr << "Seed: " << seed << std::endl;
         if (!debugflow) {
#if PROGRAM == ADVECTION || PROGRAM == PARABOLOID
            world.agentFieldVals("fields/" + subdir + "/" + picdir, FieldVars::cena); 
#elif PROGRAM == EDGE
            world.edgeInfo("fields/" + subdir + "/" + picdir);
#else
            world.gaussianDens("fields/" + subdir + "/" + picdir);
#endif
         }
      } else {
         std::cerr << "tlimit is zero; skipping simulation and final output." << std::endl;
      }
   }

   return 0;

} catch (const FieldsTooBig& e) {
   std::cerr << e.what() << std::endl;
   std::cerr << "Exception: Fields too big. Returning success from main()." << std::endl;
   return 0;
} catch (const TooManyIterationsNeeded& e) {
   std::cerr << e.what() << std::endl;
   std::cerr << "Exception: Too many iterations needed. Returning success from main()." << std::endl;
   return 0;
} catch (const InsufficientCalibration& e) {
   std::cerr << e.what() << std::endl;
   std::cerr << "Exception: Existing calibration data insufficient, refraining from creating more. Returning success from main()." << std::endl;
   return 0;
} catch (const MaxSpeedExceeded& e) {
   std::cerr << e.what() << std::endl;
   std::cerr << "Exception: Maximum speed exceeded. Returning failure from main()." << std::endl;
   return 1;
} catch (const std::exception& e) {
   std::cerr << e.what() << std::endl;
   std::cerr << "Failure due to exception. Returning failure from main()." << std::endl;
   return 1;
} catch(...) {
   std::cerr << "Failure due to exception. Returning failure from main()." << std::endl;
   return 1;
}


std::string computeCmdLine(int argc, char** argv) {
   std::ostringstream cmdLine;
   for (int i = 0; i < argc; ++i) cmdLine << argv[i] << " "; 
   return cmdLine.str();
}


double computeStalkMult(double stalkMultPrelim, bool passthrough) {
   if (stalkMultPrelim >= 0.0) { 
      return stalkMultPrelim;
   } else if (passthrough) {
      return 0.0;
   } else {
      return 1.0;
   }
}


AgentManagerParams computeAgentManagerParams(const AgentManagerParams& amParamsPrelim, bool noNoise, double diam, double mass, double eccentricity, double stalkMult) {
   AgentManagerParams amParams{noNoise ? amParamsPrelim.noNoise() : amParamsPrelim};
   const auto e2 = std::pow(eccentricity, 2.0);
   if (mass > 0.0) {
      amParams.semimajorAxis = std::sqrt(mass / M_PI) * std::pow(1.0 - e2, -0.25);
      amParams.semiminorAxis = std::sqrt(mass / M_PI) * std::pow(1.0 - e2, 0.25);
      std::cout << "Ignoring diameter. Based on requested mass of " << mass << " and eccentricity of " << eccentricity;
   } else {
      amParams.semimajorAxis = 0.001 * diam / 2.0;
      amParams.semiminorAxis = amParams.semimajorAxis * std::sqrt(1.0 - e2);
      std::cout << "Based on diameter of " << diam << " and eccentricity of " << eccentricity;
   }
   const auto minorAxis = 2.0 * amParams.semiminorAxis;
   amParams.stalkLen = stalkMult * minorAxis * std::sqrt(2.0) / Diffuser::minimumRecommendedGridCellWidthsPerMinorAxisSwimming;
   const auto oldprec = std::cout.precision(std::numeric_limits<double>::max_digits10);
   std::cout << ", semiminor and semimajor axes will be: " << amParams.semiminorAxis << ", " << amParams.semimajorAxis << std::endl;
   std::cout << "Stalk length will be " << amParams.stalkLen << std::endl;
   std::cout << "Mass will be " << amParams.mass() << std::endl;
   std::cout.precision(oldprec);
   return amParams;
}


std::tuple<double, double> computeSpatialResolution(double nCellWidthsPerDiameterPrelim, double maxCellWidth, double cellWidthPrelim, double semiminorAxis, bool passthrough) {
   const auto cellWidthExplicitlySet = cellWidthPrelim > 0.0;
   const auto h = cellWidthExplicitlySet ? cellWidthPrelim : std::min(maxCellWidth, 2.0 * semiminorAxis / nCellWidthsPerDiameterPrelim);
   const auto nCellWidthsPerDiameter = cellWidthExplicitlySet ? 2.0 * semiminorAxis / h : nCellWidthsPerDiameterPrelim;
   std::cout << "Cell widths per agent diameter: " << nCellWidthsPerDiameter << std::endl;
   if (nCellWidthsPerDiameter < std::sqrt(2.0)) 
      throw std::runtime_error("hperdiam must be at least one cell diagonal (sqrt(2)) for sufficient agent resolution");
   if (!passthrough && nCellWidthsPerDiameter < Diffuser::minimumRecommendedGridCellWidthsPerMinorAxisSwimming)
      std::cerr << "Warning: for swimming agents, hperdiam may be too low. At least 6.0 is recommended. Also check stalk length." << std::endl;
   std::cout << "Grid cell width (h) set to " << h << std::endl;
   return {h, nCellWidthsPerDiameter};
}


std::vector<double> computeDecayRates(const std::vector<double>& decayRatesPrelim, const std::vector<double>& halfLives, bool idealResults) {
   std::cerr << "Degradation rates...\n";
   std::vector<double> decayRates;

   if (idealResults) {
      std::cerr << "  Ideal results requested, so ignoring any explicit settings of degradation rates and setting to zero." << std::endl;
      decayRates.push_back(0.0);
      
   } else if (decayRatesPrelim.front() < 0.0) { 
      for (size_t i = 0; i < halfLives.size(); ++i)
         decayRates.push_back(halfMinToDecaySec(halfLives.at(i)));
   } else {
      decayRates = decayRatesPrelim;
   }

   for (auto& decayRate : decayRates) std::cerr << "  Degradation rate: " << decayRate << " per second\n";
   return decayRates;
}


std::vector<double> computeDifRates(const std::vector<double>& difRatesPrelim, const std::vector<double>& decayRates, const std::vector<double>& smoothingH) {
   std::vector<double> difRates;

   if (smoothingH.front() >= 0.0) {
      if (smoothingH.size() != decayRates.size())
         throw std::runtime_error("List of smoothing lengths and list of decay rates (or half-lives) do not match");
      std::cerr << "One or more smoothing lengths (h) explicitly set. Ignoring any requested diffusion rates, setting diffusion rates to k*h^2." << std::endl;
      for (size_t i = 0; i < decayRates.size(); ++i)
         difRates.push_back(decayRates.at(i) * std::pow(smoothingH.at(i), 2.0));

   } else {
      if (difRatesPrelim.size() != decayRates.size())
         throw std::runtime_error("List of diffusion rates and list of decay rates (or half-lives) do not match");
      for (auto& difRate : difRatesPrelim)
         difRates.push_back(difRate);
   }

   std::cerr << "Diffusion rates...\n";
   for (auto& difRate : difRates) {
      constexpr double um2PerMm2 = 1.0e6;
      std::cerr << "  Diffusion rate: " << difRate << " microns^2 per second\n";
      difRate /= um2PerMm2;
   }

   return difRates;
}


std::vector<double> computeEnvDifRates(const std::vector<double>& envDifRatesPrelim, const std::vector<double>& envDecayRates) {
   std::vector<double> envDifRates;

   if (envDifRatesPrelim.size() != envDecayRates.size())
      throw std::runtime_error("List of environmental diffusion rates and list of decay rates (or half-lives) do not match");

   if (nEnvVars > 0) {
      std::cerr << "Environmental diffusion rates...\n";
      for (const auto& difRate : envDifRatesPrelim) {
         constexpr double um2PerMm2 = 1.0e6;
         std::cerr << "  Environmental diffusion rate: " << difRate << " microns^2 per second\n";
         envDifRates.push_back(difRate / um2PerMm2);
      }
   }

   return envDifRates;
}


Substances computeSubstances(const std::vector<double>& difRates, const std::vector<double>& decayRates, const std::vector<double>& gains, const std::vector<double>& envDifRates, const std::vector<double>& envDecayRates) {
   Substances subs;
   for (size_t i = 0; i < nSPHSubs; ++i) {
      subs.at(i).difRate = repeatBack<double>(difRates, i);
      subs.at(i).decayRate = repeatBack<double>(decayRates, i);
      subs.at(i).gain = repeatBack<double>(gains, i);
   }
   for (size_t i = 0; i < nEnvVars; ++i) {
      subs.at(i + nSPHSubs).difRate = repeatBack<double>(envDifRates, i);
      subs.at(i + nSPHSubs).decayRate = repeatBack<double>(envDecayRates, i);
      subs.at(i + nSPHSubs).gain = 1.0; // Gains probably don't make sense for environmental fields
   }
   return subs;
}


std::tuple<double, double, double> computeTimeParams(double dt, double minTimestepsPerCell, double maxspeedPrelim, double tlimitPrelim, double speedupFactor, const Substances& subs, double h, bool idealResults) {
   if (dt > Diffuser::maxDt(subs, h)) 
      std::cerr << "Explicit dt (" << dt << ") greater than maximum for diffusion (" 
         << Diffuser::maxDt(subs, h) << "), so multiple subcycles will be used." << std::endl;

   double maxspeed = 0.0;
   if (idealResults) {
      std::cerr << "Agent speed irrelevant for finding ideal results. Ignoring --maxspeed." << std::endl;
      if (dt <= 0.0) {
         dt = Diffuser::maxDt(subs, h);
         std::cerr << "dt not set above zero, so we will use the maximum for diffusion, which is " << Diffuser::maxDt(subs, h) << std::endl;
      }
   } else if (maxspeedPrelim > 0.0) {
      std::cerr << "Maximum speed for simulation explicitly set to " << maxspeedPrelim << std::endl;
      maxspeed = speedupFactor * maxspeedPrelim;
      if (speedupFactor != 1.0) std::cerr << "Speedup factor is " << speedupFactor << ", so true maximum speed will be " << maxspeed << "." << std::endl;
      if (dt <= 0.0) {
         const auto dtDesired = h / (maxspeed * minTimestepsPerCell);
         dt = std::min(dtDesired, Diffuser::maxDt(subs, h));
         if (dt < dtDesired) {
            std::cerr << "Based on maximum speed and minimum timesteps per grid cell traversal, dt would be " << dtDesired << std::endl;
            std::cerr << "However, to avoid multiple diffusion substeps per timestep, setting dt to " << dt << std::endl;
         } else {
            std::cerr << "Based on maximum speed and minimum timesteps per grid cell traversal, setting dt to " << dt << std::endl;
         }
      } else {
         std::cerr << "Both --maxspeed and --dt are set, so minTimestepsPerCell will be ignored for simulation." << std::endl;
         std::cerr << "Effective minTimestepsPerCell will be " << h / (maxspeed * dt) << std::endl;
      }
   } else {
      if (dt <= 0.0) {
         dt = Diffuser::maxDt(subs, h);
         std::cerr << "Neither --maxspeed nor --dt set above zero; maximum for diffusion (" << dt << ") will be used." << std::endl;
      }
      maxspeed = h / (dt * minTimestepsPerCell);
      std::cerr << "Maximum speed not explicitly set, so based on h, dt, and minTimestepsPerCell, it will be set to " << maxspeed << std::endl;
   }

   const auto tlimit = tlimitPrelim / speedupFactor;
   return {dt, maxspeed, tlimit};
}


double findMaxAngularSpeed(double gridh, const AgentManagerParams& amParams, double dt, double proportionAllowed, bool passthrough) {
   double maxAngularSpeed = 0.0;
   if (passthrough) {
      std::cerr << "No max angular speed for simulation." << std::endl;
      maxAngularSpeed = -1.0;
   } else {
      const double cellCoveringMaxAngularSpeed = findMaxRotationWithoutCoveringCellProportion(gridh, amParams, proportionAllowed) / dt;
      // 0.2 below is an arbitrary multiplier based on informal, subjective tests.
      // 1.0 would be a more absolute limit.
      const double fluidMaxAngularSpeed = 0.2 * gridh / (dt * amParams.semimajorAxis);
      if (cellCoveringMaxAngularSpeed >= 0.0) {
         std::cerr << "cellCoveringMaxAnguarSpeed: " << cellCoveringMaxAngularSpeed << "; fluidMaxAngularSpeed: " << fluidMaxAngularSpeed << std::endl;
         maxAngularSpeed = std::min(cellCoveringMaxAngularSpeed, fluidMaxAngularSpeed);
      } else {
         maxAngularSpeed = fluidMaxAngularSpeed;
      }
   }
   std::cerr << "Max angular speed for simulation: " << maxAngularSpeed << std::endl;
   return maxAngularSpeed;
}


unsigned int computeSeed(std::string seedString) {
   unsigned int seed = 0;
   if (seedString != "") {
      std::istringstream readSeed(seedString);
      readSeed >> seed;
   } else {
      std::random_device rd;
      seed = rd();
   }
   std::cout << "Seed: " << seed << "\n";
   return seed;
}


void printSubstanceInfo(const Substances& subs, double initDens, AgentManagerParams amParams) {
   std::cerr << "Smoothing parameters (sqrt(E/k), also called h, in microns) for each SPH substance: " << std::endl;
   for (size_t i = 0; i < nSPHSubs; ++i) {
      std::cerr << "  Substance " << i << ": ";
      if (subs.at(i).decayRate > 0.0) {
         const auto smoothingH = std::sqrt(subs.at(i).difRate / subs.at(i).decayRate);
         constexpr double umPerMm = 1000.0;
         std::cerr << smoothingH * umPerMm;
         if (initDens > 0.0) {
            //const auto nAgentsPerArea = initDens / (M_PI * amParams.semimajorAxis * amParams.semiminorAxis);
            const auto nAgentsPerArea = initDens / amParams.mass();
            const auto pihh = M_PI * std::pow(smoothingH, 2.0);
            // Found with Maxima code such as (assume(a > 0), find_root(quad_qags(x * bessel_k(0, x), x, 0, a)[1] = .95, a, 1e-15, 500))^2;
            // But there's also code in Calibrator for this (rOverHIncludeingProportionOfMSF()) if we want to avoid hard-coding it.
            std::cerr << "; # agents within 40%, 68%, 95% of kernel mass: " 
               << nAgentsPerArea * 1.0 * pihh << ", "
               << nAgentsPerArea * 3.4 * pihh << ", "
               << nAgentsPerArea * 16.0 * pihh << std::endl;
         }
      } else {
         std::cerr << "Decay rate is zero; no meaningful smoothing parameter exists." << std::endl;
      }
   }
   std::cerr << std::endl;
}


void writeInfoToReproduce(const std::string& picdir, const std::string& cmdLine, unsigned int seed) {
   const auto infoFilename = picdir + ".info";
   std::ofstream infoFile{infoFilename};
   infoFile << cmdLine << "--seed " << seed << std::endl;
   dl_iterate_phdr(dlPrint, &infoFile);
#ifdef EIGEN_MACROS_H
   infoFile << "Eigen version: " 
      << EIGEN_WORLD_VERSION << "."
      << EIGEN_MAJOR_VERSION << "."
      << EIGEN_MINOR_VERSION << std::endl;
#endif
   std::system(("gcc --version >> " + infoFilename).c_str());
   infoFile.close();
   std::filesystem::create_directory("tars");
   std::system(("tar czf tars/" + picdir + ".tgz src/*.cpp src/*.h src/*.hpp GNUmakefile " + infoFilename).c_str());
   std::filesystem::remove(infoFilename);
}


SecondsPerUpdate adjustSecondsPerUpdate(const SecondsPerUpdate& secondsPerUpdatePrelim, double speedupFactor) {
   return {
      secondsPerUpdatePrelim.live / speedupFactor,
      secondsPerUpdatePrelim.pic / speedupFactor,
      secondsPerUpdatePrelim.text / speedupFactor,
   };
}


void calibrate(
      const CalibInfo& calibInfo, 
      const Substances& subs, 
      double visc, 
      const AgentManagerParams& amParams, 
      double nCellWidthsPerDiameter, 
      double minTimestepsPerCell, 
      double liveSecondsPerUpdate, 
      int pmag, 
      const SimulationParams& simParams, 
      std::default_random_engine& gen) 
{
   if (calibInfo.nocalib) return;

   // Calibrate without noise.
   Calibrator calibrator(subs, visc, amParams.noNoise(), nCellWidthsPerDiameter, minTimestepsPerCell, liveSecondsPerUpdate, pmag, simParams.passthrough, gen);

   // Create set of speeds for which to calibrate.
   std::set<double> speeds;
   if (calibInfo.onecalib >= 0.0) {
      speeds.insert(calibInfo.onecalib);
   } else {
      speeds.insert(0.0);
      const double minSpeed = 0.00001;
      //constexpr double maxSpeedBiological = 7.0;
      const double maxSpeedPractical = 0.65536;
      const double maxSpeedCalibration = maxSpeedPractical;
      for (auto s = minSpeed; s <= maxSpeedCalibration; s *= 2) speeds.insert(s);
      speeds.insert(maxSpeedCalibration);
      std::cerr << "Max speed for calibration: " << maxSpeedCalibration << std::endl;
   }
   // Read each speed into a stringstream and back out again, potentially
   // reducing precision, to ensure that all speeds will compare as equal
   // when read back in from files during calibration.  Alternatively,
   // we could just print max_digits10 to the files, but it's ugly to
   // use so many digits, and we don't need that kind of precision. C++20
   // will introduce std::format to fix this.
   std::set<double> roundedSpeeds;
   for (auto speed : speeds) {
      std::stringstream buf;
      buf << speed;
      buf >> speed;
      roundedSpeeds.insert(speed);
   }

   if (calibInfo.nonewcalib) {
      Agent::cdat = calibrator.readCalibDataFromFile(std::set<double>{});

      bool highestNeededFound = false;
      bool insufficientData = false;
      for (auto speed : roundedSpeeds) {
         if (!highestNeededFound) {
            if (Agent::cdat.find(speed) == Agent::cdat.end())
               insufficientData = true;
            if (speed > simParams.maxSpeed)
               highestNeededFound = true;
         }
      }
      if (insufficientData)
         throw InsufficientCalibration("Insufficient calibration data found; refraining from creating more.");

   } else {
      try {
         Agent::cdat = calibrator.calibrate(roundedSpeeds, calibInfo.recalib, calibInfo.onecalib >= 0.0);
      } catch (SDLQuitException& e) {
         std::cerr << "SDL quit during calibration. Exiting. Message: " << e.what() << std::endl;
         throw;
      }
   }
}


void runSimulation(
      World& world, 
      double dt, 
      double tlimit, 
      const std::string& picdir, 
      const SecondsPerUpdate& secondsPerUpdate, 
      const SubstanceVars& fixedMaxVals, 
      int pmag) 
{
   // Optional objects for graphical output.
   std::optional<Screen> screen = secondsPerUpdate.live < 0 ? std::nullopt : std::make_optional<Screen>(world.provideFields(), fixedMaxVals, pmag, picdir);
   std::optional<Picture> pic = secondsPerUpdate.pic < 0 ? std::nullopt : std::make_optional<Picture>(world.provideFields(), fixedMaxVals, "pics/" + picdir);

   // Set up timing information to inform user of time elapsed and estimated remaining.
   const auto startTime = std::chrono::steady_clock::now();
   auto lastReportedTime = startTime;
   auto t = 0.0;
   auto lastReportedSimTime = t;
   //int id = 0;
   try {
      while (tlimit < 0.0 ? world.nAgentsRemaining() > 0 : t <= tlimit) {
         if (whetherToUpdate(t, secondsPerUpdate.live, dt) && screen) screen->update();
         if (whetherToUpdate(t, secondsPerUpdate.pic, dt) && pic) pic->update();
         if (whetherToUpdate(t, secondsPerUpdate.text, dt)) {
            world.textupdate(t);
            if (world.nAgentsRemaining() > 0) world.printAgentStats();
         }

         const auto currentTime = std::chrono::steady_clock::now();
         const std::chrono::duration<double> sinceLastReported = currentTime - lastReportedTime;
         constexpr auto printInterval = 30.0;
         if (sinceLastReported.count() > printInterval) {
            const std::chrono::duration<double> elapsedSoFar = currentTime - startTime;
            std::cerr << "Minutes elapsed so far: " << elapsedSoFar.count() / 60.0 << "; estimated remaining: ";
            if (tlimit >= 0.0) {
               const auto proportionOfLimit = t / tlimit;
               const auto estimatedRemaining = elapsedSoFar.count() * (1.0 / proportionOfLimit - 1.0);
               std::cerr << estimatedRemaining / 60.0;
            } else {
               std::cerr << "unknown";
            }
            std::cerr << "; real time per simulation time: " << sinceLastReported.count() / (t - lastReportedSimTime) << " (recent), " << elapsedSoFar.count() / t << " (overall)" << std::endl;
            lastReportedTime = currentTime;
            lastReportedSimTime = t;
         }
         world.update(dt);
#if PROGRAM == GROWTH
#ifdef CUT
         constexpr auto cuttime = 43200;
         if (t >= cuttime && t < cuttime + dt)
            world.cut();
#endif
#endif
         t += dt;
      }
   } catch (SDLQuitException& e) {
      std::cerr << "SDL quit during simulation. Ending loop. Message: " << e.what() << std::endl;
   }

   const std::chrono::duration<double> elapsedTotal = std::chrono::steady_clock::now() - startTime;
   std::cerr << "Seconds elapsed total: " << elapsedTotal.count() << std::endl;
}


// Mostly a copy of runSimulation(). The two could be combined in the
// future, but all the code for "ideal results" is beside the main point
// of this simulator.
void runSimulationFDM(
      const Fparams& fparams,
      const Substances& subs,
      double targetLaplacian,
      double radius,
      double dt, 
      double tlimit, 
      const std::string& picdir, 
      const SecondsPerUpdate& secondsPerUpdate, 
      const SubstanceVars& fixedMaxVals, 
      int pmag) 
{
   Fields fields{fparams, true};
   DiffuserFDM diffuserFDM{fields, subs, targetLaplacian, radius};

   std::optional<Screen> screen = secondsPerUpdate.live < 0 ? std::nullopt : std::make_optional<Screen>(fields, fixedMaxVals, pmag, picdir);
   std::optional<Picture> pic = secondsPerUpdate.pic < 0 ? std::nullopt : std::make_optional<Picture>(fields, fixedMaxVals, "pics/" + picdir);

   const auto startTime = std::chrono::steady_clock::now();
   auto lastReportedTime = startTime;
   double t = 0.0;
   try {
      while (tlimit < 0.0 || t <= tlimit) {
         if (whetherToUpdate(t, secondsPerUpdate.live, dt) && screen) {
            screen->update();
            std::cerr << "t: " << t << "; max: " << fields.grids.at(DiffStates::one).max() << std::endl;
         }
         if (whetherToUpdate(t, secondsPerUpdate.pic, dt) && pic) pic->update();

         const auto currentTime = std::chrono::steady_clock::now();
         const std::chrono::duration<double> sinceLastReported = currentTime - lastReportedTime;
         if (sinceLastReported.count() > 30.0) {
            lastReportedTime = currentTime;
            const std::chrono::duration<double> elapsedSoFar = currentTime - startTime;
            std::cerr << "Minutes elapsed so far: " << elapsedSoFar.count() / 60.0 << "; estimated remaining: ";
            if (tlimit >= 0.0) {
               const auto proportionOfLimit = t / tlimit;
               const auto estimatedRemaining = elapsedSoFar.count() * (1.0 / proportionOfLimit - 1.0);
               std::cerr << estimatedRemaining / 60.0 << std::endl;
            } else {
               std::cerr << "unknown" << std::endl;
            }
         }

         diffuserFDM.diffuse(dt),
         t += dt;
      }
   } catch (SDLQuitException& e) {
      std::cerr << "SDL quit during simulation. Ending loop. Message: " << e.what() << std::endl;
   }

   const std::chrono::duration<double> elapsedTotal = std::chrono::steady_clock::now() - startTime;
   std::cerr << "Seconds elapsed total: " << elapsedTotal.count() << std::endl;
   fields.outputField("fields/" + picdir, DiffStates::one);
}


int dlPrint(dl_phdr_info* info, size_t, void* outfile) {
   const std::filesystem::path libpath{info->dlpi_name};
   *static_cast<std::ofstream*>(outfile) << std::filesystem::weakly_canonical(libpath) << std::endl;
   return 0;
}


template <typename T>
T repeatBack(const std::vector<T>& vec, size_t i) {
   if (i < vec.size())
      return vec.at(i);
   else
      return vec.back();
}


int subcyclesNeeded(double dt, double maxdt) {
   double dsubcycles = std::ceil(dt / maxdt);
   if (dsubcycles > std::numeric_limits<int>::max())
      throw std::runtime_error("Subcycles needed greater than max integer value.");
   const int subcycles = (int)dsubcycles;
   std::cerr << "  " << subcycles << " subcycles will be needed for each timestep.\n";
   return subcycles;
}


double halfMinToDecaySec(double half) {
   return std::log(2.0) / (60.0 * half);
}


double findMaxRotationWithoutCoveringCellProportion(double gridh, const AgentManagerParams& amParams, double proportionAllowed) {
   if (!rotationCoversCellProportion(gridh, amParams, M_PI / 2.0, proportionAllowed)) return -1.0;

   double highAngle = M_PI / 2.0;
   double lowAngle = 0.0;
   double testAngle = 0.0;
   while ((highAngle - lowAngle) / highAngle > 0.01) {
      testAngle = (highAngle + lowAngle) / 2.0;
      //std::cerr << "low: " << lowAngle << "; high: " << highAngle << "; test: " << testAngle << std::endl;
      if (rotationCoversCellProportion(gridh, amParams, testAngle, proportionAllowed))
         highAngle = testAngle;
      else
         lowAngle = testAngle;
   }
   return testAngle;
}


bool rotationCoversCellProportion(double gridh, const AgentManagerParams& amParams, double angle, double proportion) {
   constexpr auto passthrough = true;
   const auto width = std::ceil(3 * amParams.semimajorAxis / gridh) * gridh;
   Fparams fp{width, width, gridh};
   Fields fields{fp, passthrough};
   std::default_random_engine gen{}; // Dummy; won't be used because we're only calling AgentManager::rotateOnly()
   AgentManager am{amParams, AMExternalState{Point{fp.xphys / 2.0, fp.yphys / 2.0}, 0.0}, fields, Substances{}, SimulationParams{true, -1.0, -1.0}, gen};
   Grid<double> oldAgentBodies{fields.agentBodies};
   am.rotateOnly(angle);
   return (fields.agentBodies.chebyshev(oldAgentBodies) >= proportion);
}


bool growthRectangle(const Point& location, const Fparams& fp, double height) {
   const double xCenter = fp.xphys / 2.0;
   const double width = 3.0;
   const double lowerEdgeY = 1.0;
   const double yRelativeToLowerEdge = location.y - lowerEdgeY;
   return std::abs(location.x - xCenter) < width / 2.0 && yRelativeToLowerEdge > 0.0 && yRelativeToLowerEdge < height;
}


bool growxShape(const Point& location, const Fparams& fp) {
   const double xCenter = fp.xphys / 2.0;
   const double yCenter = fp.yphys / 2.0;
   const double width = 2.0;
   const double lowerMargin = 1.0;
   const bool goodDistFromBottom = (location.y > lowerMargin && location.y < lowerMargin + width);
   const bool goodDistFromSide = (location.x < fp.xphys - lowerMargin && location.x > fp.xphys - lowerMargin - width);
   const bool inColumn = std::abs(location.x - xCenter) < width / 2.0;
   const bool inRow = std::abs(location.y - yCenter) < width / 2.0;
   return (goodDistFromBottom && inColumn) || (goodDistFromSide && inRow);
}


bool circle(const Point& location, const Fparams& fp, double radius) {
   const Point center{Point{fp.xphys, fp.yphys} / 2.0};
   const Point relativeToCenter = location - center;
   return relativeToCenter.mag() < radius;
}


bool twoCircles(const Point& location, const Fparams& fp, double radius) {
   const auto lcen = Point{(1.0 / 6.0) * fp.xphys, 0.5 * fp.yphys};
   const auto rcen = Point{(5.0 / 6.0) * fp.xphys, 0.5 * fp.yphys};
   const auto ldist = (location - lcen).mag();
   const auto rdist = (location - rcen).mag();
   return ldist < radius || rdist < radius;
}
