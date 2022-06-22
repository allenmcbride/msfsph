// calibrator.h
// Allen McBride
// March 10, 2022
//
// The Calibrator class takes a set of speeds, checks for files containing
// existing calibration samples, and generates new samples until enough
// are found. Statistics are then computed, and a map with calibration
// data returned for use at runtime. After each sample, files are
// re-checked, in case another process is also computing calibration
// data at the same time. See section 3.4.2 for more on calibration.


#include "agent.h"
#include "agentmanager.h"
#include "fields.h"
#include "screen.h"
#include "substances.h"

#include <set>
#include <fstream>
#include <random>
#include <string>
#include <valarray>


// An array of calibration values, one for each substance. std::valarray
// is semantically convenient for doing arithmetic with respect to each
// substance at once.
typedef std::valarray<double> CalibValBySub;


class TooManyIterationsNeeded : public std::runtime_error {
   public: using runtime_error::runtime_error;
};


// Information for a single calibration sample.
struct CalibSample {
   CalibRecord record;
   double duration;
   double yphys;
   double dt;
   friend CalibSample operator+(const CalibSample& lhs, const CalibSample& rhs) {
      return {lhs.record + rhs.record, lhs.duration + rhs.duration, lhs.yphys + rhs.yphys, lhs.dt + rhs.dt};
   }
   friend CalibSample operator/(const CalibSample& a, const double k) {
      return {a.record / k, a.duration / k, a.yphys / k, a.dt / k};
   }
};


// For convenience in storing all the statistics we need to compute
// with respect to calibration samples.
struct SampleStats {
   CalibRecord mean;
   CalibRecord stdDev;
   CalibRecord stdErr;
   double meanDuration;
   double meanYphys;
   double meanDt;
};


class Calibrator {
   private:
      // Hard-coded parameters relating to accuracy for a given
      // calibration run.
      static constexpr double propOfOriginalToAllow = 0.01;
      static constexpr double viscosityTolerance = 0.05;
      static constexpr double pressureTolerance = 0.01;

      // Hard-coded number of samples needed for zero and nonzero speeds
      // for crawling ("Pass") and swimming ("Nopass") models.
      static constexpr size_t nsampPassNonzero = 36;
      static constexpr size_t nsampPassZero = 4 * nsampPassNonzero;
      static constexpr size_t nsampNopassNonzero = 36;
      static constexpr size_t nsampNopassZero = 4 * nsampNopassNonzero;

      // How far ahead should we look before deciding pulse is
      // sufficiently hard to detect that we have found steady state? This
      // is the same as the ratio of t2 to t1 in eq. 3.20.
      static constexpr double multOfSteadyTimeToLook = 2.0;

      // Subclass to integrate pulse, as in eq. 3.19.
      class Integrator {
         private:
            const Substances& subs;
            CalibValBySub subsDecayRates;
            CalibValBySub lastIntegral;
            CalibValBySub candidateCutoffTimes;
            std::map<double, CalibValBySub> calibIntegrals;
         public:
            Integrator(const Substances& subs);
            bool addAndCheck(CalibBySub calibBySub, double t, double dt, bool debug);
            std::array<double, nSPHSubs> provideCutoffTimes() const;
      };

      Substances subs;
      double visc;
      AgentManagerParams amParams;
      double maxTimestepsPerCell;
      std::set<double> speeds;
      double liveSecondsPerUpdate;
      int pmag;
      bool passthrough;
      std::default_random_engine& gen;
      double gridh; // Grid cell size
      std::array<int, nSPHSubs> dupLocations;

      size_t nsamplesNeeded(double speed) const { return passthrough ? nsamplesNeededPass(speed) : nsamplesNeededNopass(speed); }
      static size_t nsamplesNeededPass(double speed) { return speed == 0.0 ? nsampPassZero : nsampPassNonzero; }
      static size_t nsamplesNeededNopass(double speed) { return speed == 0.0 ? nsampNopassZero : nsampNopassNonzero; }

      double calcYphys(const Substances& subsForSpeed, double speed, double duration) const;
      AMExternalState calcStartState(double speed, Fparams fparams, double duration) const;
      std::string calcFilename(int isub, double speed = -1.0) const; // -1.0 speed means don't include speed in name
      void deleteOldData(const std::set<double>& speeds, bool onecalib) const;
      std::map<double, CalibRecord> readFile(std::string filename) const;
      void processSamples(double speed) const;
      double findDuration(const Substances& subsForSpeed, double speed) const;
      void calcDurationWithPulse(const Substances& subsForSpeed, double speed) const;
      double calibOneSpeed(double speed, const Substances& subsForSpeed, Fparams fparams, double duration, bool pulse) const;

   public:
      static constexpr double rootFindingAccuracy = 0.999;
      
      Calibrator(
            Substances subs,
            double visc,
            AgentManagerParams amParams,
            double nCellWidthsPerDiameter,
            double maxTimestepsPerCell,
            double liveSecondsPerUpdate,
            int pmag,
            bool passthrough,
            std::default_random_engine& gen
            );
      CalibData calibrate(std::set<double> speeds, bool recalib, bool onecalib) const;
      CalibData readCalibDataFromFile(const std::set<double>& speeds) const;
};
