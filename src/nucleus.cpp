// nucleus.cpp
// Allen McBride
// June 21, 2022
//
// The most important aspect of this file is that it contains
// implementations of the update() function. This is where the bulk of
// compiled Morphgen code resides. Other notes are in nucleus.h.

#include "nucleus.h"

#include "agent.h" // DEBUG
#include "agentinterface.h"
#include "agentmanager.h" // DEBUG
#include "point.h"
#include "substances.h"

#include <algorithm>
#include <chrono> // DEBUG
#include <cmath>
#include <iostream> // DEBUG
#include <random> // DEBUG
#include <tuple>
#include <utility>


int Nucleus::nextid = 0;


std::pair<double, double> toDualRail(double val);
double saturate(double val, double cap);
Point saturate(Point val, double cap);
double sigmoid(double condition, double threshold, double steepness);
double sigmoidClamped(double condition, double threshold, double steepness);


// This is to simplify setting references in field equations, given
// difference between indexing in NuclearOutput and SPHVals
double& Nucleus::fieldRef(int isub) {
   if (isub >= nDiffStates && isub < nSPHSubs) return outVals.fieldVals.at(isub - nDiffStates);
   else throw std::invalid_argument("Nucleus::fieldRef called on non-field variable");
}


#if PROGRAM == DIFFUSION
void Nucleus::update(double dt, bool debug) {
   const auto rho = ai.findVal(DiffStates::one);
   const auto gradrho = ai.findGrad(DiffStates::one);
   outVals.velocity = rho == 0.0 ? Point{} : -0.0001 * gradrho / rho;
}


#elif PROGRAM == DIFFUSION_TRANSLATION
void Nucleus::update(double dt, bool debug) {
   const auto rho = ai.findVal(DiffStates::one);
   const auto gradrho = ai.findGrad(DiffStates::one);
   const auto gravdir = ai.findGravityDir();
   const auto translation = (1.0 / 2880.0) * gravdir.rotated(M_PI / 2.0);
   const auto diffusion = rho == 0.0 ? Point{} : -0.0001 * gradrho / rho;
   outVals.velocity = translation + diffusion;
}


#elif PROGRAM == PARABOLOID
void Nucleus::update(double dt, bool debug) {
   auto& cena = fieldRef(FieldVars::cena);
   const auto laplcena = ai.findLaplacian(FieldVars::cena);
   const auto circle = ai.findCue(Cues::circle);
   if (circle)
      cena += dt * 0.001 * (laplcena + 1.0);
   else
      cena = 0.0;
}


#elif PROGRAM == ADVECTION
void Nucleus::update(double dt, bool debug) {
   const auto gravdir = ai.findGravityDir();
   outVals.velocity = 0.001 * gravdir.rotated(M_PI / 2.0);
   auto& cena = fieldRef(FieldVars::cena);
   const auto gradcena = ai.findGrad(FieldVars::cena);
#ifdef NOUPWIND
   cena += dt * outVals.velocity.dot(gradcena);
#else
   cena += dt * ai.velDotGradUpstream(FieldVars::cena, outVals.velocity);
#endif
   if (cena > 1.0) cena += dt * 0.1 * (1.0 - cena);
}


#elif PROGRAM == MIXTURE
void Nucleus::update(double dt, bool debug) {
   auto& a = outVals.diffStateVals.at(DiffStates::a);
   auto& b = outVals.diffStateVals.at(DiffStates::b);
   auto& cena = fieldRef(FieldVars::cena);
   auto& cenb = fieldRef(FieldVars::cenb);
   const auto rhoa = ai.findVal(DiffStates::a);
   const auto rhob = ai.findVal(DiffStates::b);
   const auto gradrhoa = ai.findGrad(DiffStates::a);
   const auto gradrhob = ai.findGrad(DiffStates::b);
   const auto gradcena = ai.findGrad(FieldVars::cena);
   const auto gradcenb = ai.findGrad(FieldVars::cenb);
   const auto laplcena = ai.findLaplacian(FieldVars::cena);
   const auto laplcenb = ai.findLaplacian(FieldVars::cenb);
   const auto gravdir = ai.findGravityDir();
   const bool left = ai.findCue(Cues::left);
   const bool right = ai.findCue(Cues::right);

   outVals.velocity = (cena > 0.0000001 || cenb > 0.0000001) * 0.00001 * (a - b) * gravdir.rotated(M_PI / 2);

   if (a > 0.9) {
      constexpr auto targ = 0.03;
      if (internal.timer < 0.99)
         internal.timer += dt * 0.05 * rhoa * std::exp(-200 * gradrhoa.mag()) * (1.0 - internal.timer) * (cena < 0.00001);
      if (internal.timer > 0.99) {
         cena += dt * 0.1 * (1.0 - cena);
      } else {
         cena += dt * 0.5 * laplcena;
         if (cena > 0.0000001 && rhoa <= targ) 
            outVals.velocity += 0.01 * (targ - gradrhoa.mag()) * gradcena.normed();
         if (ai.edgeDetect(DiffStates::a) > 0.5) {
            cena = 0.0;
         }
      }
      if (rhoa > targ)
         outVals.velocity += 0.05 * (targ - gradrhoa.mag()) * gradrhoa.normed();
      if (ai.edgeDetect(DiffStates::a) > 0.5)
         outVals.smoothDensityGuards.at(DiffStates::a) = 0.99;
      else 
         outVals.smoothDensityGuards.at(DiffStates::a) = 0.99;
   }
   if (b > 0.9) {
      constexpr auto targ = 0.04;
      if (internal.timer < 0.99)
         internal.timer += dt * 0.1 * rhob * std::exp(-200 * gradrhob.mag()) * (1.0 - internal.timer) * (cenb < 0.00001);
      if (internal.timer > 0.99) {
         cenb += dt * 0.1 * (1.0 - cenb);
      } else {
         cenb += dt * 0.5 * laplcenb;
         if (cenb > 0.0000001 && rhob <= targ) 
            outVals.velocity += 0.01 * (targ - gradrhob.mag()) * gradcenb.normed();
         if (ai.edgeDetect(DiffStates::b) > 0.5) {
            cenb = 0.0;
         }
      }
      if (rhob > targ)
         outVals.velocity += 0.1 * (targ - gradrhob.mag()) * gradrhob.normed();
      if (ai.edgeDetect(DiffStates::b) > 0.5)
         outVals.smoothDensityGuards.at(DiffStates::b) = 0.99;
      else 
         outVals.smoothDensityGuards.at(DiffStates::b) = 0.99;
   }

}


#elif PROGRAM == GROWX
void Nucleus::update(double dt, bool debug) {
   outVals.velocity = Point{};
   auto& a = outVals.diffStateVals.at(DiffStates::a);
   auto& b = outVals.diffStateVals.at(DiffStates::b);
   const auto rhoa = ai.findVal(DiffStates::a);
   const auto rhob = ai.findVal(DiffStates::b);
   const auto gradrhoa = ai.findGrad(DiffStates::a);
   const auto gradrhob = ai.findGrad(DiffStates::b);
   const auto gravdir = ai.findGravityDir();
   auto& meristema = fieldRef(FieldVars::meristema);
   auto& meristemb = fieldRef(FieldVars::meristemb);
   auto& colora = fieldRef(FieldVars::colora);
   auto& colorb = fieldRef(FieldVars::colorb);
   const auto gradcolora = ai.findGrad(FieldVars::colora);
   const auto gradcolorb = ai.findGrad(FieldVars::colorb);
   const auto laplmeristema = ai.findLaplacian(FieldVars::meristema);
   const auto laplmeristemb = ai.findLaplacian(FieldVars::meristemb);

   colora = 1.0;
   colorb = 1.0;

   internal.timer += dt;
   const auto afterTimer = internal.timer > 10800;
   auto velStem = Point{};
   auto velGradrho = Point{};
   if (a > 0.5) {
      const auto dir = gravdir;
      const bool edge = ai.colorEdge(FieldVars::colora) > (afterTimer ? 0.2 : 0.3);
      double relAngle = gradcolora.angle() - dir.angle();
      if (relAngle < 0.0) relAngle += 2.0 * M_PI;
      if ((relAngle < 0.4 * M_PI || relAngle > 1.6 * M_PI) && edge) {
         meristema += dt * 0.1 * (1.0 - meristema);
      } else if ((relAngle < 1.4 * M_PI && relAngle > 0.6 * M_PI) && edge) {
         meristema += dt * 0.1 * (0.0 - meristema);
      } else {
         meristema += dt * 0.02 * laplmeristema;
      }
      outVals.growthSwitches.at(DiffStates::a) = afterTimer && !edge;
      velStem = afterTimer ? -0.0001 * meristema * dir.normed() : Point{};
      velGradrho = edge ? Point{} : -0.0008 * gradrhoa;
   }
   if (b > 0.5) {
      const auto dir = gravdir.rotated(M_PI / 2.0);
      const bool edge = ai.colorEdge(FieldVars::colorb) > (afterTimer ? 0.2 : 0.3);
      double relAngle = gradcolorb.angle() - dir.angle();
      if (relAngle < 0.0) relAngle += 2.0 * M_PI;
      if ((relAngle < 0.4 * M_PI || relAngle > 1.6 * M_PI) && edge) {
         meristemb += dt * 0.1 * (1.0 - meristemb);
      } else if ((relAngle < 1.4 * M_PI && relAngle > 0.6 * M_PI) && edge) {
         meristemb += dt * 0.1 * (0.0 - meristemb);
      } else {
         meristemb += dt * 0.02 * laplmeristemb;
      }
      outVals.growthSwitches.at(DiffStates::b) = afterTimer && !edge;
      velStem = afterTimer ? -0.0001 * meristemb * dir.normed() : Point{};
      velGradrho = edge ? Point{} : -0.0008 * gradrhob;
   }

   outVals.velocity = saturate(velStem + velGradrho, 0.01);

   outVals.growthRates.at(DiffStates::a) = 0.0000;
   outVals.growthRates.at(DiffStates::b) = 0.0000;
}


#elif PROGRAM == GROWTH
void Nucleus::update(double dt, bool debug) {
   outVals.velocity = Point{};
   const auto rho = ai.findVal(DiffStates::one);
   const auto gradrho = ai.findGrad(DiffStates::one);
   const auto gravdir = ai.findGravityDir();
   auto& meristem = fieldRef(FieldVars::meristem);
   auto& color = fieldRef(FieldVars::color);
   const auto gradcolor = ai.findGrad(FieldVars::color);
   const auto laplmeristem = ai.findLaplacian(FieldVars::meristem);

   color = 1.0;

   internal.timer += dt;
   const auto afterTimer = internal.timer > 10800;

   const bool edge = ai.colorEdge(FieldVars::color) > (afterTimer ? 0.2 : 0.3);
   const double relAngle = gradcolor.angleBetween(gravdir);
   if (relAngle < 0.4 * M_PI && edge) {
      meristem += dt * 0.1 * (1.0 - meristem);
   } else if (relAngle > 0.6 * M_PI && edge) {
      meristem += dt * 0.1 * (0.0 - meristem);
   } else {
      meristem += dt * 0.02 * laplmeristem;
   }

   const auto velStem = afterTimer ? -0.0001 * meristem * gravdir.normed() : Point{};
   outVals.growthSwitches.at(DiffStates::one) = afterTimer && !edge;

   const auto velGradrho = edge ? Point{} : -0.0004 * gradrho;
   const auto velGradrhoSum = velGradrho;

   outVals.velocity = saturate(velStem + velGradrhoSum, 0.01);

   outVals.growthRates.at(DiffStates::one) = 0.0000;
}


#elif PROGRAM == PATH
void Nucleus::update(double dt, bool debug) {
   const auto rho = ai.findVal(DiffStates::one);
   const auto sk = fieldRef(FieldVars::seeker);
   const auto tr = fieldRef(FieldVars::tracer);
   const auto mk = fieldRef(FieldVars::marker);
   const auto sc = ai.findCue(Cues::source);
   const auto tg = ai.findCue(Cues::target);
   const auto gradrho = ai.findGrad(DiffStates::one);
   const auto gradsk = ai.findGrad(FieldVars::seeker);
   const auto gradtr = ai.findGrad(FieldVars::tracer);
   const auto gradmk = ai.findGrad(FieldVars::marker);
   const auto laplsk = ai.findLaplacian(FieldVars::seeker);
   const auto lapltr = ai.findLaplacian(FieldVars::tracer);
   const auto laplmk = ai.findLaplacian(FieldVars::marker);

   const auto oppositeGradientScore = -gradsk.dot(gradmk) / 0.0001;
   
   if (internal.lock < 1.0e-7) {
      outVals.velocity = -0.0001 * gradrho.normed();
      outVals.smoothDensityGuards.at(DiffStates::one) = 0.0;
      outVals.smoothVelocityGuards.at(DiffStates::one) = 0.0;
   } else {
      const auto velGradmk = gradmk.normed();

      constexpr auto lo = 0.01;
      constexpr auto hi = 0.05;

      const auto angleScore = sigmoidClamped(oppositeGradientScore, 0.5, 0.9);
      const auto edgeScore = 1.0 - sigmoidClamped(ai.edgeDetect(DiffStates::one), 0.5, 0.8);
      const double maxGradrhoOverRho = hi - (hi - lo) * std::max(edgeScore, angleScore);

      outVals.smoothDensityGuards.at(DiffStates::one) = std::exp(-maxGradrhoOverRho);
      constexpr auto slo = 0.99;
      constexpr auto shi = 0.92;
      outVals.smoothVelocityGuards.at(DiffStates::one) = shi - (shi - slo) * sigmoidClamped(mk, 0.0, 0.90);
      outVals.velocity = velGradmk;
   }

   if (sc) {
      fieldChange(FieldVars::seeker, 0.1 * (1.0 - sk), dt);
      outVals.smoothVelocityGuards.at(DiffStates::one) = 0.999;
      outVals.smoothDensityGuards.at(DiffStates::one) = std::exp(-0.003);
   } else {
      fieldChange(FieldVars::seeker, 0.027 * laplsk - 0.0002 * sk, dt);
   }

   const auto currentVel = Point{ai.speed(), 0.0}; // Agents always moving forward
   if (tg) {
      if (sk > 0.0) 
         fieldChange(FieldVars::tracer, 0.1 * (1.0 - tr), dt);
   } else {
      if (gradsk.mag() > 0.000 && sk > 0.0) {
         const double advInc = ai.velDotGradUpstream(FieldVars::tracer, -0.0007 * gradsk.normed() + currentVel);
         fieldChange(FieldVars::tracer, advInc, dt);
      }
   }

   if (tr >= internal.lock)
      internal.lock += dt * 0.1 * (tr - internal.lock); 
   else
      internal.lock -= dt * 0.000002 * internal.lock;

   if (internal.lock >= mk) {
      fieldChange(FieldVars::marker, 0.1 * (internal.lock - mk), dt);
   } else {
      fieldChange(FieldVars::marker, 0.00027 * laplmk - 0.00009 * mk, dt); // somewhere between .0001 and .0005 seeming best
   }

   if (tr > 1.0) fieldChange(FieldVars::tracer, 0.1 * (1.0 - tr), dt);
   if (mk > 1.0) fieldChange(FieldVars::marker, 0.1 * (1.0 - mk), dt);
}


#else // e.g., CALIBRATE
void Nucleus::update(double dt, bool debug) {
}


#endif


std::pair<double, double> toDualRail(double val) {
   return (val > 0.0) ? std::make_pair(val, 0.0) : std::make_pair(0.0, -val);
}


double saturate(double val, double cap) {
   if (cap == 0.0) return 0.0;
   return cap * std::tanh(val / cap);
}


Point saturate(Point val, double cap) {
   const double mag = val.mag();
   const Point saturated = mag == 0.0 ? Point{} : saturate(mag, cap) * val / mag;

   return saturated;
}


double sigmoid(double condition, double threshold, double steepness) {
   return std::tanh(steepness * (condition - threshold)) * 0.5 + 0.5;
}


// Based on https://www.desmos.com/calculator/3zhzwbfrxd
double sigmoidClamped(double condition, double threshold, double steepness) {
   if (steepness < 0.0 || steepness >= 1.0)
      throw std::invalid_argument("Steepness of " + std::to_string(steepness) + " invalid in sigmoidClamped()");
   if (threshold < 0.0 || threshold > 1.0)
      throw std::invalid_argument("Threshold of " + std::to_string(threshold) + " invalid in sigmoidClamped()");
   if (condition < 0.0) return 0.0;
   if (condition > 1.0) return 1.0;
   const auto exponent = 2.0 / (1.0 - steepness) - 1.0;
   if (condition <= threshold && threshold != 0.0) return std::pow(condition, exponent) / std::pow(threshold, exponent - 1.0);
   else return 1.0 - std::pow(1.0 - condition, exponent) / std::pow(1.0 - threshold, exponent - 1.0);
}
