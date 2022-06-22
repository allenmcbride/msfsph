// substances.h
// Allen McBride
// June 21, 2022
//
// This file contains the actual Substance struct along with
// other constants and structs that only make sense in contexts
// involving substances. This is the other main location, in addition
// to nucleus.cpp, containing code resulting from the compilation of
// Morphgen code. Specifically, this file defines the mapping between
// Morphgen field variables and physical morphogens.

#ifndef SUBSTANCES_H
#define SUBSTANCES_H

#include <array>
#include <bitset>
#include <cmath>

#define GROWTH 1
#define DIFFUSION 2
#define DIFFUSION_TRANSLATION 3
#define PATH 4
#define ADVECTION 5
#define MIXTURE 6
#define CALIBRATE 7
#define PARABOLOID 8
#define GROWX 9
#define EDGE 10

#if PROGRAM == ADVECTION || PROGRAM == DIFFUSION_TRANSLATION
#define XWRAP
#endif

constexpr bool isf = false;

// DEBUG
constexpr bool debugflow = false;
constexpr bool debugspin = false;
constexpr double debugspeed = 10.0;
constexpr int debugndots = 10;

struct Substance {
   double difRate;
   double decayRate;
   double gain;
   bool include = true;
   double timeToSteady(double propRemaining) const { return -std::log(propRemaining) / decayRate; }
   double smoothingLength() const { return std::sqrt(difRate / decayRate); }
};

// Cues are just binary environmental signals; either an agent detects it or not.
// EnvVars are passive but diffuse and decay like morphogens and have continuous levels that an agent can detect.

#if PROGRAM == GROWTH

struct Cues { enum { }; };
constexpr int nCues = 0;
struct DiffStates { enum { one }; };
constexpr int nDiffStates = 1;
struct FieldVars  { enum { meristem = nDiffStates, color }; };
constexpr int nFieldVars = 2;
constexpr int diffVarForFieldVar(int fieldVar) {
   constexpr std::array<int, nFieldVars> lookup = { DiffStates::one, DiffStates::one };
   return lookup.at(fieldVar - nDiffStates);
}
struct EnvVars    { enum { }; };
constexpr int nEnvVars = 0;

#elif PROGRAM == GROWX

struct Cues { enum { }; };
constexpr int nCues = 0;
struct DiffStates { enum { one, a, b }; };
constexpr int nDiffStates = 3;
struct FieldVars  { enum { meristema = nDiffStates, meristemb, colora, colorb }; };
constexpr int nFieldVars = 4;
constexpr int diffVarForFieldVar(int fieldVar) {
   constexpr std::array<int, nFieldVars> lookup = { DiffStates::a, DiffStates::b, DiffStates::a, DiffStates::b };
   return lookup.at(fieldVar - nDiffStates);
}
struct EnvVars    { enum { }; };
constexpr int nEnvVars = 0;

#elif PROGRAM == DIFFUSION || PROGRAM == DIFFUSION_TRANSLATION

struct Cues { enum { }; };
constexpr int nCues = 0;
struct DiffStates { enum { one }; };
constexpr int nDiffStates = 1;
struct FieldVars { enum { }; };
constexpr int nFieldVars = 0;
constexpr int diffVarForFieldVar(int fieldVar) {
   constexpr std::array<int, nFieldVars> lookup = { };
   return lookup.at(fieldVar - nDiffStates);
}
struct EnvVars    { enum { }; };
constexpr int nEnvVars = 0;

#elif PROGRAM == PATH

struct Cues { enum { source, target }; };
constexpr int nCues = 2;
struct DiffStates { enum { one }; };
constexpr int nDiffStates = 1;
struct FieldVars  { enum { seeker = nDiffStates, tracer, marker }; };
constexpr int nFieldVars = 3;
constexpr int diffVarForFieldVar(int fieldVar) {
   constexpr std::array<int, nFieldVars> lookup = { DiffStates::one, DiffStates::one, DiffStates::one };
   return lookup.at(fieldVar - nDiffStates);
}
struct EnvVars    { enum { }; };
constexpr int nEnvVars = 0;

#elif PROGRAM == ADVECTION

struct Cues { enum { }; };
constexpr int nCues = 0;
struct DiffStates { enum { one }; };
constexpr int nDiffStates = 1;
struct FieldVars  { enum { cena = nDiffStates }; };
constexpr int nFieldVars = 1;
constexpr int diffVarForFieldVar(int fieldVar) {
   constexpr std::array<int, nFieldVars> lookup = { DiffStates::one };
   return lookup.at(fieldVar - nDiffStates);
}
struct EnvVars    { enum { }; };
constexpr int nEnvVars = 0;

#elif PROGRAM == MIXTURE

struct Cues { enum { left, right }; };
constexpr int nCues = 2;
struct DiffStates { enum { one, a, b }; };
constexpr int nDiffStates = 3;
struct FieldVars  { enum { cena = nDiffStates, cenb }; };
constexpr int nFieldVars = 2;
constexpr int diffVarForFieldVar(int fieldVar) {
   constexpr std::array<int, nFieldVars> lookup = { DiffStates::a, DiffStates::b };
   return lookup.at(fieldVar - nDiffStates);
}
struct EnvVars    { enum { }; };
constexpr int nEnvVars = 0;

#elif PROGRAM == PARABOLOID

struct Cues { enum { circle }; };
constexpr int nCues = 1;
struct DiffStates { enum { one }; };
constexpr int nDiffStates = 1;
struct FieldVars  { enum { cena = nDiffStates }; };
constexpr int nFieldVars = 1;
constexpr int diffVarForFieldVar(int fieldVar) {
   constexpr std::array<int, nFieldVars> lookup = { DiffStates::one };
   return lookup.at(fieldVar - nDiffStates);
}
struct EnvVars    { enum { }; };
constexpr int nEnvVars = 0;

#elif PROGRAM == EDGE

struct Cues { enum { }; };
constexpr int nCues = 0;
struct DiffStates { enum { one }; };
constexpr int nDiffStates = 1;
struct FieldVars { enum { color = nDiffStates }; };
constexpr int nFieldVars = 1;
constexpr int diffVarForFieldVar(int fieldVar) {
   constexpr std::array<int, nFieldVars> lookup = { };
   return lookup.at(fieldVar - nDiffStates);
}
struct EnvVars    { enum { }; };
constexpr int nEnvVars = 0;

#elif PROGRAM == CALIBRATE

struct Cues { enum { }; };
constexpr int nCues = 0;
struct DiffStates { enum { one }; };
constexpr int nDiffStates = 36;
struct FieldVars { enum { }; };
constexpr int nFieldVars = 0;
constexpr int diffVarForFieldVar(int fieldVar) {
   constexpr std::array<int, nFieldVars> lookup = { };
   return lookup.at(fieldVar - nDiffStates);
}
struct EnvVars    { enum { }; };
constexpr int nEnvVars = 0;

#endif

constexpr int nSPHSubs = nDiffStates + nFieldVars;
constexpr int nSubs = nDiffStates + nFieldVars + nEnvVars;
typedef std::array<Substance, nSubs> Substances;
typedef std::array<double, nSubs> SubstanceVars; // TODO rename to SubstanceVals
typedef std::bitset<nCues> CueType; // TODO put this in fields.h?
typedef std::array<double, 3> Color; // cyan, magenta, yellow
typedef std::array<Color, nSubs + nCues> Colors;

#if PROGRAM == GROWTH
constexpr Colors defaultColors{ { {1.0, 0.0, 0.0}, {0.0, 1.0, 1.0}, {0.0, 0.0, 0.0} } };
#elif PROGRAM == GROWX
constexpr Colors defaultColors{ { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 1.0, 1.0}, {1.0, 0.0, 1.0} } };
#elif PROGRAM == DIFFUSION || PROGRAM == DIFFUSION_TRANSLATION
constexpr Colors defaultColors{ { {1.0, 1.0, 0.0} } };
#elif PROGRAM == PATH
constexpr Colors defaultColors{ { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0} } };
#elif PROGRAM == ADVECTION
constexpr Colors defaultColors{ { {1.0, 0.0, 0.0} , {0.0, 1.0, 1.0} } };
#elif PROGRAM == MIXTURE
constexpr Colors defaultColors{ { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 1.0, 1.0}, {1.0, 0.0, 1.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} } };
#elif PROGRAM == PARABOLOID
constexpr Colors defaultColors{ { {0.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.1, 0.1, 0.1} } };
#elif PROGRAM == EDGE
constexpr Colors defaultColors{ { {0.0, 0.0, 0.0}, {1.0, 1.0, 0.0} } };
#else
constexpr Colors defaultColors{ { {1.0, 1.0, 0.0} } };
#endif

#endif
