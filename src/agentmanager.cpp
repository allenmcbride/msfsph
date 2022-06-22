// agentmanager.cpp
// Allen McBride
// April 13, 2022
//
// Descriptions in agentmanager.h

#include "agentmanager.h"
#include "agent.h"
#include "grid.h"
#include "fields.h"
#include "managerinterface.h"
#include "point.h"
#include "substances.h"
#include "world.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iomanip> // DEBUG
#include <iostream>
#include <random>
#include <string>


int AgentManager::nextid = 0;


void AgentManager::updateInstrumentLoc() {
   double angleParam = 0.0;
   for (int isensor = 0; isensor < 4; ++isensor, angleParam += M_PI / 2.0) {
      sensors.at(isensor).location = surfacePoint(angleParam, amParams.stalkLen);
      prodLocations.at(isensor) = surfacePoint(angleParam + M_PI / 4.0, amParams.stalkLen);
   }
}


// "extraRadius" represents perpendicular distance from the surface
// only for the cardinal (relative to orientation) directions, or when
// eccentricity is zero.
Point AgentManager::surfacePoint(double angleParam, double extraRadius) const {
   const double effectiveSemimajor = amParams.semimajorAxis + extraRadius;
   const double effectiveSemiminor = amParams.semiminorAxis + extraRadius;
   return {
      extState.pos.x + effectiveSemimajor * std::cos(angleParam) * std::cos(extState.ori)
         - effectiveSemiminor * std::sin(angleParam) * std::sin(extState.ori),
         extState.pos.y + effectiveSemimajor * std::cos(angleParam) * std::sin(extState.ori)
            + effectiveSemiminor * std::sin(angleParam) * std::cos(extState.ori),
   }; 
}


// Calculate matrix of proportions occupied (but don't change fields yet),
// and also calculate information about collisions.  Stay agnostic about
// whether grid directions are same or opposite of physical directions.
// Must only be called while agent's own proportions are NOT asserted
// on grid.
int AgentManager::occlusion(double leeway, bool debug) {
   const auto oldPrec = debug ? std::cerr.precision(std::numeric_limits<double>::max_digits10) : std::cerr.precision();
   const auto a2 = std::pow(amParams.semimajorAxis, 2.0);
   const auto b2 = std::pow(amParams.semiminorAxis, 2.0);
   const auto e = (a2 - b2) * std::cos(2.0 * extState.ori);
   const auto ey = std::sqrt((a2 + b2 - e) / 2.0);
   const auto ex = std::sqrt((a2 + b2 + e) / 2.0);

   // These four numbers are not wrapped, so rowmin, etc. will also not
   // be wrapped.
   const auto yHigh = ey + extState.pos.y;
   const auto yLow = -ey + extState.pos.y;
   const auto xHigh = ex + extState.pos.x;
   const auto xLow = -ex + extState.pos.x;

   const auto rowmin = fields.fp.yToRow(yLow);
   const auto rowmax = fields.fp.yToRow(yHigh);
   if (rowmin > rowmax) throw std::logic_error("rowmin > rowmax");
   cornerRow = rowmin;

   const auto colmin = fields.fp.xToCol(xLow);
   const auto colmax = fields.fp.xToCol(xHigh);
   if (colmin > colmax) throw std::logic_error("colmin > colmax");
   cornerCol = colmin;

   if (debug && id == 149) std::cerr << "      xLow: " << xLow << " ... " << " cornerCol: " << cornerCol << std::endl;

   proportions.resize(rowmax - rowmin + 1);
   for (auto & row : proportions) 
      row.resize(colmax - colmin + 1);
   if (!passthrough) {
      edgesOccludedVertical.resize(rowmax - rowmin + 1);
      edgesOccludedHorizontal.resize(rowmax - rowmin);
      for (auto & row : edgesOccludedVertical) 
         row.resize(colmax - colmin);
      for (auto & row : edgesOccludedHorizontal) 
         row.resize(colmax - colmin + 1);
   }

   //Reset collision calcuation information
   int collCount = 0;

   // The arguments to findBoundIntersections must not be wrapped. The
   // maps it returns will be with respect to unwrapped coordinates.
   const auto [rowIntersections, colIntersections] = findBoundIntersections(rowmin, rowmax, colmin, colmax);
   for (int col = colmin; col <= colmax; ++col) {
      for (int row = rowmin; row <= rowmax; ++row) {
         const auto occupiedBySelf = proportionOccluded(col, row, colIntersections, rowIntersections, debug);
         proportions.at(row - cornerRow).at(col - cornerCol) = occupiedBySelf;

         if (!passthrough) {
            // row and col should be not wrapped, and Fparams should
            // return points with unwrapped coordinates.
            if (col < colmax) {
               const auto lowBound = fields.fp.boundToPt(row, col + 1).y;
               const auto highBound = fields.fp.boundToPt(row + 1, col + 1).y;
               const auto intersections = colIntersections.at(fields.fp.boundToPt(row, col + 1).x);
               edgesOccludedVertical.at(row - cornerRow).at(col - cornerCol) = findEdgeOccluded(lowBound, highBound, intersections);
            }
            if (row < rowmax) {
               const auto lowBound = fields.fp.boundToPt(row + 1, col).x;
               const auto highBound = fields.fp.boundToPt(row + 1, col + 1).x;
               const auto intersections = rowIntersections.at(fields.fp.boundToPt(row + 1, col).y);
               edgesOccludedHorizontal.at(row - cornerRow).at(col - cornerCol) = findEdgeOccluded(lowBound, highBound, intersections);
            }
         }

         // Check if this is a collision point to add to the count
         if (occupiedBySelf + fields.propOccupied(row, col) > 1.0 + leeway + std::numeric_limits<double>::epsilon()) {
            ++collCount;
            if (debug) {
               std::cerr << id << ", " << collCount << ": Collision: body-body" << std::endl;
               std::cerr << "occupiedBySelf: " << occupiedBySelf << std::endl;
               std::cerr << "propOcc: " << fields.propOccupied(row, col) << std::endl;
               std::cerr << "row, col: " << row << ", " << col << std::endl;
               std::cerr << "cornerRow, cornerCol: " << cornerRow << ", " << cornerCol << std::endl;
            }
         }
         if (!passthrough && fields.instruments->at(row, col) < 0) throw std::runtime_error("fields.instruments is negative"); // DEBUG
         if (!passthrough && fields.fullyObstructed(row, col, occupiedBySelf) && fields.instruments->at(row, col) > 0) {
            ++collCount;
            if (debug) std::cerr << id << ", " << collCount << ": Collision: this agent's body with another agent's instrument" << std::endl;
         }
      }
   }

   updateInstrumentLoc();

   // Make sure this agent's new instrument locations aren't fully
   // obstructed by another agent's body
   if (!passthrough) {
      for (const auto& sensor : sensors)
         if (fields.fullyObstructed(sensor.location)) {
            ++collCount;
            if (debug) std::cerr << id << ", " << collCount << ": Collision: this agent's sensor with another agent's body" << std::endl;
         }
      for (const auto& prodLocation : prodLocations)
         if (fields.fullyObstructed(prodLocation)) {
            ++collCount;
            if (debug) std::cerr << id << ", " << collCount << ": Collision: this agent's producer with another agent's body" << std::endl;
         }
   }

   return collCount;
}


// When boundaries are periodic, this works even when extState has
// not been wrapped, as long as the row/column parameters are also
// not wrapped.
std::pair<IntersectionsByLine, IntersectionsByLine> AgentManager::findBoundIntersections(int rowmin, int rowmax, int colmin, int colmax) const {
   IntersectionsByLine rowIntersections;
   IntersectionsByLine colIntersections;
   const double a2 = std::pow(amParams.semimajorAxis, 2.0);
   const double b2 = std::pow(amParams.semiminorAxis, 2.0);
   const double rt2ab = amParams.semimajorAxis * amParams.semiminorAxis * std::sqrt(2.0);
   const double costh = std::cos(extState.ori);
   const double sinth = std::sin(extState.ori);
   const double sincos = costh * sinth;
   const double costh2 = std::pow(costh, 2.0);
   const double sinth2 = std::pow(sinth, 2.0);
   const double cos2th = std::cos(2.0 * extState.ori);
   const double diffa2b2cos2 = (a2 - b2) * cos2th;
   const double denomX = costh2 / a2 + sinth2 / b2;
   const double denomY = costh2 / b2 + sinth2 / a2;
   for (int row = rowmin; row <= rowmax + 1; ++row) {
      const double y = fields.fp.boundToPt(row, 0).y;
      const double relY = y - extState.pos.y;
      const double e1raw = a2 + b2 - 2.0 * std::pow(relY, 2.0) - diffa2b2cos2;
      // I have forgotten why I inserted this shim instead of just using
      // e1raw. It may no longer be necessary.
      const double e1 = e1raw < 0.0 && e1raw > -std::numeric_limits<double>::epsilon() ? 0.0 : e1raw;
      if (e1 < 0.0) {
         rowIntersections[y] = {false, 0.0, 0.0};
      } else {
         const double e2 = std::sqrt(e1);
         const double e3 = -relY * sincos / a2 + relY * sincos / b2;
         const double lowX = (e3 - e2 / rt2ab) / denomX + extState.pos.x;
         const double highX = (e3 + e2 / rt2ab) / denomX + extState.pos.x;
         rowIntersections[y] = {true, lowX, highX};
      }
   }
   for (int col = colmin; col <= colmax + 1; ++col) {
      const double x = fields.fp.boundToPt(0, col).x;
      const double relX = x - extState.pos.x;
      const double e1raw = a2 + b2 - 2.0 * std::pow(relX, 2.0) + diffa2b2cos2;
      const double e1 = e1raw < 0.0 && e1raw > -std::numeric_limits<double>::epsilon() ? 0.0 : e1raw;
      if (e1 < 0.0) {
         colIntersections[x] = {false, 0.0, 0.0};
      } else {
         const double e2 = std::sqrt(e1);
         const double e3 = -relX * sincos / a2 + relX * sincos / b2;
         const double lowY = (e3 - e2 / rt2ab) / denomY + extState.pos.y;
         const double highY = (e3 + e2 / rt2ab) / denomY + extState.pos.y;
         colIntersections[x] = {true, lowY, highY};
      }
   }
   return {rowIntersections, colIntersections};
}


CornersInside AgentManager::findCornersInside(double cornerLow, double cornerHigh, const Intersections& intersections) const {
   if (!intersections.real || intersections.high <= cornerLow || intersections.low >= cornerHigh) 
      return CornersInside::neither;
   else if (intersections.high > cornerHigh && intersections.low < cornerLow) 
      return CornersInside::both;
   else if (intersections.high > cornerHigh)
      return CornersInside::high;
   else if (intersections.low < cornerLow)
      return CornersInside::low;
   else
      return CornersInside::span;
}


double AgentManager::findEdgeOccluded(double cornerLow, double cornerHigh, const Intersections& intersections) const {
   if (!intersections.real || intersections.high <= cornerLow || intersections.low >= cornerHigh) 
      return 0.0;
   else if (intersections.high >= cornerHigh && intersections.low <= cornerLow) 
      return 1.0;
   else if (intersections.high >= cornerHigh)
      return cornerHigh - intersections.low;
   else if (intersections.low <= cornerLow)
      return intersections.high - cornerLow;
   else
      return intersections.high - intersections.low;
}


// Assume agents large enough, and eccentricity low enough, that the
// following two cases are not possible:
//    * Agent crosses three or four cell faces, but occludes no corners
//    * Agent crosses two non-opposing faces, but occludes no corners
// This leaves 34 cases:
//    Seven cases that are distinct for all four directions, for 28 total cases:
//       * One corner out
//       * One corner in
//       * One side in
//       * Two corners out
//       * One side out, corner counterclockwise also out
//       * One side out, corner clockwise also out
//       * One side crossed, all corners out
//    Two cases that are distinct in two directions, for four total cases:
//       * Ellipse cuts through middle horizontally or vertically
//       * Ellipse cuts through middle diagonally
//    Two cases that are the same in all directions, for two total cases:
//       * Ellipse fully occludes cell
//       * Ellipse fully beyond cell
// The five cases consisting of one side crossed but no corners occluded,
// plus the case of the ellipse fully beyond the cell, are combined in three
// branches of the code below, for 32 total points of return.
double AgentManager::proportionOccluded(int col, int row, const IntersectionsByLine& colIntersections, const IntersectionsByLine& rowIntersections, bool debug) {
   const Point lowCorner = fields.fp.boundToPt(row, col);
   const Point highCorner = fields.fp.boundToPt(row + 1, col + 1);
   const double westX = lowCorner.x;
   const double eastX = highCorner.x;
   const double southY = lowCorner.y;
   const double northY = highCorner.y;
   const auto westIntersections = colIntersections.at(westX);
   const auto eastIntersections = colIntersections.at(eastX);
   const auto southIntersections = rowIntersections.at(southY);
   const auto northIntersections = rowIntersections.at(northY);
   const auto cornersInsideWest = findCornersInside(southY, northY, westIntersections);
   const auto cornersInsideEast = findCornersInside(southY, northY, eastIntersections);
   const auto cornersInsideSouth = findCornersInside(westX, eastX, southIntersections);
   const auto cornersInsideNorth = findCornersInside(westX, eastX, northIntersections);
   double proportion = -1.0;
   switch(cornersInsideWest) {
      case CornersInside::both:
         switch(cornersInsideEast) {
            case CornersInside::both:
               // Ellipse fully covers grid cell.
               proportion = 1.0;
               break;
            case CornersInside::neither:
               // Ellipse covers west side of cell; approximate as trapezoid.
               proportion = (southIntersections.high + northIntersections.high - 2.0 * westX) / (2.0 * fields.fp.hx);
               break;
            case CornersInside::high:
               // Ellipse covers all but southeast corner; approximate uncovered part as triangle.
               proportion = 1.0 - (eastX - southIntersections.high) * (eastIntersections.low - southY) / (2.0 * fields.fp.cellArea());
               break;
            case CornersInside::low:
               // Ellipse covers all but northeast corner; approximate uncovered part as triangle.
               proportion = 1.0 - (eastX - northIntersections.high) * (northY - eastIntersections.high) / (2.0 * fields.fp.cellArea());
               break;
            case CornersInside::span:
               // Ellipse covers all but northeast and southeast corners; approximate uncovered part as two triangles.
               proportion = 1.0 - ((eastX - northIntersections.high) * (northY - eastIntersections.high)
                     + (eastX - southIntersections.high) * (eastIntersections.low - southY)) 
                  / (2.0 * fields.fp.cellArea());
               break;
         }
         break;
      case CornersInside::neither:
         switch(cornersInsideEast) {
            case CornersInside::both:
               // Ellipse covers east side of cell; approximate as trapezoid.
               proportion = (2.0 * eastX - southIntersections.low - northIntersections.low) / (2.0 * fields.fp.hx);
               break;
            case CornersInside::neither:
               if (cornersInsideSouth == CornersInside::span && cornersInsideNorth == CornersInside::span) {
                  // Ellipse cuts vertically through middle of cell; approximate as trapezoid.
                  proportion = (southIntersections.high - southIntersections.low 
                        + northIntersections.high - northIntersections.low)
                     / (2.0 * fields.fp.hx);
               } else {
                  // Ellipse may cross north or south face, or may be completely outside grid cell;
                  // either way approximate as zero.
                  proportion = 0.0;
               }
               break;
            case CornersInside::high:
               if (cornersInsideSouth == CornersInside::span) {
                  // Ellipse covers northeast corner and covers a central part of south face;
                  // approximate uncovered part as one trapezoid to west and one triangle in southeast.
                  proportion = 1.0 - (eastX - southIntersections.high) * (eastIntersections.low - southY) / (2.0 * fields.fp.cellArea())
                     - (southIntersections.low + northIntersections.low - 2.0 * westX) / (2.0 * fields.fp.hx);
               } else {
                  // Ellipse covers northeast corner only; approximate as triangle.
                  proportion = (eastX - northIntersections.low) * (northY - eastIntersections.low) / (2.0 * fields.fp.cellArea());
               }
               break;
            case CornersInside::low:
               if (cornersInsideNorth == CornersInside::span) {
                  // Ellipse covers southeast corner and covers a central part of north face;
                  // approximate uncovered part as one trapezoid to west and one triangle in northeast.
                  proportion = 1.0 - (eastX - northIntersections.high) * (northY - eastIntersections.high) / (2.0 * fields.fp.cellArea())
                     - (southIntersections.low + northIntersections.low - 2.0 * westX) / (2.0 * fields.fp.hx);
               } else {
                  // Ellipse covers southeast corner only; approximate as triangle.
                  proportion = (eastX - southIntersections.low) * (eastIntersections.high - southY) / (2.0 * fields.fp.cellArea());
               }
               break;
            case CornersInside::span:
               // Ellipse crosses east face, but no corners. Approximate as zero.
               proportion = 0.0;
               break;
         }
         break;
      case CornersInside::high:
         switch(cornersInsideEast) {
            case CornersInside::both:
               // Ellipse covers all but southwest corner; approximate uncovered part as triangle.
               proportion = 1.0 - (southIntersections.low - westX) * (westIntersections.low - southY) / (2.0 * fields.fp.cellArea());
               break;
            case CornersInside::neither:
               if (cornersInsideSouth == CornersInside::span) {
                  // Ellipse covers northwest corner and covers a central part of south face;
                  // approximate uncovered part as one trapezoid to east and one triangle in southwest.
                  proportion = 1.0 - (southIntersections.low - westX) * (westIntersections.low - southY) / (2.0 * fields.fp.cellArea())
                     - (2.0 * eastX - southIntersections.high - northIntersections.high) / (2.0 * fields.fp.hx);
               } else {
                  // Ellipse covers northwest corner only; approximate as triangle.
                  proportion = (northIntersections.high - westX) * (northY - westIntersections.low) / (2.0 * fields.fp.cellArea());
               }
               break;
            case CornersInside::high:
               if (cornersInsideSouth == CornersInside::span) {
                  // Ellipse covers all but southeast and southwest corners; approximate uncovered part as two triangles.
                  proportion = 1.0 - ((eastX - southIntersections.high) * (eastIntersections.low - southY)
                        + (southIntersections.low - westX) * (westIntersections.low - southY)) 
                     / (2.0 * fields.fp.cellArea());
               } else {
                  // Ellipse covers north side of cell; approximate as trapezoid.
                  proportion = (2.0 * northY - eastIntersections.low - westIntersections.low) / (2.0 * fields.fp.hy);
               }
               break;
            case CornersInside::low:
               // Ellipse covers northwest and southeast corners, cutting diagonally through cell;
               // approximate uncovered parts as two triangles.
               proportion = 1.0 - ((eastX - northIntersections.high) * (northY - eastIntersections.high)
                     + (southIntersections.low - westX) * (westIntersections.low - southY)) 
                  / (2.0 * fields.fp.cellArea());
               break;
            case CornersInside::span:
               // Ellipse covers northwest corner and covers a central part of east face;
               // approximate uncovered part as one trapezoid to south and one triangle in northeast.
               proportion = 1.0 - (eastX - northIntersections.high) * (northY - eastIntersections.high) / (2.0 * fields.fp.cellArea())
                  - (eastIntersections.low + westIntersections.low - 2.0 * southY) / (2.0 * fields.fp.hy);
               break;
         }
         break;
      case CornersInside::low:
         switch(cornersInsideEast) {
            case CornersInside::both:
               // Ellipse covers all but northwest corner; approximate uncovered part as triangle.
               proportion = 1.0 - (northIntersections.low - westX) * (northY - westIntersections.high) / (2.0 * fields.fp.cellArea());
               break;
            case CornersInside::neither:
               if (cornersInsideNorth == CornersInside::span) {
                  // Ellipse covers southwest corner and covers a central part of north face;
                  // approximate uncovered part as one trapezoid to east and one triangle in northwest.
                  proportion = 1.0 - (northIntersections.low - westX) * (northY - westIntersections.high) / (2.0 * fields.fp.cellArea())
                     - (2.0 * eastX - southIntersections.high - northIntersections.high) / (2.0 * fields.fp.hx);
               } else {
                  // Ellipse covers southwest corner only; approximate as triangle.
                  proportion = (southIntersections.high - westX) * (westIntersections.high - southY) / (2.0 * fields.fp.cellArea());
               }
               break;
            case CornersInside::high:
               // Ellipse covers southwest and northeast corners, cutting diagonally through cell;
               // approximate uncovered parts as two triangles.
               proportion = 1.0 - ((eastX - southIntersections.high) * (eastIntersections.low - southY)
                     + (northIntersections.low - westX) * (northY - westIntersections.high)) 
                  / (2.0 * fields.fp.cellArea());
               break;
            case CornersInside::low:
               if (cornersInsideNorth == CornersInside::span) {
                  // Ellipse covers all but northeast and northwest corners; approximate uncovered part as two triangles.
                  proportion = 1.0 - ((eastX - northIntersections.high) * (northY - eastIntersections.high)
                        + (northIntersections.low - westX) * (northY - westIntersections.high)) 
                     / (2.0 * fields.fp.cellArea());
               } else {
                  // Ellipse covers south side of cell; approximate as trapezoid.
                  proportion = (eastIntersections.high + westIntersections.high - 2.0 * southY) / (2.0 * fields.fp.hy);
               }
               break;
            case CornersInside::span:
               // Ellipse covers southwest corner and covers a central part of east face;
               // approximate uncovered part as one trapezoid to north and one triangle in southeast.
               proportion = 1.0 - (eastX - southIntersections.high) * (eastIntersections.low - southY) / (2.0 * fields.fp.cellArea())
                  - (2.0 * northY - eastIntersections.high - westIntersections.high) / (2.0 * fields.fp.hy);
               break;
         }
         break;
      case CornersInside::span:
         switch(cornersInsideEast) {
            case CornersInside::both:
               // Ellipse covers all but northwest and southwest corners; approximate uncovered part as two triangles.
               proportion = 1.0 - ((northIntersections.low - westX) * (northY - westIntersections.high)
                     + (southIntersections.low - westX) * (westIntersections.low - southY)) 
                  / (2.0 * fields.fp.cellArea());
               break;
            case CornersInside::neither:
               // Ellipse crosses west face, but no corners. Approximate as zero.
               proportion = 0.0;
               break;
            case CornersInside::high:
               // Ellipse covers northeast corner and covers a central part of west face;
               // approximate uncovered part as one trapezoid to south and one triangle in northwest.
               proportion = 1.0 - (northIntersections.low - westX) * (northY - westIntersections.high) / (2.0 * fields.fp.cellArea())
                  - (eastIntersections.low + westIntersections.low - 2.0 * southY) / (2.0 * fields.fp.hy);
               break;
            case CornersInside::low:
               // Ellipse covers southeast corner and covers a central part of west face;
               // approximate uncovered part as one trapezoid to north and one triangle in southwest.
               proportion = 1.0 - (southIntersections.low - westX) * (westIntersections.low - southY) / (2.0 * fields.fp.cellArea())
                  - (2.0 * northY - eastIntersections.high - westIntersections.high) / (2.0 * fields.fp.hy);
               break;
            case CornersInside::span:
               // Ellipse cuts horizontally through middle of cell; approximate as trapezoid.
               proportion = (westIntersections.high - westIntersections.low 
                     + eastIntersections.high - eastIntersections.low)
                  / (2.0 * fields.fp.hy);
               break;
         }
         break;
   }
   if (proportion < 0.0 || proportion > 1.0) {
      std::cerr << static_cast<std::underlying_type<CornersInside>::type>(cornersInsideWest) << std::endl; // DEBUG
      std::cerr << static_cast<std::underlying_type<CornersInside>::type>(cornersInsideEast) << std::endl; // DEBUG
      std::cerr << static_cast<std::underlying_type<CornersInside>::type>(cornersInsideSouth) << std::endl; // DEBUG
      std::cerr << static_cast<std::underlying_type<CornersInside>::type>(cornersInsideNorth) << std::endl; // DEBUG
      std::cerr << westIntersections << std::endl; // DEBUG
      std::cerr << eastIntersections << std::endl; // DEBUG
      std::cerr << southIntersections << std::endl; // DEBUG
      std::cerr << northIntersections << std::endl; // DEBUG
      std::cerr << "lowCorner: " << lowCorner << std::endl; // DEBUG
      std::cerr << "highCorner: " << highCorner << std::endl; // DEBUG
      std::cerr << "proportion: " << proportion << std::endl; // DEBUG
      throw std::logic_error("Either no case matched in proportionOccluded, or proportion calculated outside [0.0, 1.0].");
   }
   return proportion;   
}


// Modify the agentBodies field to reflect occupancy fractions of this agent.
// Modify instruments field to reflect instrument counts.
void AgentManager::assertOccupation(bool debug) {
   if (debug && id == 149) {
      std::cerr << "149 corners: " << cornerRow << "; " << cornerCol << std::endl;
      std::cerr << "149 loc: " << extState.pos << std::endl;
      for (int row = 0; row < static_cast<int>(proportions.size()); ++row) {
         for (int col = 0; col < static_cast<int>(proportions.at(0).size()); ++col) {
            std::cerr << "   " << row + cornerRow << "," << col + cornerCol << ":" << proportions.at(row).at(col) << "     ";
         }
         std::cerr << std::endl;
      }
   }

   for (int row = 0; row < static_cast<int>(proportions.size()); ++row)
      for (int col = 0; col < static_cast<int>(proportions.at(0).size()); ++col)
         fields.agentBodies.at(row + cornerRow, col + cornerCol) += proportions.at(row).at(col);

   if (!passthrough) {
      for (const auto& sensor : sensors)
         ++fields.instruments->at(sensor.location);
      for (const auto& prodLocation : prodLocations)
         ++fields.instruments->at(prodLocation);
   }
}


// Opposite of assertOccupation.
void AgentManager::eraseOccupation(bool debug) {
   for (int row = 0; row < static_cast<int>(proportions.size()); ++row) {
      for (int col = 0; col < static_cast<int>(proportions.at(0).size()); ++col) {
         if (debug && row + cornerRow == 373 && col + cornerCol == 373) std::cerr << std::endl;
         if (debug && row + cornerRow == 373 && col + cornerCol == 373) std::cerr << "id: " << id << ";  before erase: " << fields.agentBodies.at(row + cornerRow, col + cornerCol) << " (" << cornerRow << ", " << cornerCol << "; " << proportions.size() << ", " << proportions.at(0).size() << ")" << std::endl;
         fields.agentBodies.at(row + cornerRow, col + cornerCol) -= proportions.at(row).at(col);
         if (debug && row + cornerRow == 373 && col + cornerCol == 373) std::cerr << "id: " << id << ";  after erase: " << fields.agentBodies.at(row + cornerRow, col + cornerCol) << " (" << cornerRow << ", " << cornerCol << "; " << proportions.size() << ", " << proportions.at(0).size() << ")" << std::endl;
         if (debug && row + cornerRow == 373 && col + cornerCol == 373) std::cerr << "id: " << id << ";  proportion: " << proportions.at(row).at(col) << std::endl;
         if (debug && row + cornerRow == 373 && col + cornerCol == 373) std::cerr << std::endl;
      }
   }

   if (!passthrough) {
      for (const auto& sensor : sensors)
         --fields.instruments->at(sensor.location);
      for (const auto& prodLocation : prodLocations)
         --fields.instruments->at(prodLocation);
   }
}


// This should be called after all Agent objects have Agent::update() or Agent::calibrate() called,
// and before the Diffuser update is called.
void AgentManager::assertTransientInfo() {
   for (int row = 0; row < static_cast<int>(proportions.size()); ++row) {
      for (int col = 0; col < static_cast<int>(proportions.at(0).size()); ++col) {
         const auto gridRow = row + cornerRow;
         const auto gridCol = col + cornerCol;
         const auto vectorFromCenter = fields.fp.coordsToPt(gridRow, gridCol) - extState.pos;
         const auto angleFromOri = Point::angleDiff(vectorFromCenter.angle(), extState.ori);
         const auto isFront = std::max(0.0, 1.0 - std::abs(angleFromOri) / (M_PI / 2.0));
         if (fields.front) fields.front->at(gridRow, gridCol) += isFront * proportions.at(row).at(col);

         if (!passthrough) {
            if (row < static_cast<int>(proportions.size()) - 1) {
               const auto bodyVelEdgeY = bodyVelAtPt(fields.fp.horizEdgeToPt(gridRow, gridCol)).y;
               const auto edgeOccluded = edgesOccludedHorizontal.at(row).at(col);
               fields.edgeFreepartHorizontal->at(gridRow, gridCol) -= edgeOccluded;
               fields.edgeBodyVelPropHorizontal->at(gridRow, gridCol) += edgeOccluded * bodyVelEdgeY;
            }

            if (col < static_cast<int>(proportions.at(0).size()) - 1) {
               const auto bodyVelEdgeX = bodyVelAtPt(fields.fp.vertEdgeToPt(gridRow, gridCol)).x;
               const auto edgeOccluded = edgesOccludedVertical.at(row).at(col);
               fields.edgeFreepartVertical->at(gridRow, gridCol) -= edgeOccluded;
               fields.edgeBodyVelPropVertical->at(gridRow, gridCol) += edgeOccluded * bodyVelEdgeX;
            }
         }
      }
   }
}


Point AgentManager::boundingBoxRelative() const {
   const double a2 = std::pow(amParams.semimajorAxis + amParams.stalkLen, 2.0);
   const double b2 = std::pow(amParams.semiminorAxis + amParams.stalkLen, 2.0);
   const double e = (a2 - b2) * std::cos(2.0 * extState.ori);
   return {std::sqrt((a2 + b2 + e) / 2.0), std::sqrt((a2 + b2 - e) / 2.0)};
}


void AgentManager::moveInBounds(const Point& targetLocation) {
   const Point cornerRelative{boundingBoxRelative()};
   const Point targetCornerLow{targetLocation - cornerRelative};
   const Point targetCornerHigh{targetLocation + cornerRelative};
   const Point oldPos{extState.pos}; // DEBUG

   // There is some addition and subtraction involved in determining
   // in-bounds vs.  out-of-bounds, so some shims are needed to make up
   // for the triangle inequality not always holding.
   const double xShim = std::numeric_limits<double>::epsilon() * fields.fp.xphys;
   const double yShim = std::numeric_limits<double>::epsilon() * fields.fp.yphys;

   bool xAdjustLow = false;
   bool xAdjustHigh = false;
   if (!fields.inBoundsXLow(targetCornerLow.x)) {
      extState.pos.x = cornerRelative.x + xShim;
      xAdjustLow = true;
   }
   if (!fields.inBoundsXHigh(targetCornerHigh.x)) {
      extState.pos.x = fields.fp.xphys - cornerRelative.x - xShim;
      xAdjustHigh = true;
   }
   if (!xAdjustLow && !xAdjustHigh)
      extState.pos.x = targetLocation.x;
   else if (xAdjustLow && xAdjustHigh)
      throw std::runtime_error("Agent appears to be out of bounds on both sides, in x dimension");

   bool yAdjustLow = false;
   bool yAdjustHigh = false;
   if (!fields.inBoundsYLow(targetCornerLow.y)) {
      extState.pos.y = cornerRelative.y + yShim;
      yAdjustLow = true;
   }
   if (!fields.inBoundsYHigh(targetCornerHigh.y)) {
      extState.pos.y = fields.fp.yphys - cornerRelative.y - yShim;
      yAdjustHigh = true;
   }
   if (!yAdjustLow && !yAdjustHigh)
      extState.pos.y = targetLocation.y;
   else if (yAdjustLow && yAdjustHigh)
      throw std::runtime_error("Agent appears to be out of bounds on both sides, in y dimension");

#if PROGRAM == PATH
   if (xAdjustLow || xAdjustHigh || yAdjustLow || yAdjustHigh)
      postCollisionTimer = maxSettleTime;
#endif

   if (!inBounds()) {
      const Point oldCornerLow{oldPos - cornerRelative};
      const Point oldCornerHigh{oldPos + cornerRelative};
      const Point newCornerLow{extState.pos - cornerRelative};
      const Point newCornerHigh{extState.pos + cornerRelative};
      auto oldPrec = std::cerr.precision(std::numeric_limits<double>::max_digits10);
      std::cerr << "xphys: " << fields.fp.xphys << std::endl;
      std::cerr << "yphys: " << fields.fp.yphys << std::endl;
      std::cerr << "cornerRelative: " << cornerRelative << std::endl;
      std::cerr << "target location: " << targetLocation << std::endl;
      std::cerr << "old location: " << oldPos << std::endl;
      std::cerr << "xShim: " << xShim << std::endl;
      std::cerr << "yShim: " << xShim << std::endl;
      std::cerr << "old low corner: " << oldCornerLow << std::endl;
      std::cerr << "old high corner: " << oldCornerHigh << std::endl;
      std::cerr << "new low corner: " << newCornerLow << std::endl;
      std::cerr << "new high corner: " << newCornerHigh << std::endl;
      std::cerr.precision(oldPrec);
      throw std::runtime_error("inBounds() should return true at end of moveInBounds()");
   }
}


bool AgentManager::resolveCollisions(const AMExternalState& origExtState, const Point& displacement, double oriChange, double collResDisp) {
   thread_local std::uniform_real_distribution<double> randOri(-M_PI, M_PI);
   std::uniform_real_distribution<double> randOriChange{std::min(0.0, oriChange), 
      std::nextafter(std::max(0.0, oriChange), std::numeric_limits<double>::max())};
   std::uniform_real_distribution<double> randUnit{0.0, 1.0};
   int maxTries = 10;
   int nTries = 0;
   while (occlusion() > 0 && nTries++ < maxTries) {
      const auto oriChangeTry = randOriChange(gen);
      extState.ori = origExtState.ori + oriChangeTry;
      const auto dispMagTry = std::sqrt(randUnit(gen)) * collResDisp;
      moveInBounds(origExtState.pos + Point::polarPt(dispMagTry, randOri(gen)));
   }
#if PROGRAM == PATH
   if (nTries > 0) 
      postCollisionTimer = maxSettleTime;
#endif
   return nTries <= maxTries; // False iff no non-collision point found
}



void AgentManager::sensorUpdate(double dt, bool trueValue) {
   sensorLifetime += dt;

   for (auto& sensor : sensors) {
      if (!fields.fullyObstructed(sensor.location)) {
         for (int isub = 0; isub < nSubs; ++isub) {
            const auto concentration = fields.concentration(isub, sensor.location);
            if (concentration < 0.0) {
               std::cerr << "Negative concentration found: " << concentration << std::endl;
               throw std::runtime_error("Concentrataion less than zero at sensor location: " + std::to_string(concentration));
            }
            if (amParams.sensorMeanBindingTime <= 0.0 || trueValue) {
               sensor.reading.at(isub) = (trueValue ? 1.0 : sensor.bias) * concentration;
               sensorsReady = true;
            } else {
               const auto lifetimeCorrectionFactor = -std::expm1(-sensorLifetime / amParams.sensorMeanBindingTime);
               const auto readingCoefficient = sensor.bias / (amParams.sensorMeanBindingTime * lifetimeCorrectionFactor);
               if (amParams.sensorEventRatePerConcentration > 0.0) {
                  const auto nExpectedBindings = dt * amParams.sensorEventRatePerConcentration * concentration;
                  const auto nBindings = std::poisson_distribution<NBoundType>(nExpectedBindings)(gen);
                  const auto unbindingProbabilityAlreadyBound = -std::expm1(-dt / amParams.sensorMeanBindingTime);
                  const auto unbindingProbabilityNewlyBound = -std::expm1(-dt / (2.0 * amParams.sensorMeanBindingTime));
                  const auto nUnbindingAlreadyBound = std::binomial_distribution<NBoundType>(sensor.nBound.at(isub), unbindingProbabilityAlreadyBound)(gen);
                  const auto nUnbindingNewlyBound = std::binomial_distribution<NBoundType>(nBindings, unbindingProbabilityNewlyBound)(gen);
                  sensor.nBound.at(isub) += nBindings - nUnbindingAlreadyBound - nUnbindingNewlyBound;

                  if (sensor.nBound.at(isub) < 0)
                     throw std::logic_error("It should be impossible for number of bound ligand molecules to be negative"); 

                  sensor.reading.at(isub) = readingCoefficient * sensor.nBound.at(isub) / amParams.sensorEventRatePerConcentration;
                  if (isub == DiffStates::one) 
                     readinessCount += nBindings;
                  constexpr NBoundType readinessThreshold = 100;
                  if (!sensorsReady && readinessCount >= readinessThreshold)
                     sensorsReady = true;
               } else {
                  sensor.amountBound.at(isub) += dt * concentration;
                  sensor.amountBound.at(isub) *= std::exp(-dt / amParams.sensorMeanBindingTime);
                  sensor.reading.at(isub) = readingCoefficient * sensor.amountBound.at(isub);
                  sensorsReady = true;
               }
            }
         }
      }
   }
}


// Sense concentrations at sensor locations.  Return a value of -1
// if a sensor is obstructed or failed, to signal a bad reading.
// I'm deliberately using lifetime here instead of sensorLifetime,
// so that instrument failures are consistent and not a function of
// burn-in time.  It would probably make sense to discard the idea
// of sensors being obstructed.
SensorVals AgentManager::sense(int isub) const { 
   SensorVals sensorReadings;
   for (size_t i = 0; i < sensors.size(); ++i) {
      const bool working = sensors.at(i).lifespan <= 0.0 || lifetime < sensors.at(i).lifespan;
      if (working && !fields.fullyObstructed(sensors.at(i).location))
         sensorReadings.at(i) = sensors.at(i).reading.at(isub);
      else 
         sensorReadings.at(i) = -1.0;
   }
   return sensorReadings;
}


// rate is mass of substance per unit time.  Ignore negative values to
// remain physical.  Produce nothing if secretion system has failed.
void AgentManager::produce(int isub, const double rate, const double dt) const {
   const auto productionBroken = productionLifespan > 0.0 && lifetime > productionLifespan;
   const auto waitingAfterCollision = isub >= nDiffStates && postCollisionTimer > 0.0;
   if (productionBroken || waitingAfterCollision) return;

   if (rate > 0.0) {
      for (const auto& loc : prodLocations) { 
         thread_local std::normal_distribution<double> noise(0.0, 1.0);
         fields.injectMass(isub, loc, dt * rate / 4.0 + 0.0 * noise(gen));
      }
   }
}


void AgentManager::move(double desiredSpeed, double desiredAngularVel, double collResSpeed, double dt) {
   if (propulsionLifespan > 0.0 && lifetime > propulsionLifespan) {
      // Don't move if propulsion has failed.
      vel = Point{};
      angularVel = 0.0;
   } else {

      // Calculate motion
      const double afterBiasAngularVel = desiredAngularVel + oriBias;
      const double afterBiasOriChange = dt * afterBiasAngularVel;
      const double afterBiasSpeed = desiredSpeed * speedBias;
      const Point afterBiasPropulsion = Point::polarPt(afterBiasSpeed, extState.ori + afterBiasOriChange);
      thread_local std::normal_distribution<double> motionNoiseRel(0.0, amParams.motionNoiseRelSD);
      const Point motionNoise = afterBiasSpeed * Point{motionNoiseRel(gen), motionNoiseRel(gen)};
      const Point actualPropulsion = afterBiasPropulsion + motionNoise;
      const Point bodyforceVel = 0.00 * fields.gravityDir;
      const Point displacement = dt * (actualPropulsion + bodyforceVel);
      if (maxSpeed > 0.0 && displacement.mag() - maxSpeed * dt > std::numeric_limits<double>::epsilon()) {
         std::cerr << "Effective speed requested by displacement would be " << displacement.mag() / dt << std::endl;
         throw MaxSpeedExceeded("Maximum speed exceeded by intended displacement");
      }

      //eraseOccupation(true); // DEBUG
      
      // Move agent
      const auto origExtState = extState;
      extState.ori += afterBiasOriChange;
      moveInBounds(extState.pos + displacement);

      // If agent updates are to be parallelized, the atomic section
      // needs to start here and go through the next addorSubProportions.
      // If we tried to hold off on the first addOrSubProportions call
      // to reduce the size of the atomic block, multiple agents could
      // see a potential target location as unoccupied.
      eraseOccupation(); // CORRECT LOCATION

      if (!resolveCollisions(origExtState, displacement, afterBiasOriChange, collResSpeed * dt)) {
         extState = origExtState;
         if (occlusion(std::numeric_limits<double>::epsilon()) > 0) {
            std::cerr.precision(std::numeric_limits<double>::max_digits10);std::cerr.precision(std::numeric_limits<double>::max_digits10);
            std::cerr << origExtState.pos << std::endl;
            std::cerr << origExtState.ori << std::endl;
            std::cerr << extState.pos << std::endl;
            std::cerr << extState.ori << std::endl;
            throw std::logic_error("Collision should never be found at original position (after moving).");
         }
      }

      // This should be the only place where vel is updated
      vel = (extState.pos - origExtState.pos) / dt; 
      angularVel = Point::angleDiff(extState.ori, origExtState.ori) / dt;
      assertOccupation();

#ifdef XWRAP
      // Unnecessary to simulation because of wrapping logic in grid
      // references, but simplifies final output by ensuring agent
      // positions are within domain.
      extState.pos.x = std::fmod(extState.pos.x, fields.fp.xphys);
      if (extState.pos.x < 0.0)
         extState.pos.x += fields.fp.xphys;
#endif
   }
}


void AgentManager::init() {
   if (!inBounds()) 
      throw OOBException("Position: " + std::to_string(extState.pos.x) + ", " + std::to_string(extState.pos.y));

   if (amParams.sensorEventRatePerConcentration > 0.0 && amParams.sensorMeanBindingTime <= 0.0)
      throw std::invalid_argument("If sensorEventRatePerConcentration is positive, sensorMeanBindingTime must also be positive.");

   id = nextid++;

   // NOISE
   for (auto& sensor : sensors) 
      sensor.bias = std::lognormal_distribution<double>(0.0, amParams.sensorBiasLogSD)(gen);
   speedBias = std::lognormal_distribution<double>(0.0, amParams.speedBiasLogSD)(gen);
   oriBias = std::normal_distribution<double>(0.0, amParams.oriBiasSD)(gen);

   for (auto& sensor : sensors) 
      sensor.lifespan = amParams.sensorMTTF <= 0.0 ? 0.0 : std::exponential_distribution(1.0 / amParams.sensorMTTF)(gen);
   propulsionLifespan = amParams.propulsionMTTF <= 0.0 ? 0.0 : std::exponential_distribution(1.0 / amParams.propulsionMTTF)(gen);
   productionLifespan = amParams.productionMTTF <= 0.0 ? 0.0 : std::exponential_distribution(1.0 / amParams.productionMTTF)(gen);

   if (occlusion() > 0)
      throw ObstructionException("Position: " + std::to_string(extState.pos.x) + ", " + std::to_string(extState.pos.y));
         
   assertOccupation();
}


AgentManager::AgentManager(
      const AgentManagerParams& amParams,
      const AMExternalState& extState,
      Fields& fields,
      const Substances& subs,
      const SimulationParams& simParams,
      std::default_random_engine& gen,
      NuclearOutput nucOut
      ) : 
   amParams{amParams},
   extState{extState},
   fields{fields},
   gen{gen},
   maxSpeed{simParams.maxSpeed},
   passthrough{simParams.passthrough}, 
   maxSettleTime{
      [&] {
         double maxSettleTime = 0.0;
         for (const Substance& s : subs)
            maxSettleTime = s.decayRate > 0.0 ? std::max(s.timeToSteady(0.01), maxSettleTime) : 0.0;
         return maxSettleTime;
      }()
   },
   agent{
      AgentParams{
         2.0 * (amParams.semimajorAxis + amParams.stalkLen),
         2.0 * (amParams.semiminorAxis + amParams.stalkLen),
         amParams.mass(),
         simParams.leakRate,
         simParams.maxSpeed,
         simParams.maxAngularSpeed,
         maxSettleTime,
         -std::log(propOfOriginalToAllow) * amParams.sensorMeanBindingTime, // Wait for 99% of bound ligand to be gone
         simParams.speedupFactor
      },
      subs, 
      ManagerInterface(*this), 
      nucOut
   }
{
   init();
}


AgentManager::AgentManager(
      const AgentManager& am,
      const Point& pos
      ) : 
   amParams{am.amParams}, 
   extState{pos, am.extState.ori},
   fields{am.fields},
   gen{am.gen},
   maxSpeed{am.maxSpeed},
   passthrough{am.passthrough}, 
   maxSettleTime{am.maxSettleTime},
   agent{am.agent, ManagerInterface(*this)}
{
   init();
}


// UNFINISHED
AgentManager::AgentManager(
      Fields& fields, 
      std::default_random_engine& gen, 
      std::istream& is
      ) :
   amParams{readObjectFromStream<AgentManagerParams>(is)},
   extState{readObjectFromStream<AMExternalState>(is)},
   fields{fields},
   gen{gen},
   maxSpeed{readObjectFromStream<double>(is)},
   passthrough{static_cast<bool>(is.get())},
   maxSettleTime{readObjectFromStream<double>(is)},
   agent{ManagerInterface(*this), is},
   id{readObjectFromStream<int>(is)},
   sensors{readObjectFromStream<std::array<Sensor, 4> >(is)},
   prodLocations{readObjectFromStream<std::array<Point, 4> >(is)},
   oriBias{readObjectFromStream<double>(is)},
   speedBias{readObjectFromStream<double>(is)},
   productionLifespan{readObjectFromStream<double>(is)},
   propulsionLifespan{readObjectFromStream<double>(is)},
   lifetime{readObjectFromStream<double>(is)},
   sensorLifetime{readObjectFromStream<double>(is)},
   memory{readObjectFromStream<Memory>(is)},
   sensorsReady{static_cast<bool>(is.get())},
   readinessCount{readObjectFromStream<NBoundType>(is)},
   postCollisionTimer{readObjectFromStream<double>(is)},
   vel{readObjectFromStream<Point>(is)},
   angularVel{readObjectFromStream<double>(is)},
   cornerRow{readObjectFromStream<int>(is)},
   cornerCol{readObjectFromStream<int>(is)}
{
   const auto propRows = readObjectFromStream<int>(is);
   const auto propCols = readObjectFromStream<int>(is);
   proportions.resize(propRows);
   for (auto& v : proportions) {
      v.resize(propCols);
      is.read(reinterpret_cast<char *>(v.data()), v.size() * sizeof(AgentBodyType));
   }
   if (!passthrough) {
      edgesOccludedVertical.resize(propRows);
      edgesOccludedHorizontal.resize(propRows - 1);
      for (auto& v : edgesOccludedVertical) {
         v.resize(propCols - 1);
         is.read(reinterpret_cast<char *>(v.data()), v.size() * sizeof(EdgeType));
      }
      for (auto& v : edgesOccludedHorizontal) {
         v.resize(propCols);
         is.read(reinterpret_cast<char *>(v.data()), v.size() * sizeof(EdgeType));
      }
   }
}


AgentManager::~AgentManager() {
   eraseOccupation();
}


// UNFINISHED
void AgentManager::serialize(std::ostream& os) const {
   os.write(reinterpret_cast<const char *>(&amParams), sizeof(amParams));
   os.write(reinterpret_cast<const char *>(&extState), sizeof(extState));
   os.write(reinterpret_cast<const char *>(&maxSpeed), sizeof(maxSpeed));
   os.put(static_cast<char>(passthrough));
   os.write(reinterpret_cast<const char *>(&maxSettleTime), sizeof(maxSettleTime));
   agent.serialize(os);
   os.write(reinterpret_cast<const char *>(&id), sizeof(id));
   os.write(reinterpret_cast<const char *>(&sensors), sizeof(sensors));
   os.write(reinterpret_cast<const char *>(&prodLocations), sizeof(prodLocations));
   os.write(reinterpret_cast<const char *>(&oriBias), sizeof(oriBias));
   os.write(reinterpret_cast<const char *>(&speedBias), sizeof(speedBias));
   os.write(reinterpret_cast<const char *>(&productionLifespan), sizeof(productionLifespan));
   os.write(reinterpret_cast<const char *>(&propulsionLifespan), sizeof(propulsionLifespan));
   os.write(reinterpret_cast<const char *>(&lifetime), sizeof(lifetime));
   os.write(reinterpret_cast<const char *>(&sensorLifetime), sizeof(sensorLifetime));
   os.write(reinterpret_cast<const char *>(&memory), sizeof(memory));
   os.put(static_cast<char>(sensorsReady));
   os.write(reinterpret_cast<const char *>(&readinessCount), sizeof(readinessCount));
   os.write(reinterpret_cast<const char *>(&postCollisionTimer), sizeof(postCollisionTimer));
   os.write(reinterpret_cast<const char *>(&vel), sizeof(vel));
   os.write(reinterpret_cast<const char *>(&angularVel), sizeof(angularVel));
   os.write(reinterpret_cast<const char *>(&cornerRow), sizeof(cornerRow));
   os.write(reinterpret_cast<const char *>(&cornerCol), sizeof(cornerCol));
   const auto propRows = proportions.size();
   const auto propCols = proportions.front().size();
   os.write(reinterpret_cast<const char *>(&propRows), sizeof(propRows));
   os.write(reinterpret_cast<const char *>(&propCols), sizeof(propCols));
   for (const auto& v : proportions)
      os.write(reinterpret_cast<const char *>(v.data()), v.size() * sizeof(AgentBodyType));
   if (!passthrough) {
      for (const auto& v : edgesOccludedVertical)
         os.write(reinterpret_cast<const char *>(v.data()), v.size() * sizeof(EdgeType));
      for (const auto& v : edgesOccludedHorizontal)
         os.write(reinterpret_cast<const char *>(v.data()), v.size() * sizeof(EdgeType));
   }
}


LifeEvent AgentManager::update(double dt, bool debug) { 

   lifetime += dt;
   if (postCollisionTimer > 0.0) postCollisionTimer -= dt;

   // Keep track of information for debugging, visualization, etc.
   SensorVals freeparts;
   for (size_t i = 0; i < freeparts.size(); ++i) 
      freeparts.at(i) = 1.0 - fields.propOccupied(sensors.at(i).location);
   memory.remember(Record{
         id,
         extState.pos, 
         extState.ori, 
         vel, 
         agent.provideDensGrad(), 
         agent.provideSpeed(), 
         agent.provideVel(), 
         //agent.provideAcc(), 
         agent.provideNucOut(),
         freeparts
         });

   sensorUpdate(dt);
   double lambda = sensorsReady ? agent.update(dt, debug) : 0.0;

   LifeEvent event = LifeEvent::none;
   bool success = lambda == 0.0 ? false : std::bernoulli_distribution(dt * std::abs(lambda))(gen);
   if (success) {
      if (lambda > 0.0) event = LifeEvent::reproduce;
      else event = LifeEvent::die;
   }
   return event;
}


bool AgentManager::inBounds() const {
   const Point cornerRelative{boundingBoxRelative()};
   const Point cornerLow{extState.pos - cornerRelative};
   const Point cornerHigh{extState.pos + cornerRelative};
   return fields.inBoundsXLow(cornerLow.x)
      && fields.inBoundsXHigh(cornerHigh.x)
      && fields.inBoundsYLow(cornerLow.y)
      && fields.inBoundsYHigh(cornerHigh.y);
}


// Like getting a sensor reading, but always using true concentration
// and only currently available values. Only for output (e.g. finding
// SPH smoothing lengths) and debugging.
double AgentManager::provideDensOneConc() const {
   int nUnblocked = 0;
   double sumConc = 0.0;
   for (auto& sensor : sensors) {
      if (!fields.fullyObstructed(sensor.location)) {
         sumConc += fields.concentration(DiffStates::one, sensor.location);
         ++nUnblocked;
      }
   }
   return sumConc / nUnblocked;
}


void AgentManager::Memory::remember(const Record& r) {
   if (nUsed < records.size()) ++nUsed;
   records.at(head) = r;
   head = (head + 1) % records.size();
}


// n is number of records requested, n = -1 means all available
void AgentManager::Memory::report(int n) const {
   std::cout << "Report (most recent last):\n";

   size_t nAdj = 0;
   if (n < 0) {
      nAdj = nUsed;
   } else if (static_cast<size_t>(n) > nUsed) {
      nAdj = nUsed;
      std::cerr << "Warning: " << n << " records requested; only " << nUsed << " available.\n";
   } else {
      nAdj = n;
   }

   int cur = (head + records.size() - nAdj) % records.size();
   std::cout << "begin... head: " << head << "; nAdj: " << nAdj << "; cur: " << cur << "; records.size(): " << records.size() << "\n";
   auto oldPrec = std::cout.precision(6);
   for (size_t i = 0; i < nAdj; ++i) {
      const Record& r = records.at(cur);
      std::cout 
         << "  id: " << r.id
         << "  pos: " << r.pos
         << "  ori: " << r.orientation 
         << "  absVel: " << r.absVel
         << "  densGradRel: " << r.densGradRel
         << "  speed: " << r.speed
         << "  vel: " << r.relVel
         //<< "  acc: " << r.acc
         << "\n    nucOut: ";
      for (const auto& s : r.nucOut.diffStateVals) std::cout << s << ";  ";
      for (const auto& s : r.nucOut.fieldVals) std::cout << s << ";  ";
      for (const auto& s : r.nucOut.growthSwitches) std::cout << s << ";  ";
      for (const auto& s : r.nucOut.growthRates) std::cout << s << ";  ";
      std::cout << std::endl;
      std::cout << "    Freeparts: ";
      for (auto e : r.freeparts) std::cout << e << ";  ";
      std::cout << std::endl;

      cur = (cur + 1) % records.size();
   }
   std::cout.precision(oldPrec);
   std::cout << "finish... head: " << head << "; nAdj: " << nAdj << "; cur: " << cur << "; records.size(): " << records.size() << "\n";
   std::cout << std::endl;
}
