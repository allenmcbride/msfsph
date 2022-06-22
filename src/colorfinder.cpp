// colorfinder.cpp
// Allen McBride
// March 10, 2022
//
// Comments in colorfinder.h


#include "colorfinder.h"
#include "fields.h"
#include "substances.h"

#include <cmath>
#include <limits>


double ColorFinder::updateMax(int isub) const {
   const double headroomFactor = 1.0;
   return std::max(fields.grids.at(isub).max() * headroomFactor, std::numeric_limits<double>::min());
}


double ColorFinder::colorSum(int r, int c, int channel) const {
   double total = 0.0;
   for (size_t isub = 0; isub < nSubs; ++isub)
      total += colors.at(isub).at(channel) * fields.grids.at(isub).at(r, c) / maxVals.at(isub);
   for (size_t icue = 0; icue < nCues; ++icue)
      if (fields.cues && fields.cues->at(r, c).test(icue)) total += colors.at(nSubs + icue).at(channel);
   return total * (1.0 - fields.propOccupied(r, c));
}


ColorFinder::ColorFinder(
      const Fields& fields, 
      const SubstanceVars& fixedMaxVals,
      const Colors& colors
      ) : 
   fields{fields}, 
   fixedMaxVals{fixedMaxVals},
   colors{colors}
{
   updateMaxes();
}


PixColor ColorFinder::findPixCol(int r, int c) const {
   r = fields.fp.ypixels - r - 1; // Flip image vertically so row 0 is at bottom
   PixColor col;
   
   const auto occ = fields.propOccupied(r, c);
   const auto occGray = occ - (fields.front ? fields.front->at(r, c) : 0.0);
   const double cyan = std::min(1.0,
         colorSum(r, c, 0) +
         occGray);
   const double magenta = std::min(1.0,
         colorSum(r, c, 1) + 
         occGray);
   const double yellow = std::min(1.0,
         colorSum(r, c, 2) + 
         occGray);
   col.red = 255 * (1.0 - cyan);
   col.green = 255 * (1.0 - magenta);
   col.blue = 255 * (1.0 - yellow);

   return col;
}


void ColorFinder::updateMaxes() {
   for (int isub = 0; isub < nSubs; ++isub)
      maxVals.at(isub) = fixedMaxVals.at(isub) > 0.0 ? fixedMaxVals.at(isub) : updateMax(isub);
}
