// colorfinder.h
// Allen McBride
// March 10, 2022
//
// For each field, find the pixel color to represent each grid
// cell, based on the mapping from field values to colors defined in
// substances.h. fixedMaxVals is one of the constructor arguments; it is
// an array with an entry for each field. If zero, colors representing
// that field are re-scaled at each updateMaxes() call to the max value
// for that field over space. If nonzero, colors are scaled to the
// value given.
//
// The actual specifications of colors for each substance is found
// in substances.h


#ifndef COLORFINDER_H
#define COLORFINDER_H

#include "substances.h"


class Fields;


struct PixColor {
   unsigned char red;
   unsigned char green;
   unsigned char blue;
};


class ColorFinder {
   private:
      const Fields& fields;
      const SubstanceVars fixedMaxVals;
      const Colors colors;
      SubstanceVars maxVals;
      double updateMax(int isub) const;

      // For each color channel (cyan, magenta, yellow), find appropriate
      // value based on field values and color specifications for
      // different substances.
      double colorSum(int r, int c, int channel) const;
   public:
      ColorFinder(
            const Fields& fields, 
            const SubstanceVars& fixedMaxVals,
            const Colors& colors = defaultColors 
            );
      PixColor findPixCol(int r, int c) const;

      // For each grid that is not given a fixed maximum, update the
      // value that corresponds to full color saturation based on the
      // maximum value present in the grid.
      void updateMaxes();
};

#endif
