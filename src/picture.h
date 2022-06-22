// picture.h
// Allen McBride
// March 23, 2022
//
// The Picture class manages output of a series of PNG files. It works
// similarly to the Screen class, with an owned ColorFinder object that
// keep a reference to the simulation's Fields object.


#ifndef PICTURE_H
#define PICTURE_H

#include "colorfinder.h"
#include "fparams.h"
#include "substances.h"

#include <string>

class Fields;

class Picture {
   private:

      //Data members
      const Fparams fparams;
      ColorFinder colfind;
      const std::string picdir;
      const int countInc;
      int count;


   public:

      // Constructor is analogous to Screen constructor, but with a
      // starting number and an increment for each additional picture.
      Picture(
            const Fields& fields, 
            const SubstanceVars& fixedMaxVals, 
            const std::string& picdir, 
            int countBegin = 1, 
            int countInc = 1
            );

      // Write a new picture to a successively-numbered file.
      void update();
};

#endif
