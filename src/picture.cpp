// picture.cpp
// Allen McBride
// April 20, 2022
//
// Descriptions in picture.h


#include "picture.h"

#include "colorfinder.h"
#include "fields.h"
#include "fparams.h"
#include "substances.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>


Picture::Picture(
      const Fields& fields, 
      const SubstanceVars& fixedMaxVals, 
      const std::string& picdir, 
      const int countBegin, 
      int countInc
      ) : 
   fparams{fields.fp}, 
   colfind{fields, fixedMaxVals}, 
   picdir{picdir}, 
   countInc{countInc}, 
   count{countBegin} 
{
   std::filesystem::create_directories(picdir);
}


void Picture::update() {
   colfind.updateMaxes();

   std::string basename = picdir + "/pathcuda-" + std::to_string(count);
   std::string filename = basename + ".ppm";
   std::ofstream pathpic(filename);

   pathpic << "P6\n"
      << std::to_string(fparams.xpixels)
      << " "
      << std::to_string(fparams.ypixels)
      << "\n255\n";
   for (int r = 0; r < fparams.ypixels; ++r) {
      for (int c = 0; c < fparams.xpixels; ++c) {
         PixColor col = colfind.findPixCol(r, c);
         pathpic.put(col.red);
         pathpic.put(col.green);
         pathpic.put(col.blue);
      }
   }
   pathpic.close();

   std::string command = "convert " + basename + ".ppm " + basename + ".png ; rm " + basename + ".ppm";
   std::system(command.c_str());

   count += countInc;
}
