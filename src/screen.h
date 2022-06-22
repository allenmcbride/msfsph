// screen.h
// Allen McBride
// March 23, 2022
//
// The Screen class manages real-time graphical output using SDL2.


#ifndef SCREEN_H
#define SCREEN_H

#include "colorfinder.h"
#include "fparams.h"
#include "substances.h"

#include <SDL.h>
#include <exception>
#include <string>


class Fields;


class SDLQuitException : public std::runtime_error {
   public: using runtime_error::runtime_error;
};


class Screen {

   private:

      // Data members
      const Fparams fparams;
      ColorFinder colfind;
      const int pmag;
      SDL_Window* sdlWindow;
      SDL_Renderer* sdlRenderer;
      SDL_Texture* sdlTexture;
      SDL_Surface* im;

      // If window dimensions are changed by user, reduce window size
      // as needed to keep pixels square.
      void restoreAspectRatio();


   public:

      // Construct Screen instance. fields will be passed to
      // new ColorFinder object, which will keep a reference to
      // it. fixedMaxVals is a vector with a value for each field
      // indicating either a field value for maximum color saturation
      // or a nonpositive value to signal that colors should be scaled
      // according to the maximum field value at a given time. pmag is
      // for integer scaling but is of little use now that SDL's scaling
      // is used. title is window title.
      Screen(
            const Fields& fields, 
            const SubstanceVars& fixedMaxVals, 
            int pmag, 
            std::string title
            ); // pass pmag = 0 to turn off screen output
      
      // Destructor to clean up SDL2 memory.
      ~Screen();

      // Update screen based on colors returned by colfind. 
      void update();
};

#endif
