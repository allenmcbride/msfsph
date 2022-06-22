// screen.cpp
// Allen McBride
// April 20, 2022
//
// Descriptions in screen.h

#include "screen.h"

#include "colorfinder.h"
#include "fields.h"
#include "fparams.h"
#include "substances.h"

#include <SDL.h>

#include <cmath>
#include <iostream>
#include <string>

using namespace std::literals::string_literals;


void Screen::restoreAspectRatio() {
   int width, height;
   SDL_GetWindowSize(sdlWindow, &width, &height);
   const auto correctAspectRatio = static_cast<double>(fparams.xpixels) / fparams.ypixels;
   const auto actualAspectRatio = static_cast<double>(width) / height;
   if (actualAspectRatio > correctAspectRatio) {
      const auto correctWidth = static_cast<int>(std::round(height * correctAspectRatio));
      if (width != correctWidth) {
         std::cerr << "Changing window size from " << width << " to " << correctWidth << std::endl;
         std::cerr << "Keeping window height of " << height << std::endl;
         SDL_SetWindowSize(sdlWindow, correctWidth, height);
      }
   } else {
      const auto correctHeight = static_cast<int>(std::round(width / correctAspectRatio));
      if (height != correctHeight) {
         std::cerr << "Changing window height from " << height << " to " << correctHeight << std::endl;
         std::cerr << "Keeping window width of " << width << std::endl;
         SDL_SetWindowSize(sdlWindow, width, correctHeight);
      }
   }
}


// Construct with an internal ColorFinder object that's initialized just
// once with given Fields object
Screen::Screen(
      const Fields& fields, 
      const SubstanceVars& fixedMaxVals,
      int pmag, 
      std::string title
      ) : 
   fparams{fields.fp},
   colfind{fields, fixedMaxVals},
   pmag{pmag} 
{
   if (pmag <= 0) 
      throw std::runtime_error("Screen called with pmag <= 0: " + std::to_string(pmag));

   if (SDL_Init(SDL_INIT_VIDEO) < 0)
      throw std::runtime_error("Error in SDL_Init: "s + SDL_GetError());

   if (SDL_CreateWindowAndRenderer(
         fparams.xpixels * pmag,
         fparams.ypixels * pmag,
         SDL_WINDOW_RESIZABLE,
         &sdlWindow,
         &sdlRenderer)
         < 0)
      throw std::runtime_error("Error in SDL_CreateWindowAndRenderer: "s + SDL_GetError());

   int topBorder, leftBorder, bottomBorder, rightBorder;
   SDL_GetWindowBordersSize(sdlWindow, &topBorder, &leftBorder, &bottomBorder, &rightBorder);
   int winWidth, winHeight;
   SDL_GetWindowSize(sdlWindow, &winWidth, &winHeight);
   const int winTotalWidth = winWidth + leftBorder + rightBorder;
   const int winTotalHeight = winHeight + topBorder + bottomBorder;
   SDL_Rect usableScreenRect;
   SDL_GetDisplayUsableBounds(SDL_GetWindowDisplayIndex(sdlWindow), &usableScreenRect);
   if (winTotalWidth > usableScreenRect.x)
      SDL_SetWindowSize(sdlWindow, usableScreenRect.w, winHeight);
   if (winTotalHeight > usableScreenRect.y)
      SDL_SetWindowSize(sdlWindow, winWidth, usableScreenRect.h);
   restoreAspectRatio();

   SDL_SetWindowTitle(sdlWindow, title.c_str());
         
   if (SDL_SetRenderDrawColor(sdlRenderer, 0, 0, 0, 255) < 0)
      throw std::runtime_error("Error in SDL_SetRenderDrawColor: "s + SDL_GetError());
   if (SDL_RenderClear(sdlRenderer) < 0)
      throw std::runtime_error("Error in SDL_RenderClear: "s + SDL_GetError());
   //SDL_RenderPresent(sdlRenderer);

   SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "0");

   sdlTexture = SDL_CreateTexture(
         sdlRenderer,
         SDL_PIXELFORMAT_ARGB8888,
         SDL_TEXTUREACCESS_STREAMING,
         fparams.xpixels,
         fparams.ypixels);
   if (!sdlTexture)
      throw std::runtime_error("Error in SDL_CreateTexture: "s + SDL_GetError());

   im = SDL_CreateRGBSurface(0, fparams.xpixels, fparams.ypixels, 32, 0, 0, 0, 0);
   if (!im)
      throw std::runtime_error("Error in SDL_CreateRGBSurface: "s + SDL_GetError());
}


// Destroy according to SDL rules. The SDL_SetVideoMode() memory will be freed on SDL_Quit(), but not
// the SDL_CreateRGBSurface() memory.
Screen::~Screen() {
   SDL_FreeSurface(im);
   SDL_Quit();
}


void Screen::update() {
   SDL_Event event;
   while (SDL_PollEvent(&event))
     if (event.type == SDL_QUIT)
      throw SDLQuitException("SDL_QUIT received");

   colfind.updateMaxes();

#pragma omp parallel for
   for (int r = 0; r < fparams.ypixels; ++r) {
      for (int c = 0; c < fparams.xpixels; ++c) {
         const PixColor col = colfind.findPixCol(r, c);
         const Uint32 rgb = SDL_MapRGB(im->format, col.red, col.green, col.blue); 
         ((Uint32*)(im->pixels)) [(int) (r * fparams.xpixels) + c] = rgb;
      }
   }

   restoreAspectRatio();
   if (SDL_UpdateTexture(sdlTexture, NULL, im->pixels, fparams.xpixels * sizeof(Uint32)) < 0)
      throw std::runtime_error("Error in SDL_UpdateTexture: "s + SDL_GetError());
   if (SDL_RenderClear(sdlRenderer) < 0)
      throw std::runtime_error("Error in SDL_RenderClear: "s + SDL_GetError());
   if (SDL_RenderCopy(sdlRenderer, sdlTexture, NULL, NULL) < 0)
      throw std::runtime_error("Error in SDL_RenderCopy: "s + SDL_GetError());
   SDL_RenderPresent(sdlRenderer);
}
