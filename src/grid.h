// grid.h
// Allen McBride
// March 23, 2022
//
// This file contains two similar classes, Grid and EdgeGrid, each
// templated for grids of a given type. The main purpose of these classes
// is to encapsulate accesses by row and column. If XWRAP is defined,
// column accesses are wrapped for a cylindrical topology. The main
// purpose of making EdgeGrid a separate class is so that a grid of
// values for the edges of Grid cells can share the same Fparams object
// as the corresponding Grid.


#ifndef GRID_H
#define GRID_H

#include "fparams.h"
#include "point.h"
#include "serialization.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>


template <typename T>
class Grid {

   private:

      // Data members (size information and actual grid data)
      const Fparams fp;
      std::vector<T> data;
      
      // Check that sizes of two Grids are the same
      bool checkFpMatch(const Grid& other) const { 
         return fp.xpixels == other.fp.xpixels && fp.ypixels == other.fp.ypixels; 
      }


   public:

      // Construct Grid of given dimensions, filling with given value.
      Grid(const Fparams& fp, const T val) : fp{fp}, data(fp.xpixels * fp.ypixels, val) {}

      // Unfinished deserializing constructor
      Grid(std::istream& is) :
         fp{is}
      {
         data.resize(fp.xpixels * fp.ypixels);
         is.read(reinterpret_cast<char *>(data.data()), data.size() * sizeof(T));
      }

      // Fill Grid with given value
      void fill(const T val) { std::fill(data.begin(), data.end(), val); }

      // Swap data with another Grid
      void swap(Grid& other) { 
         if (checkFpMatch(other)) data.swap(other.data);
         else throw std::runtime_error("Cannot swap grids with different dimensions");
      }

      // Copy another Grid's data into this one.
      void copy(const Grid& other) {
         if (checkFpMatch(other)) data = other.data;
         else throw std::runtime_error("Cannot copy from grid with different dimensions");
      }

      // Reference to grid cell corresponding to given row, column
#ifdef XWRAP
      T& at(int row, int col) { return data[row * fp.xpixels + (col % fp.xpixels + fp.xpixels) % fp.xpixels]; }
      const T& at(int row, int col) const { return data[row * fp.xpixels + (col % fp.xpixels + fp.xpixels) % fp.xpixels]; }
#else
      T& at(int row, int col) { return data[row * fp.xpixels + col]; }
      const T& at(int row, int col) const { return data[row * fp.xpixels + col]; }
#endif

      // Reference to grid cell corresponding to point
      T& at(const Point& p) { return at(fp.yToRow(p.y), fp.xToCol(p.x)); }
      const T& at(const Point& p) const { return at(fp.yToRow(p.y), fp.xToCol(p.x)); }

      // Sum of data, value of min and max elements
      T sum() const { return std::accumulate(data.begin(), data.end(), T{}); }
      T min() const { return *std::min_element(data.begin(), data.end()); }
      T max() const { return *std::max_element(data.begin(), data.end()); }

      // Calculate Chebyshev distance with another Grid (absolute
      // difference between the pair of corresponding cells with maximum
      // absolute difference).
      T chebyshev(const Grid<T>& other) const {
         if (other.data.size() != data.size())
            throw std::logic_error("Chebyshev called on two grids of different size.");
         T max{};
         for (size_t i = 0; i < data.size(); ++i)
            max = std::max(max, std::abs(data.at(i) - other.data.at(i)));
         return max;
      }

      // Print Grid to file of given name, unformatted
      void print(std::string fn) { 
         std::ofstream f(fn);
         for (const auto& e : data) f << e << std::endl;
      }

      // Print Grid to given stream, formatted
      friend std::ostream& operator<<(std::ostream& s, const Grid<T>& a) {
         //s << "ysize: " << a.fp.ypixels << "; xsize: " << a.fp.xpixels << std::endl;
         for (int r = 0; r < a.fp.ypixels; ++r) {
            for (int c = 0; c < a.fp.xpixels; ++c)
               s << std::setprecision(3) << a.at(r, c) << " ";
            s << std::endl;
         }
         return s;
      }

      // Unfinished serializer
      void serialize(std::ostream& os) const {
         os.write(reinterpret_cast<const char *>(&fp), sizeof(fp));
         os.write(reinterpret_cast<const char *>(data.data()), data.size() * sizeof(T));
      }
};


template <typename T>
class EdgeGrid {

   private:

      // Data members (size and actual data)
      const int xsize;
      const int ysize;
      std::vector<T> data;


   public:

      // Construct EdgeGrid of given size, filling it with given
      // value. "Horizontal" specifies whether the edges stored are
      // those parallel with the x-axis (true) or y-axis (false).
      EdgeGrid(const Fparams& fp, const T val, bool horizontal) : 
         xsize{fp.xpixels - (horizontal ? 0 : 1)},
         ysize{fp.ypixels - (horizontal ? 1 : 0)},
         data(xsize * ysize, val)
         {}
      
      // Unfinished deserializing constructor
      EdgeGrid(std::istream& is) :
         xsize{readObjectFromStream<int>(is)},
         ysize{readObjectFromStream<int>(is)}
      {
         data.resize(xsize * ysize);
         is.read(reinterpret_cast<char *>(data.data()), data.size() * sizeof(T));
      }
      
      // Fill EdgeGrid with given value.
      void fill(const T val) { std::fill(data.begin(), data.end(), val); }

      // Return reference to grid cell for given row, column.
      T& at(int row, int col) { return data[row * xsize + col]; }
      const T& at(int row, int col) const { return data[row * xsize + col]; }

      // Unfinished serializer
      void serialize(std::ostream& os) const {
         os.write(reinterpret_cast<const char *>(&xsize), sizeof(xsize));
         os.write(reinterpret_cast<const char *>(&ysize), sizeof(ysize));
         os.write(reinterpret_cast<const char *>(data.data()), data.size() * sizeof(T));
      }
};

#endif
