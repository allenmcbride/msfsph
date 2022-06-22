// serialization.h
// Allen McBride
// April 20, 2022
//
// Unfinished. Reinterpret bytes from stream as object of desired type.


#ifndef SERIALIZATION_H
#define SERIALIZATION_H

//TODO consider using decltype for the template parameter every time
//this is called

template <typename T>
T readObjectFromStream(std::istream& is) {
   T x;
   is.read(reinterpret_cast<char *>(&x), sizeof(x));
   return x;
}

#endif
