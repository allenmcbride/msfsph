// point.h
// Allen McBride
// June 21, 2022
//
// Point is a simple struct to encapsulate 2D vectors along with the
// various operations performed on them. Its style could be improved
// based on the links below.

#ifndef POINT_H
#define POINT_H

#include <algorithm>
#include <cmath>
#include <iostream>

// TODO might help: https://codereview.stackexchange.com/questions/26608/review-of-2d-vector-class
struct Point {
   double x;
   double y;
   double mag() const { return std::hypot(x, y); }
   Point normed() const { return isZero() ? Point{0.0, 0.0} : *this / mag(); }
   Point rotated(const double th) const { return Point{x * std::cos(th) - y * std::sin(th), x * std::sin(th) + y * std::cos(th)}; }
   double angle() const { return isZero() ? 0.0 : std::atan2(y, x); }
   bool isZero() const { return x == 0.0 && y == 0.0; }
   double dot(const Point& rhs) const { return x * rhs.x + y * rhs.y; }
   Point projection(const Point& onto) const { return (dot(onto) / onto.dot(onto)) * onto; }
   double scalarProjection(const Point& onto) const { return dot(onto) / onto.mag(); }
   double angleDiff(const Point& rhs) const { return angleDiff(angle(), rhs.angle()); }
   double angleBetween(const Point& rhs) const { return isZero() || rhs.isZero() ? 0.0 : std::acos(std::clamp(dot(rhs) / (mag() * rhs.mag()), -1.0, 1.0)); }

   // These are quick-and-dirty; canonical forms require compound versions and move constructors
   // TODO https://stackoverflow.com/questions/4581961/c-how-to-overload-operator?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
   Point& operator+=(const Point& a) { x += a.x; y += a.y; return *this; }
   Point& operator-=(const Point& a) { x -= a.x; y -= a.y; return *this; }
   Point& operator*=(const double k) { x *= k; y *= k; return *this; }
   friend Point operator-(const Point& p) { return Point{-p.x, -p.y}; }
   friend Point operator+(const Point& lhs, const Point& rhs) { return Point{lhs.x + rhs.x, lhs.y + rhs.y}; }
   friend Point operator-(const Point& lhs, const Point& rhs) { return Point{lhs.x - rhs.x, lhs.y - rhs.y}; }
   friend Point operator*(const double k, const Point& a) { return Point{k * a.x, k * a.y}; }
   friend Point operator/(const Point& a, const double k) { return Point{a.x / k, a.y / k}; }
   friend bool operator==(const Point& lhs, const Point& rhs) { return lhs.x == rhs.x && lhs.y == rhs.y; }
   friend bool operator!=(const Point& lhs, const Point& rhs) { return lhs.x != rhs.x || lhs.y != rhs.y; }
   friend std::ostream& operator<<(std::ostream& s, const Point& a) { s << "<" << a.x << "," << a.y << ">"; return s; }

   static Point polarPt(double r, double th) { return r * Point{std::cos(th), std::sin(th)}; }
   static double angleDiff(double angle1, double angle2) { return std::remainder(angle1 - angle2, 2.0 * M_PI); }
};

#endif
