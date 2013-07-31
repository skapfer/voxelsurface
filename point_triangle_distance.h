#ifndef POINT_TRIANGLE_DISTANCE 
#define POINT_TRIANGLE_DISTANCE 
#include <ostream>
#include <utility>

namespace point_triangle_distance
{
    template <typename POINT>
    std::pair <double, POINT> closest_point_on_triangle (const POINT &a, const POINT *tri);
}

#endif /* POINT_TRIANGLE_DISTANCE */
