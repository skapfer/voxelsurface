#include "triangulation_search_structure.h"
#include "point_triangle_distance.h"
#include "util.h"
#include <cmath>
#include <algorithm>
#include <cassert>
#include <CGAL/centroid.h>

namespace triangulation_search_structure
{

typedef K::Plane_3 Plane_3;

SearchStructure::SearchStructure ()
{
    max_dist_from_centroid_ = 0.;
    num_tri_ = 0;
}

SearchStructure::~SearchStructure ()
{
}

void SearchStructure::add_triangle (const Point_3 *x)
{
    // determine min. angle. if one angle is smaller than 5 degrees,
    // we don't consider this triangle for determining orientations.
    Vector_3 a = x[1]-x[0];
    Vector_3 b = x[2]-x[1];
    Vector_3 c = x[0]-x[2];
    double len_a = CGAL::to_double (a.squared_length ());
    double len_b = CGAL::to_double (b.squared_length ());
    double len_c = CGAL::to_double (c.squared_length ());
    double ang_a = std::acos (-CGAL::to_double (b*c) / std::sqrt (len_b*len_c));
    double ang_b = std::acos (-CGAL::to_double (c*a) / std::sqrt (len_c*len_a));
    double ang_c = std::acos (-CGAL::to_double (a*b) / std::sqrt (len_a*len_b));
    Point_3 centroid = CGAL::centroid (x, x+3);
    double ang_error = ang_a+ang_b+ang_c - M_PI;
    max_dist_from_centroid_ = util::max4 (
        max_dist_from_centroid_,
        CGAL::squared_distance (x[0], centroid),
        CGAL::squared_distance (x[1], centroid),
        CGAL::squared_distance (x[2], centroid));

    bool dont_consider_triangle = false;
    if (! (fabs (ang_error) < 0.01))
        // symptom of weird stuff happening
        std::cerr << "[SearchStructure::add_triangle] WARNING: triangle "
                  << num_tri_ << " is degenerate and will be ignored for "
                  << "orienting.  interior angles: " 
                  << ang_a << " " << ang_b << " " << ang_c << "\n",
        dont_consider_triangle = true;
    if (! (ang_a > M_PI/180.*5.))
        dont_consider_triangle = true;
    if (! (ang_b > M_PI/180.*5.))
        dont_consider_triangle = true;
    if (! (ang_c > M_PI/180.*5.))
        dont_consider_triangle = true;

    // store away necessary triangle data
    Tri_ *final = tri_alloc_.allocate ();
    final->near_degenerate = dont_consider_triangle;
    assert (dont_consider_triangle ||
        !Plane_3 (x[0], x[1], x[2]).is_degenerate ());
    std::copy (x, x+3, final->vertices.begin ());
    final->centroid = centroid;

    // store vertices
    tree_.insert (&final->centroid);
    ++num_tri_;
}

Vector_3 SearchStructure::Tri_::normal () const
{
    return CGAL::normal (vertices[0], vertices[1], vertices[2]);
}

double SearchStructure::Tri_::distance_squared_to_point (const Point_3 &query) const
{
    return point_triangle_distance::closest_point_on_triangle (query, &vertices[0]).first;
}

int SearchStructure::Tri_::halfspace_of_point (const Point_3 &query) const
{
    if (near_degenerate)
        return LookupResult::INDETERMINATE;
    
    const Vector_3 nrml = this->normal ();
    return util::sgn (CGAL::to_double (nrml * (query - vertices[0])));
}

// about 20 triangles are traversed for each voxel.
unsigned long long total_triangles_tested = 0,
                   total_lookups = 0,
                   max_triangles_in_one_lookup = 0;

LookupResult SearchStructure::lookup_point (const Point_3 &query) const
{
    const Point_3 *x = &query;
    Neighbor_search s (tree_, TreeTraits::Point_d (x));
    Neighbor_search::iterator it = s.begin ();
    LookupResult found;
    found.distance_squared = INFINITY;
    found.halfspace = LookupResult::INDETERMINATE;
    double last_vertex_dist;
    ++total_lookups;
    unsigned long long tri_in_this_lookup = 0;
    assert (it != s.end ());
    do {
        const Tri_ *t = tri_alloc_.find_object (&*(*it).first);
        last_vertex_dist = CGAL::to_double ((*it).second);

        double dist_to_tri = t->distance_squared_to_point (query);
        if (found.distance_squared > dist_to_tri) {
            found.distance_squared = dist_to_tri;
            int new_halfspace = t->halfspace_of_point (query);
            found.halfspace = new_halfspace;
            ++it;
            continue;
        } else {
            ++it;
            ++total_triangles_tested;
            ++tri_in_this_lookup;
        }
    } while (found.distance_squared + 1.05*max_dist_from_centroid_ > last_vertex_dist && it != s.end ());

    max_triangles_in_one_lookup = std::max (max_triangles_in_one_lookup, tri_in_this_lookup);
    return found;
}

}
