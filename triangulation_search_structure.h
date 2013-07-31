#ifndef X45bca305598010bed8798032d285a9860e83c459
#define X45bca305598010bed8798032d285a9860e83c459

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <boost/array.hpp>
#include <vector>
#include "aligned_allocator.h"

namespace triangulation_search_structure {

typedef CGAL::Simple_cartesian <double> K;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;

struct TreeTraits {
    typedef CGAL::Search_traits_3 <K> mother;
    typedef mother::FT FT;
    typedef mother::Cartesian_const_iterator_d Cartesian_const_iterator_d;

    struct Point_d {

        Point_d (const Point_3 *p) {
            ptr_ = p;
        }

        const Point_3 &operator* () const {
            return *ptr_;
        }

        Cartesian_const_iterator_d cartesian_begin () const {
            return mother::Construct_cartesian_const_iterator_d ()(*ptr_);
        }

        Cartesian_const_iterator_d cartesian_end () const {
            return mother::Construct_cartesian_const_iterator_d ()(*ptr_, 0);
        }

    private:
        Point_d ();
        const Point_3 *ptr_;
    };

    struct Construct_cartesian_const_iterator_d {
        Cartesian_const_iterator_d operator() (const Point_d &p) const {
            return p.cartesian_begin ();
        }
        Cartesian_const_iterator_d operator() (const Point_d &p, int) const {
            return p.cartesian_end ();
        }
    };
};

struct LookupResult {
    double distance_squared;
    enum { POSITIVE_HALFSPACE = 1, NEGATIVE_HALFSPACE = -1, INDETERMINATE = 0 };
    int halfspace;
};

class SearchStructure {
public:
    SearchStructure ();
    ~SearchStructure ();
    
    void add_triangle (const Point_3 *x);
    void clear ();
    LookupResult lookup_point (const Point_3 &) const;

private:
    typedef CGAL::Orthogonal_incremental_neighbor_search <TreeTraits> Neighbor_search;
    // FIXME do we really need this to be mutable?
    mutable Neighbor_search::Tree tree_;
    class Tri_ {
    public:
        boost::array <Point_3, 3> vertices;
        Point_3 centroid;
        // flag to tell that tri is skinny and won't be used for orienting
        bool near_degenerate;

        Vector_3 normal () const;
        double distance_squared_to_point (const Point_3 &query) const;
        int halfspace_of_point (const Point_3 &query) const;
    };
    AlignedAllocator <Tri_> tri_alloc_;
    int num_tri_;
    double max_dist_from_centroid_;
};

extern
unsigned long long total_triangles_tested,
                   total_lookups,
                   max_triangles_in_one_lookup;

}

#endif // X45bca305598010bed8798032d285a9860e83c459
