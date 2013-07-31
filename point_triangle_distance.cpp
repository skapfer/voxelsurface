#include "point_triangle_distance.h"
#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Simple_cartesian.h>
#include "gauss_elim.h"
#include <algorithm>

namespace point_triangle_distance {

    // helper computing the inner product on CGAL Points
    template <typename POINT>
    static inline
    double inner (const POINT &a, const POINT &b) {
        return CGAL::to_double (a.x()*b.x() + a.y()*b.y() + a.z()*b.z());
    }

    // return the point closest to a in the subspace spanned by b[0]...b[NUM-1]
    // if the point is not part of the convex hull of b[0]...b[NUM-1],
    // then a tuple of NANs is returned.
    template <int NUM, typename POINT>
    static
    POINT closest_point_in_aff_subspace (const POINT &a, const POINT *b)
    {
        GaussSolver <NUM+1> gs;

        for (int i = 0; i != NUM; ++i) {
            for (int j = 0; j <= i; ++j) {
                double ip = inner (b[i], b[j]);
                gs.set_coeff (i, j, ip);
                gs.set_coeff (j, i, ip);
            }

            gs.set_coeff (i, NUM, 1.);
            gs.set_coeff (NUM, i, 1.);
            gs.set_right (i, inner (b[i], a));
        }

        gs.set_right (NUM, 1.);
        gs.set_coeff (NUM, NUM, 0.);
        gs.solve ();

        double x = 0., y = 0., z = 0.;

        for (int i = 0; i != NUM; ++i) {
            double co = gs.get_solution (i);
            int is_in_interior = (co >= 0. && co <= 1.);
            co *= is_in_interior;
            co /= is_in_interior;
            //std::cerr << co << std::endl;
            x += CGAL::to_double (b[i].x()) * co;
            y += CGAL::to_double (b[i].y()) * co;
            z += CGAL::to_double (b[i].z()) * co;
        }

        return POINT (x, y, z);
    }

    template <typename POINT>
    static inline
    void consider (std::pair <double, POINT> *ret, const POINT &a, const POINT &new_p)
    {
        double dist = CGAL::to_double (CGAL::squared_distance (a, new_p));
        if (dist < ret->first)
            ret->first = dist, ret->second = new_p;
    }

    template <typename POINT>
    std::pair <double, POINT> closest_point_on_triangle (const POINT &a, const POINT *tri)
    {
        std::pair <double, POINT> ret;
        ret.first = INFINITY;
        consider (&ret, a, tri[0]);
        consider (&ret, a, tri[1]);
        consider (&ret, a, tri[2]);

        POINT new_p = closest_point_in_aff_subspace <3> (a, tri);
        consider (&ret, a, new_p);

        POINT line[2] = { tri[0], tri[1] };
        new_p = closest_point_in_aff_subspace <2> (a, line);
        consider (&ret, a, new_p);
        line[1] = tri[2];
        new_p = closest_point_in_aff_subspace <2> (a, line);
        consider (&ret, a, new_p);
        line[0] = tri[1];
        new_p = closest_point_in_aff_subspace <2> (a, line);
        consider (&ret, a, new_p);

        return ret;
    }

    typedef CGAL::Simple_cartesian <double>::Point_3 cart_point_t;
    template std::pair <double, cart_point_t> closest_point_on_triangle (const cart_point_t &a, const cart_point_t *tri);
    typedef CGAL::Simple_homogeneous <double>::Point_3 homo_point_t;
    template std::pair <double, homo_point_t> closest_point_on_triangle (const homo_point_t &a, const homo_point_t *tri);
}
