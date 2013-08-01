#include "distance_map.h"
#include "util.h"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <algorithm>

namespace distance_map_utils {

    template <typename DISTMAPTYPE>
    void check_distance_map (const DISTMAPTYPE &dmap)
    {
        typedef typename DISTMAPTYPE::VoxelNeighbors vox_neighbors_type;
        bool found_inconsistencies = false;

        // simple consistency check, uses the triangle inequality
        ITERATE_OVER_VOXEL_FIELD (i,j,k, dmap)
        {
            using util::field_width;
            vox_neighbors_type n (dmap, i,j,k);
            for (typename vox_neighbors_type::key_type it = n.begin (); it != n.end (); ++it)
            {
                double it_offset = dmap.voxel_squared_distance (
                    it.i () - i, it.j () - j, it.k () - k);
                if (! (sqrt (it_offset) + sqrt (it[dmap]) >= sqrt (dmap[i][j][k])*0.995))
                {
                    if (!found_inconsistencies)
                        std::cerr << "inconsistencies in distance map:\n";
                    found_inconsistencies = true;
                    std::cerr << field_width (4, i) << " "
                              << field_width (4, j) << " "
                              << field_width (4, k) << " "
                              << field_width (10, sqrt (it_offset)) << " +  "
                              << field_width (10, sqrt (it[dmap])) << " >= "
                              << field_width (10, sqrt (dmap[i][j][k])) << " viol = "
                              << field_width (10, sqrt (it_offset) + sqrt (it[dmap]) - sqrt (dmap[i][j][k])) << "\n";
                }

                if (! (sqrt (it[dmap])*0.995 <= sqrt (it_offset) + sqrt (dmap[i][j][k])))
                {
                    if (!found_inconsistencies)
                        std::cerr << "inconsistencies in distance map:\n";
                    found_inconsistencies = true;
                    std::cerr << field_width (4, i) << " "
                              << field_width (4, j) << " "
                              << field_width (4, k) << " "
                              << field_width (10, sqrt (it[dmap])) << " <= "
                              << field_width (10, sqrt (it_offset)) << " +  "
                              << field_width (10, sqrt (dmap[i][j][k])) << " viol = "
                              << field_width (10, sqrt (it[dmap]) - (sqrt (it_offset) + sqrt (dmap[i][j][k]))) << "\n";
                }

            }
        }

        if (found_inconsistencies)
            std::cerr << "inconsistencies found, aborting" << util::ABORT;
    }

    template <typename DISTMAPTYPE, typename HALFSPACEMAP>
    bool is_orientable      (const DISTMAPTYPE &dmap,
                             const HALFSPACEMAP &hsmap)
    {
        typedef typename DISTMAPTYPE::VoxelNeighbors vox_neighbors_type;
        double vox_ld = util::max3 (dmap.voxel_diagonal ()[0],
            dmap.voxel_diagonal ()[1], dmap.voxel_diagonal ()[2]);

        ITERATE_OVER_VOXEL_FIELD (i,j,k, dmap)
        {
            vox_neighbors_type n (dmap, i,j,k);
            if (dmap[i][j][k] > vox_ld)
                for (typename vox_neighbors_type::key_type it = n.begin (); it != n.end (); ++it)
                    if (hsmap[i][j][k] != it[hsmap])
                        return false;
        }

        return true;
    }

    template <typename DISTMAP>
    double improve_quantile_guess (const DISTMAP &dmap, double min_, double max_,
                                   double prec, double quantile, int remaining_iter = 50)
    {
        double mid_ = .5 * (max_ + min_);
        util::PercentageBelow <double> a (mid_);
        a = for_each_voxel (dmap, a);

        if (remaining_iter == 0) {
            std::cerr << "[distance_map_utils::improve_quantile_guess] Could "
                         "not achieve prescribed volume fraction.  Actual vf = "
                      << a.value () << std::endl;
        }

        if (max_ - min_ < prec)
            return mid_;

        if (a.value () < quantile)
            return improve_quantile_guess (dmap, mid_, max_, prec, quantile);
        else
            return improve_quantile_guess (dmap, min_, mid_, prec, quantile);
    }

    template <typename DISTMAP>
    double find_quantile (const DISTMAP &dmap, double quantile)
    {
        /*
        std::vector <typename DISTMAP::value_type> tmp;
        tmp.reserve (size_t (dmap.nvox[0]) * dmap.nvox[1] * dmap.nvox[2]);
        ITERATE_OVER_VOXEL_FIELD (i,j,k, dmap)
            tmp.push_back (dmap[i][j][k]);
        size_t n = tmp.size () * quantile;
        std::nth_element (tmp.begin (), tmp.begin () + n, tmp.end ());
        return tmp[n];
        */
        util::MinMax <double> mm;
        mm = for_each_voxel (dmap, mm);
        return improve_quantile_guess (dmap, mm.min (), mm.max (),
                                       (mm.max () - mm.min ()) * 1e-6, quantile);
    }

    template <typename PHASEMAP, typename DISTMAP>
    void segment_by_quantile (PHASEMAP *phmap, const DISTMAP &dmap, double quantile)
    {
        double q = find_quantile (dmap, quantile);
        std::cerr << "quant = " << q << "\n";
        ITERATE_OVER_VOXEL_FIELD (i,j,k, dmap)
            (*phmap)[i][j][k] = (dmap[i][j][k] < q);
    }

    template <typename PHASEMAP, typename DISTMAP>
    void segment_by_threshold (PHASEMAP *phmap, const DISTMAP &dmap, double thr)
    {
        ITERATE_OVER_VOXEL_FIELD (i,j,k, dmap)
            (*phmap)[i][j][k] = (dmap[i][j][k] < thr);
    }

    template <typename PHASEMAP>
    void write_phasemap_to_bin_file (std::ostream &os, const PHASEMAP &phmap)
    {
        os << phmap.nvox[2] << " " << phmap.nvox[1] << " " << phmap.nvox[0] << "\n";
        for (int i = 0; i != phmap.nvox[0]; ++i)
        for (int j = 0; j != phmap.nvox[1]; ++j) {
            for (int k = 0; k != phmap.nvox[2]; ++k)
                os << char (phmap[i][j][k] + '0');
            os << "\n";
        }

        os.flush ();
    }

    template <typename PHASEMAP>
    typename volume_fractions_type <PHASEMAP>::type
    compute_volume_fractions (const PHASEMAP &phmap)
    {
        typename volume_fractions_type <PHASEMAP>::type ret;
        ITERATE_OVER_VOXEL_FIELD (i,j,k, phmap)
            ++ret[phmap[i][j][k]];
        return ret;
    }

    template <typename TYPE>
    const TYPE &no_char (const TYPE &v) { return v; }

    int no_char (const char &v) { return v; }
    int no_char (const signed char &v) { return v; }
    int no_char (const unsigned char &v) { return v; }

    // write the contents of a std::map to an iostream
    // (for debugging, not in public interface)
    template <typename STDMAP>
    void dump_a_map (std::ostream &os, const STDMAP &map)
    {
        typename STDMAP::const_iterator it = map.begin ();
        typename STDMAP::const_iterator end_ = map.end ();
        if (it == end_) {
            os << "{ }";
        } else {
            --end_;
            os << "{ ";
            for (; it != end_; ++it)
                os << no_char (it->first) << " -> " << it->second << ", ";
            os << no_char (it->first) << " -> " << it->second << "}";
        }
    }

    // write the contents of a volume fractions object to an iostream
    // (for debugging)
    template <typename VOLFRACTYPE>
    void dump_volume_fractions (std::ostream &os,
                                const VOLFRACTYPE &vf)
    {
        dump_a_map (os, vf);
    }

    template <typename DISTANCMAP>
    void write_distance_map_to_ai_file (string filename, const DISTANCMAP &dmap)
    {
        std::ofstream os (filename.c_str ());
        if (!os)
            throw std::runtime_error ("Unable to open " + filename + " for writing.");
        for (int i = 0; i != dmap.nvox[0]; ++i)
        for (int j = 0; j != dmap.nvox[1]; ++j) {
            for (int k = 0; k != dmap.nvox[2]; ++k)
                os << dmap[i][j][k] << " ";
            os << "\n";
        }
        if (!os)
            throw std::runtime_error ("Error writing data to " + filename + ".");
        os.close ();
    }

    template <typename PHASEMAP>
    void write_phasemap_to_bin_file (string filename, const PHASEMAP &phmap)
    {
        std::ofstream os (filename.c_str ());
        if (!os)
            throw std::runtime_error ("Unable to open " + filename + " for writing.");
        write_phasemap_to_bin_file (os, phmap);
        if (!os)
            throw std::runtime_error ("Error writing data to " + filename + ".");
        os.close ();
    }
}
