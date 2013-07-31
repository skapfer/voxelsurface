#ifndef X56764b0bd70296236efadcbd3a45b8f4b14a96d1 
#define X56764b0bd70296236efadcbd3a45b8f4b14a96d1 
#include "cubic_voxel_field.h"
#include <string>

namespace distance_map_utils {

    using std::string;

    template <typename GEOMETRY>
    struct distance_map_type
    {
        typedef cubic_voxel_field::CubicVoxelField <GEOMETRY, double> type;
    };

    template <typename GEOMETRY>
    struct halfspace_map_type
    {
        typedef cubic_voxel_field::CubicVoxelField <GEOMETRY, int8_t> type;
    };
    
    template <typename DISTMAPTYPE>
    void check_distance_map (const DISTMAPTYPE &dmap);

    // return whether the halfspacemap given looks consistent;
    // in particular whether -1's and 1's are adjacent.
    template <typename DISTMAPTYPE, typename HALFSPACEMAP>
    bool is_orientable      (const DISTMAPTYPE &dmap,
                             const HALFSPACEMAP &hsmap);

    template <typename PHASEMAP, typename DISTMAP>
    void segment_by_quantile (PHASEMAP *phmap, const DISTMAP &dmap, double quantile = 0.5);

    template <typename PHASEMAP, typename DISTMAP>
    void segment_by_threshold (PHASEMAP *phmap, const DISTMAP &dmap, double threshold = 0.);

    template <typename PHASEMAP, typename TRANSTBL>
    void translate_phase_indices (PHASEMAP *phmap, TRANSTBL tbl)
    {
        ITERATE_OVER_VOXEL_FIELD (i,j,k, *phmap)
        {
            assert ((*phmap)[i][j][k] >= 0);
            (*phmap)[i][j][k] = tbl[(*phmap)[i][j][k]];
        }
    }

    template <typename PHASEMAP>
    struct volume_fractions_type {
        typedef
        std::map <typename PHASEMAP::value_type, typename PHASEMAP::count_type>
        type;

        typedef typename type::iterator iterator;
        typedef typename type::const_iterator const_iterator;
    };

    template <typename PHASEMAP>
    typename volume_fractions_type <PHASEMAP>::type
    compute_volume_fractions (const PHASEMAP &);

    template <typename VOLFRACTYPE>
    void dump_volume_fractions (std::ostream &os,
                                const VOLFRACTYPE &vf);

    template <typename PHASEMAP>
    void write_phasemap_to_bin_file (string filename, const PHASEMAP &);

    template <typename DISTANCEMAP>
    void write_distance_map_to_ai_file (string, const DISTANCEMAP &);
}

#endif /* X56764b0bd70296236efadcbd3a45b8f4b14a96d1 */
