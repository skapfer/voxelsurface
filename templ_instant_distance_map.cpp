#include "pointtype.h"
#include "distance_map_impl.h"

namespace distance_map_utils
{
    template
    void check_distance_map (const distance_map_utils::distance_map_type <Geometry>::type &dmap);

    template
    bool is_orientable (const DistanceMap &dmap, const HalfspaceMap &hsmap);

    template
    void segment_by_quantile (PhaseMap *, const DistanceMap &, double quantile);

    template
    void segment_by_threshold (PhaseMap *, const DistanceMap &, double);

    typedef cubic_voxel_field::CubicVoxelField <Geometry, int16_t> int16_field;
    typedef cubic_voxel_field::CubicVoxelField <Geometry, int8_t> int8_field;

    template
    volume_fractions_type <int8_field>::type
    compute_volume_fractions (const int8_field &phmap);
    
    template
    volume_fractions_type <int16_field>::type
    compute_volume_fractions (const int16_field &phmap);
    
    template
    void dump_volume_fractions (std::ostream &os,
                                const volume_fractions_type <int8_field>::type &vf);
    
    template
    void dump_volume_fractions (std::ostream &os,
                                const volume_fractions_type <int16_field>::type &vf);

    template
    void write_phasemap_to_bin_file (string filename, const int8_field &phmap);

    template
    void write_phasemap_to_bin_file (string filename, const int16_field &phmap);

    template
    void write_distance_map_to_ai_file (string filename, const DistanceMap &dmap);
}
