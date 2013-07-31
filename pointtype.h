// FIXME streamline this
#include "triangulation_search_structure.h"
#include "cubic_voxel_field.h"
#include "distance_map.h"

using triangulation_search_structure::Point_3;
using triangulation_search_structure::Vector_3;

typedef cubic_voxel_field::Geometry <Point_3, Vector_3> Geometry;
typedef distance_map_utils::distance_map_type <Geometry>::type DistanceMap;
typedef distance_map_utils::halfspace_map_type <Geometry>::type HalfspaceMap;
typedef cubic_voxel_field::CubicVoxelField <Geometry, int8_t> PhaseMap;
