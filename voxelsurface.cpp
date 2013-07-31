#include <boost/program_options.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>
#include "pointtype.h"
#include "progressbar.h"
#include "readpoly.h"
#include "util.h"
#include <fstream>

using triangulation_search_structure::SearchStructure;
namespace po = boost::program_options;

class SearchStructureAdapter : public TriangleSink {
public:
    SearchStructureAdapter (SearchStructure *x) {
        sink_ = x;
    }

    virtual void insert_triangle (const double *v0,
        const double *v1, const double *v2)
    {
        Point_3 p[3] = { Point_3 (v0[0], v0[1], v0[2]),
            Point_3 (v1[0], v1[1], v1[2]), Point_3 (v2[0], v2[1], v2[2]) };
        sink_->add_triangle (p);
    }

private:
    SearchStructure *sink_;
};

static
Point_3 parse_point_3 (std::string str)
{
    {
        std::vector <std::string> tokens;
        using namespace boost::algorithm;
        split (tokens, str, is_any_of (","), token_compress_on);
        if (tokens.size () != 3u)
            goto error;
        double parsed_floats[3];
        std::transform (tokens.begin (), tokens.end (), parsed_floats,
            &boost::lexical_cast <double, std::string>);
        return Point_3 (parsed_floats[0], parsed_floats[1], parsed_floats[2]);
    }
error:
    throw std::runtime_error (str + " is not a valid floating-point triplet.");
}

static
std::vector <std::string> split_at_commas (std::string str)
{
    std::vector <std::string> tokens;
    using namespace boost::algorithm;
    split (tokens, str, is_any_of (","), token_compress_off);
    return tokens;
}

template <typename DISTMAP, typename HALFSPACEMAP>
void heal_unoriented_voxels (const DISTMAP &dmap, HALFSPACEMAP *orimap)
{
    bool another_iteration;
    do {
        another_iteration = false;
        double furthest_unoriented_voxel = 0.;
        ITERATE_OVER_VOXEL_FIELD(i,j,k,dmap)
        {
            if ((*orimap)[i][j][k] != 0) continue;
            int a = cubic_voxel_field::for_each_neighbor (*orimap, i,j,k, util::summer (0));
            if (a > 5 || a < -5) {
                a = util::sgn (a);
                (*orimap)[i][j][k] = a;
                //std::cerr << "fixed voxel: " << dmap[i][j][k] << "\n";
                another_iteration = true;
            } else {
                if (dmap[i][j][k] > furthest_unoriented_voxel)
                    furthest_unoriented_voxel = dmap[i][j][k];
                //std::cerr << "unfixable voxel: " << dmap[i][j][k] << "\n";
            }
        }
        if (another_iteration)
            std::cerr << "\n    furthest indeterminate voxel has distance "
                      << furthest_unoriented_voxel << "\n";
    } while (another_iteration);
}


int main (int argc, char **argv) {

    po::options_description opt_descr ("Allowed options");
    std::string operation_mode;
    std::string filled_output;
    std::string single_input_file, multi_input_file, output_file;
    std::string inflation_series_control;
    std::string parallel_series_control;
    std::string origin, upper_corner;
    std::string quantile_filename;
    int nvoxx, nvoxy = -1, nvoxz = -1;
    bool aperiodic = false;
    double single_inflate_volfrac = .5;
    opt_descr.add_options ()
        ("mode",
            po::value <std::string> (&operation_mode)->default_value ("default"),
            "Mode of operation")
        ("help",
            "Show help message")
        ("filled-output",
            po::value <std::string> (&filled_output),
            "Output file for labyrinth filling mode")
        ("parallel-series-control",
            po::value <std::string> (&parallel_series_control),
            "Map volfrac --> filename")
        ("inflated-output",
            po::value <std::string> (&output_file),
            "Output file for surface inflation mode")
        ("inflation-volfrac",
            po::value <double> (&single_inflate_volfrac)->default_value (.5),
            "Volume fraction for surface inflation mode")
        ("inflation-series-control",
            po::value <std::string> (&inflation_series_control),
            "Map volfrac --> filename")
        ("quantile-output",
            po::value <std::string> (&quantile_filename),
            "File to write distance quantiles to")
        ("surface",
            po::value <std::string> (&single_input_file)->default_value ("input.poly"),
            "Input file (triangulated surface)")
        ("multi-surface",
            po::value <std::string> (&multi_input_file)->default_value ("<unset>"),
            "Input files (multiple triangulated surfaces, separated by commas)")
        ("origin",
            po::value <std::string> (&origin)->default_value ("0,0,0"),
            "Origin of the discretization box")
        ("upper-corner",
            po::value <std::string> (&upper_corner)->default_value ("1,1,1"),
            "1,1,1 point of the discretization box")
        ("aperiodic",
            po::value <bool> (&aperiodic)->default_value (false),
            "whether the discretization box is a periodic unit cell")
        ("discret",
            po::value <int> (&nvoxx)->default_value (256),
            "Discretization in x direction")
        ;

    po::variables_map vm;
    po::store (po::parse_command_line (argc, argv, opt_descr), vm);
    po::notify (vm); 

    if (vm.count ("help"))
    {
        std::cerr << opt_descr;
        return 0;
    }

    // determine input files
    if (single_input_file != "input.poly" && multi_input_file != "<unset>")
    {
        std::cerr << "Only one of --surface and --multi-surface may "
            "be given at the command line.\n";
        return 1;
    }
    std::vector <std::string> input_files;
    if (multi_input_file != "<unset>")
        input_files = split_at_commas (multi_input_file);
    else
        input_files.push_back (single_input_file);

    // set up geometry
    Geometry g;
    g.origin = parse_point_3 (origin);
    g.upper_corner = parse_point_3 (upper_corner);
    g.is_periodic = !aperiodic;

    if (nvoxx < 1 || nvoxx > 4096)
        std::cerr << "nvoxx makes no sense: " << nvoxx << "\n" << util::ABORT;
    if (nvoxy != -1 || nvoxz != -1)
        std::cerr << "Setting discretization-y and discretization-z presently unsupported.\n"
                  << util::ABORT;
    if (nvoxy == -1)
        nvoxy = nvoxx * (g.upper_corner - g.origin)[1] / (g.upper_corner - g.origin)[0];
    if (nvoxz == -1)
        nvoxz = nvoxx * (g.upper_corner - g.origin)[2] / (g.upper_corner - g.origin)[0];
    g.nvox[0] = nvoxx;
    g.nvox[1] = nvoxy;
    g.nvox[2] = nvoxz;
    std::cerr << "geometry = " << g << std::endl;

    if (operation_mode == "default")
    {
        bool got_something_to_do = false;

        if (output_file != "")
            std::cerr << "todo += inflated surface\n", got_something_to_do = true;
        if (inflation_series_control != "")
            std::cerr << "todo += inflation series\n", got_something_to_do = true;
        if (filled_output != "")
            std::cerr << "todo += filled labyrinth\n", got_something_to_do = true;
        if (parallel_series_control != "")
            std::cerr << "todo += parallel body series\n", got_something_to_do = true;
        if (quantile_filename != "")
            std::cerr << "todo += distance quantiles\n", got_something_to_do = true;

        if (!got_something_to_do)
            std::cerr << "you have not told me what to do. refusing to compute nothing.\n"
                      << util::ABORT;

        SearchStructure tss;
        for (unsigned i = 0; i != input_files.size (); ++i)
        {
            std::cerr << "input file[" << i << "] = " << input_files[i] << "\n";
            SearchStructureAdapter tssa (&tss);
            PolyFileSink *pfs = poly_file_sink_from_triangle_sink (&tssa);
            std::ifstream is (input_files[i].c_str ());
            parse_poly_file (pfs, is);
            delete pfs;
        }

        DistanceMap dmap (g, NAN);
        HalfspaceMap omap (g, 0);
        std::cerr << "computing distance map ";
        ProgressBar pbar (g.num_voxels ());
        ITERATE_OVER_VOXEL_FIELD (i,j,k, g)
        {
            ++pbar;
            triangulation_search_structure::LookupResult l
                = tss.lookup_point (g.voxel_center (i, j, k));
            dmap[i][j][k] = l.distance_squared;
            omap[i][j][k] = l.halfspace;
        }
        std::cerr << "checking distance map ";
        distance_map_utils::check_distance_map (dmap);
        std::cerr << "\nortienting distance map\n    ";
        distance_map_utils::dump_volume_fractions (std::cerr,
            distance_map_utils::compute_volume_fractions (omap));
        heal_unoriented_voxels (dmap, &omap);
        bool io;
        std::cerr << "    surface is "
            << ((io = distance_map_utils::is_orientable (dmap, omap)) ? "" : "not ")
            << "orientable\n";

        std::cerr << "    ";
        distance_map_utils::dump_volume_fractions (std::cerr,
            distance_map_utils::compute_volume_fractions (omap));
        std::cerr << "\n";

        if (!io)
        {
            // check if we were asked to do something where we need an orientable surface.
            if (parallel_series_control != "" || filled_output != "")
                std::cerr << "sorry, I could not orient the surface.  aborting...\n"
                          << util::ABORT;
        }

        std::vector <std::string> fns;
        std::vector <double> vfs;
        if (inflation_series_control != "")
        {
            std::ifstream ctrl (inflation_series_control.c_str ());
            for (;;) 
            {
                std::string fn;
                double vf;
                if (! (ctrl >> vf >> std::ws))
                    break;
                std::getline (ctrl, fn);
                if (!ctrl)
                    break;
                vfs.push_back (vf);
                fns.push_back (fn);
            }
        }
        if (output_file != "")
        {
            vfs.push_back (single_inflate_volfrac);
            fns.push_back (output_file);
        }
        assert (vfs.size () == fns.size ());
        PhaseMap phases (dmap);
        for (size_t i = 0; i != fns.size (); ++i)
        {
            std::cerr << "[main] segmenting: vf = " << vfs[i] << "; fn = "
                      << fns[i] << "\n";
            distance_map_utils::segment_by_quantile (&phases, dmap, vfs[i]);
            int transl[2] = { 2, 1 };
            distance_map_utils::translate_phase_indices (&phases, transl);
            distance_map_utils::write_phasemap_to_bin_file (fns[i], phases);
        }

        fns.clear ();
        vfs.clear ();

        // make the signed distance map
        ITERATE_OVER_VOXEL_FIELD (i,j,k, g)
            dmap[i][j][k] *= omap[i][j][k];

        if (filled_output != "")
        {
            std::cerr << "[main] producing filled; fn = " << filled_output << "\n";
            distance_map_utils::segment_by_threshold (&phases, dmap, 0);
            int transl[2] = { 2, 1 };
            distance_map_utils::translate_phase_indices (&phases, transl);
            distance_map_utils::write_phasemap_to_bin_file (filled_output, phases);
            filled_output = "";
        }

        if (parallel_series_control != "")
        {
            std::ifstream ctrl (parallel_series_control.c_str ());
            for (;;) 
            {
                std::string fn;
                double vf;
                if (! (ctrl >> vf >> std::ws))
                    break;
                std::getline (ctrl, fn);
                if (!ctrl)
                    break;
                vfs.push_back (vf);
                fns.push_back (fn);
            }
        }
        if (filled_output != "")
        {
            std::cerr << "filling is currently broken.\n";
            std::abort ();
            vfs.push_back (0.);
            fns.push_back (filled_output);
        }

        assert (vfs.size () == fns.size ());
        for (size_t i = 0; i != fns.size (); ++i)
        {
            if (!io)
                std::cerr << "internal inconsistency: surface is not oriented\n"
                          << util::ABORT;
            std::cerr << "[main] segmenting parallel body: vf = " << vfs[i] << "; fn = "
                      << fns[i] << "\n";
            distance_map_utils::segment_by_quantile (&phases, dmap, vfs[i]);
            int transl[2] = { 2, 1 };
            distance_map_utils::translate_phase_indices (&phases, transl);
            distance_map_utils::write_phasemap_to_bin_file (fns[i], phases);
        }

        if (quantile_filename != "")
        {
            // finally, make the distance quantiles.
            // strip the sign again from the distance map
            ITERATE_OVER_VOXEL_FIELD(i,j,k,dmap)
                dmap[i][j][k] = fabs (dmap[i][j][k]);
            // linearize the array
            size_t num_voxels = dmap.num_voxels ();
            boost::array <int, 3> new_shape = {{ num_voxels, 1, 1 }};
            dmap.reshape (new_shape);
            // sort it
            std::sort (&dmap[0][0][0], num_voxels + &dmap[0][0][0]);
            // read off quantiles
            std::ofstream of (quantile_filename.c_str ());
            size_t is = num_voxels / 100;
            double voxvol = dmap.voxel_volume ();
            for (size_t i = is; i < num_voxels; i += is)
            {
                double de = dmap[i][0][0] - dmap[i-is][0][0];
                double sf_area = is * voxvol / de;
                of << std::setw (15) << i * voxvol
                   << std::setw (15) << dmap[i][0][0]
                   << std::setw (15) << sf_area << "\n";
            }

            if (!of)
                std::cerr << "Error writing " << quantile_filename << "\n";
            of.close ();
            if (!of)
                std::cerr << "Error writing " << quantile_filename << "\n";
        }
    }
    else
    {
        std::cerr << "Invalid mode " << operation_mode << "\n";
        std::abort ();
    }

    std::cerr << "[DEBUG] average tri's considered per voxel: "
              << double (triangulation_search_structure::total_triangles_tested)
                    / triangulation_search_structure::total_lookups << "\n";
    std::cerr << "[DEBUG] max. tri's in one lookup: "
              << triangulation_search_structure::max_triangles_in_one_lookup
              << "\n";
    return 0;
}
