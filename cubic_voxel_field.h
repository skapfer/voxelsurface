#ifndef CUBIC_VOXEL_FIELD_INCLUDED 
#define CUBIC_VOXEL_FIELD_INCLUDED 
#include <stdexcept>    
#include <boost/multi_array.hpp>

#define ITERATE_OVER_VOXEL_FIELD(i,j,k,dmap) \
    for (int i = 0; i != (dmap).nvox[0]; ++i) \
    for (int j = 0; j != (dmap).nvox[1]; ++j) \
    for (int k = 0; k != (dmap).nvox[2]; ++k)


namespace cubic_voxel_field {

class GeometryBase {
protected:
    GeometryBase () {
        is_periodic = true;
    }

    static double sq (double x) { return x*x; }

public:
    typedef size_t count_type;
    int nvox[3];
    bool is_periodic;

    void reduce_voxel_coordinates (int *i, int *j, int *k) const
    {
        assert (is_periodic);
        *i += nvox[0]; *i %= nvox[0];
        *j += nvox[1]; *j %= nvox[1];
        *k += nvox[2]; *k %= nvox[2];
    }

    count_type num_voxels () const
    {
        count_type ret = nvox[0] * nvox[1] * nvox[2];
        assert (ret > 0);
        return ret;
    }

    friend std::ostream &operator<< (std::ostream &os, const GeometryBase &g)
    {
        os << "(nvox = " << g.nvox[0] << ", " << g.nvox[1] << ", "
           << g.nvox[2] << (g.is_periodic ? "; periodic)" : "; aperiodic)");
        return os;
    }

    class VoxelNeighbors;

    class VoxelNeighborsKey {
    public:
        template <typename FIELD>
        typename FIELD::value_type &operator[] (FIELD &f) const
        {
            return f[is_[ci_]][js_[cj_]][ks_[ck_]];
        }

        template <typename FIELD>
        typename FIELD::value_type const &operator[] (const FIELD &f) const
        {
            return f[is_[ci_]][js_[cj_]][ks_[ck_]];
        }

        VoxelNeighborsKey &operator++ ()
        {
            if (++ck_ == nk_)
            {
                ck_ = 0;
                if (++cj_ == nj_)
                {
                    cj_ = 0;
                    ++ci_;
                }
            }
            return *this;
        }

        int i () const { return is_[ci_]; }
        int j () const { return js_[cj_]; }
        int k () const { return ks_[ck_]; }

        friend bool operator!= (const VoxelNeighborsKey &a, const VoxelNeighborsKey &b)
        {
            return a.ci_ != b.ci_ || a.cj_ != b.cj_ || a.ck_ != b.ck_;
        }

    private:
        // construct the 'end' key.
        VoxelNeighborsKey (const GeometryBase &g, int i, int j, int k)
        {
            int ni_;
            if (g.is_periodic) {
                int mins[3]  = { i-1, j-1, k-1 };
                int mids[3]  = { i,   j,   k   };
                int maxes[3] = { i+1, j+1, k+1 };
                g.reduce_voxel_coordinates (&mins[0], &mins[1], &mins[2]);
                g.reduce_voxel_coordinates (&maxes[0], &maxes[1], &maxes[2]);
                ni_ = nj_ = nk_ = 3;
                is_[0] = mins[0]; is_[1] = mids[0]; is_[2] = maxes[0];
                js_[0] = mins[1]; js_[1] = mids[1]; js_[2] = maxes[1];
                ks_[0] = mins[2]; ks_[1] = mids[2]; ks_[2] = maxes[2];
            } else {
                ni_ = nj_ = nk_ = 0;
                for (int ii = std::max (0, i-1); ii < std::min (g.nvox[0], i+2); ++ii)
                    is_[ni_++] = ii;
                for (int jj = std::max (0, j-1); jj < std::min (g.nvox[1], j+2); ++jj)
                    js_[nj_++] = jj;
                for (int kk = std::max (0, k-1); kk < std::min (g.nvox[2], k+2); ++kk)
                    ks_[nk_++] = kk;
            }

            cj_ = ck_ = 0;
            ci_ = ni_;
        }

        int is_[3], js_[3], ks_[3];
        int nj_, nk_;
        int ci_, cj_, ck_;

        friend class VoxelNeighbors;
    };

    class VoxelNeighbors {
    public:
        VoxelNeighbors (const GeometryBase &g, int i, int j, int k) :
            end_ (g, i, j, k)
        {
        }

        typedef VoxelNeighborsKey key_type;

        key_type begin () const {
            key_type ret (end_);
            ret.ci_ = 0;
            return ret;
        }

        const key_type &end () const {
            return end_;
        }

    private:
        key_type end_;
    };

    void check_validity () {
        if (nvox[0] <= 0 || nvox[1] <= 0 || nvox[2] <= 0)
            throw std::runtime_error ("Invalid instance of GeometryBase, "
                "voxel number is zero or negative.");
    }
};

template <typename POINT, typename VECTOR = POINT>
class Geometry : public GeometryBase {
public:
    POINT origin, upper_corner;

    POINT voxel_center (int i, int j, int k) const
    {
        // let's hope VECTOR is not a integer datatype
        assert (VECTOR (1, 0, 0) + VECTOR (.5, 0, 0) != VECTOR (1, 0, 0));

        VECTOR diag = upper_corner - origin;
        VECTOR off = VECTOR (
            diag[0] * (.5 + i) / nvox[0],
            diag[1] * (.5 + j) / nvox[1],
            diag[2] * (.5 + k) / nvox[2]);
        return origin + off;
    }

    VECTOR voxel_diagonal () const
    {
        VECTOR diag = upper_corner - origin;
        return VECTOR (diag[0] / nvox[0], diag[1] / nvox[1], diag[2] / nvox[2]);
    }

    double voxel_squared_distance (int di, int dj, int dk) const
    {
        const VECTOR diag = this->voxel_diagonal ();
        return sq (di*diag[0]) + sq (dj*diag[1]) + sq (dk*diag[2]);
    }

    double voxel_volume () const
    {
        const VECTOR diag = this->voxel_diagonal ();
        return (diag[0]) * (diag[1]) * (diag[2]);
    }
};

template <typename GEOMETRY, typename VALUE>
class CubicVoxelField : public GEOMETRY, public boost::multi_array <VALUE, 3> {
public:
    typedef VALUE value_type;
    typedef GEOMETRY geometry_type;

    CubicVoxelField (const geometry_type &geo, const value_type &v = value_type ())
        : geometry_type (geo),
          boost::multi_array <value_type, 3> (
              boost::extents[geo.nvox[0]][geo.nvox[1]][geo.nvox[2]])
    {
        GEOMETRY::check_validity ();
        for (int i = 0; i != geometry_type::nvox[0]; ++i)
        for (int j = 0; j != geometry_type::nvox[1]; ++j)
        for (int k = 0; k != geometry_type::nvox[2]; ++k)
            (*this)[i][j][k] = v;
    }

    const value_type wrapped (int i, int j, int k) const
    {
        geometry_type::reduce_voxel_coordinates (&i, &j, &k);
        return (*this)[i][j][k];
    }
};

template <typename FIELD, typename PRED>
PRED for_each_voxel (FIELD &field, PRED pred)
{
    ITERATE_OVER_VOXEL_FIELD (i,j,k, field)
        pred (field[i][j][k]);
    return pred;
}

template <typename FIELD, typename PRED>
PRED for_each_neighbor (FIELD &field, int i, int j, int k, PRED p)
{
    typename FIELD::VoxelNeighbors n (field, i,j,k);
    typename FIELD::VoxelNeighbors::key_type it = n.begin ();
    for (; it != n.end (); ++it)
        p (it[field]);
    return p;
}

}

#endif // CUBIC_VOXEL_FIELD_INCLUDED
