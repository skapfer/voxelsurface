#ifndef Z2bbd892249eace8fbd7f0f61729e703cab3bc8e9 
#define Z2bbd892249eace8fbd7f0f61729e703cab3bc8e9 
#include <limits>
#include <algorithm>
#include <cstdlib>
#include <ostream>

namespace util {

    template <typename TYPE>
    struct FieldWidth {
        FieldWidth (int width, const TYPE &val)
            : w_ (width), val_ (val)
        {
        }

        friend std::ostream &operator<< (std::ostream &os, const FieldWidth <TYPE> &v)
        {
            os.width (v.w_);
            return os << v.val_ << ' ';
        }

    private:
        int w_;
        const TYPE &val_;
    };

    template <typename TYPE>
    const FieldWidth <TYPE> field_width (int width, const TYPE &val)
    {
        return FieldWidth <TYPE> (width, val);
    }

    inline static
    std::ostream &ABORT (std::ostream &os) {
        os << std::endl;
        std::abort ();
        return os;
    }

    inline static
    int sgn (double x)
    {
        if (x == 0)
            return 0;
        else if (x > 0)
            return 1;
        else if (x < 0)
            return -1;
        else
            return -42;
    }

    template <typename TYPE>
    inline
    TYPE sq (const TYPE &x)
    {
        return x*x;
    }

    // NaN-aware max functions
    // normally, return the max of two arguments;
    // if either of a,b is NaN, return that.
    template <typename TYPE>
    inline
    const TYPE &max2 (const TYPE &a, const TYPE &b)
    {
        if (a > b)
            return a;
        if (b > a)
            return b;
        // you would think they should be equal here.
        // BUT one of them could be NaN.
        // find the NaN, and return that.
        if (! (a==a))
            return a;
        else
            return b;
    }

    template <typename TYPE>
    inline
    const TYPE &max3 (const TYPE &a, const TYPE &b, const TYPE &c)
    {
        return max2 (a, max2 (b, c));
    }

    template <typename TYPE>
    inline
    const TYPE &max4 (const TYPE &a, const TYPE &b,
                      const TYPE &c, const TYPE &d)
    {
        return max2 (max2 (a, b), max2 (c, d));
    }

    template <typename TYPE>
    class SummerType {
    public:
        SummerType (const TYPE &initvalue) : val_ (initvalue) {}
        operator const TYPE & () const { return val_; }
        void operator() (const TYPE &a) { val_ += a; }
    private:
        TYPE val_;
    };

    // helper function to avoid having to specify template argument to 
    // constructor
    template <typename TYPE>
    SummerType <TYPE> summer (const TYPE &initvalue)
    {
        return SummerType <TYPE> (initvalue);
    }

    template <typename TYPE>
    class MinMax {
    public:
        MinMax () {
            min_ = std::numeric_limits <TYPE>::max ();
            max_ = std::numeric_limits <TYPE>::min ();
        }

        const TYPE &min () { return min_; }
        const TYPE &max () { return max_; }

        void operator() (const TYPE &a) {
            min_ = std::min (a, min_);
            max_ = std::max (a, max_);
        }

    private:
        TYPE min_, max_;
    };

    template <typename TYPE, typename COUNTER_TYPE = unsigned long>
    struct PercentageBelow {
    public:
        PercentageBelow (TYPE thresh) : thresh_ (thresh), below_ (0u), totalcount_ (0u) {}

        double value () const {
            return double (below_) / totalcount_;
        }

        void operator() (TYPE a) {
            ++totalcount_;
            if (a < thresh_)
                ++below_;
        }

    private:
        TYPE thresh_;
        COUNTER_TYPE below_, totalcount_;
    };
}
 
#endif /* Z2bbd892249eace8fbd7f0f61729e703cab3bc8e9 */
