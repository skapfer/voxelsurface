#ifndef Ce7faa0f942e72552b012e8376911a719806028bb 
#define Ce7faa0f942e72552b012e8376911a719806028bb 

#include <iostream>

class ProgressBar
{
public:
    ProgressBar (size_t job_size)
        : job_size_ (job_size), done_ (0u), next_step_ (1)
    {
        std::cerr << "0% ";
    }

    void operator++ ()
    {
        *this += 1u;
    }

    void operator+= (size_t increment)
    {
        done_ += increment;
        if (double (done_) / job_size_ >= next_step_ / 10.)
        {
            std::cerr << next_step_++ << "0% ";
            if (next_step_ == 11)
                std::cerr << "\n";
        }
    }

private:
    size_t job_size_, done_;
    int next_step_;
};
    
#endif /* Ce7faa0f942e72552b012e8376911a719806028bb */
