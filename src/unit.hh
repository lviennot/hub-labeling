#pragma once

#define CHECK(x)                                                            \
    do { if (!(x)) {                                                        \
            std::cerr <<"CHECK failed: "<< #x <<"\n"                        \
                      <<" at: " << __FILE__ <<":" << __LINE__ << "\n"       \
                      << " in function: " << __func__ << "\n"               \
                      << std::flush;                                        \
            std::abort();                                                   \
        } } while (0)
