#ifndef LOGGING_HH
#define LOGGING_HH

#include <sys/resource.h> // getrusage
#include <sys/time.h> // gettimeofday
#include <string.h> // strncmp
#include <iostream>
#include <chrono>
#include <thread>
#include <atomic>

class logging {

private:

    static long long int mem_usage_kb(){ // in kB
        FILE* file = fopen("/proc/self/status", "r");
        if (file) {
            long long int result = -1;
            char line[128];

            while (fgets(line, 128, file) != NULL){
                if (strncmp(line, "VmRSS:", 6) == 0){
                    // parse line :
                    int i = strlen(line);
                    const char* p = line;
                    while (*p <'0' || *p > '9') p++;
                    line[i-3] = '\0'; // assumes line ends in " kB"
                    result = atoll(p);
                    break;
                }
            }
            fclose(file);
            return result;
        } else {
            struct rusage usage;
            getrusage(RUSAGE_SELF, &usage);
            return usage.ru_maxrss / 1000;
        }
    }

    static double time_s() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return tv.tv_sec + tv.tv_usec * 1e-6;
    }
    
    std::atomic<double> t_now;
    std::atomic<long long int> mem_now;
    double t_end;
    std::thread update_now;
    
    double t_init, t_last, t_prog;
    std::string prefix;
    
public:
    logging(std::string pref = "") : prefix(pref) {
        t_now = time_s();
        mem_now = mem_usage_kb();
        t_end = t_now + 50 * 365 * 24 * 3600;
        update_now = std::thread([this](){
                std::chrono::milliseconds ms10(10);
                while (t_now.load(std::memory_order_acquire) < t_end) {
                    std::this_thread::sleep_for(ms10);
                    t_now.store(time_s(), std::memory_order_release);
                    mem_now.store(mem_usage_kb(), std::memory_order_release);
                }
            });
        
        t_init = t_now;
        t_last = 0.;
        t_prog = 1.;
    }

    ~logging() {
        t_end = t_init;
        update_now.join();
    }

    double lap() {
        t_now = time_s();
        return t_now;
    }

    bool progress(float fact = 1.0) {
        if (t_now.load(std::memory_order_acquire) >= t_last + fact*t_prog) {
            t_last = t_now.load(std::memory_order_acquire);
            t_prog *= 1.2;
            return true;
        }
        return false;
    }

    std::ostream & cerr(double t_lap = 0.) {
        bool incr = true;
        if (t_lap == 0.) { t_lap = t_init; incr = false; }
        std::cerr << prefix << (prefix.size() > 0 ? " " : "")
                  << (incr ? "+" : "")
                  << (t_now.load(std::memory_order_acquire) - t_lap) <<"s "
                  << mem_now.load(std::memory_order_acquire) / 1000 << "m "
                  <<": ";
        return std::cerr;
    }

};

#endif // LOGGING_HH
