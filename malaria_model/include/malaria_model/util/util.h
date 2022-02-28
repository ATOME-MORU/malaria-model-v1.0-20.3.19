#ifndef UTIL_H
#define UTIL_H


#define DEBUG

#ifdef DEBUG
#define DEBUG_MSG(str) do { std::cout << str << std::endl; } while( false )
#else
#define DEBUG_MSG(str) do { } while ( false )
#endif



#include <string>
#include <iomanip> //std::put_time
#include <ctime>
#include <sstream> // std::ostringstream

#include "rapidjson/fwd.h"

namespace util {

const char kPathSeparator =
#ifdef _WIN32
                            '\\';
#else
                            '/';
#endif


class InputParser {
    
    private:
        std::vector <std::string> tokens;

    public:
        InputParser(int &argc, char **argv);
        
        const std::string& getCmdOption(const std::string &option) const;
        
        bool cmdOptionExists(const std::string &option) const;

};

// template <class T>
// class Printables {
//     public:

//     Printables<T>();
    
//     void print_data_table(const T* DataManager, int total, int num_rows, int num_columns_in_full);
        
// };

template <class T>
void print_data_table(const T* DataManager, int total, int num_rows, int num_columns_in_full);


std::string get_output_prefix();

inline bool checkFileExists(const std::string& name);
inline bool createDirectory(const std::string& directory_name);

rapidjson::Document get_json_from_file(const std::string);
bool validate_json_against_schema(const rapidjson::Document* doc_schema, const rapidjson::Document* doc_input);


template<typename T>
class Enum_iterator {
public:
    class Iterator{
    public:
        Iterator(int value) :
            m_value(value)
            {}

        T operator*(void) const {
           return (T)m_value;
        }

        void operator++(void) {
           ++m_value;
        }

        bool operator!=(Iterator rhs) {
           return m_value != rhs.m_value;
        }
        
    private:
        int m_value;
    };
};
template<typename T>
typename Enum_iterator<T>::Iterator begin(Enum_iterator<T>) {
   return typename Enum_iterator<T>::Iterator((int)T::First);
}
template<typename T>
typename Enum_iterator<T>::Iterator end(Enum_iterator<T>) {
   return typename Enum_iterator<T>::Iterator(((int)T::Last) + 1);
}



//https://stackoverflow.com/questions/17432502/how-can-i-measure-cpu-time-and-wall-clock-time-on-both-linux-windows
//  Windows
// #ifdef _WIN32
// #include <Windows.h>
// double get_wall_time(){
//     LARGE_INTEGER time,freq;
//     if (!QueryPerformanceFrequency(&freq)){
//         //  Handle error
//         return 0;
//     }
//     if (!QueryPerformanceCounter(&time)){
//         //  Handle error
//         return 0;
//     }
//     return (double)time.QuadPart / freq.QuadPart;
// }
// double get_cpu_time(){
//     FILETIME a,b,c,d;
//     if (GetProcessTimes(GetCurrentProcess(),&a,&b,&c,&d) != 0){
//         //  Returns total user time.
//         //  Can be tweaked to include kernel times as well.
//         return
//             (double)(d.dwLowDateTime |
//             ((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
//     }else{
//         //  Handle error
//         return 0;
//     }
// }

// //  Posix/Linux
// #else
// #include <time.h>
// #include <sys/time.h>
double get_wall_time();
// double get_cpu_time(){
//     return (double)clock() / CLOCKS_PER_SEC;
// }
// #endif

double get_cpu_time();

}

#endif