#include <algorithm> //std::find
#include <iostream>
// #include <vector> // not if -std=c++11

#include <iterator>

#include <math.h> //ceil
#include <algorithm> //max

#ifdef WIN32
#include <direct.h> // not tested on WIN32
#else
#include <sys/stat.h> //stat in checkFileExists, mkdir in createDirectory
#endif

#include <fstream>


#include <time.h>
#include <sys/time.h>

#include <vector>

#include "rapidjson/document.h"
#include "rapidjson/error/en.h"
#include "rapidjson/schema.h"

#include "rapidjson/stringbuffer.h"
#include "rapidjson/prettywriter.h"


#include "util/util.h"
#include "human/human.h"
#include "human/blood_system.h"
#include "village/village.h"

namespace util {


InputParser::InputParser(int &argc, char **argv) {

    for (int i=1; i < argc; i++){
        this->tokens.push_back(std::string(argv[i]));
    }

}

const std::string& InputParser::getCmdOption(const std::string &option) const {

    std::vector<std::string>::const_iterator itr;
    itr =  std::find(this->tokens.begin(), this->tokens.end(), option);

    if (itr != this->tokens.end() && ++itr != this->tokens.end()){
        return *itr;
    }

    static const std::string empty_string("");
    return empty_string;

}

bool InputParser::cmdOptionExists(const std::string &option) const {

    return std::find(this->tokens.begin(), this->tokens.end(), option)
           != this->tokens.end();

}

rapidjson::Document get_json_from_file(const std::string file_name) {

    std::ifstream ifs(file_name);
    if (!ifs) {
        std::cerr << "Error reading file (" << file_name << ")" << std::endl;
        exit(EXIT_FAILURE);
    }

    rapidjson::Document rj_doc;
    std::string file_contents( (std::istreambuf_iterator<char>(ifs)),
                                std::istreambuf_iterator<char>()
                             );
    ifs.close();

    rj_doc.Parse(file_contents.c_str());
    if (rj_doc.HasParseError()){
        std::cerr << "File " << file_name << " is not a valid JSON\n";
        std::cerr << "Error at offset " << static_cast<unsigned>(rj_doc.GetErrorOffset());
        std::cerr << ": " << rapidjson::GetParseError_En(rj_doc.GetParseError()) << std::endl;
        exit(EXIT_FAILURE);
    }

    return rj_doc;

}

bool validate_json_against_schema(const rapidjson::Document* doc_schema, const rapidjson::Document* doc_input) {
    rapidjson::SchemaDocument rj_schema_doc(*doc_schema);
    rapidjson::SchemaValidator rj_validator(rj_schema_doc);

    if(!doc_input->Accept(rj_validator)){

        rapidjson::StringBuffer rj_sb;
        rj_validator.GetInvalidSchemaPointer().StringifyUriFragment(rj_sb);
        printf("Invalid schema: %s\n", rj_sb.GetString());
        printf("Invalid keyword: %s\n", rj_validator.GetInvalidSchemaKeyword());
        rj_sb.Clear();
        rj_validator.GetInvalidDocumentPointer().StringifyUriFragment(rj_sb);
        printf("Invalid document: %s\n", rj_sb.GetString());

        rj_sb.Clear();
        rapidjson::PrettyWriter<rapidjson::StringBuffer> w(rj_sb);
        rj_validator.GetError().Accept(w);
        fprintf(stderr, "Error report:\n%s\n", rj_sb.GetString());

        // exit(EXIT_FAILURE);
        return false;

    }
    return true;
}

// template <class T>
// Printables<T>::Printables(){}

// template <class T>
// void Printables<T>::print_data_table(const T* DataManager, int total, int num_rows, int num_columns_in_full) {
//     int c;
//     int c_max = total/num_rows;
//     for (int dd = 0; dd < total; dd++){
//         c = dd%c_max;
//         if (c < num_columns_in_full) {
//             std::cout << "|";
//             DataManager->print_one(dd);
//         } else if (dd == total - 1) {
//             std::cout << "| .. |";
//             DataManager->print_one(dd);
//             std::cout << "|" << std::endl;
//         } else if (c == c_max-1) {
//             std::cout << "| .. | .. |";
//             DataManager->print_one(dd);
//             std::cout << "|" << std::endl;
//         }
//     }
// }

// template class Printables<human::HumanManager>;

//TODO: add function to print random samples
template <class T>
void print_data_table(const T* DataManager, int total, int num_rows, int num_columns_in_full) {
    int items_per_row = ceil(static_cast<float>(total)/num_rows);

    num_columns_in_full = std::min(num_columns_in_full, items_per_row);

    for (int rr = 0; rr < num_rows; rr++) {
        for (int cc = 0; cc < num_columns_in_full; cc++) {
            std::cout << "|";
            // int to_print = std::min(rr*items_per_row+cc, total);
            DataManager->print_one(std::min(rr*items_per_row+cc, total-1));
        }
        std::cout << "\n";
        // if (rr == num_rows-1) {
        //     std::cout << "| .. |";
        //     DataManager->print_one(total-1);
        //     std::cout << "|" << std::endl;
        // } else {
        //     std::cout << "| .. | .. |";
        //     DataManager->print_one((rr+1)*items_per_row-1);
        //     std::cout << "|" << std::endl;
        // }
    }
}

template void print_data_table<village::VillageManager>(const village::VillageManager* DataManager, int total, int num_rows, int num_columns_in_full);
template void print_data_table<human::HumanManager>(const human::HumanManager* DataManager, int total, int num_rows, int num_columns_in_full);
template void print_data_table<human::BloodSystemManager>(const human::BloodSystemManager* DataManager, int total, int num_rows, int num_columns_in_full);

std::string get_output_prefix(){
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y%m%d%H%M%S");
    return oss.str();
}

// Non-library filesystem functions

inline bool checkFileExists(const std::string& name) { //file_name better

    struct stat buffer; //#include <sys/stat.h>
    return (stat (name.c_str(), &buffer) == 0);

}

inline bool createDirectory(const std::string& directory_name) {
#ifdef WIN32
    if (_mkdir(directory_name.c_str())==-1){ // not tested this
#else
    if (mkdir(directory_name.c_str(), 0755)==-1){
#endif
        if( errno == EEXIST ) {
            //directory already exist which is OK
        } else {
            std::cout << "error while creating output directory ("<< directory_name <<")"
            // << strerror_s(errno) << std::endl;
            << errno << std::endl;
            return false;
        }
    }
    return true;
}





double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}


} //namespace util