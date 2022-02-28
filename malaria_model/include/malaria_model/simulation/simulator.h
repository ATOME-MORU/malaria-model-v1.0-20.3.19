#ifndef SIMULATOR_H
#define SIMULATOR_H

namespace simulation{

// const std::string kJson_schema_config_file_name = "schema.config.json"; // to remove

class Simulator {

public:

    const rapidjson::Document config;
    // const bool config_checked;
    const std::string output_directory;
    const std::string output_prefix;

    const int num_time_steps = 0;
    const bool kIf_batch_mode;

    village::VillageManager vll_mgr;
    intervention::Mda mda;
    intervention::Fmda fmda;
    intervention::Tmda tmda;
    intervention::Vc vc;
    human::HumanManager hmn_mgr;
    human::BloodSystemManager bsys_mgr;
    simulation::ReportManager rpt_mgr;

    village::MosquitoManager* msq_mgr;

    village::MosquitoManagerIbmRa* msq_mgr_ibm;

    //////////////////////////////////////////////////////////////////////
    //// Initialisation

    Simulator(
        const std::string& configuration_file_name,
        const std::string output_directory_name,
        const std::string& output_prefix_name
    );

    ~Simulator();

    static bool validate_config(const std::string& configuration_file_name);
    
    //// - Villages
    //// - MDA
    //// - Humans
    //// - Bsys
    //// - Reporter
    void init();

    //// Initialisation - Village Domain
    //// The followings are initialised
    //// - action/property (config dependencies)
    //// - Data Space Allocation
    //// - Read File (village.data_file)
    ////    0 Pop
    ////    1 Longitude
    ////    2 Latitude
    ////    3 Transmission Coefficient
    ////    4 Number of Malaria Posts
    ////    5 Opening day (time step) of Malaria Post
    ////    6 Migrant Percentage
    ////    7 Artemisinin Resistance Percentage
    //// - Village Distances
    //// - Static Mobility Network
    //// - dev. Treatment Rate (treatment)
    //// - Mosquito Seasonality (mosquito.seasonality.data_file)

    //// Initialisation - Human Domain
    //// The followings are initialised
    //// - action/property (config dependencies)
    //// - Data Space Allocation
    //// - Age (human.age_weight_vector_file)
    //// - Gender (human.male_percentage)
    //// - ITN (human.itn_probability)

    //// Initialisation - Blood Domain
    //// The followings are initialised
    //// - action/property (config dependencies)
    //// - Data Space Allocation


    //////////////////////////////////////////////////////////////////////
    //// Execution & Finish

    void run();
    void close();

    //////////////////////////////////////////////////////////////////////
    //// Output

    void output_mda();

};

}
#endif
