#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include "pcap2kitti.h"

int main( int argc, char* argv[] )
{
    // Open pcap2kitti that retrieve from Sensor
/*    const boost::asio::ip::address address = boost::asio::ip::address::from_string( "192.168.1.21" );
    const unsigned short port = 2368;
    velodyne::VLP16Capture capture( address, port );*/
    //velodyne::HDL32ECapture capture( address, port );

    // Open pcap2kitti that retrieve from PCAP

    std::string file, path, filename, pcapFile;

    // Retrieve the (non-option) argument:
    if ( (argc <= 1) || (argv[argc-1] == NULL) )
        std::cerr << "Please provide input file...(example: /../../Velodyne.pcap)" << std::endl;
    else pcapFile = argv[argc - 1];

    std::vector<std::string> filevec;
    std::istringstream ss(pcapFile);
    while(getline(ss, file, '/' )){
        filevec.push_back(file);
    }
    for(int i = 0; i < filevec.size() - 1; ++i){
        path = path + filevec[i] + "/";
    }
    file = filevec.back();

    auto pos = file.find(".pcap");
    if(pos != std::string::npos)
        filename = file.substr(0, pos);
    std::cout << "filename " << filename << std::endl;


    //velodyne::VLP16Capture capture( pcapFile );
    velodyne::HDL32ECapture capture( pcapFile );

    if( !capture.isOpen() ){
        std::cerr << "Can't open pcap2kitti." << std::endl;
        return -1;
    }

    int index = 0;
    while( capture.isRun() ){
        std::vector<velodyne::Laser> lasers; // Capture One Rotation Data
        capture >> lasers;
        if( lasers.empty() ){
            continue;
        }

        std::ostringstream dirname;
        dirname << std::setw(10) << std::setfill('0') << std::to_string(index);
        std::string outdir = path + "/" + dirname.str() + ".bin";
        std::fstream outfile(outdir);

        outfile.open(outdir, std::fstream::out | std::fstream::trunc | std::fstream::binary); // opens the file
        if( !outfile ) { // file couldn't be opened
            std::cerr << "Error: file could not be opened" << std::endl;
            exit(1);
        }


        for( const velodyne::Laser& laser : lasers ){ // Access to Laser Data
            std::vector<double> rowlaser;
            const double azimuth = laser.azimuth; // Laser Azimuth ( degree )
            const double vertical = laser.vertical; // Laser Vertical ( degree )
            const unsigned short distance = laser.distance; // Laser Distance ( centimeter )
            const unsigned int intensity = static_cast<unsigned int>( laser.intensity ); // Laser Intensity
            const unsigned int id = static_cast<unsigned int>( laser.id ); // Laser ID ( VLP-16 : 0 - 15, HDL-32E : 0 - 31 )
            float x = static_cast<float>( ( distance * std::cos( vertical ) ) * std::sin( azimuth ) );
            float y = static_cast<float>( ( distance * std::cos( vertical ) ) * std::cos( azimuth ) );
            float z = static_cast<float>( ( distance * std::sin( vertical ) ) );
            const long long timestamp = laser.time; // Laser TimeStamp ( microseconds )

            if(x != 0.0f && y != 0.0f && z != 0.0f ) {
                rowlaser.push_back(x);
                rowlaser.push_back(y);
                rowlaser.push_back(z);
                rowlaser.push_back(intensity);
                rowlaser.push_back(timestamp);
                rowlaser.push_back(id);
                for (const auto &e : rowlaser) outfile << e << "\t";
                outfile << "\n";
            }

            // Laser TimeStamp ( MM/DD/YY hh:mm:ss.ffffff )
/*            std::time_t time = static_cast<int>( laser.time / 1000000 );
            const tm* localtime = std::localtime( &time );
            std::cout << "localtime : " << std::put_time( localtime, "%c" ) << "." << static_cast<int>( laser.time % 1000000 ) << "\n";*/

        }

        outfile.close();
        dirname.clear();
        ++ index;
    }
    return 0;
}
