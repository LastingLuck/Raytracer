#include "EasyBMP.h"
#include "parse.h"
#include "raytrace.h"
#include <chrono>

/**
 * Compile command:
 * $ g++ -g -o raytrace main.cpp parse.cpp raytrace.cpp pixel.cpp image.cpp EasyBMP.cpp
 * 
 * can also just call make
 */

void displayElapsed(double sec);

int main( int argc, char* argv[] ){
	Image *img = NULL;
    Scene* data = new Scene();
	bool did_output = false;
    char* scenename;
    //Variables for timing
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsed;
	// first argument is program name
	argv++, argc--;
    
	// no argument case
	if (argc == 0) {
		printf("No Scene File Specified\n");
        std::exit(-1);
	}
    else if (argc == 1) {
        scenename = argv[0];
    }
    else {
        printf("Too many arguments. Defaulting to first one\n");
        scenename = argv[0];
    }

    //Parse the scene file
    start = std::chrono::system_clock::now();
    if (data->parseScene(scenename) < 0) {
        printf("Scene Parse Failed\n");
        std::exit(-1);
    }
    end = std::chrono::system_clock::now();
    elapsed = end - start;
    printf("Scene parsed in ");
    displayElapsed(elapsed.count());

    #ifdef DEBUG
    printf("Printing\n");
    printScene(data);
    #endif
    
    //Ray trace
    RayTrace t = RayTrace();
    start = std::chrono::system_clock::now();
    img = t.rayTrace(*data);
    end = std::chrono::system_clock::now();
    elapsed = end - start;
    printf("Scene traced in ");
    displayElapsed(elapsed.count());

    //Write image
    std::string out = data->getImage().getFileName().substr(0, data->getImage().getFileName().find(".bmp")+4);
    delete data;
    int len = out.length();
    char* output = new char[len+1];
    strncpy(output, out.c_str(), len);
    output[len] = 0;
    img->Write(output);
    delete img;
    did_output = true;
    #ifdef DEBUG
    printf("Outputed image to file: %s\n", output);
    #endif
    delete[] output;
	if (!did_output) {
		fprintf( stderr, "Did not output image\n" );
	}
    
	return 0;
}

void displayElapsed(double sec) {
	int days, hours, himutes;

	days = sec / 86400.0;
	sec -= days * 86400.0;

	hours = sec / 3600.0;
	sec -= hours * 3600.0;

	minutes = sec / 60.0;
	sec -= minutes * 60.0;

	printf("%dd %dh %dm %fs\n", days, hours, minutes, sec);
}