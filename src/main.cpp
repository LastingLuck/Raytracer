#include "EasyBMP.h"
#include "parse.h"
#include "raytrace.h"

/**
 * Compile command:
 * $ g++ -g -o raytrace main.cpp parse.cpp raytrace.cpp pixel.cpp image.cpp EasyBMP.cpp
 * 
 * can also just call make
 */

int main( int argc, char* argv[] ){
	Image *img = NULL;
    Scene* data = new Scene();
	bool did_output = false;
    char* scenename;
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
    if (data->parseScene(scenename) < 0) {
        printf("Scene Parse Failed\n");
        std::exit(-1);
    }
    #ifdef DEBUG
    printf("Printing\n");
    printScene(data);
    #endif
    
    //Ray trace
    RayTrace t = RayTrace();
    img = t.rayTrace(*data);
    //Write image
    std::string out = data->getImage().getFileName().substr(0, data->getImage().getFileName().find(".bmp")+4);
    int len = out.length();
    char* output = new char[len+1];
    strncpy(output, out.c_str(), len);
    output[len] = 0;
    img->Write(output);
    did_output = true;
    delete[] output;
    #ifdef DEBUG
    printf("Outputed image to file: %s\n", output);
    #endif
    
	if (!did_output) {
		fprintf( stderr, "Did not output image\n" );
	}
    
	delete img;
    delete data;
    
	return 0;
}
