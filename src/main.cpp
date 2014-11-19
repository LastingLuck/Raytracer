#include "EasyBMP.h"
#include "parse.h"
#include "raytrace.h"

/**
 * Compile command:
 * $ g++ -g -o raytrace main.cpp parse.cpp raytrace.cpp pixel.cpp image.cpp EasyBMP.cpp
 * 
 * can also just call make
 */

SceneData init() {
    SceneData d;
    //Camera
    d.camera.position = (Vector){0, 0, 0};
    d.camera.direction = (Vector){0, 0, 1};
    d.camera.cameraUp = (Vector){0, 1, 0};
    d.camera.cameraHa = 45.0 * (M_PI / 180.0);
    //Vertex/Normal
    d.vertexes.reserve(0);
    d.normals.reserve(0);
    //Image
    d.file.filename = "raytraced.bmp";
    d.file.width = 640;
    d.file.height = 480;
    //Light
    d.ambient = (Light){(Vector){0, 0, 0}, (Vector){0, 0, 0}, (Vector){0, 0, 0}, 0, 0};\
    //Misc
    d.BGColor = (Vector){0, 0, 0};
    d.depth = 5;
    d.eyeray = PERSP;
    d.objNum = 0;
    d.bvhthresh = 500;
    d.bvhdepth = 5;
    d.bvh = 0;
    d.sampleNum = 1;
    #ifdef DEBUG
    printf("Setup Done\n");
    #endif
    return d;
}

int main( int argc, char* argv[] ){
	Image *img = NULL;
    SceneData data = init();
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
    if (parseScene(scenename, &data) < 0) {
        printf("Scene Parse Failed\n");
        std::exit(-1);
    }
    #ifdef DEBUG
    printf("Printing\n");
    printScene(&data);
    #endif
    
    //Ray trace
    img = RayTrace(&data);
    //Write image
    std::string out = data.file.filename.substr(0, data.file.filename.find(".bmp")+4);
    int len = out.length();
    char* output = new char[len+1];
    strncpy(output, out.c_str(), len);
    output[len] = 0;
    img->Write(output);
    did_output = true;
    #ifdef DEBUG
    printf("Outputed image to file: %s\n", output);
    #endif
    
	if (!did_output) {
		fprintf( stderr, "Did not output image\n" );
	}
    
	delete img;
    
	return 0;
}
