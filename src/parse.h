#ifndef PARSE_H_
#define PARSE_H_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <array>

enum ProjType {
    ORTHO,
    PERSP
};

typedef struct f3d {
	float x, y, z;
} Vector;

typedef struct material {
    Vector ambientColor;
    Vector diffuseColor;
    Vector specularColor;
    double cosPow; //Phong cosine power for specular highlights
    Vector transmissiveColor;
    double ior;
} Material;

typedef struct sphere {
    Vector position;
    double radius;
    Material mat;
} Sphere;

typedef struct light {
    Vector color;        //Color of the light (All)
    Vector position;   //Position of the Light (Point, Spot)
    Vector direction;  //Direction the light is pointing (Directional, Spot)
    float angle1;       //Used for Spot Light
    float angle2;       //Used for Spot Light
} Light;

typedef struct cam {
    Vector position;   //Camera Position
    Vector direction;  //Camera Viewing Direction
    Vector cameraUp;   //Up vector
    double cameraHa;    //Half of the height angle of the viewing frustrum
} Camera;

typedef struct image {
    std::string filename;     //Output filename
    int width;          //Width of image
    int height;         //Height of image
} File;

typedef struct triangle {
	Material mat;
    Vector vertices[3];
    Vector normals[6]; //First 3 are the first normal for each point, last 3 are the second normal
    bool ntri;
} Triangle;

typedef struct flatplane {
	Material mat;
	Vector point; //A point on the plane
	Vector normal;
	Vector inormal; //Inverse of normal
} Plane;

typedef struct flatrectangle {
	Material mat;
	Vector points[4];
	Vector normal;
	Vector inormal;
} Rectangle;

typedef struct aabb {
	//std::array<Vector, 8> vertices; //(TL TR BL[min] BR) (TL TR[max] BL BR)
	//std::array<Vector, 6> normals; //F T R Bo L Ba
	Vector max;
	Vector min;
	//Vector max; //Contains the max x, y, z values of the bounding box
	//Vector min; //Contains the min x, y, z values of the bounding box
	struct aabb* sub1;
	struct aabb* sub2;
	struct aabb* parent;
	
	int objNum;
	std::vector<Sphere> spheres;
	std::vector<Triangle> triangles;
    std::vector<Plane> planes;
    //std::vector<Rectangle> rectangles;
} Box;

typedef struct data {
    Camera camera;   
    std::vector<Sphere> spheres; //List of spheres to put in the scene
    std::vector<Vector> vertexes;
    std::vector<Vector> normals;
    std::vector<Triangle> triangles;
    std::vector<Plane> planes;
    File file;
    Vector BGColor;      //RGB value of the background
    std::vector<Light> directional;
    std::vector<Light> point;
    std::vector<Light> spot;
    Light ambient;
    int depth;
    enum ProjType eyeray; //Orthographic or Perspective
    
    int bvhthresh;
    int bvhdepth;
    int objNum;
    Box* bvh;
    
    int sampleNum;
} SceneData;

//Parses the scene text file and fills out a SceneData struct
int parseScene(char* file, SceneData *scene);

//Print functions used to print the contents of SceneData to test it
void printScene(SceneData *scn);
void printCamera(Camera camera);
void printSphere(Sphere sph);
void printTriangle(Triangle tri);
void printPlane(Plane pl);
void printRect(Rectangle rec);
void printMaterial(Material mat);
void printLight(Light l);
void printImage(File img);
void printVector(Vector fVec);
/*
class Vector {
    public:
        Vector();
        Vector(float xval, float yval, float zval);
        
        float x;
        float y;
        float z;
        
        //Performs a dot product operation on the 2 vectors
        //float dot(Vector u, Vector v);
        float operator*(const Vector&);
        //Performs a cross prooduct
        Vector cross(Vector u, Vector v);
        //Multiple vector
        //Vector multiply(Vector u, float c);
        Vector operator*(const float& c);
        //Adds 2 vectors
        //Vector add(Vector u, Vector v);
        Vector operator+(const Vector& u);
        //Vector add(Vector u, float c);
        Vector operator+(const float& c);
        //Subtracts 2 vectors
        //Vector sub(Vector u, Vector v);
        Vector operator-(const Vector& u);
        Vector operator-(const float& c);
        //Returns the length of the vector
        float length();
        //Finds the length but does not take the square root
        float lengthsq();
        //Normalizes the vector
        Vector norm();
};

class Material {
    public:
        Material();
        Material();
};

class Sphere {
    
};

class Triangle {
    public:
        Triangle();
        Triangle(Vector v1, Vector v2, Vector v3);
        Triangle(Vector v1, Vector v2, Vector v3, Vector n1, Vector n2, Vector n3);
        
        Vector vertices[3];
        Vector normals[6];
        bool ntri;
        
        Vector getNormal();
};

class Light {
    
};

class Camera {
    public:
        Camera();
        Camera(Vector p, Vector d, Vector u, double ha);
    
        Vector position;
        Vector direction;
        Vector up;
        double halfAngle;
};

class Scene {
    public:
        //Camera
        Camera camera;
        enum ProjType proj;
        //Image
        std::string filename;
        int width, height;
        //Sphere
        Sphere* spheres;
        int sphereNum;
        //Triangle
        Vector* vertexPool;
        int vertNum;
        Vector* normalPool;
        int normNum;
        Triangle* triangles;
        int triNum;
        //Background
        Vector bgColor;
        //depth
        int maxDepth;
        //Lights
        Light* directional;
        Light* point;
        Light* spot;
        Light ambient;
        int lightNum[3];
        
        double getPlaneDist();
};
*/
#endif
