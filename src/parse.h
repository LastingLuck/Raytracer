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
/*
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
*/
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

class Vector {
    public:
        Vector();
        Vector(float xval, float yval, float zval);
        
        float x;
        float y;
        float z;
        
        //Performs a dot product operation on the 2 vectors
        static float dot(const Vector& u, const Vector& v);
        //Performs a cross prooduct
        static Vector cross(const Vector& u, const Vector& v);
        //Finds the length of the distance between the 2 vectors (same as magnitude(u-v))
        static float length(const Vector& u, const Vector& v);
        //Square version of length()
        static float lengthSq(const Vector& u, const Vector& v);
        
        //Returns the magnitude of the vector
        float magnitude();
        //Finds the magnitude but does not take the square root
        float magnitudeSq();
        //Normalizes the vector
        Vector norm();
        
        //Multiply vector
        Vector operator*(const float& c);
        Vector& operator*=(const Vector& u);
        //Adds 2 vectors
        Vector operator+(const Vector& u);
        Vector operator+(const float& c);
        Vector& operator+=(const Vector& u);
        Vector& operator+=(const float& c);
        //Subtracts 2 vectors
        Vector operator-(const Vector& u);
        Vector operator-(const float& c);
        Vector& operator-=(const Vector& u);
        Vector& operator-=(const float& c);
};

class Material {
    public:
        Material();
        Material(Vector ambient, Vector diffuse, Vector specular);
        Material(Vector amb, Vector dif, Vector spec, Vector trans, float pow, float iref);
        
        void setAmbient(Vector amb);
        void setDiffuse(Vector dif);
        void setSpecular(Vector spec);
        void setTransmissive(Vector trns);
        void setCosPower(float pow);
        void setIndexRefract(float index);
        
        Vector getAmbient() const { return ambientColor; }
        Vector getDiffuse() const { return diffuseColor; }
        Vector getSpecular() const { return specularColor; }
        Vector getTransmissive() const { return transmissiveColor; }
        float getCosPower() const { return cosPow; }
        float getIndexRefract() const { return ior; }
    private:
        Vector ambientColor;
        Vector diffuseColor;
        Vector specularColor;
        
        float cosPow; //Phong cosine power for specular highlights
        Vector transmissiveColor;
        float ior;
};

class Sphere {
    public:
        Sphere();
        Sphere(Vector pos, float rad, Material material);
        
        void setPosition(const Vector& pos);
        void setRadius(const float& rad);
        void setMaterial(Material material);
        
        Vector getPosition() const { return position; }
        float getRadius() const { return radius; }
        float getMaterial() const { return mat; }
    private:
        Vector position;
        float radius;
        Material mat;
};

class Triangle {
    public:
        Triangle();
        Triangle(Vector v1, Vector v2, Vector v3);
        Triangle(Vector v1, Vector v2, Vector v3, Vector n1, Vector n2, Vector n3);
        
        void setVertex(const Vector& vert, int index);
        void setVertices(const Vector[3] verts);
        void setNormal(const Vector& norm);
        void setNormal(const Vector& norm, int index);
        void setNormals(const Vector[3] norms);
        
        Vector getNormal();
        Vector getNormal(int index);
        Vector getVertex(int index);
        Vector[3] getVertices();
    private:
        Vector vertices[3];
        Vector normals[3];
        bool ntri;
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

#endif
