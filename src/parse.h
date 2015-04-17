#ifndef PARSE_H_
#define PARSE_H_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <array>
//#include "raytrace.h"


enum ProjType {
    ORTHO,
    PERSP
};

enum LightType {
    AMBIENT,
    DIRECTIONAL,
    POINT, 
    SPOT
};

class Vector {
    public:
        Vector();
        Vector(float val);
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
        static Vector zero() { return Vector(); }
        static Vector one() { return Vector(1); }
        static Vector up() { return Vector(0, 1, 0); }
        static Vector down() { return Vector(0, -1, 0); }
        static Vector left() { return Vector(-1, 0, 0); }
        static Vector right() { return Vector(1, 0, 0); }
        static Vector forward() { return Vector(0, 0, 1); }
        static Vector backward() { return Vector(0, 0, -1); }
        
        //Returns the magnitude of the vector
        float magnitude();
        //Finds the magnitude but does not take the square root
        float magnitudeSq();
        //Normalizes the vector
        Vector norm();
        
        //Multiply vector
        Vector operator*(const float& c) const;
        Vector& operator*=(const Vector& u);
        //Adds 2 vectors
        Vector operator+(const Vector& u) const;
        Vector operator+(const float& c) const;
        Vector& operator+=(const Vector& u);
        Vector& operator+=(const float& c);
        //Subtracts 2 vectors
        Vector operator-(const Vector& u) const;
        Vector operator-(const float& c) const;
        Vector& operator-=(const Vector& u);
        Vector& operator-=(const float& c);
        //Divides each entry in the vector by a constant
        Vector operator/(const float& c) const;
        Vector& operator/=(const float& c);
};

class Material {
    public:
        Material();
        Material(const Vector& ambient, const Vector& diffuse, const Vector& specular);
        Material(const Vector& amb, const Vector& dif, const Vector& spec, const Vector& trans, float pow, float iref);
        
        void setAmbient(const Vector& amb) { ambientColor = amb; }
        void setDiffuse(const Vector& dif) { diffuseColor = dif; }
        void setSpecular(const Vector& spec) { specularColor = spec; }
        void setTransmissive(const Vector& trns) { transmissiveColor = trns; }
        void setCosPower(float pow) { cosPow = pow; }
        void setIndexRefract(float index) { ior = index; }
        
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
        Sphere(const Vector& pos, float rad, const Material& material);
        
        void setPosition(const Vector& pos) { position = pos; }
        void setRadius(const float& rad) { radius = rad; }
        void setMaterial(Material material) { mat = material; }
        
        Vector getPosition() const { return position; }
        float getRadius() const { return radius; }
        Material getMaterial() const { return mat; }
    private:
        Vector position;
        float radius;
        Material mat;
};

class Triangle {
    public:
        Triangle();
        Triangle(const Vector& v1, const Vector& v2, const Vector& v3);
        Triangle(const Vector& v1, const Vector& v2, const Vector& v3, const Material& m);
        Triangle(const Vector& v1, const Vector& v2, const Vector& v3, const Vector& n1, const Vector& n2, const Vector& n3);
        Triangle(const Vector& v1, const Vector& v2, const Vector& v3, const Vector& n1, const Vector& n2, const Vector& n3, const Material& m);
        
        void setMaterial(const Material& material) { mat = material; }
        void setVertex(const Vector& vert, int index) { vertices[clamp(index, 0, 2)] = vert; }
        void setVertices(const Vector verts[3]);
        void setNormal(const Vector& norm) { normals[0] = normals[1] = normals[2] = norm; ntri = false;}
        void setNormal(const Vector& norm, int index) { normals[clamp(index, 0, 2)] = norm; ntri = true; }
        void setNormals(const Vector norms[3]);
        
        Material getMaterial() const { return mat; }
        Vector getNormal() const { return normals[0]; }
        Vector getNormal(int index) const { return normals[clamp(index, 0, 2)]; }
        Vector getVertex(int index) const { return vertices[clamp(index, 0, 2)]; }
        std::vector<Vector> getVertices() const { return std::vector<Vector>(vertices, vertices+3); }
        std::vector<Vector> getNormals() const { return std::vector<Vector>(normals, normals+3); }
        bool isNormal() { return ntri; } //< Returns wether the triangle has seperate normals for each point
        Vector& operator[](const int index);
    private:
        int clamp(int num, int min, int max) const { return (num < min) ? min : (num > max) ? max : num; }
        
        Vector vertices[3];
        Vector normals[3];
        bool ntri;
        Material mat;
};

class Plane {
    public:
        Plane();
        Vector normal;
        Vector point;
    private:
        
};

class Light {
    public:
        Light();
        Light(const Vector& lightColor);
        Light(const Vector& lightColor, const Vector& lightPosDir, enum LightType type);
        Light(const Vector& lightColor, const Vector& lightPos, const Vector& lightDir, float spotAngle, float maxAngle);
        
        //these do not affect the type of light (i.e. changing direction of a spot
        //light has no effect)
        void setColor(const Vector& newColor) { color = newColor; }
        void setPosition(const Vector& newPos) { position = newPos; }
        void setDirection(const Vector& newDir) { direction = newDir; }
        void setSpotAngle(float angle) { angle1 = angle; }
        void setMaxAngle(float angle) { angle2 = angle; }
        void setType(enum LightType type) { ltype = type; }
        
        Vector getColor() const { return color; }
        Vector getPosition() const { return position; }
        Vector getDirection() const { return direction; }
        float getSpotAngle() const { return angle1; }
        float getMaxAngle() const { return angle2; }
        enum LightType getType() const { return ltype; }
        
    private:
        Vector color;         //Color of the light (All)
        Vector position;      //Position of the Light (Point, Spot)
        Vector direction;     //Direction the light is pointing (Directional, Spot)
        float angle1;         //Angle where it acts as a point light (Spot)
        float angle2;         //Max angle that the spot light reaches (Spot)
        enum LightType ltype; //Light type (0-ambient, 1-directional, 2-point, 3-spot)
};

class Camera {
    public:
        Camera();
        Camera(const Vector& pos, const Vector& dir, const Vector& upDir, float fov);
        
        void setPosition(const Vector& pos) { position = pos; }
        void setDirection(const Vector& dir) { direction = dir; }
        void setUpDirection(const Vector& upDir) { up = upDir; }
        void setFOV(float angle) { halfAngle = angle / 2.0f; }
        void setHalfFOV(float angle) { halfAngle = angle; }
        
        Vector getPosition() const { return position; }
        Vector getDirection() const { return direction; }
        Vector getUpDirection() const { return up; }
        float getFOV() const { return 2.0f * halfAngle; }
        float getHalfFOV() const { return halfAngle; }
    private:
        Vector position;
        Vector direction;
        Vector up;
        float halfAngle;
};

namespace rt { //Image conflicts with EasyBMP class of same name.
class Image {
    public:
        Image();
        Image(const std::string& name);
        Image(const std::string& name, int imageWidth, int imageHeight);
        
        void setFileName(const std::string& name) { filename = name; }
        void setWidth(int newWidth) { width = newWidth; }
        void setHeight(int newHeight) { height = newHeight; }
        
        std::string getFileName() const { return filename; }
        int getWidth() const { return width; }
        int getHeight() const { return height; }
    private:
        std::string filename;     //Output filename
        int width;                //Width of image
        int height;               //Height of image
};
} //end rt namespace

class AABB {
    public:
        AABB();
        AABB(const Vector& minimum, const Vector& maximum);
        AABB(AABB* p);
        AABB(const std::vector<Sphere> sphs, const std::vector<Triangle> tris);
        
        void setMin(const Vector& m) { min = m; }
        void setMax(const Vector& m) { max = m; }
        void setLeftChild(AABB* child) { left = child; }
        void setRightChild(AABB* child) { right = child; }
        void setParent(AABB* p) { parent = p; }
        void addSphere(const Sphere& sph) { spheres.push_back(sph); }
        void addTriangle(const Triangle& tri) { triangles.push_back(tri); }
        void setSpheres(const std::vector<Sphere> sphs) { spheres = sphs; }
        void setTriangles(const std::vector<Triangle> tris) { triangles = tris; }
        
        Vector getMin() const { return min; }
        Vector getMax() const { return max; }
        std::vector<Sphere> getSpheres() const { return spheres; }
        std::vector<Triangle> getTriangles() const { return triangles; }
        bool isLeaf() const { return !left && !right; }
        AABB* getLeftChild() const { return left; }
        AABB* getRightChild() const { return right; }
        AABB* getParent() const { return parent; }
        int getObjectNum() const { return spheres.size()+triangles.size(); }
        
        bool isInBox(const Sphere& sph) const;
        bool isInBox(const Triangle& tri) const;
        bool isInBox(const Plane& pln) const;
    private:
        std::vector<Sphere> spheres;
        std::vector<Triangle> triangles;
        Vector min;
        Vector max;
        //May not actually be left and right in space.
        AABB* left;
        AABB* right;
        AABB* parent;
        
        void findMinMax(float x0, float x1, float x2, float& min, float& max) const;
};

class Scene {
    public:
        Scene();
        Scene(enum ProjType view);
        
        void init(); //< (re)initialize everything to default values
        //Parses the scene text file and fills out a SceneData struct
        int parseScene(char* file);
        float getPlaneDist() const;
        void clearPools();
        
        void setCamera(const Camera& cam) { camera = cam; }
        void setProjType(enum ProjType view) { proj = view; }
        void setImage(const rt::Image& image) { img = image; }
        void addSphere(const Sphere& sphere) { spheres.push_back(sphere); objNum++; }
        void addVertex(const Vector& vert) { vertexPool.push_back(vert); }
        void addNormal(const Vector& norm) { normalPool.push_back(norm); }
        void addTriangle(const Triangle& tri) { triangles.push_back(tri); objNum++; }
        void setBackground(const Vector& bg) { bgColor = bg; }
        void setMaxDepth(int maximumDepth) { maxDepth = maximumDepth; }
        void addLight(const Light& dirLight) { lights.push_back(dirLight); }
        void setAmbientLight(const Light& ambLight) { ambient = ambLight; }
        void setBVHDepth(int depth) { bvhDepth = depth; }
        void setBVHThreshold(int threshold) { bvhThresh = threshold; }
        void setBVHRoot(AABB* bvhRoot) { bvh = bvhRoot; }
        void setSampleRate(int rate) { sampleNum = rate; }
        
        Camera getCamera() const { return camera; }
        enum ProjType getProjType() const { return proj; }
        rt::Image getImage() const { return img; }
        std::vector<Sphere> getSpheres() const { return spheres; }
        std::vector<Vector> getVertices() const { return vertexPool; }
        std::vector<Vector> getNormals() const { return normalPool; }
        std::vector<Triangle> getTriangles() const { return triangles; }
        Vector getBGColor() const { return bgColor; }
        int getDepth() const { return maxDepth; }
        std::vector<Light> getLights() const { return lights; }
        Light getAmbLight() const { return ambient; }
        int getNumObjects() const { return objNum; }
        int getBVHDepth() const { return bvhDepth; }
        int getBVHThreshold() const { return bvhThresh; }
        bool useBvh() const { return useBVH; }
        AABB* getBVHRoot() const { return bvh; }
        int getSampleRate() const { return sampleNum; }
        
    private:
        //Camera
        Camera camera;
        enum ProjType proj;
        //Image
        rt::Image img;
        //Sphere
        std::vector<Sphere> spheres;
        //Triangle
        std::vector<Vector> vertexPool;
        std::vector<Vector> normalPool;
        std::vector<Triangle> triangles;
        //Background
        Vector bgColor;
        //depth
        int maxDepth;
        //Lights
        std::vector<Light> lights;
        //std::vector<Light> directional;
        //std::vector<Light> point;
        //std::vector<Light> spot;
        Light ambient;
        //Misc
        int objNum;
        //BVH
        int bvhDepth;
        int bvhThresh;
        bool useBVH;
        AABB* bvh;
        //Super Sample
        int sampleNum;
};

//Print functions used to print the contents of SceneData to test it
void printScene(Scene *scn);
void printCamera(Camera camera);
void printSphere(Sphere sph);
void printTriangle(Triangle tri);
void printMaterial(Material mat);
void printLight(Light l);
void printImage(rt::Image img);
void printVector(Vector fVec);
void printBVH(AABB* b, int depth);

#endif
