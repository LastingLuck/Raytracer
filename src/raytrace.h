#ifndef RAYTRACE_H_
#define RAYTRACE_H_

#include "EasyBMP.h"
#include "image.h"
#include "pixel.h"
#include "parse.h"

class Ray {
    public:
        Ray();
        Ray(const Vector& pos, const Vector& dir);
        
        void setPosition(const Vector& pos) { p = pos; }
        void setDirection(const Vector& dir) { d = dir; }
        
        Vector getPosition() const { return p; }
        Vector getDirection() const { return d; }
    private:
        Vector p;
        Vector d;
};

class Intersect {
    //Normal of intersect point
    Vector normal;
	//Intersected ray
	Ray r;
	//t value of intersection (Distancesq)
	double dist;
	//t replaced with point of intersection
	Vector p;
	//Intersection material
	Material m;
};

class RayTrace {
    public:
        RayTrace();
        
        //Performs a ray trace and returns an image
        Image* rayTrace(const Scene& scn);
        void forceBVH(bool force) { forcebvh = force; }
    private:
        bool forcebvh;
        
        Vector ave(Vector v[], int num);
        //Recursive ray structure
        Vector evaluateRayTree(const Scene& scn, const Ray& ray, int depth) const;
        //Returns a ray calculated from the pixel position and the viewing angle
        Ray getRay(int x, int y, int w, int h, const Vector& p1, const Vector& p2, 
                   double pd, const Camera& c, ProjType proj, bool sample) const;
        //Returns a pixel with the background color
        Vector getColor(const Intersect& i, const Scene& scn, int depth) const;
        //Distance to the viewing plane
        double getPlaneDist(double angle, int h) const;
        //Finds 2 extreme points of the viewing plane
        void getExtremePoints(const Camera& c, double d, double w, double h, Vector& p1, Vector& p2) const;
        //Tests for intersections
        Intersect* intersect(const Ray& trace, const Scene& scn, double tmin, double tmax);
        Intersect* intersectSpheres(const Ray& trace, const std::vector<Sphere>& objList, double dmin, double dmax);
        Intersect* intersectTriangle(const Ray& trace, const std::vector<Triangle>& triList, double dmin, double dmax);
};

//Intersect* intersectPlane(Ray* trace, std::vector<Plane> planeList, double dmin, double dmax);
//Intersect* intersectRectangle(Ray* trace, std::vector<Rectangle> rectList, double dmin, double dmax);

/**
 * Acceleration Structure
 */
void printBVH(Box* b, int depth);
void makeBVH(SceneData* scn, Box* box, int depth);
Intersect* intersectBVH(Ray* trace, Box* bvh, double dmin, double dmax);
//void findBoundingVerts(SceneData* scn, std::array<Vector, 8> (&verts));
void findBoundingVerts(SceneData* scn, Vector& bl, Vector& tr);
//char findLongestAxis(std::array<Vector, 8> vrts);
char findLongestAxis(Vector vmin, Vector vmax);
bool isSphereInBox(Box* b, Sphere* sph);
bool isTriangleInBox(Box* b, Triangle* t);
bool isPlaneInBox(Box* b, Plane* pln);
bool intersectRayAABB(Ray* trace, Box* bvh, double dmin, double dmax);
void findMinMax(float x0, float x1, float x2, float& min, float& max);

#endif
