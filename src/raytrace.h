#ifndef RAYTRACE_H_
#define RAYTRACE_H_

#include "EasyBMP.h"
#include "image.h"
#include "pixel.h"
#include "parse.h"

typedef struct rayline {
    Vector p;
    Vector d;
} Ray;

//Type of intersection
enum insectType {
    SPHERE,
    TRIANGLE
};

typedef struct intersection {
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
} Intersect;

//Performs a ray trace and returns an image
Image* RayTrace(SceneData* scn);
//Recursive ray structure
Vector evaluateRayTree(SceneData* scn, Ray* ray, int depth, bool useBVH);
//Tests for intersections
Intersect* intersect(Ray* trace, SceneData* scn, double tmin, double tmax, bool useBVH);
Intersect* intersectSpheres(Ray* trace, std::vector<Sphere> objList, double dmin, double dmax);
Intersect* intersectTriangle(Ray* trace, std::vector<Triangle> triList, double dmin, double dmax);
Intersect* intersectPlane(Ray* trace, std::vector<Plane> planeList, double dmin, double dmax);
Intersect* intersectRectangle(Ray* trace, std::vector<Rectangle> rectList, double dmin, double dmax);
//Returns a ray calculated from the pixel position and the viewing angle
Ray getRay(int x, int y, int w, int h, Vector p1, Vector p2, double pd, Camera* c, ProjType proj, bool sample);
//Returns a pixel with the background color
Vector getColor(Intersect* i, SceneData* scn, int depth, bool useBVH);
//Distance to the viewing plane
double getPlaneDist(double angle, int h);
//Finds 2 extreme points of the viewing plane
void getExtremePoints(Camera* c, double d, double w, double h, Vector& p1, Vector& p2);

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
/**
 *  Vector Operations
 */
//Performs a dot product operation on the 2 vectors
double dot(Vector u, Vector v);
//Performs a cross prooduct
Vector cross(Vector u, Vector v);
//Multiple vector
Vector multiply(Vector u, float c);
Vector div(Vector u, float c);
//Adds 2 vectors
Vector add(Vector u, Vector v);
Vector add(Vector u, float c);
//Subtracts 2 vectors
Vector sub(Vector u, Vector v);
Vector sub(Vector u, float c);
//Returns the length of the vector
double length(Vector u);
//Finds the length but does not take the square root
double lengthsq(Vector u);
//Normalizes the vector
Vector norm(Vector u);

Vector ave(Vector v[], int num);

#endif
