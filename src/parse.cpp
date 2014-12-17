#include "parse.h"
#include <cstdio>
#include <cstring>
#include <fstream>
#include <stdexcept>

#define MAX_LINE 100

/********************
 * Scene
 ********************/
Scene::Scene() {
    init();
}

Scene::Scene(enum ProjType view) {
    init();
    proj = view;
}

void Scene::init() {
    //Camera
    camera = Camera((Vector){0, 0, 0}, (Vector){0, 0, 1}, (Vector){0, 1, 0}, 45.0*(M_PI / 180.0));
    //Image
    img = rt::Image("raytraced.bmp", 640, 480);
    //Light
    ambient = Light(Vector(0, 0, 0));
    //Misc
    bgColor = Vector(0, 0, 0);
    maxDepth = 5;
    proj = PERSP;
    objNum = 0;
    //BVH
    bvhThresh = 500;
    bvhDepth = 5;
    useBVH = false;
    //Super sampling
    sampleNum = 1;
    #ifdef DEBUG
    printf("Setup Done\n");
    #endif
}

int Scene::parseScene(char* file) {
    char trash[256];
    std::ifstream scn;
    scn.open(file);
    if(scn.fail()) {
        printf("File could not be opened: %s\n", file);
        return -1;
    }
    Material curMat(Vector::zero(), Vector::one(), Vector::zero(), Vector::zero(), 5, 1);
    std::string str;
    while(scn >> str) {
        if(str == "" || str[0] == '#') {
            scn.getline(trash, 256);
            continue;
        }
        else if(str == "camera") {
            float x, y, z, dx, dy, dz, ux, uy, uz, ha;
            scn >> x >> y >> z >> dx >> dy >> dz >> ux >> uy >> uz >> ha;
            camera.setPosition(Vector(x, y, z));
            camera.setDirection(Vector(dx, dy, dz));
            camera.setUpDirection(Vector(ux, uy, uz));
            camera.setHalfFOV(ha);
        }
        else if(str == "film_resolution") {
            int w, h;
            scn >> w >> h;
            img.setWidth(w);
            img.setHeight(w);
        }
        else if(str == "output_image") {
            std::string n;
            scn >> n;
            img.setFileName(n);
        }
        else if(str == "sphere") {
            //std::string tmp;
            //scn >> tmp;
            float x, y, z, r;
            scn >> x >> y >> z >> r;
            //printf("Sphere at (%s)\n", tmp.c_str());
            addSphere(Sphere(Vector(x, y, z), r, curMat));
        }
        else if(str == "background") {
            float r, g, b;
            scn >> r >> g >> b;
            setBackground(Vector(r, g, b));
        }
        else if(str == "material") {
            float ar, ag, ab, dr, dg, db, sr, sg, sb, ns, tr, tg, tb, ior;
            scn >> ar >> ag >> ab >> dr >> dg >> db >> sr >> sg >> sb >> ns 
                >> tr >> tg >> tb >> ior;
            curMat.setAmbient(Vector(ar, ag, ab));
            curMat.setDiffuse(Vector(dr, dg, db));
            curMat.setSpecular(Vector(sr, sg, sb));
            curMat.setTransmissive(Vector(tr, tg, tb));
            curMat.setCosPower(ns);
            curMat.setIndexRefract(ior);
        }
        else if(str == "directional_light") {
            float r, g, b, x, y, z;
            scn >> r >> g >> b >> x >> y >> z;
            addDirectionalLight(Light(Vector(r, g, b), Vector(x, y, z), DIRECTIONAL));
        }
        else if(str == "point_light") {
            float r, g, b, x, y, z;
            scn >> r >> g >> b >> x >> y >> z;
            addSpotLight(Light(Vector(r, g, b), Vector(x, y, z), POINT));
        }
        else if(str == "spot_light") {
            float r, g, b, px, py, pz, dx, dy, dz, angle1, angle2;
            scn >> r >> g >> b >> px >> py >> pz >> dx >> dy >> dz >> angle1 >> angle2;
            addSpotLight(Light(Vector(r, g, b), Vector(px, py, pz), Vector(dx, dy, dz), angle1, angle2));
        }
        else if(str == "ambient_light") {
            float r, g, b;
            scn >> r >> g >> b;
            setAmbientLight(Light(Vector(r, g, b)));
        }
        else if(str == "max_depth") {
            int n;
            scn >> n;
            setMaxDepth(n);
        }
        else if(str == "ray_type") {
            std::string type;
            scn >> type;
            if(type == "perspective") {
                setProjType(PERSP);
            }
            else if(type == "orthographic") {
                setProjType(ORTHO);
            }
            else {
                printf("Unknown Projection Type: %s\n", type.c_str());
				printf("Defaulting to Perspective\n");
                setProjType(PERSP);
            }
        }
        else if(str == "vertex") {
            float x, y, z;
            scn >> x >> y >> z;
            addVertex(Vector(x, y, z));
        }
        else if(str == "normal") {
            float x, y, z;
            scn >> x >> y >> z;
            addNormal(Vector(x, y, z).norm());
        }
        else if(str == "triangle") {
            int v1, v2, v3;
            scn >> v1 >> v2 >> v3;
            std::vector<Vector> vl = getVertices();
            int vn = vl.size();
            if(v1 >= vn || v1 < 0) {
                printf("Error: Specified vertex (%d) in triangle (%d) does not exist\n", v1, (int)getTriangles().size());
                return -1;
            }
            if(v2 >= vn || v2 < 0) {
                printf("Error: Specified vertex (%d) in triangle (%d) does not exist\n", v2, (int)getTriangles().size());
                return -1;
            }
            if(v3 >= vn || v3 < 0) {
                printf("Error: Specified vertex (%d) in triangle (%d) does not exist\n", v3, (int)getTriangles().size());
                return -1;
            }
            addTriangle(Triangle(vl[v1], vl[v2], vl[v3]));
        }
        else if(str == "normal_triangle") {
            int v1, v2, v3, n1, n2, n3;
            scn >> v1 >> v2 >> v3 >> n1 >> n2 >> n3;
            std::vector<Vector> vl = getVertices();
            int vn = vl.size();
            std::vector<Vector> nl = getNormals();
            int nn = nl.size();
            if(v1 >= vn || v1 < 0) {
                printf("Error: Specified vertex (%d) in triangle (%d) does not exist\n", v1, (int)getTriangles().size());
                return -1;
            }
            if(v2 >= vn || v2 < 0) {
                printf("Error: Specified vertex (%d) in triangle (%d) does not exist\n", v2, (int)getTriangles().size());
                return -1;
            }
            if(v3 >= vn || v3 < 0) {
                printf("Error: Specified vertex (%d) in triangle (%d) does not exist\n", v3, (int)getTriangles().size());
                return -1;
            }
            if(n1 >= nn || n1 < 0) {
                printf("Error: Specified normal (%d) in triangle (%d) does not exist\n", n1, (int)getTriangles().size());
                return -1;
            }
            if(n2 >= nn || n2 < 0) {
                printf("Error: Specified normal (%d) in triangle (%d) does not exist\n", n2, (int)getTriangles().size());
                return -1;
            }
            if(n3 >= nn || n3 < 0) {
                printf("Error: Specified normal (%d) in triangle (%d) does not exist\n", n3, (int)getTriangles().size());
                return -1;
            }
            addTriangle(Triangle(vl[v1], vl[v2], vl[v3], nl[n1], nl[n2], nl[n3]));
        }
        else if(str == "bvh_threshold") {
            int thresh;
            scn >> thresh;
            setBVHThreshold(thresh);
        }
        else if(str == "bvh_depth") {
            int depth;
            scn >> depth;
            setBVHDepth(depth);
        }
        else if(str == "sample_rate") {
            int rate;
            scn >> rate;
            setSampleRate(rate);
        }
    }
    clearPools();
    scn.close();
    return 1;
}

//Returns the distance to the viewing plane based on the angle from the
//viewing point and the height of the plane
float Scene::getPlaneDist() const {
    return img.getHeight() / (2.0 * tan(camera.getHalfFOV()));
}

void Scene::clearPools() {
    vertexPool.clear();
    normalPool.clear();
}

/********************
 * Vector
 ********************/
Vector::Vector() {
    x = y = z = 0;
}

Vector::Vector(float val) {
    x = y = z = val;
}

Vector::Vector(float xval, float yval, float zval) {
    x = xval;
    y = yval;
    z = zval;
}

float Vector::dot(const Vector& u, const Vector& v) {
    return (u.x * v.x) + (u.y * v.y) + (u.z * v.z);
}

Vector Vector::cross(const Vector& u, const Vector& v) {
    return Vector(u.y*v.z-u.z*v.y, u.z*v.x-u.x*v.z, u.x*v.y-u.y*v.x);
}

float Vector::length(const Vector& u, const Vector& v) {
    Vector dist = u - v;
    return sqrt((dist.x * dist.x) + (dist.y * dist.y) + (dist.z * dist.z));
}

float Vector::lengthSq(const Vector& u, const Vector& v) {
    const Vector dist = u - v;
    return (dist.x * dist.x) + (dist.y * dist.y) + (dist.z * dist.z);
}

float Vector::magnitude() {
    return sqrt((x * x) + (y * y) + (z * z));
}

float Vector::magnitudeSq() {
    return (x * x) + (y * y) + (z * z);
}

Vector Vector::norm() {
    float mag = magnitude();
    return Vector(x/mag, y/mag, z/mag);
}

Vector Vector::operator*(const float& c) const {
    return Vector(x*c, y*c, z*c);
}

Vector& Vector::operator*=(const Vector& u) {
    this->x *= u.x;
    this->y *= u.y;
    this->z *= u.z;
    return *this;
}

Vector Vector::operator+(const Vector& u) const {
    return Vector(x+u.x, y+u.y, z+u.z);
}

Vector Vector::operator+(const float& c) const {
    return Vector(x+c, y+c, z+c);
}

Vector& Vector::operator+=(const Vector& u) {
    this->x += u.x;
    this->y += u.y;
    this->z += u.z;
    return *this;
}

Vector& Vector::operator+=(const float& c) {
    this->x += c;
    this->y += c;
    this->z += c;
    return *this;
}

Vector Vector::operator-(const Vector& u) const {
    return Vector(x-u.x, y-u.y, z-u.z);
}

Vector Vector::operator-(const float& c) const {
    return Vector(x-c, y-c, z-c);
}

Vector& Vector::operator-=(const Vector& u) {
    this->x -= u.x;
    this->y -= u.y;
    this->z -= u.z;
    return *this;
}

Vector& Vector::operator-=(const float& c) {
    this->x -= c;
    this->y -= c;
    this->z -= c;
    return *this;
}

Vector Vector::operator/(const float& c) const {
    return Vector(x/c, y/c, z/c);
}

Vector& Vector::operator/=(const float& c) {
    this->x /= c;
    this->y /= c;
    this->z /= c;
    return *this;
}

/********************
 * Material
 ********************/
Material::Material() {
    ambientColor = Vector();
    diffuseColor = Vector(1);
    specularColor = Vector();
    cosPow = 5;
    transmissiveColor = Vector();
    ior = 1;
}
 
Material::Material(const Vector& ambient, const Vector& diffuse, const Vector& specular) {
    ambientColor = ambient;
    diffuseColor = diffuse;
    specularColor = specular;
    cosPow = 5;
    transmissiveColor = Vector();
    ior = 1;
}

Material::Material(const Vector& amb, const Vector& dif, const Vector& spec, const Vector& trans, float pow, float iref) {
    ambientColor = amb;
    diffuseColor = dif;
    specularColor = spec;
    cosPow = pow;
    transmissiveColor = trans;
    ior = iref;
}

/********************
 * Sphere
 ********************/
Sphere::Sphere() {
    position = Vector();
    radius = 1.0f;
    mat = Material();
}

Sphere::Sphere(const Vector& pos, float rad, const Material& material) {
    position = pos;
    radius = rad;
    mat = material;
}

/********************
 * Triangle
 ********************/
Triangle::Triangle() {
    vertices[0] = vertices[1] = vertices[2] = Vector();
    normals[0] = normals[1] = normals[2] = Vector(0, 1, 0);
    mat = Material();
    ntri = false;
}

Triangle::Triangle(const Vector& v1, const Vector& v2, const Vector& v3) {
    vertices[0] = v1;
    vertices[1] = v2;
    vertices[2] = v3;
    normals[0] = normals[1] = normals[2] = Vector::cross(v1, v2).norm();
    mat = Material();
    ntri = false;
}

Triangle::Triangle(const Vector& v1, const Vector& v2, const Vector& v3, const Vector& n1, const Vector& n2, const Vector& n3) {
    vertices[0] = v1;
    vertices[1] = v2;
    vertices[2] = v3;
    normals[0] = n1;
    normals[1] = n2;
    normals[2] = n3;
    mat = Material();
    ntri = true;
}

Vector& Triangle::operator[](const int index) {
    if(index < 0 || index > 2) {
        throw std::out_of_range("Triangle index out of range");
    }
    return vertices[index];
}

/********************
 * Light
 ********************/
Light::Light() {
    color = Vector();
    position = Vector();
    direction = Vector(0, 1, 0);
    angle1 = M_PI_4;
    angle2 = M_PI_2;
}

Light::Light(const Vector& lightColor) {
    color = lightColor;
    position = Vector();
    direction = Vector(0, 1, 0);
    angle1 = M_PI_4;
    angle2 = M_PI_2;
}

Light::Light(const Vector& lightColor, const Vector& lightPosDir, enum LightType type) {
    color = lightColor;
    if(type == POINT) {
        position = lightPosDir;
        direction = Vector(0, 1, 0);
    }
    else {
        position = Vector();
        direction = lightPosDir;
    }
    angle1 = M_PI_4;
    angle2 = M_PI_2;
}

Light::Light(const Vector& lightColor, const Vector& lightPos, const Vector& lightDir, float spotAngle, float maxAngle) {
    color = lightColor;
    position = lightPos;
    direction = lightDir;
    angle1 = spotAngle;
    angle2 = maxAngle;
}

/********************
 * Camera
 ********************/
Camera::Camera() {
    position = Vector();
    direction = Vector(1, 0, 0);
    up = Vector(0, 1, 0);
    halfAngle = 45.0f;
}

Camera::Camera(const Vector& pos, const Vector& dir, const Vector& upDir, float fov) {
    position = pos;
    direction = dir;
    up = upDir;
    halfAngle = fov;
}

/********************
 * Image
 ********************/
rt::Image::Image() {
    filename = "raytraced.bmp";
    width =  640;
    height = 480;
}

rt::Image::Image(const std::string& name) {
    filename = name;
    width = 640;
    height = 480;
}

rt::Image::Image(const std::string& name, int imageWidth, int imageHeight) {
    filename = name;
    width = imageWidth;
    height = imageHeight;
}

/********************
 * AABB
 ********************/
AABB::AABB() {
    left = right = parent = 0;
}

AABB::AABB(AABB* p) {
    left = right = 0;
    parent = p;
}

AABB::AABB(std::vector<Sphere> sphs, std::vector<Triangle> tris) {
    spheres = sphs;
    triangles = tris;
    left = right = parent = 0;
}



/********************
 * Printing
 ********************/
void printScene(Scene *scn) {
    printf("\nCamera:\n");
    printCamera(scn->getCamera());
    printf("\n");
    
    printf("Spheres:\n");
    int num = scn->getSpheres().size();
    if(num != 0) {
        for (int i = 0; i < num; i++) {
            printSphere(scn->getSpheres()[i]);
        }
    }
    else {
        printf("No Spheres\n");
    }
    printf("\n");
    
    printf("Image:\n");
    printImage(scn->getImage());
    printf("\n");
    
    printf("Background Color: ");
    printVector(scn->getBGColor());
    printf("\n");
    
    printf("Directional Lights:\n");
    num = scn->getDirLights().size();
    if (num > 0) {
        for (int i = 0; i < num; i++) {
            printLight(scn->getDirLights()[i]);
        }
    }
    else {
        printf("No Directional Lights\n");
    }
    printf("\n");
    
    printf("Point Lights:\n");
    num = scn->getPointLights().size();
    if (num > 0) {
        for (int i = 0; i < num; i++) {
            printLight(scn->getPointLights()[i]);
        }
    }
    else {
        printf("No Point Lights\n");
    }
    printf("\n");
    
    printf("Spot Lights:\n");
    num = scn->getSpotLights().size();
    if (num > 0) {
        for (int i = 0; i < num; i++) {
            printLight(scn->getSpotLights()[i]);
        }
    }
    else {
        printf("No Spot Lights\n");
    }
    printf("\n");
    
    printf("Ambient Light:\n");
    printLight(scn->getAmbLight());
    printf("\n");
    
    printf("Depth: %d\n", scn->getDepth());
    
    if(scn->getProjType() == PERSP) {
        printf("Projection Type: Perspective\n\n");
    }
    else {
        printf("Projection Type: Orthographic\n\n");
    }
    
    printf("Triangles:\n");
    num = scn->getTriangles().size();
    if(num != 0) {
		for(int i = 0; i < num; i++) {
			printTriangle(scn->getTriangles()[i]);
		}	
	}
	else {
		printf("No Triangles\n");
	}
	printf("\n");
	/*
	printf("Planes:\n");
	num = scn->planes.size();
	if(num != 0) {
		for(int i = 0; i < num; i++) {
			printPlane(scn->planes[i]);
		}
	}
	else {
		printf("No Planes\n");
	}
    printf("\n");
    */
}

void printCamera(Camera camera) {
    printf("Position: ");
    printVector(camera.getPosition());
    printf("Direction: ");
    printVector(camera.getDirection());
    printf("Up Vector: ");
    printVector(camera.getUpDirection());
    printf("Half Angle: %f\n", camera.getHalfFOV());
}

void printTriangle(Triangle tri) {
    printf("Vertex1: ");
    printVector(tri.getVertex(0));
    printf("Vertex2: ");
    printVector(tri.getVertex(1));
    printf("Vertex3: ");
    printVector(tri.getVertex(2));
    printf("Normal: ");
    printVector(tri.getNormal(0));
    if(tri.isNormal()) {
        printf("Normal2: ");
        printVector(tri.getNormal(1));
        printf("Normal3: ");
        printVector(tri.getNormal(2));
    }
    printf("Material:\n");
    printMaterial(tri.getMaterial());
}

void printSphere(Sphere sph) {
    printf("Position: ");
    printVector(sph.getPosition());
    printf("Radius: %f\n", sph.getRadius());
    printf("Material:\n");
    printMaterial(sph.getMaterial());
}

void printMaterial(Material mat) {
    printf("Ambient Color: ");
    printVector(mat.getAmbient());
    printf("Diffuse Color: ");
    printVector(mat.getDiffuse());
    printf("Specular Color: ");
    printVector(mat.getSpecular());
    printf("Cosin Power: %f\n", mat.getCosPower());
    printf("Transmissive Color: ");
    printVector(mat.getTransmissive());
    printf("ior: %f\n", mat.getIndexRefract());
}

void printLight(Light l) {
    printf("Color: ");
    printVector(l.getColor());
    printf("Position: ");
    printVector(l.getPosition());
    printf("Direction: ");
    printVector(l.getDirection());
    printf("Angle: (%f %f)\n", l.getSpotAngle(), l.getMaxAngle());
}

void printImage(rt::Image img) {
    printf("File Name: %s\n", img.getFileName().c_str());
    printf("Width: %d\n", img.getWidth());
    printf("Height: %d\n", img.getHeight());
}

void printVector(Vector fVec) {
    printf("(%f %f %f)\n", fVec.x, fVec.y, fVec.z);
}
