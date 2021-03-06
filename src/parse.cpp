#include "parse.h"
#include <cstdio>
#include <cstring>
#include <fstream>
#include <stdexcept>

#define DEGTORAD 0.01745329251

/*======================== X-tests ========================*/

#define AXISTEST_X01(a, b, fa, fb)			   \
	p0 = a*v0.y - b*v0.z;			       	   \
	p2 = a*v2.y - b*v2.z;			       	   \
        if(p0<p2) {fmin=p0; fmax=p2;} else {fmin=p2; fmax=p0;} \
	rad = fa * boxhalfsize.y + fb * boxhalfsize.z;   \
	if(fmin>rad || fmax<-rad) return false;

#define AXISTEST_X2(a, b, fa, fb)			   \
	p0 = a*v0.y - b*v0.z;			           \
	p1 = a*v1.y - b*v1.z;			       	   \
        if(p0<p1) {fmin=p0; fmax=p1;} else {fmin=p1; fmax=p0;} \
	rad = fa * boxhalfsize.y + fb * boxhalfsize.z;   \
	if(fmin>rad || fmax<-rad) return false;

/*======================== Y-tests ========================*/

#define AXISTEST_Y02(a, b, fa, fb)			   \
	p0 = -a*v0.x + b*v0.z;		      	   \
	p2 = -a*v2.x + b*v2.z;	       	       	   \
        if(p0<p2) {fmin=p0; fmax=p2;} else {fmin=p2; fmax=p0;} \
	rad = fa * boxhalfsize.x + fb * boxhalfsize.z;   \
	if(fmin>rad || fmax<-rad) return false;

#define AXISTEST_Y1(a, b, fa, fb)			   \
	p0 = -a*v0.x + b*v0.z;		      	   \
	p1 = -a*v1.x + b*v1.z;	     	       	   \
        if(p0<p1) {fmin=p0; fmax=p1;} else {fmin=p1; fmax=p0;} \
	rad = fa * boxhalfsize.x + fb * boxhalfsize.z;   \
	if(fmin>rad || fmax<-rad) return false;

/*======================== Z-tests ========================*/

#define AXISTEST_Z12(a, b, fa, fb)			   \
	p1 = a*v1.x - b*v1.y;			           \
	p2 = a*v2.x - b*v2.y;			       	   \
        if(p2<p1) {fmin=p2; fmax=p1;} else {fmin=p1; fmax=p2;} \
	rad = fa * boxhalfsize.x + fb * boxhalfsize.y;   \
	if(fmin>rad || fmax<-rad) return false;

#define AXISTEST_Z0(a, b, fa, fb)			   \
	p0 = a*v0.x - b*v0.y;				   \
	p1 = a*v1.x - b*v1.y;			           \
    if(p0<p1) {fmin=p0; fmax=p1;} else {fmin=p1; fmax=p0;} \
	rad = fa * boxhalfsize.x + fb * boxhalfsize.y;   \
	if(fmin>rad || fmax<-rad) return false;

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
    camera = Camera(Vector(), Vector(0, 0, 1), Vector(0, 1, 0), 45.0*DEGTORAD);
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
            camera.setHalfFOV(ha*DEGTORAD);
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
            addLight(Light(Vector(r, g, b), Vector(x, y, z), DIRECTIONAL));
        }
        else if(str == "point_light") {
            float r, g, b, x, y, z;
            scn >> r >> g >> b >> x >> y >> z;
            addLight(Light(Vector(r, g, b), Vector(x, y, z), POINT));
        }
        else if(str == "spot_light") {
            float r, g, b, px, py, pz, dx, dy, dz, angle1, angle2;
            scn >> r >> g >> b >> px >> py >> pz >> dx >> dy >> dz >> angle1 >> angle2;
            addLight(Light(Vector(r, g, b), Vector(px, py, pz), Vector(dx, dy, dz), angle1, angle2));
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
            addTriangle(Triangle(vl[v1], vl[v2], vl[v3], curMat));
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
            addTriangle(Triangle(vl[v1], vl[v2], vl[v3], nl[n1], nl[n2], nl[n3], curMat));
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
    normals[0] = normals[1] = normals[2] = Vector::cross(v2-v1, v3-v1).norm();
    mat = Material();
    ntri = false;
}

Triangle::Triangle(const Vector& v1, const Vector& v2, const Vector& v3, const Material& m) {
    vertices[0] = v1;
    vertices[1] = v2;
    vertices[2] = v3;
    normals[0] = normals[1] = normals[2] = Vector::cross(v2-v1, v3-v1).norm();
    mat = m;
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

Triangle::Triangle(const Vector& v1, const Vector& v2, const Vector& v3, const Vector& n1, const Vector& n2, const Vector& n3, const Material& m) {
    vertices[0] = v1;
    vertices[1] = v2;
    vertices[2] = v3;
    normals[0] = n1;
    normals[1] = n2;
    normals[2] = n3;
    mat = m;
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
    ltype = AMBIENT;
}

Light::Light(const Vector& lightColor) {
    color = lightColor;
    position = Vector();
    direction = Vector(0, 1, 0);
    angle1 = M_PI_4;
    angle2 = M_PI_2;
    ltype = AMBIENT;
}

Light::Light(const Vector& lightColor, const Vector& lightPosDir, enum LightType type) {
    color = lightColor;
    if(type == POINT) {
        position = lightPosDir;
        direction = Vector(0, 1, 0);
    }
    else if(type == DIRECTIONAL) {
        position = Vector();
        direction = lightPosDir;
    }
    ltype = type;
    angle1 = M_PI_4;
    angle2 = M_PI_2;
}

Light::Light(const Vector& lightColor, const Vector& lightPos, const Vector& lightDir, float spotAngle, float maxAngle) {
    color = lightColor;
    position = lightPos;
    direction = lightDir;
    angle1 = spotAngle;
    angle2 = maxAngle;
    ltype = SPOT;
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
 * Plane
 ********************/

Plane::Plane() {
    normal = Vector();
    point = Vector();
}

/********************
 * AABB
 ********************/
AABB::AABB() {
    left = right = parent = 0;
    min = max = Vector();
}

AABB::AABB(const Vector& minimum, const Vector& maximum) {
    left = right = parent = 0;
    min = minimum;
    max = maximum;
}

AABB::AABB(AABB* p) {
    left = right = 0;
    parent = p;
    //parent = new AABB(p.getSpheres(), p.getTriangles());
    //parent->setMin(p.getMin());
    //parent->setMax(p.getMax());
    //parent->setParent(p.getParent());
    //parent->setLeftChild(p.getLeftChild());
    //parent->setRightChild(p.getRightChild());
}

AABB::AABB(const std::vector<Sphere> sphs, const std::vector<Triangle> tris) {
    spheres = sphs;
    triangles = tris;
    left = right = parent = 0;
}

//http://www.mrtc.mdh.se/projects/3Dgraphics/paperF.pdf
//On Faster Sphere-Box Overlap Testing
bool AABB::isInBox(const Sphere& sph) const {
    double d = 0, e, rad = sph.getRadius();
	Vector spos = sph.getPosition();
	//X
	if((e = spos.x - min.x) < 0) {
		if(e < -rad) {
			return false;
		}
		d += e*e;
	}
	else if((e = spos.x - max.x) > 0) {
		if(e > rad) {
			return false;
		}
		d += e*e;
	}
	//Y
	if((e = spos.y - min.y) < 0) {
		if(e < -rad) {
			return false;
		}
		d += e*e;
	}
	else if((e = spos.y - max.y) > 0) {
		if(e > rad) {
			return false;
		}
		d += e*e;
	}
	//Z
	if((e = spos.z - min.z) < 0) {
		if(e < -rad) {
			return false;
		}
		d += e*e;
	}
	else if((e = spos.z - max.z) > 0) {
		if(e > rad) {
			return false;
		}
		d += e*e;
	}
	
	return (d <= (rad*rad));
}

//http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/pubs/tribox.pdf
//Fast 3D Triangle-Box Overlap Testing
bool AABB::isInBox(const Triangle& tri) const {
	Vector v0, v1, v2, e0, e1, e2, normal;
    float fex, fey, fez, fmin, fmax, p0, p1, p2, rad;
    Vector boxhalfsize = (max - min) / 2.0;
    Vector bc = min + boxhalfsize;
    std::vector<Vector> verts = tri.getVertices();
    v0 = verts[0] - bc;
    v1 = verts[1] - bc;
    v2 = verts[2] - bc;
    
    e0 = v1 - v0;
    e1 = v2 - v1;
    e2 = v0 - v2;
    
    //Bullet 3
    fex = fabs(e0.x);
    fey = fabs(e0.y);
    fez = fabs(e0.z);
    AXISTEST_X01(e0.z, e0.y, fez, fey);
    AXISTEST_Y02(e0.z, e0.x, fez, fex);
    AXISTEST_Z12(e0.y, e0.x, fey, fex);
    
    fex = fabs(e1.x);
    fey = fabs(e1.y);
    fez = fabs(e1.z);
    AXISTEST_X01(e1.z, e1.y, fez, fey);
    AXISTEST_Y02(e1.z, e1.x, fez, fex);
    AXISTEST_Z0(e1.y, e1.x, fey, fex);
    
    fex = fabs(e2.x);
    fey = fabs(e2.y);
    fez = fabs(e2.z);
    AXISTEST_X2(e2.z, e2.y, fez, fey);
    AXISTEST_Y1(e2.z, e2.x, fez, fex);
    AXISTEST_Z12(e2.y, e2.x, fey, fex);
    
    //Bullet 1
    findMinMax(v0.x, v1.x, v2.x, fmin, fmax);
    if(fmin > boxhalfsize.x || fmax < -boxhalfsize.x) {
        return false;
    }
    
    findMinMax(v0.y, v1.y, v2.y, fmin, fmax);
    if(fmin > boxhalfsize.y || fmax < -boxhalfsize.y) {
        return false;
    }
    
    findMinMax(v0.z, v1.z, v2.z, fmin, fmax);
    if(fmin > boxhalfsize.z || fmax < -boxhalfsize.z) {
        return false;
    }
    
    //Bullet 2
    normal = Vector::cross(e0, e1);
    Plane pl;
    pl.point = v0;
    pl.normal = normal;
    if(!isInBox(pl)) {
        return false;
    }
    
    return true;
}

//Also From Triangle-Box Testing
bool AABB::isInBox(const Plane& pln) const {
	Vector normal = pln.normal, vert = pln.point, maxbox = max;
	Vector vmin, vmax; 
    float v;
    
    v = vert.x;
    if(normal.x>0.0f) {
      vmin.x = -maxbox.x - v;
      vmax.x = maxbox.x - v;
    }
    else {
      vmin.x = maxbox.x - v;
      vmax.x = -maxbox.x - v;
    }
    
    if(normal.y>0.0f) {
      vmin.y = -maxbox.y - v;
      vmax.y = maxbox.y - v;
    }
    else {
      vmin.y = maxbox.y - v;
      vmax.y = -maxbox.y - v;
    }
    
    if(normal.z>0.0f) {
      vmin.z = -maxbox.z - v;
      vmax.z = maxbox.z - v;
    }
    else {
      vmin.z = maxbox.z - v;
      vmax.z = -maxbox.z - v;
    }
    
    if(Vector::dot(normal, vmin) > 0.0f) {
        return false;
    }
    if(Vector::dot(normal, vmax) >= 0.0f) {
        return true;
    }
    return false;
}

void AABB::findMinMax(float x0, float x1, float x2, float& min, float& max) const {
	min = max = x0;
	if(x1<min) min=x1;
	if(x1>max) max=x1;
	if(x2<min) min=x2;
	if(x2>max) max=x2;
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
    
    printf("Lights:\n");
    num = scn->getLights().size();
    if (num > 0) {
        for (int i = 0; i < num; i++) {
            printLight(scn->getLights()[i]);
        }
    }
    else {
        printf("No Lights\n");
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
    printf("Type: ");
    switch(l.getType()) {
        case DIRECTIONAL:
            printf("Directional\n");
            printf("Direction: ");
            printVector(l.getDirection());
            break;
        case POINT:
            printf("Point\n");
            printf("Position: ");
            printVector(l.getPosition());
            break;
        case SPOT:
            printf("Spot\n");
            printf("Position: ");
            printVector(l.getPosition());
            printf("Direction: ");
            printVector(l.getDirection());
            printf("Angle: (%f %f)\n", l.getSpotAngle(), l.getMaxAngle());
            break;
        default:
            printf("Ambient\n");
    }
    printf("Color: ");
    printVector(l.getColor());
}

void printImage(rt::Image img) {
    printf("File Name: %s\n", img.getFileName().c_str());
    printf("Width: %d\n", img.getWidth());
    printf("Height: %d\n", img.getHeight());
}

void printVector(Vector fVec) {
    printf("(%f %f %f)\n", fVec.x, fVec.y, fVec.z);
}
