#include "raytrace.h"
#include <cmath>
#include <limits>
#include <cstdlib>
#include <ctime>

Ray::Ray() {
    p = Vector(); //< init to (0, 0, 0)
    d = Vector();
}

Ray::Ray(const Vector& pos, const Vector& dir) {
    p = pos;
    d = dir;
}

Intersect::Intersect() {
    normal = Vector();
    r = Ray();
    dist = 0;
    p = Vector();
    m = Material();
}

Intersect::Intersect(const Vector& norm, const Ray& ray, double distsq, const Vector& point, const Material& mat) {
    normal = norm;
    r = ray;
    dist = distsq;
    p = point;
    m = mat;
}

RayTrace::RayTrace() {
    std::srand(std::time(0));
    forcebvh = false;
}

RayTrace::RayTrace(bool bvhforce) {
    std::srand(std::time(0));
    forcebvh = bvhforce;
}

//Performs a raytrace on the Scene scn with projection type proj
Image* RayTrace::rayTrace(Scene& scn) {
	#ifdef DEBUG
    printf("Beginning Trace\n");
    #endif
    bool useBVH = true;
    //#ifdef DEBUG
    
    AABB* rootbox = new AABB();
    if(useBVH) {
        scn.setBVHRoot(rootbox);
        bvh->make(scn, rootbox, 0);
        printBVH(rootbox, 0);   
    }
    
    //#endif
    /*
    if(scn.getNumObjects() > scn.getBVHThreshold()) {
		AABB* rootbox = new AABB();
		rootbox->setMin(botleft);
		rootbox->setMax(topright);
		useBVH = true;
		makeBVH(scn, rootbox, 0);
	}
	*/
	ProjType proj = scn.getProjType();
    int width = scn.getImage().getWidth(), height = scn.getImage().getHeight();
    Image* dest = new Image(width, height);
    Camera c = scn.getCamera();
    Vector p1, p2;
    //Distance to the viewing plane
    double planeDist = scn.getPlaneDist();
    getExtremePoints(c, planeDist, width, height, p1, p2);
    #ifdef DEBUG
    printf("P1: (%f %f %f)\n", p1.x, p1.y, p1.z);
    printf("P2: (%f %f %f)\n", p2.x, p2.y, p2.z);
    printf("Plane Dist: %f\n", planeDist);
    #endif
    int totalpix = width * height;
    int tenper = .1 * totalpix;
    int count = 0;
    int sampleNm = scn.getSampleRate();
    bool sample = false;
    if(sampleNm > 1) {
		sample = true;
	}
    Pixel pix;
    #ifdef PARALLEL
    #pragma omp parallel for collapse(2) schedule(static, 125)
    #endif
    for(int x = 0; x < width; x++) {
        for(int y = 0; y < height; y++) {
			if((y + (x * height)) % tenper == 0) {
				printf("%d%%\n", count);
				#ifdef PARALLEL
				#pragma omp atomic
				#endif
				count += 10;
			}
            //Get ray from the viewpoint to (x, y) point
            if(sample) {
                Vector* color = new Vector[sampleNm];
				for(int j = 0; j < sampleNm; j++) {
					Ray trace = getRay(x, y, width, height, p1, p2, planeDist, c, proj, sample);
					color[j] = evaluateRayTree(scn, trace, 0, useBVH);
				}
				Vector col = ave(color, sampleNm);
				pix.SetClamp(col.x*255.0, col.y*255.0, col.z*255.0);
				dest->GetPixel(x, y) = pix;
                delete[] color;
			}
            else {
				Ray trace = getRay(x, y, width, height, p1, p2, planeDist, c, proj, sample);
				Vector col = evaluateRayTree(scn, trace, 0, useBVH);
				pix.SetClamp(col.x*255.0, col.y*255.0, col.z*255.0);
				dest->GetPixel(x, y) = pix;
			}
        }
        //printf("Row %d\n", x);
    }
    printf("100%%\n");
    delete rootbox;
    return dest;
}

Vector RayTrace::evaluateRayTree(const Scene& scn, const Ray& ray, int depth, bool useBVH) const {
    if(depth > scn.getDepth()) {
        return scn.getBGColor();
        //return Vector();
    }
    Intersect* in;
    bool hit = false;
    //Go through each object and check for intersection
    if((in = intersect(ray, scn, .001, std::numeric_limits<double>::infinity(), useBVH))) {
        hit = true;
    }
    if(hit) {
        Vector pixcol = getColor(in, scn, depth, useBVH);
        delete in;
        return pixcol;
    }
    else {
        return scn.getBGColor();
    }
}

//Tests if there is an intersection between a ray and objects
Intersect* RayTrace::intersect(const Ray& trace, const Scene& scn, double dmin, double dmax, bool useBVH) const {
	Intersect* min = 0;					 
	if(!useBVH) {
        //printf("Intersect Regular\n");
		Intersect* sph = intersectSpheres(trace, scn.getSpheres(), dmin, dmax);
		Intersect* tri = intersectTriangle(trace, scn.getTriangles(), dmin, dmax);
		//Intersect* pln = intersectPlane(trace, scn->planes, dmin, dmax);
		if(sph != 0) {
			min = new Intersect(sph->normal, sph->r, sph->dist, sph->p, sph->m);
			delete sph;
			sph = 0;
		}
		if(tri != 0) {
			if(min == 0) {
				min = new Intersect(tri->normal, tri->r, tri->dist, tri->p, tri->m);
				delete tri;
				tri = 0;
			}
			else if(min->dist > tri->dist) {
				min->normal = tri->normal;
				min->r = tri->r;
				min->dist = tri->dist;
				min->p = tri->p;
				min->m = tri->m;
				delete tri;
				tri = 0;
			}
		}
        /*
		if(pln != 0) {
			if(min == 0) {
				//min = pln;
				min = new Intersect(pln->normal, pln->r, pln->dist, pln->p, pln->m);
				delete pln;
				pln = 0;
			}
			else if(min->dist > pln->dist) {
				//delete min;
				//min = pln;
				min->normal = pln->normal;
				min->r = pln->r;
				min->dist = pln->dist;
				min->p = pln->p;
				min->m = pln->m;
				delete pln;
				pln = 0;
			}
		}
        */
	}
	else {
        //printf("Intersecting BVH\n");
		min = intersectBVH(trace, scn.getBVHRoot(), dmin, dmax);
	}
	return min;
}

Intersect* RayTrace::intersectSpheres(const Ray& trace, const std::vector<Sphere>& objList, double dmin, double dmax) const {
    Intersect* in = 0;
    int sphnum = objList.size();
    if(sphnum == 0) {
		return 0;
	}
    Vector posit = trace.getPosition();
    Vector direc = trace.getDirection();
	double dd = Vector::dot(direc, direc);
	for(int i = 0; i < sphnum; i++) {
		Sphere obj = objList[i];
		Vector ec = posit - obj.getPosition();
		double dec = Vector::dot(direc, ec);
		
		double det = (dec*dec) - dd*(Vector::dot(ec, ec) - (obj.getRadius()*obj.getRadius()));
		if (det < 0.0) {
			continue;
		}
		
		double sqrtdet = sqrt(det);
		double t1 = (-dec + sqrtdet) / dd;
		double t2 = (-dec - sqrtdet) / dd;
		//t is the number of rays it takes to get to the object. Make that into
		//a number to check bounds with (distance to intersect)
		Vector p1 = direc * t1;
		Vector p2 = direc * t2;
		double d1 = p1.magnitudeSq();
		double d2 = p2.magnitudeSq();
		if(t1 < t2) {
			if(t1 > 0 && d1 > dmin && d1 < dmax) {
				Vector point = posit + p1;
				if(in == 0) {
					in = new Intersect((point - obj.getPosition()).norm(), trace, d1, point, obj.getMaterial());
				}
				else if(d1 < in->dist) {
					in->r = trace;
					in->dist = d1;
					in->m = obj.getMaterial();
					in->p = point;
					in->normal = (in->p - obj.getPosition()).norm();
				}
			}
			else if(t2 > 0 && d2 > dmin && d2 < dmax) {
				Vector point = posit + p2;
				if(in == 0) {
					in = new Intersect((point - obj.getPosition()).norm(), trace, d2, point, obj.getMaterial());
				}
				else if(d2 < in->dist) {
					in->r = trace;
					in->dist = d2;
					in->m = obj.getMaterial();
					in->p = point;
					in->normal = (in->p - obj.getPosition()).norm();
				}
			}
		}
		else {
			if(t2 > 0 && d2 > dmin && d2 < dmax) {
				Vector point = posit + p2;
				if(in == 0) {
					in = new Intersect((point - obj.getPosition()).norm(), trace, d2, point, obj.getMaterial());
				}
				else if(d2 < in->dist) {
					in->r = trace;
					in->dist = d2;
					in->m = obj.getMaterial();
					in->p = point;
					in->normal = (in->p - obj.getPosition()).norm();
				}
			}
			else if(t1 > 0 && d1 > dmin && d1 < dmax) {
				Vector point = posit + p1;
				if(in == 0) {
					in = new Intersect((point-obj.getPosition()).norm(), trace, d1, point, obj.getMaterial());
				}
				else if(d1 < in->dist) {
					in->r = trace;
					in->dist = d1;
					in->m = obj.getMaterial();
					in->p = point;
					in->normal = (point - obj.getPosition()).norm();
				}
			}
		}
	}
    return in;
}

Intersect* RayTrace::intersectTriangle(const Ray& trace, const std::vector<Triangle>& triList, double dmin, double dmax) const {
    Intersect* in = 0;
    int trinum = triList.size();
    if(trinum == 0) {
		return 0;
	}
    Triangle t;
    Vector rd = trace.getDirection();
    Vector rp = trace.getPosition();
    for(int i = 0; i < trinum; i++) {
		t = triList[i];
		Vector ae = t.getVertex(0) - rp;
		Vector ab = t.getVertex(0) - t.getVertex(1);
		Vector ac = t.getVertex(0) - t.getVertex(2);
		double detA = ab.x*(ac.y*rd.z-rd.y*ac.z) + ab.y*(ac.z*rd.x-rd.z*ac.x) + ab.z*(ac.x*rd.y-ac.y*rd.x);
		double dett = ac.z*(ab.x*ae.y-ab.y*ae.x) + ac.y*(ae.x*ab.z-ab.x*ae.z) + ac.x*(ab.y*ae.z-ab.z*ae.y);
		double tval = dett / detA;
		if(tval > 0) {
			continue;
		}
		double detb = ae.x*(ac.y*rd.z-rd.y*ac.z) + ae.y*(ac.z*rd.x-rd.z*ac.x) + ae.z*(ac.x*rd.y-ac.y*rd.x);
		double detg = rd.z*(ab.x*ae.y-ae.x*ab.y) + rd.y*(ab.z*ae.x-ae.z*ab.x) + rd.x*(ab.y*ae.z-ae.y*ab.z);
		double beta = detb / detA;
		double gamma = detg / detA;
		if(beta >= 0 && beta <= 1 && gamma >= 0  && gamma <= 1 && beta + gamma <= 1) {
            Vector point = t.getVertex(0) + 
                    (((t.getVertex(1) - t.getVertex(0)) * beta) + ((t.getVertex(2) - t.getVertex(0)) * gamma));
			double dist = Vector::lengthSq(point, rp);
			if(dist < dmin || dist > dmax) {
				continue;
			}
			if(in == 0) {
				Vector normal;
				if(t.isNormal()) {
					//Linearly Interpolate
					//multiply each vertex norm by the corresponding barycentric coordinant
					//add and renormalize
					double alpha = 1.0 - (beta + gamma);
					double ang1 = Vector::dot(t.getNormal(0), rd);
					double ang2 = Vector::dot(t.getNormal(1), rd);
					double ang3 = Vector::dot(t.getNormal(2), rd);
					Vector norm1 = (ang1 > 0) ? t.getNormal(0)*-1.0f : t.getNormal(0);
					Vector norm2 = (ang2 > 0) ? t.getNormal(1)*-1.0f : t.getNormal(1);
					Vector norm3 = (ang3 > 0) ? t.getNormal(2)*-1.0f : t.getNormal(2);
                    normal = ((norm1 * beta) + (norm2 * gamma) + (norm3 * alpha)).norm();
				}
				else {
					double ang = Vector::dot(t.getNormal(0), rd);
					if(ang > 0) {
						normal = t.getNormal(0) * -1.0f;
					}
					else {
						normal = t.getNormal(0);
					}
				}
				in = new Intersect(normal, trace, dist, point, t.getMaterial());
			}
			else if(dist < in->dist) {
				Vector normal;
				if(t.isNormal()) {
					//Linearly Interpolate
					double alpha = 1.0 - (beta + gamma);
					double ang1 = Vector::dot(t.getNormal(0), rd);
					double ang2 = Vector::dot(t.getNormal(1), rd);
					double ang3 = Vector::dot(t.getNormal(2), rd);
					Vector norm1 = (ang1 > 0) ? t.getNormal(0)*-1.0f : t.getNormal(0);
					Vector norm2 = (ang2 > 0) ? t.getNormal(1)*-1.0f : t.getNormal(1);
					Vector norm3 = (ang3 > 0) ? t.getNormal(2)*-1.0f : t.getNormal(2);
                    normal = ((norm1 * beta) + (norm2 * gamma) + (norm3 * alpha)).norm();
				}
				else {
					double ang = Vector::dot(t.getNormal(0), rd);
					if(ang > 0) {
						normal = t.getNormal(0) * -1.0f;
					}
					else {
						normal = t.getNormal(0);
					}
				}
				in->r = trace;
				in->dist = dist;
				in->p = point;
				in->m = t.getMaterial();
				in->normal = normal;
			}
		}
	}
    
    return in;
}
/*
Intersect* intersectPlane(Ray* trace, std::vector<Plane> planeList, double dmin, double dmax) {
	Intersect* in = 0;
	int planenum = planeList.size();
	if(planenum == 0) {
		return 0;
	}
	Plane pl;
	Vector rp = trace->p;
	Vector rd = trace->d;
	for(int i = 0; i < planenum; i++) {
		pl = planeList[i];
		Vector normal = pl.normal;
		float dni;
		if((dni = dot(normal, rd)) > 0) {
			normal = pl.inormal;
			dni = dot(normal, rd);
		}
		if(dni != 0) {
			float t = dot(normal, sub(pl.point, rp)) / dni;
			if(t > 0) {
				Vector ipoint = add(rp, multiply(rd, t));
				double dist = lengthsq(sub(ipoint, rp));
				if(dist > dmin && dist < dmax) {
					if(in == 0) {
						in = new Intersect;
						in->normal = normal;
						in->r = *trace;
						in->dist = dist;
						in->p = ipoint;
						in->m = pl.mat;
					}
					else if(dist < in->dist) {
						in->normal = normal;
						in->r = *trace;
						in->dist = dist;
						in->p = ipoint;
						in->m = pl.mat;
					}
				}
			}
		}
	}
	return in;
}
*/
//Returns a ray from the viewpoint to the pixel in position (x, y)
Ray RayTrace::getRay(int x, int y, int w, int h, const Vector& p1, const Vector& p2, double pd, const Camera& c, ProjType proj, bool sample) const {
    Ray e;
    Vector pos = c.getPosition();
    Vector right = (Vector::cross(c.getDirection(), c.getUpDirection())).norm();
    Vector down = (c.getUpDirection() * -1.0f).norm();
    double r1, r2;
    if(sample) {
		r1 = (double)std::rand() / RAND_MAX;
		r2 = (double)std::rand() / RAND_MAX;
	}
	else {
		r1 = r2 = .5;
	}
    if (proj == PERSP) {
        Vector dir;
        //Vector from p1 to p2
        Vector ptop = p2 - p1;
        //Project that vector onto the right vector to get upper right point
        Vector ur = right * Vector::dot(ptop, right);
        ur = ur * ((double)x/w + (r1/w));
        //Project ptop onto down vector to get bottom left point
        Vector bl = down * Vector::dot(ptop, down);
        bl = bl * ((double)y/h + (r2/h));
        //Add x percent of first projection and y percent of second projection
        dir = ur + bl;
        dir = dir + p1;
        e = Ray(pos, dir.norm());
    }
    else {
        Vector p;
        //Vector from p1 to p2
        Vector ptop = p2 - p1;
        //Project that vector onto the right vector to get upper right point
        Vector ur = right * Vector::dot(ptop, right);
        ur = ur * ((double)x/w + (r1/w));
        //Project ptop onto down vector to get bottom left point
        Vector bl = down * Vector::dot(ptop, down);
        bl = bl * ((double)y/h + (r2/h));
        //Add x percent of first projection and y percent of second projection
        p = ur + bl;
        p = p + p1;
        e = Ray(p, c.getDirection());
    }
    return e;
}

//Returns a pixel with the background color
Vector RayTrace::getColor(const Intersect* i, const Scene& scn, int depth, bool useBVH) const {
	std::vector<Light> d = scn.getDirLights();
	std::vector<Light> p = scn.getPointLights();
	std::vector<Light> s = scn.getSpotLights();
	Light a = scn.getAmbLight();
	int lnum[3] = {(signed)d.size(), (signed)p.size(), (signed)s.size()};
	Material m = i->m;
	Vector dif = m.getDiffuse();
    Vector spec = m.getSpecular();
	float rval, gval, bval;
	//Get the base color(La)
	Vector aColor = m.getAmbient();
	rval = aColor.x;
	gval = aColor.y;
	bval = aColor.z;
	//Figure out diffuse light using Lambertian shading (Ld)
	//Point that the ray intersected the sphere
	Vector point = i->p;
	Vector direct = i->r.getDirection();
	//Normal norm(P-C)
	Vector normal = i->normal;
	//Multiply material ambient color by the light ambient color
	rval *= a.getColor().x;
	gval *= a.getColor().y;
	bval *= a.getColor().z;
	//Loop through each lightsource
	Ray shadow;
    shadow.setPosition(point);
	Vector I;
	Light L;
    Intersect* in = 0;
	//Directional
	for(int j = 0; j < lnum[0]; j++) {
		L = d[j];
		//Lambert
		//If the ray from the point in the direction of the light intersects
		//anything, skip the rest
		I = L.getDirection() * -1.0;
		I = I.norm();
		shadow.setDirection(I);
		//Use really small number in intersect so that it doesn't intersect itself
		if((in = intersect(shadow, scn, .001, std::numeric_limits<double>::infinity(), useBVH))) {
            delete in;
			continue;
		}
		double dotni = Vector::dot(normal, I);
		double cosAlpha = (dotni > 0) ? dotni : 0;
		Vector illum = L.getColor();
		rval += dif.x * illum.x * cosAlpha;
		gval += dif.y * illum.y * cosAlpha;
		bval += dif.z * illum.z * cosAlpha;
		//Phong
		Vector ref = (normal * (2.0*dotni)) - I;
		Vector V = direct * -1.0f;
        double powval = Vector::dot(V.norm(), ref.norm());
        double pspec;
        if(powval < 0.0) {
            pspec = 0.0;
        }
        else {
            pspec = pow(powval, m.getCosPower());
        }
		rval += spec.x * pspec * illum.x;
		gval += spec.y * pspec * illum.y;
		bval += spec.z * pspec * illum.z;
	}
	//Point
	for(int k = 0; k < lnum[1]; k++) {
		L = p[k];
		I = L.getPosition() - point;
		//Lambert
		double dist = I.magnitudeSq();
		I = I.norm();
		shadow.setDirection(I);
		if((in = intersect(shadow, scn, .001, dist, useBVH))) {
            delete in;
			continue;
		}
		Vector illum = L.getColor() * (1.0/dist);
		double dotni = Vector::dot(normal, I);
		double cosAlpha = (dotni > 0) ? dotni : 0;
		rval += dif.x * illum.x * cosAlpha;
		gval += dif.y * illum.y * cosAlpha;
		bval += dif.z * illum.z * cosAlpha;
		//Phong
		Vector ref = (normal * (2.0*dotni)) - I;
		Vector V = direct * -1.0;
		double powval = Vector::dot(V.norm(), ref.norm());
        double pspec;
        if(powval < 0.0) {
            pspec = 0.0;
        }
        else {
            pspec = pow(powval, m.getCosPower());
        }
		rval += spec.x * pspec * illum.x;
		gval += spec.y * pspec * illum.y;
		bval += spec.z * pspec * illum.z;
	}
	//Spot
	for(int l = 0; l < lnum[2]; l++) {
		L = s[l];
		I = L.getPosition() - point;
		//Lambert
		double dist = I.magnitudeSq();
		I = I.norm();
		shadow.setDirection(I);
		if((in = intersect(shadow, scn, .001, dist, useBVH))) {
			delete in;
            continue;
		}
		Vector illum;
		double dotni = Vector::dot(normal, I);
		double cosAlpha = (dotni > 0) ? dotni : 0;
		double alpha = acos(cosAlpha);
		//Acts like a point light
		if(alpha < L.getSpotAngle()) {
			illum = L.getColor() * (1.0/dist);
		}
		//Greater than angle2, light contributes nothing
		else if(alpha > L.getMaxAngle()) {
			illum = Vector();
		}
		//Linearly interpolate
		else {
			//Get amount alpha is between the 2 angles (angle1=1, angle2=0)
			double amt = 1 - ((alpha - L.getSpotAngle()) / (L.getMaxAngle() - L.getSpotAngle()));
			//Multiply light by amount
			illum = L.getColor() * ((1.0/dist)*amt);
		}
		rval += dif.x * illum.x * cosAlpha;
		gval += dif.y * illum.y * cosAlpha;
		bval += dif.z * illum.z * cosAlpha;
		//Phong
		Vector ref = (normal * (2.0*dotni)) - I;
		Vector V = direct * -1.0f;
		double powval = Vector::dot(V.norm(), ref.norm());
        double pspec;
        if(powval < 0.0) {
            pspec = 0.0;
        }
        else {
            pspec = pow(powval, m.getCosPower());
        }
		rval += spec.x * pspec * illum.x;
		gval += spec.y * pspec * illum.y;
		bval += spec.z * pspec * illum.z;
	}
    Vector refColor;
    Vector rdir = direct;
	Vector irdir = rdir * -1.0f;
    ///Reflect Ray
    if(spec.x != 0 && spec.y != 0 && spec.z != 0) {
		Ray reflection;
		reflection.setPosition(point);
		reflection.setDirection((normal * (2.0f*Vector::dot(normal, irdir))) - irdir);
		refColor = evaluateRayTree(scn, reflection, depth+1, useBVH);
		rval += spec.x * refColor.x;
		gval += spec.y * refColor.y;
		bval += spec.z * refColor.z;
	}
    ///Refract
    Vector trans = m.getTransmissive();
    if(trans.x != 0 && trans.y != 0 && trans.z != 0) {
		Ray refraction;
		refraction.setPosition(point);
		double ior;
		double dni = Vector::dot(irdir, normal);
		if(dni <= 0) { //going into object
			ior = m.getIndexRefract();
		}
		else {
			ior = 1.0 / m.getIndexRefract();
		}
		// (nr*dot(N,I)-sqrt(1-nr^2(1-dot(N,I)^2)))*N - nr*I
		// nr*I + (nr*dot(I,N)-sqrt(1-(nr*nr)*(1-dot(I,N)^2)))*N
		//Total internal refraction if sqrt < 0
		double tir = 1.0 - (ior*ior) * (1.0 - (dni*dni));
		if(tir >= 0) {
			Vector refdir;
			if(dni >= 0) {
                refdir = (normal * (ior*dni-sqrt(tir))) - (irdir * ior);
			}
			else {
                refdir = (normal * (ior*dni+sqrt(tir))) - (irdir * ior);
			}
			refraction.setDirection(refdir.norm());
			refColor = evaluateRayTree(scn, refraction, depth+1, useBVH);
			rval += trans.x * refColor.x;
			gval += trans.y * refColor.y;
			bval += trans.z * refColor.z;
		}
	}
    
    return Vector(rval, gval, bval);
}

//Gets the extreme points of top left and bottom right
void RayTrace::getExtremePoints(const Camera& c, double d, double w, double h, Vector& p1, Vector& p2) const {
	//p1 is the top left, p2 is the bottom right
	Vector p0 = c.getPosition();
	Vector up = c.getUpDirection();
	Vector right = Vector::cross(c.getDirection(), up);
	//Vector from the viewpoint to the center
	p0 = p0 + (c.getDirection() * d);
	//Add on the offset for the 2 extreme points (-x, -y) and (+x, +y)
	p1 = p0 + (right * (w/2.0)); //< Add x part
	p2 = p0 - (right * (w/2.0));
	p1 = p1 + (up * (h/2.0)); //< Add y part
	p2 = p2 - (up * (h/2.0));
}

/**
 * Acceleration Structure (BVH)
 */
void BVH::make(Scene& scn, AABB* box, int depth) {
	//If depth is 0, start by making first box
		//Initialize values so that box starts with every object
	if(depth == 0) {
		//box = new AABB(scn->spheres, scn->triangles); //< Already made before this call
        Vector botleft, topright;
        findBoundingVerts(scn, botleft, topright);
        box->setMin(botleft);
        box->setMax(topright);
        box->setSpheres(scn.getSpheres());
        box->setTriangles(scn.getTriangles());
		//box->parent = 0;
		//box->sub1 = 0;
		//box->sub2 = 0;
		//box->spheres = scn->spheres;
		//box->triangles = scn->triangles;
		//box->planes = scn->planes;
	}
    if(box->getObjectNum() < 2) {
        return;
    }
	//Find Longest dimention and split along that
		//Loop through each object in previous box
		//Figure out which box they belong in, add it
	Vector bmin = box->getMin();
	Vector bmax = box->getMax();
	//char axis = findLongestAxis(vrts);
	char axis = findLongestAxis(bmin, bmax);
	Vector w = bmax - bmin;
	if(depth == 1) {
		//printf("Axis: %c\n", axis);
		//printf("W: (%f %f %f)\n", w.x, w.y, w.z);
	}
	AABB* s1 = new AABB(box);
	//s1->parent = box;
	//s1->sub1 = 0;
	//s1->sub2 = 0;
	//s1->objNum = 0;
	AABB* s2 = new AABB(box);
	//s2->parent = box;
	//s2->sub1 = 0;
	//s2->sub2 = 0;
	//s2->objNum = 0;
	float hp;
	switch(axis) {
		case 0: //X
			hp = w.x / 2.0;
			///Sub Box 1 (Left)
			s1->setMin(bmin);
			s1->setMax(Vector(bmax.x-hp, bmax.y, bmax.z));
			///Sub Box 2 (Right)
			s2->setMin(Vector(bmin.x+hp, bmin.y, bmin.z));
			s2->setMax(bmax);
			break;
		case 1: //Y
			hp = w.y / 2.0;
			///Sub Box 1 (Top)
			s1->setMin(Vector(bmin.x, bmin.y+hp, bmin.z));
			s1->setMax(bmax);
			///Sub Box 2 (Bottom)
			s2->setMin(bmin);
			s2->setMax(Vector(bmax.x, bmax.y-hp, bmax.z));
			break;
		case 2: //Z
			hp = w.z / 2.0;
			///Sub Box 1 (Front)
			s1->setMin(bmin);
			s1->setMax(Vector(bmax.x, bmax.y, bmax.z-hp));
			///Sub Box 2 (Back)
			s2->setMin(Vector(bmin.x, bmin.y, bmin.z+hp));
			s2->setMax(bmax);
			break;
	}
	//Find which box objects are in 
	//bool noitems1 = true;
	//bool noitems2 = true;
    std::vector<Sphere> sphs = scn.getSpheres();
	int num = sphs.size();
	for(int i = 0; i < num; i++) {
		Sphere sphere = sphs[i];
		if(s1->isInBox(sphere)) {
			s1->addSphere(sphere);
			//noitems1 = false;
		}
		if(s2->isInBox(sphere)) {
			s2->addSphere(sphere);
			//noitems2 = false;
		}
	}
    std::vector<Triangle> tris = scn.getTriangles();
	num = tris.size();
	for(int i = 0; i < num; i++) {
		Triangle triangle = tris[i];
		if(s1->isInBox(triangle)) {
			s1->addTriangle(triangle);
			//noitems1 = false;
		}
		if(s2->isInBox(triangle)) {
			s2->addTriangle(triangle);
			//noitems2 = false;
		}
	}
    /*
	num = scn->planes.size();
	for(int i = 0; i < num; i++) {
		Plane plane = scn->planes[i];
		if(isPlaneInBox(s1, &plane)) {
			s1->planes.push_back(plane);
			s1->objNum++;
			//noitems1 = false;
		}
		if(isPlaneInBox(s2, &plane)) {
			s2->planes.push_back(plane);
			s2->objNum++;
			//noitems2 = false;
		}
	}
	*/
	box->setLeftChild(s1);
	box->setRightChild(s2);
	//If depth > bvhdepth
	//return. Don't do recursive calls
	if(depth == scn.getBVHDepth()) {
		return;
	}
	//Recursively call this function with the 2 new subboxes
	//Don't recurse on boxes with 0-2 items in it
	if(s1->getObjectNum() > 2) {
		make(scn, s1, depth+1);
	}
	if(s2->getObjectNum() > 2) {
		make(scn, s2, depth+1);
	}
}

Intersect* RayTrace::intersectBVH(const Ray& trace, const AABB* bvh, double dmin, double dmax) const {
	Intersect* in = 0;
	//Vector min = bhv->min, max - bvh->max;
	//Vector orig = trace->p, dir = trace->d;
	//Check if intersects bvh
		//If subboxes are null, check intersection of objects
		//Else recurse down non-null boxes
	if(intersectRayAABB(trace, bvh, dmin, dmax)) {
		//Leaf. Check intersections
		if(bvh->isLeaf()) {
			//printf("Leaf Box\n");
			Intersect* tmpSph = intersectSpheres(trace, bvh->getSpheres(), dmin, dmax);
			Intersect* tmpTri = intersectTriangle(trace, bvh->getTriangles(), dmin, dmax);
			if(tmpSph != 0) {
				in = new Intersect(tmpSph->normal, tmpSph->r, tmpSph->dist, tmpSph->p, tmpSph->m);
                delete tmpSph;
                tmpSph = 0;
			}
            if(tmpTri != 0) {
                if(in == 0) {
                    in = new Intersect(tmpTri->normal, tmpTri->r, tmpTri->dist, tmpTri->p, tmpTri->m);
                    delete tmpTri;
                    tmpTri = 0;
                }
                else if(in->dist > tmpTri->dist) {
                    in->normal = tmpTri->normal;
                    in->r = tmpTri->r;
                    in->dist = tmpTri->dist;
                    in->p = tmpTri->p;
                    in->m = tmpTri->m;
                    delete tmpTri;
                    tmpTri = 0;
                }
            }
            /*
			tmp = intersectPlane(trace, bvh->planes, dmin, dmax);
			if(in == 0) {
				in = tmp;
			}
			else if(tmp != 0 && in->dist > tmp->dist) {
				in->normal = tmp->normal;
				in->r = tmp->r;
				in->dist = tmp->dist;
				in->p = tmp->p;
				in->m = tmp->m;
			}
            */
		}
		//Recurse Down the subboxes
		else {
			Intersect* s1 = 0, *s2 = 0;
			//printf("Box1\n");
			if(bvh->getLeftChild() != 0 && bvh->getLeftChild()->getObjectNum() != 0) {
				s1 = intersectBVH(trace, bvh->getLeftChild(), dmin, dmax);
			}
			//printf("Box2\n");
			if(bvh->getRightChild() != 0 && bvh->getRightChild()->getObjectNum() != 0) {
				s2 = intersectBVH(trace, bvh->getRightChild(), dmin, dmax);
			}
            if(s1 != 0) {
                in = new Intersect(s1->normal, s1->r, s1->dist, s1->p, s1->m);
                delete s1;
                s1 = 0;
            }
            if(s2 != 0) {
                if(in == 0) {
                    in = new Intersect(s2->normal, s2->r, s2->dist, s2->p, s2->m);
                    delete s2;
                    s2 = 0;
                }
                else if(in->dist > s2->dist) {
                    in->normal = s2->normal;
					in->r = s2->r;
					in->dist = s2->dist;
					in->p = s2->p;
					in->m = s2->m;
                    delete s2;
                    s2 = 0;
                }
            }
		}
	}
	return in;
}

void BVH::findBoundingVerts(const Scene& scn, Vector& bl, Vector& tr) {
	Vector max, min;
	bool init = false;
	Sphere s;
	std::vector<Sphere> sph = scn.getSpheres();
	int num = sph.size();
	for(int i = 0; i < num; i++) {
		s = sph[i];
		Vector pnt = s.getPosition() + s.getRadius();
		if(!init) {
			max.x = pnt.x;
			max.y = pnt.y;
			max.z = pnt.z;
			//printf("Max Init: (%f %f %f)\n", pnt.x, pnt.y, pnt.z);
		}
		else {
			if(max.x < pnt.x) {
				max.x = pnt.x;
			}
			if(max.y < pnt.y) {
				max.y = pnt.y;
			}
			if(max.z < pnt.z) {
				max.z = pnt.z;
			}
		}
		pnt = s.getPosition() - s.getRadius();
		if(!init) {
			min.x = pnt.x;
			min.y = pnt.y;
			min.z = pnt.z;
			init = true;
			//printf("Min Init: (%f %f %f)\n", pnt.x, pnt.y, pnt.z);
		}
		else {
			if(min.x > pnt.x) {
				min.x = pnt.x;
			}
			if(min.y > pnt.y) {
				min.y = pnt.y;
			}
			if(min.z > pnt.z) {
				min.z = pnt.z;
			}
		}
	}
	
	Triangle t;
	std::vector<Triangle> tri = scn.getTriangles();
	num = tri.size();
	for(int i = 0; i < num; i++) {
		t = tri[i];
		for(int j = 0; j < 3; j++) {
			Vector v = t.getVertex(j);
			if(!init) {
				max = v;
				min = v;
				init = true;
			}
			else {
				if(max.x < v.x) {
					max.x = v.x;
				}
				if(max.y < v.y) {
					max.y = v.y;
				}
				if(max.z < v.z) {
					max.z = v.z;
				}
				if(min.x > v.x) {
					min.x = v.x;
				}
				if(min.y > v.y) {
					min.y = v.y;
				}
				if(min.z > v.z) {
					min.z = v.z;
				}
			}
		}	
	}
	bl = min;
	tr = max;
    
	//Rectangle r;
	//std::vector<Rectangle> rec = scn->rectangles;
	//num = rec.size();
	//for(int i = 0; i < num; i++) {
		//r = rec[i];
	//}
}

//std::array<Vector, 8> vrts
char BVH::findLongestAxis(const Vector& vmin, const Vector& vmax) {
	Vector w = vmax - vmin;
	w = (Vector){(float)fabs(w.x), (float)fabs(w.y), (float)fabs(w.z)};
	if(w.x >= w.y && w.x >= w.z) {
		return 0;
	}
	if(w.y >= w.x && w.y >= w.z) {
		return 1;
	}
	if(w.z >= w.x && w.z >= w.y) {
		return 2;
	}
	return -1;
}

//An Efficient and Robust Ray-Box Intersection Algorithm
bool RayTrace::intersectRayAABB(const Ray& trace, const AABB* box, double dmin, double dmax) const {
	Vector max = box->getMax(), min = box->getMin();
	Vector orig = trace.getPosition(), dir = trace.getDirection();
	if((orig.x >= min.x && orig.x <= max.x) && (orig.y >= min.y && orig.y <= max.y) && 
       (orig.z >= min.z && orig.z <= max.z)) {
        return true;
    }
    float tmin, tmax, tymin, tymax, tzmin, tzmax, div;
    ///X
    div = 1 / dir.x;
    if(div >= 0) {
        tmin = (min.x - orig.x) * div;
        tmax = (max.x - orig.x) * div;
    }
    else {
        tmin = (max.x - orig.x) * div;
        tmax = (min.x - orig.x) * div;
    }
    ///Y
    div = 1 / dir.y;
    if(div >= 0) {
        tymin = (min.y - orig.y) * div;
        tymax = (max.y - orig.y) * div;
    }
    else {
        tymin = (max.y - orig.y) * div;
        tymax = (min.y - orig.y) * div;
    }
    if((tmin > tymax) || (tymin > tmax)) {
        return false;
    }
    if(tymin > tmin) {
        tmin = tymin;
    }
    if(tymax < tmax) {
        tmax = tymax;
    }
    ///Z
    div = 1 / dir.z;
    if(div >= 0) {
        tzmin = (min.z - orig.z) * div;
        tzmax = (max.z - orig.z) * div;
    }
    else {
        tzmin = (max.z - orig.z) * div;
        tzmax = (min.z - orig.z) * div;
    }
    if((tmin > tzmax) || (tzmin > tmax)) {
        return false;
    }
    if(tzmin > tmin) {
        tmin = tzmin;
    }
    if(tzmax < tmax) {
        tmax = tzmax;
    }
	double dist;
    Vector distv;
    if(tmin >= 0) {
        distv = orig + (dir * tmin);
        dist = distv.magnitudeSq();
    }
    else if(tmax >= 0) {
        distv = orig - (orig + (dir * tmax));
        dist = distv.magnitudeSq();
    }
    return ((dist < dmax) && (dist > dmin));
}

/**
* Vector Operations
*/

Vector RayTrace::ave(Vector v[], int num) {
	if(num == 1) {
		return v[0];
	}
	Vector u = v[0];
	for(int i = 1; i < num; i++) {
		u = u + v[i];
	}
	u = u / num;
	return u;
}

void printBVH(AABB* b, int depth) {
	for(int i = 0; i < depth; i++) {
		printf("--");
	}
    Vector min = b->getMin();
    Vector max = b->getMax();
    std::vector<Sphere> sphs = b->getSpheres();
	printf("Min(%f %f %f), Max(%f %f %f)\n", min.x, min.y, min.z, max.x, max.y, max.z);
	for(int j = 0; j < (signed)sphs.size(); j++) {
		for(int i = 0; i < depth; i++) {
			printf("  ");
		}
		printf("Sphere:\n");
		for(int i = 0; i < depth; i++) {
			printf("  ");
		}
		printf("Pos: (%f %f %f)\n", sphs[j].getPosition().x, sphs[j].getPosition().y, sphs[j].getPosition().z);
		for(int i = 0; i < depth; i++) {
			printf("  ");
		}
		printf("Rad: %f\n", sphs[j].getRadius());
	}
    std::vector<Triangle> t = b->getTriangles();
	for(int j = 0; j < (signed)t.size(); j++) {
		for(int i = 0; i < depth; i++) {
			printf("  ");
		}
		printf("Triangle:\n");
		for(int i = 0; i < depth; i++) {
			printf("  ");
		}
		printf("V1: (%f %f %f)\n", t[j].getVertices()[0].x, t[j].getVertices()[0].y, t[j].getVertices()[0].z);
		for(int i = 0; i < depth; i++) {
			printf("  ");
		}
		printf("V2: (%f %f %f)\n", t[j].getVertices()[1].x, t[j].getVertices()[1].y, t[j].getVertices()[1].z);
		for(int i = 0; i < depth; i++) {
			printf("  ");
		}
		printf("V3: (%f %f %f)\n", t[j].getVertices()[2].x, t[j].getVertices()[2].y, t[j].getVertices()[2].z);
	}
    /*
	for(int j = 0; j < (signed)b->planes.size(); j++) {
		std::vector<Plane> p = b->planes;
		for(int i = 0; i < depth; i++) {
			printf("  ");
		}
		printf("Plane:\n");
		for(int i = 0; i < depth; i++) {
			printf("  ");
		}
		printf("Point: (%f %f %f)\n", p[j].point.x, p[j].point.y, p[j].point.z);
		for(int i = 0; i < depth; i++) {
			printf("  ");
		}
		printf("Normal: (%f %f %f)\n", p[j].normal.x, p[j].normal.y, p[j].normal.z);
	}
	*/
	if(b->getLeftChild() != 0) {
		printBVH(b->getLeftChild(), depth+1);
	}
	if(b->getRightChild() != 0) {
		printBVH(b->getRightChild(), depth+1);
	}
}

