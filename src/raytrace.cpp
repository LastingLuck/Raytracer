#include "raytrace.h"
#include <cmath>
#include <limits>
#include <cstdlib>
#include <ctime>

/*======================== X-tests ========================*/

#define AXISTEST_X01(a, b, fa, fb)			   \
	p0 = a*v0.y - b*v0.z;			       	   \
	p2 = a*v2.y - b*v2.z;			       	   \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
	rad = fa * boxhalfsize.y + fb * boxhalfsize.z;   \
	if(min>rad || max<-rad) return false;

#define AXISTEST_X2(a, b, fa, fb)			   \
	p0 = a*v0.y - b*v0.z;			           \
	p1 = a*v1.y - b*v1.z;			       	   \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * boxhalfsize.y + fb * boxhalfsize.z;   \
	if(min>rad || max<-rad) return false;

/*======================== Y-tests ========================*/

#define AXISTEST_Y02(a, b, fa, fb)			   \
	p0 = -a*v0.x + b*v0.z;		      	   \
	p2 = -a*v2.x + b*v2.z;	       	       	   \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
	rad = fa * boxhalfsize.x + fb * boxhalfsize.z;   \
	if(min>rad || max<-rad) return false;

#define AXISTEST_Y1(a, b, fa, fb)			   \
	p0 = -a*v0.x + b*v0.z;		      	   \
	p1 = -a*v1.x + b*v1.z;	     	       	   \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * boxhalfsize.x + fb * boxhalfsize.z;   \
	if(min>rad || max<-rad) return false;

/*======================== Z-tests ========================*/

#define AXISTEST_Z12(a, b, fa, fb)			   \
	p1 = a*v1.x - b*v1.y;			           \
	p2 = a*v2.x - b*v2.y;			       	   \
        if(p2<p1) {min=p2; max=p1;} else {min=p1; max=p2;} \
	rad = fa * boxhalfsize.x + fb * boxhalfsize.y;   \
	if(min>rad || max<-rad) return false;

#define AXISTEST_Z0(a, b, fa, fb)			   \
	p0 = a*v0.x - b*v0.y;				   \
	p1 = a*v1.x - b*v1.y;			           \
    if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * boxhalfsize.x + fb * boxhalfsize.y;   \
	if(min>rad || max<-rad) return false;


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
Image* RayTrace::rayTrace(const Scene& scn) {
	#ifdef DEBUG
    printf("Beginning Trace\n");
    #endif
    bool useBVH = false;
    //#ifdef DEBUG
    //Vector botleft, topright;
    //findBoundingVerts(scn, botleft, topright);
    //printf("Bounding Box:\n");
    //printf("Bottom Left Close: (%f %f %f)\n", botleft.x, botleft.y, botleft.z);
    //printf("Top Right Far: (%f %f %f)\n", topright.x, topright.y, topright.z);
    //printf("\n");
    //Box* rootbox = new Box;
	//rootbox->min = add(botleft, -.5);
	//rootbox->max = add(topright, .5);
	//scn->bvh = rootbox;
	//makeBVH(scn, scn->bvh, 0);
	//printBVH(scn->bvh, 0);
	//useBVH = true;
    //#endif
    /*
    if(scc->objNum > scn->bvhthresh) {
		Box* rootbox = new Box;
		rootbox->min = botleft;
		rootbox->max = topright;
		useBVH = true;
		makeBVH(scn, scn->bvh, 0);
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
                //for(int k = 0; k < sampleNm; k++) {
                ////    printf("Color %d: (%f %f %f)\n", k, color.x, color.y,
                //        trace.p.z, trace.d.x, trace.d.y, trace.d.z);
                //}
                //printf("Ave: (%f %f %f)");
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
    }
    printf("100%%\n");
    return dest;
}

Vector RayTrace::evaluateRayTree(const Scene& scn, const Ray& ray, int depth, bool useBVH) const {
    if(depth > scn.getDepth()) {
        //return scn->BGColor;
        return Vector();
    }
    Intersect* in;
    //printf("Dir2: (%f %f %f)\n", ray->d.x, ray->d.y, ray->d.z);
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
		if(min != 0) {
			//printf("Dir3: (%f %f %f)\n", min->r.d.x, min->r.d.y, min->r.d.z);
		}
	}
	else {
		//min = intersectBVH(trace, scn.getBVH(), dmin, dmax);
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
			//Vector point = add(t.vertices[0], add(multiply(sub(t.vertices[1], t.vertices[0]),beta), 
			//					multiply(sub(t.vertices[2], t.vertices[0]),gamma)));
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
					//normal = norm(add(multiply(norm1, beta), add(multiply(norm2, gamma), multiply(norm3, alpha))));
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
					//normal = norm(add(multiply(norm1, beta), add(multiply(norm2, gamma), multiply(norm3, alpha))));
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
    //printf("Dir1: (%f %f %f)\n", e.d.x, e.d.y, e.d.z);
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
	//Vector point = add(r.p, multiply(r.d, i->t));
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
			// 1 - (alpha-angle1 / angle2-angle1)
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
        // = sub(multiply(normal, 2.0*dot(normal, irdir)), irdir);
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
			//normal = multiply(normal, -1);
		}
		// (nr*dot(N,I)-sqrt(1-nr^2(1-dot(N,I)^2)))*N - nr*I
		// nr*I + (nr*dot(I,N)-sqrt(1-(nr*nr)*(1-dot(I,N)^2)))*N
		//Total internal refraction if sqrt < 0
		double tir = 1.0 - (ior*ior) * (1.0 - (dni*dni));
		if(tir >= 0) {
			//Vector refdir = add(multiply(rdir, ior), multiply(normal, (ior*dni)-sqrt(tir)));
			//Vector refdir = sub(multiply(normal, ior*dni-sqrt(tir)), multiply(irdir, ior));
			Vector refdir;
			if(dni >= 0) {
				//refdir = sub(multiply(normal, ior*dni-sqrt(tir)), multiply(irdir, ior));
                refdir = (normal * (ior*dni-sqrt(tir))) - (irdir * ior);
			}
			else {
				//refdir = sub(multiply(normal, ior*dni+sqrt(tir)), multiply(irdir, ior));
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
/*
void BVH::make(Scene& scn, AABB* box, int depth) {
	//If depth is 0, start by making first box
		//Initialize values so that box starts with every object
	if(depth == 0) {
		//box = new Box;
		box->parent = 0;
		box->sub1 = 0;
		box->sub2 = 0;
		//findBoundingVerts(scn, box->min, box->max);
		box->spheres = scn->spheres;
		box->triangles = scn->triangles;
		box->planes = scn->planes;
	}
	//Else find Longest dimention and split along that
		//Loop through each object in previous box
		//Figure out which box they belong in, add it
	Vector bmin = box->min;
	Vector bmax = box->max;
	//char axis = findLongestAxis(vrts);
	char axis = findLongestAxis(bmin, bmax);
	Vector w = sub(bmax, bmin);
	if(depth == 1) {
		//printf("Axis: %c\n", axis);
		//printf("W: (%f %f %f)\n", w.x, w.y, w.z);
	}
	Box* s1 = new Box;
	s1->parent = box;
	s1->sub1 = 0;
	s1->sub2 = 0;
	s1->objNum = 0;
	Box* s2 = new Box;
	s2->parent = box;
	s2->sub1 = 0;
	s2->sub2 = 0;
	s2->objNum = 0;
	float hp;
	switch(axis) {
		case 0: //X
			hp = w.x / 2.0;
			///Sub Box 1 (Left)
			s1->min = bmin;
			s1->max = (Vector){bmax.x-hp, bmax.y, bmax.z};
			///Sub Box 2 (Right)
			s2->min = (Vector){bmin.x+hp, bmin.y, bmin.z};
			s2->max = bmax;
			break;
		case 1: //Y
			hp = w.y / 2.0;
			///Sub Box 1 (Top)
			s1->min = (Vector){bmin.x, bmin.y+hp, bmin.z};
			s1->max = bmax;
			///Sub Box 2 (Bottom)
			s2->min = bmin;
			s2->max = (Vector){bmax.x, bmax.y-hp, bmax.z};
			break;
		case 2: //Z
			hp = w.z / 2.0;
			///Sub Box 1 (Front)
			s1->min = bmin;
			s1->max = (Vector){bmax.x, bmax.y, bmax.z-hp};
			///Sub Box 2 (Back)
			s2->min = (Vector){bmin.x, bmin.y, bmin.z+hp};
			s2->max = bmax;
			break;
	}
	//Find which box objects are int
	//bool noitems1 = true;
	//bool noitems2 = true;
	int num = scn->spheres.size();
	for(int i = 0; i < num; i++) {
		Sphere sphere = scn->spheres[i];
		if(isSphereInBox(s1, &sphere)) {
			s1->spheres.push_back(sphere);
			s1->objNum++;
			//noitems1 = false;
		}
		if(isSphereInBox(s2, &sphere)) {
			s2->spheres.push_back(sphere);
			s2->objNum++;
			//noitems2 = false;
		}
	}
	num = scn->triangles.size();
	for(int i = 0; i < num; i++) {
		Triangle triangle = scn->triangles[i];
		if(isTriangleInBox(s1, &triangle)) {
			s1->triangles.push_back(triangle);
			s1->objNum++;
			//noitems1 = false;
		}
		if(isTriangleInBox(s2, &triangle)) {
			s2->triangles.push_back(triangle);
			s2->objNum++;
			//noitems2 = false;
		}
	}
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
	
	box->sub1 = s1;
	box->sub2 = s2;
	//If depth > bvhdepth
	//return. Don't do recursive calls
	if(depth == scn->bvhdepth) {
		return;
	}
	//Recursively call this function with the 2 new subboxes
	//Don't recurse on boxes with 0-1 items in it
	if(s1->objNum > 1) {
		makeBVH(scn, box->sub1, depth+1);
	}
	if(s2->objNum > 1) {
		makeBVH(scn, box->sub2, depth+1);
	}
}
 
Intersect* intersectBVH(Ray* trace, Box* bvh, double dmin, double dmax) {
	Intersect* in = 0;
	//Vector min = bhv->min, max - bvh->max;
	//Vector orig = trace->p, dir = trace->d;
	//Check if intersects bvh
		//If subboxes are null, check intersection of objects
		//Else recurse down non-null boxes
	if(intersectRayAABB(trace, bvh, dmin, dmax)) {
		//Leaf. Check intersections
		if((bvh->sub1 == 0) && (bvh->sub2 == 0)) {
			//printf("Leaf Box\n");
			Intersect* tmp = 0;
			tmp = intersectSpheres(trace, bvh->spheres, dmin, dmax);
			in = tmp;
			tmp = intersectTriangle(trace, bvh->triangles, dmin, dmax);
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
		}
		//Recurse Down the subboxes
		else {
			Intersect* s1 = 0, *s2 = 0;
			//printf("Box1\n");
			if(bvh->sub1 != 0 && bvh->sub1->objNum != 0) {
				s1 = intersectBVH(trace, bvh->sub1, dmin, dmax);
			}
			//printf("Box2\n");
			if(bvh->sub2 != 0 && bvh->sub2->objNum != 0) {
				s2 = intersectBVH(trace, bvh->sub2, dmin, dmax);
			}
			in = s1;
			if(in == 0) {
				in = s2;
			}
			else if(s2 != 0) {
				if(in->dist > s2->dist) {
					in->normal = s2->normal;
					in->r = s2->r;
					in->dist = s2->dist;
					in->p = s2->p;
					in->m = s2->m;
				}
			}
		}
	}
	return in;
}

//std::array<Vector, 8> (&verts)
void findBoundingVerts(SceneData* scn, Vector& bl, Vector& tr) {
	Vector max, min;
	bool init = false;
	Sphere s;
	std::vector<Sphere> sph = scn->spheres;
	int num = sph.size();
	for(int i = 0; i < num; i++) {
		s = sph[i];
		Vector pnt = add(s.position, s.radius);
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
		pnt = sub(s.position, s.radius);
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
	std::vector<Triangle> tri = scn->triangles;
	num = tri.size();
	for(int i = 0; i < num; i++) {
		t = tri[i];
		for(int j = 0; j < 3; j++) {
			Vector v = t.vertices[j];
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
    
	Rectangle r;
	std::vector<Rectangle> rec = scn->rectangles;
	num = rec.size();
	for(int i = 0; i < num; i++) {
		r = rec[i];
	}
}

//std::array<Vector, 8> vrts
char findLongestAxis(Vector vmin, Vector vmax) {
	Vector w = sub(vmax, vmin);
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
bool intersectRayAABB(Ray* trace, Box* box, double dmin, double dmax) {
	Vector max = box->max, min = box->min;
	Vector orig = trace->p, dir = trace->d;
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
        distv = add(orig, multiply(dir, tmin));
        dist = lengthsq(distv);
    }
    else if(tmax >= 0) {
        distv = sub(orig, add(orig, multiply(dir, tmax)));
        dist = lengthsq(distv);
    }
    return ((dist < dmax) && (dist > dmin));
}

//http://www.mrtc.mdh.se/projects/3Dgraphics/paperF.pdf
//On Faster Sphere-Box Overlap Testing
bool isSphereInBox(Box* b, Sphere* sph) {
	double d = 0, e, rad = sph->radius;
	Vector min = b->min, max = b->max;
	Vector spos = sph->position;
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
bool isTriangleInBox(Box* b, Triangle* t) {
	Vector bmax = b->max, bmin = b->min;
	Vector v0, v1, v2, e0, e1, e2, normal;
    float fex, fey, fez, min, max, p0, p1, p2, rad;
    Vector boxhalfsize = div(sub(bmax, bmin), 2.0);
    Vector bc = add(bmin, boxhalfsize);
    v0 = sub(t->vertices[0], bc);
    v1 = sub(t->vertices[1], bc);
    v2 = sub(t->vertices[2], bc);
    
    e0 = sub(v1, v0);
    e1 = sub(v2, v1);
    e2 = sub(v0, v2);
    
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
    findMinMax(v0.x, v1.x, v2.x, min, max);
    if(min > boxhalfsize.x || max < -boxhalfsize.x) {
        return false;
    }
    
    findMinMax(v0.y, v1.y, v2.y, min, max);
    if(min > boxhalfsize.y || max < -boxhalfsize.y) {
        return false;
    }
    
    findMinMax(v0.z, v1.z, v2.z, min, max);
    if(min > boxhalfsize.z || max < -boxhalfsize.z) {
        return false;
    }
    
    //Bullet 2
    normal = cross(e0, e1);
    Plane* pl = new Plane;
    pl->point = v0;
    pl->normal = normal;
    if(!isPlaneInBox(b, pl)) {
        return false;
    }
    
    return true;
}

void findMinMax(float x0, float x1, float x2, float& min, float& max) {
	min = max = x0;
	if(x1<min) min=x1;
	if(x1>max) max=x1;
	if(x2<min) min=x2;
	if(x2>max) max=x2;
}

//Also From Triangle-Box Testing
bool isPlaneInBox(Box* b, Plane* pln) {
	Vector normal = pln->normal, vert = pln->point, maxbox = b->max;
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
    
    if(dot(normal, vmin) > 0) {
        return false;
    }
    if(dot(normal, vmax) >= 0) {
        return true;
    }
    return false;
}
*/
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
/*
void printBVH(Box* b, int depth) {
	for(int i = 0; i < depth; i++) {
		printf("--");
	}
	printf("Min(%f %f %f), Max(%f %f %f)\n", b->min.x, b->min.y, b->min.z, b->max.x, b->max.y, b->max.z);
	for(int j = 0; j < (signed)b->spheres.size(); j++) {
		std::vector<Sphere> s = b->spheres;
		for(int i = 0; i < depth; i++) {
			printf("  ");
		}
		printf("Sphere:\n");
		for(int i = 0; i < depth; i++) {
			printf("  ");
		}
		printf("Pos: (%f %f %f)\n", s[j].position.x, s[j].position.y, s[j].position.z);
		for(int i = 0; i < depth; i++) {
			printf("  ");
		}
		printf("Rad: %f\n", s[j].radius);
	}
	for(int j = 0; j < (signed)b->triangles.size(); j++) {
		std::vector<Triangle> t = b->triangles;
		for(int i = 0; i < depth; i++) {
			printf("  ");
		}
		printf("Triangle:\n");
		for(int i = 0; i < depth; i++) {
			printf("  ");
		}
		printf("V1: (%f %f %f)\n", t[j].vertices[0].x, t[j].vertices[0].y, t[j].vertices[0].z);
		for(int i = 0; i < depth; i++) {
			printf("  ");
		}
		printf("V2: (%f %f %f)\n", t[j].vertices[1].x, t[j].vertices[1].y, t[j].vertices[1].z);
		for(int i = 0; i < depth; i++) {
			printf("  ");
		}
		printf("V3: (%f %f %f)\n", t[j].vertices[2].x, t[j].vertices[2].y, t[j].vertices[2].z);
	}
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
	
	if(b->sub1 != 0) {
		printBVH(b->sub1, depth+1);
	}
	if(b->sub2 != 0) {
		printBVH(b->sub2, depth+1);
	}
}
*/
