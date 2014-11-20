#include "parse.h"
#include <cstdio>
#include <cstring>
#include <fstream>

#define MAX_LINE 100

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
            float x, y, z, r;
            scn >> x >> y >> z >> r;
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
    scn.close();
    return 1;
}

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
