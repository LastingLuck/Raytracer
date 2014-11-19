#include "parse.h"
#include "raytrace.h"
#include <cstdio>
#include <cstring>

#define MAX_LINE 100

Scene::Scene() {
    init();
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
    bvhthresh = 500;
    bvhdepth = 5;
    useBVH = false;
    //Super sampling
    sampleNum = 1;
    #ifdef DEBUG
    printf("Setup Done\n");
    #endif
}

int Scene::parseScene(char* file) {
    FILE* scn;
    if ((scn = fopen(file, "r")) == 0) {
        printf("File could not be opened: %s\n", file);
        return -1;
    }
    Material curMat = (Material){(Vector){0, 0, 0}, (Vector){1, 1, 1}, (Vector){0, 0, 0}, 5, (Vector){0, 0, 0}};
    char* line = new char[MAX_LINE];
    while (std::fgets(line, MAX_LINE, scn)) {
        //Parse the line
        if (line == 0 || *line == '#') {
            continue;
        }
        char* tok = strtok(line, " ");
        if (strcmp(tok, "camera") == 0) {
            tok = strtok(0, " ");
            float x = strtod(tok, 0);
            tok = strtok(0, " ");
            float y = strtod(tok, 0);
            tok = strtok(0, " ");
            float z = strtod(tok, 0);
            scene->camera.position = (Vector){x, y, z};
            tok = strtok(0, " ");
            x = strtod(tok, 0);
            tok = strtok(0, " ");
            y = strtod(tok, 0);
            tok = strtok(0, " ");
            z = strtod(tok, 0);
            scene->camera.direction = (Vector){x, y, z};
            tok = strtok(0, " ");
            x = strtod(tok, 0);
            tok = strtok(0, " ");
            y = strtod(tok, 0);
            tok = strtok(0, " ");
            z = strtod(tok, 0);
            scene->camera.cameraUp = (Vector){x, y, z};
            tok = strtok(0, " ");
            //Convert degrees to radians
            scene->camera.cameraHa = strtod(tok, 0) * (M_PI / 180.0);
        }
        else if (strcmp(tok, "film_resolution") == 0) {
            tok = strtok(0, " ");
            int x = strtol(tok, 0, 10);
            scene->file.width = x;
            tok = strtok(0, " ");
            x = strtol(tok, 0, 10);
            scene->file.height = x;
        }
        else if (strcmp(tok, "output_image") == 0) {
            tok = strtok(0, " ");
            scene->file.filename = std::string(tok);
        }
        else if (strcmp(tok, "sphere") == 0) {
            tok = strtok(0, " ");
            float x = strtod(tok, 0);
            tok = strtok(0, " ");
            float y = strtod(tok, 0);
            tok = strtok(0, " ");
            float z = strtod(tok, 0);
            tok = strtok(0, " ");
            float r = strtod(tok, 0);
            Sphere s = (Sphere){(Vector){x, y, z}, r, curMat};
            //Add the new sphere with the current material
            scene->spheres.push_back(s);
            scene->objNum++;
        }
        else if (strcmp(tok, "background") == 0) {
            tok = strtok(0, " ");
            float r = strtod(tok, 0);
            tok = strtok(0, " ");
            float g = strtod(tok, 0);
            tok = strtok(0, " ");
            float b = strtod(tok, 0);
            scene->BGColor = (Vector){r, g, b};
        }
        else if (strcmp(tok, "material") == 0) {
            tok = strtok(0, " ");
            float r = strtod(tok, 0);
            tok = strtok(0, " ");
            float g = strtod(tok, 0);
            tok = strtok(0, " ");
            float b = strtod(tok, 0);
            curMat.ambientColor = (Vector){r, g, b};
            tok = strtok(0, " ");
            r = strtod(tok, 0);
            tok = strtok(0, " ");
            g = strtod(tok, 0);
            tok = strtok(0, " ");
            b = strtod(tok, 0);
            curMat.diffuseColor = (Vector){r, g, b};
            tok = strtok(0, " ");
            r = strtod(tok, 0);
            tok = strtok(0, " ");
            g = strtod(tok, 0);
            tok = strtok(0, " ");
            b = strtod(tok, 0);
            curMat.specularColor = (Vector){r, g, b};
            tok = strtok(0, " ");
            curMat.cosPow = strtod(tok, 0);
            tok = strtok(0, " ");
            r = strtod(tok, 0);
            tok = strtok(0, " ");
            g = strtod(tok, 0);
            tok = strtok(0, " ");
            b = strtod(tok, 0);
            curMat.transmissiveColor = (Vector){r, g, b};
            tok = strtok(0, " ");
            curMat.ior = strtod(tok, 0);
        }
        else if (strcmp(tok, "directional_light") == 0) {
            tok = strtok(0, " ");
            float r = strtod(tok, 0);
            tok = strtok(0, " ");
            float g = strtod(tok, 0);
            tok = strtok(0, " ");
            float b = strtod(tok, 0);
            tok = strtok(0, " ");
            float x = strtod(tok, 0);
            tok = strtok(0, " ");
            float y = strtod(tok, 0);
            tok = strtok(0, " ");
            float z = strtod(tok, 0);
            Light d = (Light){(Vector){r, g, b}, (Vector){0, 0, 0}, (Vector){x, y, z}, 0, 0};
            scene->directional.push_back(d);
        }
        else if (strcmp(tok, "point_light") == 0) {
            tok = strtok(0, " ");
            float r = strtod(tok, 0);
            tok = strtok(0, " ");
            float g = strtod(tok, 0);
            tok = strtok(0, " ");
            float b = strtod(tok, 0);
            tok = strtok(0, " ");
            float x = strtod(tok, 0);
            tok = strtok(0, " ");
            float y = strtod(tok, 0);
            tok = strtok(0, " ");
            float z = strtod(tok, 0);
            Light p = (Light){(Vector){r, g, b}, (Vector){x, y, z}, (Vector){0, 0, 0}, 0, 0};
            scene->point.push_back(p);
        }
        else if (strcmp(tok, "spot_light") == 0) {
            tok = strtok(0, " ");
            float x = strtod(tok, 0);
            tok = strtok(0, " ");
            float y = strtod(tok, 0);
            tok = strtok(0, " ");
            float z = strtod(tok, 0);
            Vector intense = (Vector){x, y, z};
            tok = strtok(0, " ");
            x = strtod(tok, 0);
            tok = strtok(0, " ");
            y = strtod(tok, 0);
            tok = strtok(0, " ");
            z = strtod(tok, 0);
            Vector pos = (Vector){x, y, z};
            tok = strtok(0, " ");
            x = strtod(tok, 0);
            tok = strtok(0, " ");
            y = strtod(tok, 0);
            tok = strtok(0, " ");
            z = strtod(tok, 0);
            Vector dir = (Vector){x, y, z};
            tok = strtok(0, " ");
            //Convert x and y to radians
            x = strtod(tok, 0) * (M_PI / 180.0);
            tok = strtok(0, " ");
            y = strtod(tok, 0) * (M_PI / 180.0);
            Light s = (Light){intense, pos, dir, x, y};
            scene->spot.push_back(s);
        }
        else if (strcmp(tok, "ambient_light") == 0) {
            tok = strtok(0, " ");
            float r = strtod(tok, 0);
            tok = strtok(0, " ");
            float g = strtod(tok, 0);
            tok = strtok(0, " ");
            float b = strtod(tok, 0);
            scene->ambient = (Light){(Vector){r, g, b}, (Vector){0, 0, 0}, (Vector){0, 0, 0}, 0, 0};
        }
        else if (strcmp(tok, "max_depth") == 0) {
            tok = strtok(0, " ");
            scene->depth = strtol(tok, 0, 10);
        }
        else if(strcmp(tok, "ray_type") == 0) {
			tok = strtok(0, " ");
			if(strcmp(tok, "perspective") == 0) {
				scene->eyeray = PERSP;
			}
			else if(strcmp(tok, "orthographic") == 0) {
				scene->eyeray = ORTHO;
			}
			else {
				printf("Unknown Projection Type: %s\n", tok);
				printf("Defaulting to Perspective\n");
				scene->eyeray = PERSP;
			}
		}
        else if(strcmp(tok, "max_vertices") == 0) {
            tok = strtok(0, " ");
            int max = strtol(tok, 0, 10);
            scene->vertexes.reserve(max);
        }
        else if(strcmp(tok, "max_normals") == 0) {
            tok = strtok(0, " ");
            int max = strtol(tok, 0, 10);
            scene->normals.reserve(max);
        }
        else if(strcmp(tok, "vertex") == 0) {
            if(scene->vertexes.capacity() == 0) {
                printf("Error: max_vertices must be specified before vertexes\n");
                return -1;
            }
            tok = strtok(0, " ");
            float x = strtod(tok, 0);
            tok = strtok(0, " ");
            float y = strtod(tok, 0);
            tok = strtok(0, " ");
            float z = strtod(tok, 0);
            scene->vertexes.push_back((Vector){x, y, z});
        }
        else if(strcmp(tok, "normal") == 0) {
            if(scene->normals.capacity() == 0) {
                printf("Error: max_normals must be specified before normals\n");
                return -1;
            }
            tok = strtok(0, " ");
            float x = strtod(tok, 0);
            tok = strtok(0, " ");
            float y = strtod(tok, 0);
            tok = strtok(0, " ");
            float z = strtod(tok, 0);
            scene->normals.push_back((Vector){x, y, z});
        }
        else if(strcmp(tok, "triangle") == 0) {
            tok = strtok(0, " ");
            int one = strtol(tok, 0, 10);
            if(one > (signed)scene->vertexes.size()-1) {
                printf("Error: Specified vertex (%d) does not exist\n", one);
                return -1;
            }
            tok = strtok(0, " ");
            int two = strtol(tok, 0, 10);
            if(two > (signed)scene->vertexes.size()-1) {
                printf("Error: Specified vertex (%d) does not exist\n", two);
                return -1;
            }
            tok = strtok(0, " ");
            int three = strtol(tok, 0, 10);
            if(three > (signed)scene->vertexes.size()-1) {
                printf("Error: Specified vertex (%d) does not exist\n", three);
                return -1;
            }
            Triangle t;
            t.vertices[0] = scene->vertexes[one];
            t.vertices[1] = scene->vertexes[two];
            t.vertices[2] = scene->vertexes[three];
            t.normals[0] = norm(cross(sub(t.vertices[1], t.vertices[0]), sub(t.vertices[2], t.vertices[0])));
            t.normals[2] = t.normals[1] = t.normals[0];
            t.normals[3] = t.normals[4] = t.normals[5] = multiply(t.normals[0], -1);
            t.ntri = false;
            t.mat = curMat;
            scene->triangles.push_back(t);
            scene->objNum++;
        }
        else if(strcmp(tok, "normal_triangle") == 0) {
            tok = strtok(0, " ");
            int x = strtol(tok, 0, 10);
            if(x > (signed)scene->vertexes.size()-1) {
                printf("Error: Specified vertex (%d) does not exist\n", x);
                return -1;
            }
            tok = strtok(0, " ");
            int y = strtol(tok, 0, 10);
            if(y > (signed)scene->vertexes.size()-1) {
                printf("Error: Specified vertex (%d) does not exist\n", y);
                return -1;
            }
            tok = strtok(0, " ");
            int z = strtol(tok, 0, 10);
            if(z > (signed)scene->vertexes.size()-1) {
                printf("Error: Specified vertex (%d) does not exist\n", z);
                return -1;
            }
            tok = strtok(0, " ");
            int nx = strtol(tok, 0, 10);
            if(nx > (signed)scene->normals.size()-1) {
                printf("Error: Specified normal (%d) does not exist\n", nx);
                return -1;
            }
            tok = strtok(0, " ");
            int ny = strtol(tok, 0, 10);
            if(ny > (signed)scene->normals.size()-1) {
                printf("Error: Specified normal (%d) does not exist\n", ny);
                return -1;
            }
            tok = strtok(0, " ");
            int nz = strtol(tok, 0, 10);
            if(nz > (signed)scene->normals.size()-1) {
                printf("Error: Specified normal (%d) does not exist\n", nz);
                return -1;
            }
            Triangle t;
            t.vertices[0] = scene->vertexes[x];
            t.vertices[1] = scene->vertexes[y];
            t.vertices[2] = scene->vertexes[z];
            t.normals[0] = scene->normals[nx];
            t.normals[1] = scene->normals[ny];
            t.normals[2] = scene->normals[nz];
            t.normals[3] = multiply(t.normals[0], -1);
            t.normals[4] = multiply(t.normals[1], -1);
            t.normals[5] = multiply(t.normals[2], -1);
            t.ntri = true;
            t.mat = curMat;
            scene->triangles.push_back(t);
            scene->objNum++;
        }
        else if(strcmp(tok, "plane") == 0) {
			tok = strtok(0, " ");
			float x = strtod(tok, 0);
			tok = strtok(0, " ");
			float y = strtod(tok, 0);
			tok = strtok(0, " ");
			float z = strtod(tok, 0);
			tok = strtok(0, " ");
			float nx = strtod(tok, 0);
			tok = strtok(0, " ");
			float ny = strtod(tok, 0);
			tok = strtok(0, " ");
			float nz = strtod(tok, 0);
			Plane pl;
			pl.point = (Vector){x, y, z};
			pl.normal = norm((Vector){nx, ny, nz});
			pl.inormal = multiply(pl.normal, -1);
			pl.mat = curMat;
			scene->planes.push_back(pl);
			scene->objNum++;
		}
		else if(strcmp(tok, "bvh_threshold") == 0) {
			tok = strtok(0, " ");
			int thresh = strtol(tok, 0, 10);
			scene->bvhthresh = thresh;
		}
		else if(strcmp(tok, "bvh_depth") == 0) {
			tok = strtok(0, " ");
			int depth = strtol(tok, 0, 10);
			scene->bvhdepth = depth;
		}
		else if(strcmp(tok, "sample_rate") == 0) {
			tok = strtok(0, " ");
			int sample = strtol(tok, 0, 10);
			scene->sampleNum = sample;
		}
    }
    fclose(scn);
    return 1;
}

void printScene(Scene *scn) {
    printf("\nCamera:\n");
    printCamera(scn->camera);
    printf("\n");
    
    printf("Spheres:\n");
    int num = scn->spheres.size();
    if(num != 0) {
        for (int i = 0; i < num; i++) {
            printSphere(scn->spheres[i]);
        }
    }
    else {
        printf("No Spheres\n");
    }
    printf("\n");
    
    printf("Image:\n");
    printImage(scn->file);
    printf("\n");
    
    printf("Background Color: ");
    printVector(scn->BGColor);
    printf("\n");
    
    printf("Directional Lights:\n");
    num = scn->directional.size();
    if (num > 0) {
        for (int i = 0; i < num; i++) {
            printLight(scn->directional[i]);
        }
    }
    else {
        printf("No Directional Lights\n");
    }
    printf("\n");
    
    printf("Point Lights:\n");
    num = scn->point.size();
    if (num > 0) {
        for (int i = 0; i < num; i++) {
            printLight(scn->point[i]);
        }
    }
    else {
        printf("No Point Lights\n");
    }
    printf("\n");
    
    printf("Spot Lights:\n");
    num = scn->spot.size();
    if (num > 0) {
        for (int i = 0; i < num; i++) {
            printLight(scn->spot[i]);
        }
    }
    else {
        printf("No Spot Lights\n");
    }
    printf("\n");
    
    printf("Ambient Light:\n");
    printLight(scn->ambient);
    printf("\n");
    
    printf("Depth: %d\n", scn->depth);
    
    if(scn->eyeray == PERSP) {
        printf("Projection Type: Perspective\n\n");
    }
    else {
        printf("Projection Type: Orthographic\n\n");
    }
    
    printf("Triangles:\n");
    num = scn->triangles.size();
    if(num != 0) {
		for(int i = 0; i < num; i++) {
			printTriangle(scn->triangles[i]);
		}	
	}
	else {
		printf("No Triangles\n");
	}
	printf("\n");
	
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
}

void printCamera(Camera camera) {
    printf("Position: ");
    printVector(camera.position);
    printf("Direction: ");
    printVector(camera.direction);
    printf("Up Vector: ");
    printVector(camera.cameraUp);
    printf("Half Angle: %f\n", camera.cameraHa);
}

void printTriangle(Triangle tri) {
    printf("Vertex1: ");
    printVector(tri.vertices[0]);
    printf("Vertex2: ");
    printVector(tri.vertices[1]);
    printf("Vertex3: ");
    printVector(tri.vertices[2]);
    printf("Normal: ");
    printVector(tri.normals[0]);
    if(tri.ntri) {
        printf("Normal2: ");
        printVector(tri.normals[1]);
        printf("Normal3: ");
        printVector(tri.normals[2]);
    }
    printf("RNormal: ");
    printVector(tri.normals[3]);
    if(tri.ntri) {
        printf("RNormal2: ");
        printVector(tri.normals[4]);
        printf("RNormal3: ");
        printVector(tri.normals[5]);
    }
    printf("Material:\n");
    printMaterial(tri.mat);
}

void printSphere(Sphere sph) {
    printf("Position: ");
    printVector(sph.position);
    printf("Radius: %f\n", sph.radius);
    printf("Material:\n");
    printMaterial(sph.mat);
}

void printPlane(Plane pl) {
	printf("Point: ");
	printVector(pl.point);
	printf("Normal: ");
	printVector(pl.normal);
	printf("Inverse Normal: ");
	printVector(pl.inormal);
	printf("Material:\n");
	printMaterial(pl.mat);
}

void printRect(Rectangle rec) {
	printf("Point 1: ");
	printVector(rec.points[0]);
	printf("Point 2: ");
	printVector(rec.points[1]);
	printf("Point 3: ");
	printVector(rec.points[2]);
	printf("Point 4: ");
	printVector(rec.points[3]);
	printf("Normal: ");
	printVector(rec.normal);
	printf("Inverse Normal: ");
	printVector(rec.inormal);
	printf("Material:\n");
	printMaterial(rec.mat);
}

void printMaterial(Material mat) {
    printf("Ambient Color: ");
    printVector(mat.ambientColor);
    printf("Diffuse Color: ");
    printVector(mat.diffuseColor);
    printf("Specular Color: ");
    printVector(mat.specularColor);
    printf("Cosin Power: %f\n", mat.cosPow);
    printf("Transmissive Color: ");
    printVector(mat.transmissiveColor);
    printf("ior: %f\n", mat.ior);
}

void printLight(Light l) {
    printf("Color: ");
    printVector(l.color);
    printf("Position: ");
    printVector(l.position);
    printf("Direction: ");
    printVector(l.direction);
    printf("Angle: (%f %f)\n", l.angle1, l.angle2);
}

void printImage(File img) {
    printf("File Name: %s\n", img.filename.c_str());
    printf("Width: %d\n", img.width);
    printf("Height: %d\n", img.height);
}

void printVector(Vector fVec) {
    printf("(%f %f %f)\n", fVec.x, fVec.y, fVec.z);
}
