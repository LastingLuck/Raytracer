# README #

This is a simple C++ implementation of a ray tracer that outputs images in bmp format. This project makes use of the EasyBMP library by Paul Macklin. All credit for that goes to him.

More info and images can be found [here](https://sites.google.com/site/alexdahl5607/ray-tracer-pt-2). The next feature I'm working on is implementing a Bounding Volume Hierarchy structure to accelerate it further. The structure is there, however it is still incomplete. There are currently issues with triangles not testing inside of a box when fully encompassed.

**Supports**:

* Spheres
* Triangles
* Triangles with independent vertex normals
* Directional Lights
* Point Lights
* Spot Lights
* Shadows
* Reflection
* Refraction
* Parallelization
* Jittered Supersampling

To run this, just call the make command in the terminal and run the command:

$ ./raytrace scene_file

Run make with the option PARALLEL=1 to turn on parallelization.

Run make with the option DEBUG=1 to turn on debugging output. (Warning: This could be a lot in info)


Format for the scene file is a regular text file. Each line is a single component and the first word specifies how the rest of the line is interpreted. The rest of the line is a space-seperated list of parameters. Lines starting with a '#' are comments and will be skipped. 

**Commands**:

camera px py pz dx dy dz ux uy uz ha
 
* (px py pz) is the postition of the camera
 
* (dx dy dz) is the camera's viewing direction
 
* (ux uy uz) is the camera's up vector
 
* ha is one-half of the height angle of the viewing frustum
    
* Default is 0 0 0 0 0 1 0 1 0 45

film_resolution width height
 
* (width, height) is the resolution of the image
    
* Default is 640x480

output_image filename
 
* filename is the name of the outputted bmp file

* Default is raytraced.bmp

sphere x y z r
 
* (x y z) is the position
 
* r is the radius of the sphere

background r g b
 
* (r g b) is the color of the background. If there are no intersections, this is the color of the pixel
    
* Default is (0 0 0)

Note: The color values for material are assumed to be between 0 and 1. Values are not clamped.

material ar ag ab dr dg db sr sg sb ns tr tg tb ior
 
* (ar ag ab) is the ambient color of the object
 
* (dr dg db) is the diffuse color of the object
 
* (sr sg sb) is the specular color of the object
 
* ns is the Phong cosine power. The higher the number the larger the lighing dropoff
 
* tr, tg, tb, ior are not implemented in this version so ignore them
    
* Default is 0 0 0 1 1 1 0 0 0 5 0 0 0 1

directional_light r g b x y z
 
* (r g b) is the color of the light
 
* (x y x) is the direction the light shines. The light shines everywhere in this direction

point_light r g b x y z
 
* (r g b) is the color of the light
 
* (x y z) is the position of the light. Light falls off proportional to the distance squared

spot_light r g b px py pz dx dy dz angle1 angle2
 
* (r g b) is the color of the light
 
* (px py pz) is the position of the light
 
* (dx dy dz) is the direction of the light
 
* angle1 and angle2 determine how the light falls off. Points at an angle less than angle1 have the light act as a point light. Points 
at an angle greater than angle2 don't contribute to the light. Inbetween the light falls off linearly

ambient_light r g b
 
* (r g b) is the color of the light that is applied to everything in the scene
    
* Default is (0 0 0)

max_depth n
 
* n Maximum recursion depth for the rays.
    
* Default is 5

Note: Vertexes must be specified before it can be indexed by a triangle. This is due to the scene parser interpreting it line by line.

vertex x y z
 
* (x, y, z) - Position of the vertex

normal x y z
 
* (x, y, z) - Direction of the normal

triangle v1 v2 v3
 
* v1 v2 v3 - The index of the vertices specified before the triangle. Starts from 0.

normal_triangle v1 v2 v3 n1 n2 n3
 
* v1 v2 v3 - Same as the triangle
 
* n1 n2 n3 - Similar to the triangle except for the normals specified before the triangle

sample_rate n
 
* n - the number of times to sample per pixel. Default is 1(no additional sampling)