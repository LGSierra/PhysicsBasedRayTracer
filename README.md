# PhysicsBasedRayTracer!

![BaseRenderPBR](https://github.com/user-attachments/assets/1b45bb78-1798-4129-a5d0-dd9fd6266da3)
![Metallic shiny](https://github.com/user-attachments/assets/f39632ae-bced-4b9b-9250-67c0e21c4f86)
![TriangleReflectionRenderPBR](https://github.com/user-attachments/assets/f3f770db-aeae-4c00-b704-f34f8c86f173)
![MetalSpheres](https://github.com/user-attachments/assets/8d93da4e-6089-48d3-8f49-0485eff46a10)

## Welcome to a simple implementation for PBR Ray Tracing!

### Getting Started
This code compiles with C++11. To run this program:

- Run this command in your terminal inside the folder you downloaded all the files:
```
  g++ -o ray_tracer mainRayTracerPBR.cpp -std=c++11
```
This will compile the code and create "ray_tracer.exe"

- To run the program. execute this command:
```
./ray_tracer.exe
```

Otherwise, there are C++ code runner extensions in Visual Studio Code that could also work.

## How to change the image to be rendered?

You can change the parameters in ```headerRayTracerPBR.h``` to change the image being rendered. These parameters represent light, objects, and their respective properties.

### Lights

You can change the lights being rendered in the scene.
```
// Scene constants
Light scene_lights[] = {
    Light("ambient", 0.2, Vector3D(), Vector3D()), 
    Light("point", 0.6, Vector3D(2, 1, 0), Vector3D()),
    Light("directional", 0.2, Vector3D(), Vector3D(1, 4, 4))};
```
Parameters for a light object:
```Light( type of light, intensity, position of light(vec3), direction of light(vec3));```

The parameters for their respective type of light:
- ambient (intensity)
- point (intensity, position)
- directional (intensity, direction)

All lights have a intensity parameter, but ambient light only has intensity, while point light also has position, and directional light has direction. 


### Objects

This is the base image:
```
Object *scene_objects[] = {
    new Triangle(
        Vector3D(0, 0, 5),
        Vector3D(1, 2, 4),
        Vector3D(-1, 2, 4),
        Vector3D(1.0, 0.0, 0.0), // (255,0,0)
        0.0,                     // metallic
        0.06,                    // roughness
        0.5                      // reflective
        ),
    new Sphere( // red sphere
        Vector3D(0, -1, 3),
        1,
        Vector3D(1.0, 0.0, 0.0), // (255,0,0)
        0.0,                     // metallic
        0.06,                    // roughness
        0.2                      // reflective
        ),
    new Sphere( // blue sphere
        Vector3D(2, 0, 4),
        1,
        Vector3D(0.0, 0.0, 1.0), // (0,0,255)
        0.0,                     // metallic
        0.06,                    // roughness
        0.3                      // reflective
        ),
    new Sphere( // green sphere
        Vector3D(-2, 0, 4),
        1,
        Vector3D(0.0, 1.0, 0.0), // (0,255,0)
        0.0,                     // metallic
        0.41,                    // roughness
        0.4                      // reflective
        ),
    new Sphere( // yellow sphere
        Vector3D(0, -5001, 0),
        5000,
        Vector3D(1.0, 1.0, 0.0), // (255,255,0)
        0.0,                     // metallic
        0.045,                   // roughness
        0.5                      // reflective
        )};
```
The two types of objects supported are triangles and spheres:
```
new Triangle(
        Vector3D(0, 0, 5),
        Vector3D(1, 2, 4),
        Vector3D(-1, 2, 4),
        Vector3D(1.0, 0.0, 0.0), // (255,0,0)
        0.0,                     // metallic
        0.06,                    // roughness
        0.5                      // reflective
        ),
    new Sphere( // red sphere
        Vector3D(0, -1, 3),
        1,
        Vector3D(1.0, 0.0, 0.0), // (255,0,0)
        0.0,                     // metallic
        0.06,                    // roughness
        0.2                      // reflective
        )
```
The triangle object needs three points for its position. The sphere object needs one point for the position of its center, and the radius of the sphere. The rest of the parameters are:
- albedo (color): Expressed as a double, maximum value of 1.0 for all three numbers of "RGB". Note, the values accepted are not literal RGB values, but rather the RGB values divided by 255. As we can see from the color red, it's RGB is (255,0,0) but since 255/255 = 1, the values are (1.0,0,0).
- metallic: Affects the albedo value. Anything higher than 0.4 will be considered a metallic object, and will model such appearance.
- roughness: Affects diffuse reflection. Higher value means rougher surface, lower values mean smoother surfaces.
- reflective: Affects specular reflection. Higher value means a more reflective surface.


