#include <iostream>
#include <fstream>
#include <ostream>
#include <string>
#include <limits>
#define _USE_MATH_DEFINES
#include <cmath>
using namespace std;

class Vector3D
{
public:
    Vector3D(double x, double y, double z);
    Vector3D();
    ~Vector3D();

    double getX() const;
    double getY() const;
    double getZ() const;

    Vector3D normalize() const;
    double dotProduct(const Vector3D &vec) const;
    Vector3D crossProduct(const Vector3D &vec) const;

    double length() const;

    friend ostream &operator<<(ostream &, const Vector3D &);
    friend Vector3D operator+(const Vector3D &a, const Vector3D &b);
    friend Vector3D operator-(const Vector3D &a, const Vector3D &b);
    friend Vector3D operator*(const Vector3D &a, double s);
    friend Vector3D operator*(double s, const Vector3D &a);
    friend Vector3D operator*(const Vector3D &a, const Vector3D &b); // element-wise multiplication
    friend Vector3D operator/(const Vector3D &a, double s);

    // unary minus:
    friend Vector3D operator-(const Vector3D &a);

private:
    double v[3];
};

class Light
{
public:
    Light();
    ~Light();
    Light(string type, double intensity, const Vector3D &position, const Vector3D &direction);
    string getType() { return type; }
    double getIntensity() { return intensity; }
    Vector3D getPosition() { return position; }
    Vector3D getDirection() { return direction; }

private:
    string type;
    double intensity;
    Vector3D position;
    Vector3D direction;
};

class Object
{
public:
    ~Object();
    virtual bool intersect(const Vector3D &rayOrigin, const Vector3D &rayDirection, double &t) const = 0;
    virtual Vector3D getNormalAt(const Vector3D &point) const = 0;

    Vector3D albedo;  // Base color
    float metallic;   // 0.0 = dielectric, 1.0 = metal
    float roughness;  // Surface roughness (0.0 = smooth, 1.0 = rough)
    float reflective; // Reflectivity (unchanged)
};

class Triangle : public Object
{
public:
    Triangle();
    ~Triangle();
    Triangle(const Vector3D &v0, const Vector3D &v1, const Vector3D &v2, const Vector3D &albedo, float metallic, float roughness, float reflective);

    virtual bool intersect(const Vector3D &rayOrigin, const Vector3D &rayDirection, double &t) const override;
    virtual Vector3D getNormalAt(const Vector3D &point) const override;

    Vector3D v0, v1, v2; // Triangle vertices
};

class Sphere : public Object
{
public:
    Sphere();
    ~Sphere();
    Sphere(const Vector3D &center, double radius, const Vector3D &albedo, float metallic, float roughness, float reflective);

    virtual bool intersect(const Vector3D &rayOrigin, const Vector3D &rayDirection, double &t) const override;
    virtual Vector3D getNormalAt(const Vector3D &point) const override;

    Vector3D center;
    double radius;
};

// Function prototypes for intersection
bool IntersectRaySphere(
    const Vector3D &rayOrigin,
    const Vector3D &rayDirection,
    const Sphere &sphere,
    double &t);

bool IntersectRayTriangle(
    const Vector3D &rayOrigin,
    const Vector3D &rayDirection,
    const Triangle &triangle,
    double &intersectionDistance);

// Scene constants
Light scene_lights[] = {
    Light("ambient", 0.2, Vector3D(), Vector3D()),
    Light("point", 0.6, Vector3D(2, 1, 0), Vector3D()),
    Light("directional", 0.2, Vector3D(), Vector3D(1, 4, 4))};

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

const int Cw = 300;    // Canvas width
const int Ch = 300;    // Canvas height
const double Vw = 1.0; // Viewport width
const double Vh = 1.0; // Viewport height
const double D = 1.0;  // Distance from camera to viewport
const int MAX_RECURSION_DEPTH = 5;
// const int BACKGROUND_COLOR[3] = {135, 206, 235};
const int BACKGROUND_COLOR[3] = {0, 0, 0};
