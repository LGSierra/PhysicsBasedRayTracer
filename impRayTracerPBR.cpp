#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm> // For std::min
#include "headerRayTracerPBR.h"
// #include <corecrt_math_defines.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
using namespace std;

// Object base class
Object::~Object() {}

// Vector3D implementation
Vector3D::~Vector3D() {}
Vector3D::Vector3D(double x, double y, double z)
{
    v[0] = x;
    v[1] = y;
    v[2] = z;
}
Vector3D::Vector3D() : v{0, 0, 0} {}
double Vector3D::getX() const { return v[0]; }
double Vector3D::getY() const { return v[1]; }
double Vector3D::getZ() const { return v[2]; }
Vector3D Vector3D::normalize() const
{
    double mag = length();
    if (mag == 0)
        return Vector3D(0, 0, 0);
    return *this / mag;
}
double Vector3D::dotProduct(const Vector3D &vec) const
{
    return v[0] * vec.v[0] + v[1] * vec.v[1] + v[2] * vec.v[2];
}
Vector3D Vector3D::crossProduct(const Vector3D &vec) const
{
    return Vector3D(
        v[1] * vec.v[2] - v[2] * vec.v[1],
        v[2] * vec.v[0] - v[0] * vec.v[2],
        v[0] * vec.v[1] - v[1] * vec.v[0]);
}
double Vector3D::length() const
{
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

// Light class
Light::~Light() {}
Light::Light(string type, double intensity, const Vector3D &position, const Vector3D &direction)
{
    this->type = type;
    this->intensity = intensity;
    this->position = position;
    this->direction = direction;
}

// Triangle class
Triangle::~Triangle() {}
Triangle::Triangle(const Vector3D &v0, const Vector3D &v1, const Vector3D &v2, const Vector3D &albedo, float metallic, float roughness, float reflective)
{
    this->v0 = v0;
    this->v1 = v1;
    this->v2 = v2;
    this->albedo = albedo;
    this->metallic = metallic;
    this->roughness = roughness;
    this->reflective = reflective;
}
bool Triangle::intersect(const Vector3D &rayOrigin, const Vector3D &rayDirection, double &t) const
{
    return IntersectRayTriangle(rayOrigin, rayDirection, *this, t);
}
Vector3D Triangle::getNormalAt(const Vector3D &point) const
{
    Vector3D edge1 = v1 - v0;
    Vector3D edge2 = v2 - v0;
    return edge1.crossProduct(edge2).normalize();
}

// Sphere class
Sphere::~Sphere() {}
Sphere::Sphere(const Vector3D &center, double radius, const Vector3D &albedo, float metallic, float roughness, float reflective)
{
    this->center = center;
    this->radius = radius;
    this->albedo = albedo;
    this->metallic = metallic;
    this->roughness = roughness;
    this->reflective = reflective;
}
bool Sphere::intersect(const Vector3D &rayOrigin, const Vector3D &rayDirection, double &t) const
{
    return IntersectRaySphere(rayOrigin, rayDirection, *this, t);
}
Vector3D Sphere::getNormalAt(const Vector3D &point) const
{
    return (point - center).normalize();
}

// Operator Overloads for Vector3D
Vector3D operator+(const Vector3D &a, const Vector3D &b) { return Vector3D(a.getX() + b.getX(), a.getY() + b.getY(), a.getZ() + b.getZ()); }
Vector3D operator-(const Vector3D &a, const Vector3D &b) { return Vector3D(a.getX() - b.getX(), a.getY() - b.getY(), a.getZ() - b.getZ()); }
Vector3D operator*(const Vector3D &a, double s) { return Vector3D(a.getX() * s, a.getY() * s, a.getZ() * s); }
Vector3D operator*(double s, const Vector3D &a) { return Vector3D(a.getX() * s, a.getY() * s, a.getZ() * s); }
Vector3D operator*(const Vector3D &a, const Vector3D &b)
{
    return Vector3D(a.getX() * b.getX(), a.getY() * b.getY(), a.getZ() * b.getZ());
}

Vector3D operator/(const Vector3D &a, double s) { return Vector3D(a.getX() / s, a.getY() / s, a.getZ() / s); }

Vector3D operator-(const Vector3D &a) { return Vector3D(-a.getX(), -a.getY(), -a.getZ()); }

// For main file
Vector3D CanvasToViewport(double x, double y)
{
    return Vector3D(x * Vw / Cw, y * Vh / Ch, D);
}

void createPPM(ofstream &file, int w, int h, int image[][Ch][3])
{
    file << "P3" << endl;
    file << w << " " << h << endl;
    file << 255 << endl;

    for (int j = 0; j < h; j++)
    {
        for (int i = 0; i < w; i++)
        {
            int r = image[i][j][0];
            int g = image[i][j][1];
            int b = image[i][j][2];
            file << r << " " << g << " " << b << " ";
        }
        file << endl;
    }
}

// Utility: ReflectRay
Vector3D ReflectRay(const Vector3D &R, const Vector3D &N)
{
    return R - 2 * N * R.dotProduct(N);
}

// IntersectRaySphere
bool IntersectRaySphere(
    const Vector3D &rayOrigin,
    const Vector3D &rayDirection,
    const Sphere &sphere,
    double &t)
{
    Vector3D CO = rayOrigin - sphere.center;

    double a = rayDirection.dotProduct(rayDirection);
    double b = 2 * CO.dotProduct(rayDirection);
    double c = CO.dotProduct(CO) - sphere.radius * sphere.radius;

    double discriminant = b * b - 4 * a * c;

    if (discriminant < 0)
    {
        // the ray doesn't intersect the sphere
        return false;
    }

    // find the nearest t
    double sqrt_discriminant = sqrt(discriminant);
    double t1 = (-b + sqrt_discriminant) / (2 * a);
    double t2 = (-b - sqrt_discriminant) / (2 * a);

    // select the smallest positive t
    if (t1 >= 0 && t2 >= 0)
        t = std::min(t1, t2);
    else if (t1 >= 0)
        t = t1;
    else if (t2 >= 0)
        t = t2;
    else
        return false; // t1 and t2 are negative

    return true;
}

bool IntersectRayTriangle(
    const Vector3D &rayOrigin,
    const Vector3D &rayDirection,
    const Triangle &triangle,
    double &intersectionDistance)
{
    // calculate triangle normal
    Vector3D edgeVector1 = triangle.v1 - triangle.v0;
    Vector3D edgeVector2 = triangle.v2 - triangle.v0;
    Vector3D triangleNormal = edgeVector1.crossProduct(edgeVector2).normalize();

    // calculate denominator of the plane equation
    double normalDotRayDirection = triangleNormal.dotProduct(rayDirection);

    // check if the ray is parallel to the triangle plane
    if (normalDotRayDirection == 0)
    {
        // the ray is parallel to the plane
        return false;
    }

    // calculate the distance from ray origin to the intersection point
    double d = -triangleNormal.dotProduct(triangle.v0);
    intersectionDistance = -(triangleNormal.dotProduct(rayOrigin) + d) / normalDotRayDirection;

    if (intersectionDistance < 0)
    {
        // the triangle is behind the ray origin, so it gets discarded
        return false;
    }

    // calculate the intersection point
    Vector3D intersectionPoint = rayOrigin + rayDirection * intersectionDistance;

    // inside-out test to see if the point is inside the triangle
    // gets repeated on all triangle edges
    // point is outside of the triangle if the dot product of the normal and crossProduct of edge is negative
    Vector3D edgeCrossProduct;

    // edge 0
    Vector3D edge0 = triangle.v1 - triangle.v0;
    Vector3D vectorToPoint0 = intersectionPoint - triangle.v0;
    edgeCrossProduct = edge0.crossProduct(vectorToPoint0);
    if (triangleNormal.dotProduct(edgeCrossProduct) < 0)
        return false; // intersection point is outside the triangle

    // edge 1
    Vector3D edge1 = triangle.v2 - triangle.v1;
    Vector3D vectorToPoint1 = intersectionPoint - triangle.v1;
    edgeCrossProduct = edge1.crossProduct(vectorToPoint1);
    if (triangleNormal.dotProduct(edgeCrossProduct) < 0)
        return false; // Intersection point is outside the triangle

    // edge 2
    Vector3D edge2 = triangle.v0 - triangle.v2;
    Vector3D vectorToPoint2 = intersectionPoint - triangle.v2;
    edgeCrossProduct = edge2.crossProduct(vectorToPoint2);
    if (triangleNormal.dotProduct(edgeCrossProduct) < 0)
        return false; // intersection point is outside the triangle

    // the intersection point is inside the triangle
    return true;
}

// Utility: Closest Intersection
const Object *closestIntersection(
    const Vector3D &rayOrigin,
    const Vector3D &rayDirection,
    double t_min,
    double t_max,
    double &closest_t)
{
    closest_t = std::numeric_limits<double>::infinity();
    const Object *closest_object = nullptr;

    int num_objects = sizeof(scene_objects) / sizeof(scene_objects[0]);

    for (int i = 0; i < num_objects; i++)
    {
        double t;
        if (scene_objects[i]->intersect(rayOrigin, rayDirection, t))
        {
            if (t_min < t && t < t_max && t < closest_t)
            {
                closest_t = t;
                closest_object = scene_objects[i];
            }
        }
    }

    return closest_object;
}

// PBR Utility Functions
Vector3D FresnelSchlick(float cosTheta, const Vector3D &F0)
{
    return F0 + (Vector3D(1.0, 1.0, 1.0) - F0) * pow(1.0 - cosTheta, 5.0);
}
float DistributionGGX(const Vector3D &N, const Vector3D &H, float roughness)
{
    float a = roughness * roughness;
    float a2 = a * a;
    float nDotH = max((float)N.dotProduct(H), 0.0f);
    float denom = (nDotH * nDotH * (a2 - 1.0) + 1.0);
    return a2 / (M_PI * denom * denom);
}
float GeometrySchlickGGX(float nDotV, float roughness)
{
    float r = roughness + 1.0;
    float k = (r * r) / 8.0;
    return nDotV / (nDotV * (1.0 - k) + k);
}

// ComputeLighting with PBR
Vector3D ComputeLighting(const Vector3D &point, const Vector3D &normal, const Vector3D &viewDirection, const Object &object)
{
    Vector3D result(0.0, 0.0, 0.0);
    Vector3D F0 = object.metallic > 0.0f ? object.albedo : Vector3D(0.04, 0.04, 0.04);

    for (size_t i = 0; i < sizeof(scene_lights) / sizeof(scene_lights[0]); i++)
    {
        Light &light = scene_lights[i];
        Vector3D lightDirection = (light.getType() == "point")
                                      ? (light.getPosition() - point).normalize()
                                      : light.getDirection().normalize();

        // Shadow check: determine t_max based on light type.
        double t_max = (light.getType() == "point") ? 1.0 : std::numeric_limits<double>::infinity();
        double shadow_t;
        const Object *shadow_object = closestIntersection(point, lightDirection, 0.001, t_max, shadow_t);
        if (shadow_object != NULL)
        {
            continue; // point is in shadow; skip this light
        }

        Vector3D H = (viewDirection + lightDirection).normalize();
        Vector3D F = FresnelSchlick(max((float)viewDirection.dotProduct(H), 0.0f), F0);

        Vector3D kd = (Vector3D(1.0f, 1.0f, 1.0f) - F) * (1.0f - object.metallic);

        float nDotL = max((float)normal.dotProduct(lightDirection), 0.0f);
        Vector3D diffuse = kd * object.albedo * nDotL;

        float D = DistributionGGX(normal, H, object.roughness);
        float G = GeometrySchlickGGX(normal.dotProduct(viewDirection), object.roughness) *
                  GeometrySchlickGGX(normal.dotProduct(lightDirection), object.roughness);
        Vector3D specular = (D * G * F) / max(4.0f * (float)normal.dotProduct(viewDirection) * nDotL, 0.001f);

        result = result + light.getIntensity() * (diffuse + specular) * nDotL;
    }

    return result;
}

// TraceRay with PBR
void TraceRay(const Vector3D &rayOrigin, const Vector3D &rayDirection, double t_min, double t_max, int recursion_depth, int color[3])
{
    double closest_t;
    const Object *closest_object = closestIntersection(rayOrigin, rayDirection, t_min, t_max, closest_t);

    // If no intersection, just return background in sRGB since no blending will occur
    if (!closest_object)
    {
        color[0] = BACKGROUND_COLOR[0];
        color[1] = BACKGROUND_COLOR[1];
        color[2] = BACKGROUND_COLOR[2];
        return;
    }

    // Compute intersection and normal
    Vector3D P = rayOrigin + rayDirection * closest_t;
    Vector3D N = closest_object->getNormalAt(P);
    Vector3D V = (-rayDirection).normalize();

    // Compute lighting in linear space
    Vector3D lighting = ComputeLighting(P, N, V, *closest_object);

    // We'll store the final linear color here
    Vector3D finalColor = lighting; // linear color

    // If reflective, compute reflection
    if (recursion_depth > 0 && closest_object->reflective > 0.0f)
    {
        Vector3D R = ReflectRay(rayDirection, N).normalize();
        int reflected_color[3];
        TraceRay(P, R, 0.001, t_max, recursion_depth - 1, reflected_color);

        // The reflected_color you get here is gamma-corrected [0–255].
        // Convert it back to linear before blending:
        double invGamma = 2.2;
        Vector3D reflectedLinear(
            pow((reflected_color[0] / 255.0), invGamma),
            pow((reflected_color[1] / 255.0), invGamma),
            pow((reflected_color[2] / 255.0), invGamma));

        // Blend in linear space
        double rfl = closest_object->reflective;
        finalColor = finalColor * (1.0 - rfl) + reflectedLinear * rfl;
    }

    // Now apply gamma correction once at the end
    double gamma = 1.0 / 2.2;
    double lr = pow(finalColor.getX(), gamma);
    double lg = pow(finalColor.getY(), gamma);
    double lb = pow(finalColor.getZ(), gamma);

    color[0] = std::min((int)(lr * 255), 255);
    color[1] = std::min((int)(lg * 255), 255);
    color[2] = std::min((int)(lb * 255), 255);
}
