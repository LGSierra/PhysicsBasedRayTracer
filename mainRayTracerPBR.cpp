#include <iostream>
#include <fstream>
#include <limits>
#include <cmath>
#include <chrono>
#include "impRayTracerPBR.cpp" // Includes all logic for ray tracing
using namespace std;

int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    // Image buffer
    int image[Cw][Ch][3];

    // Camera setup
    Vector3D cameraPosition(0, 0, 0);

    // Generate the image pixel by pixel
    for (int x = -Cw / 2; x < Cw / 2; x++)
    {
        for (int y = -Ch / 2; y < Ch / 2; y++)
        {
            Vector3D rayDirection = CanvasToViewport(x, y);

            int color[3];
            TraceRay(cameraPosition, rayDirection, 1.0, std::numeric_limits<double>::infinity(), MAX_RECURSION_DEPTH, color);

            int i = x + Cw / 2;     // Convert canvas x to image buffer x
            int j = Ch / 2 - y - 1; // Convert canvas y to image buffer y

            // Store computed color in the buffer
            image[i][j][0] = color[0];
            image[i][j][1] = color[1];
            image[i][j][2] = color[2];
        }
    }

    // Write the image to a PPM file
    ofstream myFile("render.ppm");
    createPPM(myFile, Cw, Ch, image);
    myFile.close();
    cout << "Image file created! Please check \"render.ppm\"." << endl;

    // Cleanup dynamically allocated objects
    for (size_t i = 0; i < sizeof(scene_objects) / sizeof(scene_objects[0]); i++)
    {
        delete scene_objects[i];
    }

    // Record the end time
    auto end = std::chrono::high_resolution_clock::now();

    // Compute the duration in a suitable time unit (e.g. milliseconds)
    std::chrono::duration<double, std::milli> elapsed = end - start;

    // Print the elapsed time
    std::cout << "Program execution time: " << elapsed.count() << " ms\n";

    return 0;
}
