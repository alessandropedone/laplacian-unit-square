/**
 * @file vtk.hpp
 * @brief Utilities for reading and writing 2D grid data in VTK STRUCTURED_GRID format.
 *
 * This header provides functions to export a flattened 2D grid to a VTK file for visualization,
 * and to read grid data and coordinates from a VTK file in STRUCTURED_GRID format.
 */

#ifndef VTK_HPP
#define VTK_HPP

#include <vector>
#include <fstream>
#include <iomanip>
#include <string>

namespace vtk
{
    /**
     * @brief Writes a 2D grid to a VTK file in STRUCTURED_GRID format.
     *
     * This function exports a flattened 2D grid of size n x n to a VTK file,
     * allowing visualization with tools that support the VTK format.
     *
     * @param grid      Flattened 2D grid data of size n x n (row-major order).
     * @param n         The dimension of the grid (number of rows and columns).
     * @param filename  Output VTK file name (default is "output.vtk").
     */
    void write(const std::vector<double> &grid, int n, const std::string &filename = "output.vtk")
    {
        std::cout << "Writing VTK file: " << filename << std::endl;
        std::ofstream vtkFile(filename);
        vtkFile << "# vtk DataFile Version 3.0\n";
        vtkFile << "vtk output\n";
        vtkFile << "ASCII\n";
        vtkFile << "DATASET STRUCTURED_GRID\n";
        vtkFile << "DIMENSIONS " << n << " " << n << " 1\n";
        vtkFile << "POINTS " << n * n << " float\n";
        vtkFile << std::setprecision(8);
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                vtkFile << static_cast<double>(i) / n << " " << static_cast<double>(j) / n << " 0\n";
            }
        }
        vtkFile << "\n\n";
        vtkFile << "POINT_DATA " << n * n << "\n";
        vtkFile << "SCALARS values float\n";
        vtkFile << "LOOKUP_TABLE default\n";
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                vtkFile << grid[i * n + j] << "\n";
            }
        }
        vtkFile.close();
    }
    
    /**
     * @brief Reads a VTK file and extracts grid data and coordinates.
     *
     * This function reads a VTK file in STRUCTURED_GRID format and extracts
     * the grid values and their corresponding coordinates.
     *
     * @param filename  Input VTK file name.
     * @param grid      Output vector to store the grid values.
     * @param coords    Output vector to store the coordinates of the points.
     */
    void read(const std::string &filename, std::vector<double> &grid, std::vector<std::pair<double, double>> &coords)
    {
        std::ifstream vtkIn("output.vtk");
        std::string line;

        // Skip header lines until "POINTS"
        while (std::getline(vtkIn, line))
        {
            if (line.find("POINTS") == 0)
                break;
        }

        // Read coordinates
        for (size_t idx = 0; idx < grid.size(); ++idx)
        {
            double x, y, z;
            vtkIn >> x >> y >> z;
            coords[idx] = {x, y};
        }

        // Skip to "SCALARS" and "LOOKUP_TABLE"
        while (std::getline(vtkIn, line))
        {
            if (line.find("LOOKUP_TABLE") == 0)
                break;
        }

        // Read grid values
        for (size_t idx = 0; idx < grid.size(); ++idx)
        {
            vtkIn >> grid[idx];
        }
        vtkIn.close();
    }
}

#endif // VTK_HPP