#ifndef PLOT_HPP
#define PLOT_HPP

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <cmath>
#include <filesystem>
#include <algorithm>
#include <iomanip>

namespace fs = std::filesystem;

struct DataRow
{
    int n;
    double h;
    double serial;
    double omp;
    double mpi;
    double hybrid;
    double l2_error;
};

class CSVReader
{
public:
    static std::vector<DataRow> readCSV(const std::string &filename)
    {
        std::vector<DataRow> data;
        std::ifstream file(filename);
        std::string line;

        if (!file.is_open())
        {
            std::cerr << "Error opening file: " << filename << std::endl;
            return data;
        }

        // Skip header line
        std::getline(file, line);

        while (std::getline(file, line))
        {
            std::stringstream ss(line);
            std::string cell;
            DataRow row;

            std::getline(ss, cell, ',');
            row.n = std::stoi(cell);

            std::getline(ss, cell, ',');
            row.serial = std::stod(cell);

            std::getline(ss, cell, ',');
            row.omp = std::stod(cell);

            std::getline(ss, cell, ',');
            row.mpi = std::stod(cell);

            std::getline(ss, cell, ',');
            row.hybrid = std::stod(cell);

            // Skip 3 columns (h, serial, omp columns are handled separately)
            for (int i = 0; i < 3; ++i) {
                std::getline(ss, cell, ',');
            }
            std::getline(ss, cell, ',');
            row.l2_error = std::stod(cell);

            row.h = 1.0 / (row.n - 1);
            data.push_back(row);
        }

        file.close();
        return data;
    }
};

class Plotter
{
private:
    static int plot_counter;

public:
    // Write data to file for plotting
    static void writeDataFile(const std::string &filename,
                              const std::vector<double> &x,
                              const std::vector<std::vector<double>> &y_data,
                              const std::vector<std::string> &labels)
    {
        std::ofstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "Error creating data file: " << filename << std::endl;
            return;
        }

        // Write header
        file << "# X";
        for (const auto &label : labels)
        {
            file << "\t" << label;
        }
        file << std::endl;

        // Write data
        for (size_t i = 0; i < x.size(); ++i)
        {
            file << std::scientific << std::setprecision(6) << x[i];
            for (const auto &y_series : y_data)
            {
                if (i < y_series.size())
                {
                    file << "\t" << std::scientific << std::setprecision(6) << y_series[i];
                }
            }
            file << std::endl;
        }
        file.close();
    }

    // Create gnuplot script
    static void createGnuplotScript(const std::string &script_name,
                                    const std::string &data_file,
                                    const std::string &title,
                                    const std::string &xlabel,
                                    const std::string &ylabel,
                                    const std::vector<std::string> &labels,
                                    bool logx = true, bool logy = true,
                                    const std::string &output_file = "")
    {
        std::ofstream script(script_name);
        if (!script.is_open())
        {
            std::cerr << "Error creating gnuplot script: " << script_name << std::endl;
            return;
        }

        script << "set terminal png enhanced font 'Arial,12' size 800,600" << std::endl;
        if (!output_file.empty())
        {
            script << "set output '" << output_file << "'" << std::endl;
        }
        else
        {
            script << "set output '" << script_name.substr(0, script_name.find_last_of('.')) << ".png'" << std::endl;
        }

        script << "set title '" << title << "'" << std::endl;
        script << "set xlabel '" << xlabel << "'" << std::endl;
        script << "set ylabel '" << ylabel << "'" << std::endl;
        script << "set grid" << std::endl;
        script << "set key outside right" << std::endl;

        if (logx)
            script << "set logscale x 2" << std::endl;
        if (logy)
            script << "set logscale y 2" << std::endl;

        script << "plot ";
        for (size_t i = 0; i < labels.size(); ++i)
        {
            if (i > 0)
                script << ", ";
            script << "'" << data_file << "' using 1:" << (i + 2)
                   << " with linespoints title '" << labels[i] << "'";
        }
        script << std::endl;

        script.close();
    }

    // Console output for quick data inspection
    static void printDataSummary(const std::string &title,
                                 const std::vector<double> &x,
                                 const std::vector<std::vector<double>> &y_data,
                                 const std::vector<std::string> &labels)
    {
        std::cout << "\n=== " << title << " ===" << std::endl;
        std::cout << std::fixed << std::setprecision(4);

        // Print header
        std::cout << std::setw(12) << " ";
        for (const auto &label : labels)
        {
            std::cout << std::setw(12) << label;
        }
        std::cout << std::endl;

        // Print separator
        std::cout << std::string(12 + labels.size() * 12, '-') << std::endl;

        // Print data (first few and last few rows)
        size_t rows_to_show = std::min(size_t(5), x.size());
        for (size_t i = 0; i < rows_to_show; ++i)
        {
            std::cout << std::setw(12) << x[i];
            for (const auto &y_series : y_data)
            {
                if (i < y_series.size())
                {
                    std::cout << std::setw(12) << y_series[i];
                }
            }
            std::cout << std::endl;
        }

        if (x.size() > 10)
        {
            std::cout << std::setw(12) << "..." << std::endl;
            for (size_t i = x.size() - rows_to_show; i < x.size(); ++i)
            {
                std::cout << std::setw(12) << x[i];
                for (const auto &y_series : y_data)
                {
                    if (i < y_series.size())
                    {
                        std::cout << std::setw(12) << y_series[i];
                    }
                }
                std::cout << std::endl;
            }
        }
        std::cout << std::endl;
    }
};

int Plotter::plot_counter = 0;

void gridSizeTest(const std::string &filename)
{
    std::vector<DataRow> data = CSVReader::readCSV("test/data/" + filename);

    if (data.empty())
    {
        std::cerr << "No data loaded from file: " << filename << std::endl;
        return;
    }

    // Create plots directory if it doesn't exist
    fs::create_directories("test/plots");

    // Extract data for plotting
    std::vector<double> n_values, h_values;
    std::vector<double> serial_times, omp_times, mpi_times, hybrid_times, l2_errors;

    for (const auto &row : data)
    {
        n_values.push_back(static_cast<double>(row.n));
        h_values.push_back(row.h);
        serial_times.push_back(row.serial);
        omp_times.push_back(row.omp);
        mpi_times.push_back(row.mpi);
        hybrid_times.push_back(row.hybrid);
        l2_errors.push_back(row.l2_error);
    }

    // Plot 1: Timing vs n (with log2 scaling on both axes)
    std::vector<std::vector<double>> timing_data = {serial_times, omp_times, mpi_times, hybrid_times};
    std::vector<std::string> timing_labels = {"Serial", "OMP", "MPI", "Hybrid"};

    Plotter::writeDataFile("test/plots/timing_vs_n.dat", n_values, timing_data, timing_labels);
    Plotter::createGnuplotScript("test/plots/timing_vs_n.gp", "test/plots/timing_vs_n.dat",
                                       "Timing vs Grid Size (n)",
                                       "n", "Time (s)", timing_labels, true, true, "test/plots/timing_vs_n.png");
    Plotter::printDataSummary("Timing vs Grid Size (n)", n_values, timing_data, timing_labels);

    // Plot 2: L2 Error vs n
    std::vector<std::vector<double>> error_data = {l2_errors};
    std::vector<std::string> error_labels = {"L2 Error"};

    Plotter::writeDataFile("test/plots/l2error_vs_n.dat", n_values, error_data, error_labels);
    Plotter::createGnuplotScript("test/plots/l2error_vs_n.gp", "test/plots/l2error_vs_n.dat",
                                       "L2 Error vs Grid Size (n)",
                                       "n", "L2 Error", error_labels, true, true, "test/plots/l2error_vs_n.png");

    // Plot 3: Timing vs h (with log2 scaling on both axes)
    Plotter::writeDataFile("test/plots/timing_vs_h.dat", h_values, timing_data, timing_labels);
    Plotter::createGnuplotScript("test/plots/timing_vs_h.gp", "test/plots/timing_vs_h.dat",
                                       "Timing vs Grid Spacing (h)",
                                       "h = 1/(n-1)", "Time (s)", timing_labels, true, true, "test/plots/timing_vs_h.png");
    Plotter::printDataSummary("Timing vs Grid Spacing (h)", h_values, timing_data, timing_labels);

    // Plot 4: L2 Error vs h
    Plotter::writeDataFile("test/plots/l2error_vs_h.dat", h_values, error_data, error_labels);
    Plotter::createGnuplotScript("test/plots/l2error_vs_h.gp", "test/plots/l2error_vs_h.dat",
                                       "L2 Error vs Grid Spacing (h)",
                                       "h = 1/(n-1)", "L2 Error", error_labels, true, true, "test/plots/l2error_vs_h.png");
}

void scalabilityTest()
{
    std::vector<int> processes = {1, 2, 4, 8};
    std::vector<int> grid_sizes = {56, 60, 64};

    std::vector<double> proc_double;
    for (int p : processes)
    {
        proc_double.push_back(static_cast<double>(p));
    }

    std::vector<std::vector<double>> scalability_data;
    std::vector<std::string> scalability_labels;

    for (int n : grid_sizes)
    {
        std::vector<double> timings;

        // Read all CSV files in data directory
        for (const auto &entry : fs::directory_iterator("test/data"))
        {
            if (entry.path().extension() == ".csv")
            {
                std::vector<DataRow> data = CSVReader::readCSV(entry.path().string());

                // Find row with matching n value
                for (const auto &row : data)
                {
                    if (row.n == n)
                    {
                        timings.push_back(row.hybrid);
                        break;
                    }
                }
            }
        }

        if (timings.size() == processes.size())
        {
            scalability_data.push_back(timings);
            scalability_labels.push_back("n=" + std::to_string(n));
        }
    }

    Plotter::writeDataFile("test/plots/scalability.dat", proc_double, scalability_data, scalability_labels);
    Plotter::createGnuplotScript("test/plots/scalability.gp", "test/plots/scalability.dat",
                                       "Scalability Test",
                                       "Number of Processes", "Time (s)",
                                       scalability_labels, true, false);
    Plotter::printDataSummary("Scalability Test", proc_double, scalability_data, scalability_labels);
}

#endif // End of PLOT_HPP guard