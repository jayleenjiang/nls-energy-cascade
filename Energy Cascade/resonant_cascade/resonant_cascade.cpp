#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <set>
#include <random>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class ResonantSetSimulator
{
private:
    struct Coordinate
    {
        vector<int> z; // z1, ..., z_{N-1} where each zi ∈ {0,1,i,1+i}
        int generation;

        bool operator<(const Coordinate &other) const
        {
            return z < other.z;
        }
    };

    struct NuclearFamily
    {
        Coordinate parent1, parent2;
        Coordinate child_over, child_under; // over-achiever and under-achiever
        double angle;                       // θ(F) for this family

        complex<double> p1_freq, p2_freq; // Placement in frequency space
        complex<double> co_freq, cu_freq;
    };

    int N;                               // Number of generations
    double R;                            // Base frequency scale
    vector<set<Coordinate>> generations; // Λ_1, ..., Λ_N
    vector<NuclearFamily> families;
    map<Coordinate, complex<double>> placement; // f: Λ → ℂ
    map<Coordinate, int> coord_to_index;        // Map coordinates to mode indices

    // Energy and mass at each mode
    VectorXd energy;
    VectorXd mass;
    VectorXcd amplitudes; // Complex amplitudes c_k

    // Parameters
    double gamma = 0.1;
    double T1 = 10.0;
    double T3 = 0.01;

    // Convert coordinate element to complex number
    complex<double> elem_to_complex(int e)
    {
        if (e == 0)
            return complex<double>(0, 0);
        if (e == 1)
            return complex<double>(1, 0);
        if (e == 2)
            return complex<double>(0, 1); // i
        if (e == 3)
            return complex<double>(1, 1); // 1+i
        return complex<double>(0, 0);
    }

public:
    ResonantSetSimulator(int n_gen, double r_scale) : N(n_gen), R(r_scale)
    {
        constructResonantSet();
    }

    void constructResonantSet()
    {
        // Build combinatorial model: Λ_j = S_2^{j-1} × S_1^{N-j}
        // where S_1 = {1, i} and S_2 = {0, 1+i}
        generations.resize(N + 1); // 1-indexed

        for (int j = 1; j <= N; j++)
        {
            generateGeneration(j);
        }

        // Create nuclear families and compute placement
        createNuclearFamilies();
        computePlacement();

        // Initialize mode arrays
        int total_modes = 0;
        for (int j = 1; j <= N; j++)
        {
            total_modes += generations[j].size();
        }

        energy.resize(total_modes);
        mass.resize(total_modes);
        amplitudes.resize(total_modes);

        // Map coordinates to indices
        int idx = 0;
        for (int j = 1; j <= N; j++)
        {
            for (const auto &coord : generations[j])
            {
                coord_to_index[coord] = idx++;
            }
        }
    }

    void generateGeneration(int j)
    {
        // Generate all (N-1)-tuples for generation j
        vector<int> tuple(N - 1);
        generateTuplesRecursive(j, 0, tuple, generations[j]);
    }

    void generateTuplesRecursive(int gen, int pos, vector<int> &tuple, set<Coordinate> &result)
    {
        if (pos == N - 1)
        {
            Coordinate coord;
            coord.z = tuple;
            coord.generation = gen;
            result.insert(coord);
            return;
        }

        vector<int> choices;
        if (pos < gen - 1)
        {
            choices = {0, 3}; // S_2 = {0, 1+i}
        }
        else
        {
            choices = {1, 2}; // S_1 = {1, i}
        }

        for (int c : choices)
        {
            tuple[pos] = c;
            generateTuplesRecursive(gen, pos + 1, tuple, result);
        }
    }

    void createNuclearFamilies()
    {
        // For each generation pair (j, j+1), create families
        for (int j = 1; j < N; j++)
        {
            // Each family is determined by fixing z_1,...,z_{j-1} and z_{j+1},...,z_{N-1}
            // and varying z_j over S = {0,1,i,1+i}

            map<vector<int>, vector<Coordinate>> groups;

            // Group elements by their fixed parts
            for (const auto &elem : generations[j])
            {
                vector<int> key;
                for (int k = 0; k < N - 1; k++)
                {
                    if (k != j - 1)
                        key.push_back(elem.z[k]);
                }
                groups[key].push_back(elem);
            }

            for (const auto &elem : generations[j + 1])
            {
                vector<int> key;
                for (int k = 0; k < N - 1; k++)
                {
                    if (k != j - 1)
                        key.push_back(elem.z[k]);
                }
                groups[key].push_back(elem);
            }

            // Create families from groups
            for (auto &[key, members] : groups)
            {
                if (members.size() == 4)
                {
                    NuclearFamily fam;

                    // Identify parents (generation j) and children (generation j+1)
                    for (const auto &m : members)
                    {
                        if (m.generation == j)
                        {
                            if (m.z[j - 1] == 1)
                                fam.parent1 = m;
                            else if (m.z[j - 1] == 2)
                                fam.parent2 = m; // i
                        }
                        else
                        {
                            if (m.z[j - 1] == 0)
                                fam.child_under = m;
                            else if (m.z[j - 1] == 3)
                                fam.child_over = m; // 1+i
                        }
                    }

                    // Assign angle - perturb from π/2 for better cascade
                    fam.angle = M_PI / 2 + 0.1 * (drand48() - 0.5);

                    families.push_back(fam);
                }
            }
        }
    }

    void computePlacement()
    {
        // Initial placement for generation 1
        for (const auto &coord : generations[1])
        {
            complex<double> prod(R, 0);
            for (int k = 0; k < N - 1; k++)
            {
                prod *= elem_to_complex(coord.z[k]);
            }
            placement[coord] = prod;
        }

        // Recursive placement for later generations
        for (int j = 2; j <= N; j++)
        {
            for (auto &fam : families)
            {
                if (fam.parent1.generation == j - 1)
                {
                    complex<double> p1 = placement[fam.parent1];
                    complex<double> p2 = placement[fam.parent2];

                    complex<double> phase = exp(complex<double>(0, fam.angle));

                    // Place children according to rectangle rule
                    placement[fam.child_over] = 0.5 * ((1.0 + phase) * p1 + (1.0 - phase) * p2);
                    placement[fam.child_under] = 0.5 * ((1.0 + phase) * p1 - (1.0 - phase) * p2);

                    fam.p1_freq = p1;
                    fam.p2_freq = p2;
                    fam.co_freq = placement[fam.child_over];
                    fam.cu_freq = placement[fam.child_under];
                }
            }
        }
    }

    void initializeEnergy(double initial_energy)
    {
        // Put all energy in first generation
        energy.setZero();
        mass.setZero();
        amplitudes.setZero();

        int modes_in_gen1 = generations[1].size();
        for (const auto &coord : generations[1])
        {
            int idx = coord_to_index[coord];
            energy(idx) = initial_energy / modes_in_gen1;
            mass(idx) = sqrt(energy(idx));
            amplitudes(idx) = sqrt(mass(idx)) * exp(complex<double>(0, 2 * M_PI * drand48()));
        }
    }

    void evolveResonantDynamics(double dt, double T_final, const string &output_file)
    {
        ofstream out(output_file);
        out << "time,generation,total_energy,max_energy\n";

        int steps = T_final / dt;
        int output_freq = steps / 100; // Output 100 data points

        for (int step = 0; step <= steps; step++)
        {
            double t = step * dt;

            // Output statistics
            if (step % output_freq == 0)
            {
                for (int j = 1; j <= N; j++)
                {
                    double gen_energy = 0, max_energy = 0;
                    for (const auto &coord : generations[j])
                    {
                        int idx = coord_to_index[coord];
                        gen_energy += energy(idx);
                        max_energy = max(max_energy, energy(idx));
                    }
                    out << t << "," << j << "," << gen_energy << "," << max_energy << "\n";
                }
            }

            // Evolve using resonant interactions within nuclear families
            VectorXcd new_amplitudes = amplitudes;

            for (const auto &fam : families)
            {
                int i1 = coord_to_index[fam.parent1];
                int i2 = coord_to_index[fam.parent2];
                int io = coord_to_index[fam.child_over];
                int iu = coord_to_index[fam.child_under];

                // Resonant four-wave interaction
                complex<double> c1 = amplitudes(i1);
                complex<double> c2 = amplitudes(i2);
                complex<double> co = amplitudes(io);
                complex<double> cu = amplitudes(iu);

                // Energy transfer rates (simplified model)
                double transfer_rate = 0.01;

                // Parents lose energy to children
                new_amplitudes(i1) += dt * transfer_rate * conj(c2) * co * cu;
                new_amplitudes(i2) += dt * transfer_rate * conj(c1) * co * cu;

                // Over-achiever gets most energy
                new_amplitudes(io) += dt * transfer_rate * 0.9 * c1 * c2 * conj(cu);

                // Under-achiever gets little energy
                new_amplitudes(iu) += dt * transfer_rate * 0.1 * c1 * c2 * conj(co);
            }

            // Add dissipation at boundaries
            for (const auto &coord : generations[1])
            {
                int idx = coord_to_index[coord];
                new_amplitudes(idx) *= (1 - gamma * dt);
                // Add forcing
                new_amplitudes(idx) += sqrt(dt * gamma * T1) *
                                       complex<double>(randn(), randn()) / sqrt(2.0);
            }

            for (const auto &coord : generations[N])
            {
                int idx = coord_to_index[coord];
                new_amplitudes(idx) *= (1 - gamma * dt);
                // Add dissipation
                new_amplitudes(idx) += sqrt(dt * gamma * T3) *
                                       complex<double>(randn(), randn()) / sqrt(2.0);
            }

            amplitudes = new_amplitudes;

            // Update energy and mass
            for (int i = 0; i < amplitudes.size(); i++)
            {
                mass(i) = abs(amplitudes(i));
                energy(i) = mass(i) * mass(i);
            }
        }

        out.close();
    }

    void analyzeEnergyDistribution(const string &output_file)
    {
        ofstream out(output_file);
        out << "generation,mode_index,frequency_magnitude,energy\n";

        for (int j = 1; j <= N; j++)
        {
            for (const auto &coord : generations[j])
            {
                int idx = coord_to_index[coord];
                double freq_mag = abs(placement[coord]);
                out << j << "," << idx << "," << freq_mag << "," << energy(idx) << "\n";
            }
        }
        out.close();
    }

    
    double computeCascadeEfficiency(int s = 2)
    {
        // Compare weighted norms between early and late generations
        double early_norm = 0, late_norm = 0;

        for (const auto &coord : generations[min(3, N)])
        {
            int idx = coord_to_index[coord];
            double freq_mag = abs(placement[coord]);
            early_norm += pow(freq_mag, 2 * s) * energy(idx);
        }

        for (const auto &coord : generations[max(1, N - 2)])
        {
            int idx = coord_to_index[coord];
            double freq_mag = abs(placement[coord]);
            late_norm += pow(freq_mag, 2 * s) * energy(idx);
        }

        return late_norm / (early_norm + 1e-10);
    }

private:
    double randn()
    {
        static random_device rd;
        static mt19937 gen(rd());
        static normal_distribution<> dis(0, 1);
        return dis(gen);
    }
};

int main()
{
    // Parameters
    int N = 7;      // Number of generations
    double R = 100; // Base frequency scale
    double T_final = 1000.0;
    double dt = 0.01;

    // Create and run simulation
    ResonantSetSimulator sim(N, R);
    sim.initializeEnergy(1.0);

    cout << "Running resonant cascade simulation..." << endl;
    sim.evolveResonantDynamics(dt, T_final, "cascade_evolution.csv");

    cout << "Analyzing energy distribution..." << endl;
    sim.analyzeEnergyDistribution("energy_distribution.csv");

    double efficiency = sim.computeCascadeEfficiency();
    cout << "Cascade efficiency factor: " << efficiency << endl;
    cout << "Theoretical factor for N=" << N << ": " << pow(2, N - 5) << endl;

    return 0;
}