#include <iostream>
#include <complex>
#include <vector>
#include <map>
#include <boost/rational.hpp>
#include <cmath>
#include <iomanip>
#include <fstream> 

using namespace std;
using Rational = boost::rational<long long>; 
using RationalComplex = complex<Rational>;

// elements in Sigma
struct SigmaElement {
    vector<RationalComplex> coords;
    int generation;
    int index;
};

// a nuclear family
struct NuclearFamily {
    int parent1_idx, parent2_idx;  // F_1 and F_i
    int child1_idx, child2_idx;    // F_0 and F_{1+i}
};

class FrequencyEmbedding {
private:
    int N; // number of generations
    vector<vector<SigmaElement>> sigma; // all elements in generation j
    vector<map<int, RationalComplex>> f; // maps element index to complex frequency
    
    RationalComplex S1[2] = { RationalComplex(1, 0), RationalComplex(0, 1) }; // {1, i}
    RationalComplex S2[2] = { RationalComplex(0, 0), RationalComplex(1, 1) }; // {0, 1+i}
 
    
public:
    FrequencyEmbedding(int n) : N(n), sigma(n+1), f(n+1) {}
    
    int phi(const RationalComplex& z) {
        if (z == S2[0] || z == S1[0]) return 0; // 0 or 1
        if (z == S2[1] || z == S1[1]) return 1; // 1+i or i
        return -1;
    }
    
    int computeIndex(const vector<RationalComplex>& coords) {
        int idx = 0;
        if (N <= 1) return 0;
        int power = 1 << (N - 2);

        //  treat coordinates as binary digits to compute unique index
        for (const auto& z : coords) {
            idx += phi(z) * power;
            power >>= 1;
        }
        return idx;
    }
    


    void generateSigmaJ(int j) {
        if (j < 1 || j > N || N <= 1) return;
        vector<RationalComplex> current(N-1); // N-1 coordinates per element
        generateSigmaJRecursive(j, 0, current);
    }
    
    // generate all 2^(N-1) possible combinations
    void generateSigmaJRecursive(int j, int pos, vector<RationalComplex>& current) {
        if (pos == N-1) {
            // Base case: filled all N-1 positions
            sigma[j].push_back({current, j, computeIndex(current)});
            return;
        }

        // Choose set based on position
        RationalComplex* set = (pos < j-1) ? S2 : S1;

        // try both elements from the chosen set
        for (int i = 0; i < 2; i++) {
            current[pos] = set[i];
            generateSigmaJRecursive(j, pos+1, current);
        }
    }
    
    // first generation 
    void initializeF1(int R) {
        for (const auto& elem : sigma[1]) {
            RationalComplex product(R, 0);
            for (const auto& z : elem.coords) product *= z;
            f[1][elem.index] = product;
        }
    }

    
    NuclearFamily findFamily(const SigmaElement& child) {
        NuclearFamily family;
        if (child.generation <= 1) return family;

        int j = child.generation - 1; // parent generation

        // find parents by changing the (j-1)th coordinate to S1 values
        vector<RationalComplex> parent_coords = child.coords;
        parent_coords[j-1] = S1[0]; // set to 1
        family.parent1_idx = computeIndex(parent_coords);

        parent_coords[j-1] = S1[1]; // set to i
        family.parent2_idx = computeIndex(parent_coords);
        
        // find sibling by flipping the (j-1)th coordinate
        vector<RationalComplex> child_coords = child.coords;
        if (child.coords[j-1] == S2[0]) { // if 0
            family.child1_idx = child.index;
            child_coords[j-1] = S2[1]; // sibling has 1+i
            family.child2_idx = computeIndex(child_coords);
        } else {  // if 1+i
            family.child2_idx = child.index;
            child_coords[j-1] = S2[0]; // sibling has 0
            family.child1_idx = computeIndex(child_coords);
        }
        return family;
    }
    
    pair<Rational, Rational> getPythagoreanAngle(int choice) {
        switch (choice % 4) {
            case 0: return {Rational(3, 5), Rational(4, 5)};
            case 1: return {Rational(5, 13), Rational(12, 13)};
            case 2: return {Rational(8, 17), Rational(15, 17)};
            case 3: return {Rational(7, 25), Rational(24, 25)};
            default: return {Rational(1, 1), Rational(0, 1)};
        }
    }

    void computeNextGeneration(int j) {
    if (j < 1 || j >= N) return;
    
    // map each parent pair to a unique * angle
    map<pair<int,int>, int> family_angle_map;  
    int angle_choice = 0;
    
    for (const auto& child : sigma[j+1]) {
        NuclearFamily family = findFamily(child);
        
        // create a unique key for this family 
        auto parent_pair = make_pair(family.parent1_idx, family.parent2_idx);
        
        // Only assign a new angle if we haven't seen this family before
        if (family_angle_map.find(parent_pair) == family_angle_map.end()) {
            family_angle_map[parent_pair] = angle_choice++;
        }
        //  Get the angle for this family 
        auto [cos_theta, sin_theta] = getPythagoreanAngle(family_angle_map[parent_pair]);
        RationalComplex exp_theta(cos_theta, sin_theta);
        
        // get parent frenquencies
        RationalComplex f_1 = f[j][family.parent1_idx];
        RationalComplex f_i = f[j][family.parent2_idx];
        
        // compute rectangles 
        RationalComplex coeff_A = (RationalComplex(1,0) + exp_theta) / Rational(2);
        RationalComplex coeff_B = (RationalComplex(1,0) - exp_theta) / Rational(2);

        if (child.coords[j-1] == S2[0]) { 
            f[j+1][child.index] = coeff_B * f_1 + coeff_A * f_i;
        } else {
            f[j+1][child.index] = coeff_A * f_1 + coeff_B * f_i;
        }
    }
}
    
    void buildEmbedding(int R) {
        for (int j = 1; j <= N; j++) generateSigmaJ(j); 
        initializeF1(R);
        for (int j = 1; j < N; j++) computeNextGeneration(j);
    }

    void exportFamiliesForPlotting(int parent_gen, int R) {
        if (parent_gen < 1 || parent_gen >= N) return;
        int child_gen = parent_gen + 1;
        string filename = "rectangles_N" + to_string(N) + "_R" + to_string(R) + "_f" + to_string(parent_gen) + "_to_f" + to_string(child_gen) + ".csv";
        ofstream outfile(filename);

        if (!outfile.is_open()) {
            cerr << "Error: Could not open file " << filename << " for writing." << endl;
            return;
        }
        outfile << "family_id,point_type,real,imag" << endl;

        int family_id = 0;
        for (const auto& child : sigma[child_gen]) {
            if (child.coords[parent_gen - 1] != S2[0]) continue;
            NuclearFamily family = findFamily(child);
            RationalComplex p1 = f[parent_gen][family.parent1_idx];
            RationalComplex p2 = f[parent_gen][family.parent2_idx];
            RationalComplex c1 = f[child_gen][family.child1_idx];
            RationalComplex c2 = f[child_gen][family.child2_idx];

            outfile << family_id << ",parent," << p1.real() << "," << p1.imag() << endl;
            outfile << family_id << ",parent," << p2.real() << "," << p2.imag() << endl;
            outfile << family_id << ",child," << c1.real() << "," << c1.imag() << endl;
            outfile << family_id << ",child," << c2.real() << "," << c2.imag() << endl;
            family_id++;
        }
        outfile.close();
        cout << "Exported data to " << filename << endl;
    }
};

int main() {
    int N, R;
    cout << "number of generations: ";
    cin >> N;
    cout << "scaling factor: ";
    cin >> R;

    if (N < 2) {
        cout << "N must be at least 2." << endl;
        return 1;
    }

    FrequencyEmbedding embedding(N);
    embedding.buildEmbedding(R);

    for (int j = 1; j < N; ++j) {
        embedding.exportFamiliesForPlotting(j, R);
    }
    
    return 0;
}

