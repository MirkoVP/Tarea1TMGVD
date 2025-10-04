//calibrar countsketch con ground truth
//compilación: g++ -std=c++17 -O3 EvaluacionTS.cpp -o evTS
//ejecución: ./evCS
//nota: ajustar tamaño del sketch, valor de k y ruta a ground truth manualmente

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <chrono>
#include <iomanip>
#include "FastaKmerRunner.hpp"
#include "TowerSketch.hpp"

struct GTEntry {
    std::string kmer;
    int64_t count;
};
std::vector<GTEntry> read_ground_truth(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) {
        throw std::runtime_error("No se pudo abrir " + filename);
    }
    std::vector<GTEntry> data;
    std::string line;
    bool header = true;
    while (std::getline(in, line)) {
        if (header) { header = false; continue; }
        std::stringstream ss(line);
        std::string field;
        std::vector<std::string> cols;
        while (std::getline(ss, field, ',')) cols.push_back(field);
        if (cols.size() < 5) continue;
        GTEntry e;
        e.kmer  = cols[2];
        e.count = std::stoll(cols[3]);
        data.push_back(e);
    }
    return data;
}

int main() {
    //tamaño tower sketch y valor de k
    const int d= 10;
    const int base_w = 4'000'000;
    const int L = 6; //niveles
    const size_t K = 31;

    TowerSketch_CMCU ts(d, base_w, L, 0xCAFEBABEu);
    auto files = list_fna_files("Genomas");
    if (files.empty()) {
        std::cerr << "[INFO] No se encontraron .fna en ./Genomas\n";
        return 0;
    }

    std::cout << "[INFO] Archivos: " << files.size() << "\n";
    std::cout << "[INFO] d=" << d << " base_w=" << base_w << " L=" << L << " K=" << K << "\n";

    uint64_t total_kmers = 0;
    auto t0 = std::chrono::high_resolution_clock::now();

    for (const auto& f : files) {
        std::cout << "[INFO] Procesando: " << f << "\n";
        for_each_kmer_in_fasta(f, K, [&](const std::string& kmer) {
            ts.insert(kmer);
            ++total_kmers;
        });
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    double secs = std::chrono::duration<double>(t1 - t0).count();

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "[STATS] k-mers procesados: " << total_kmers << "\n";
    std::cout << "[STATS] Tiempo (s):        " << secs << "\n";
    std::cout << "[STATS] k-mers/s:          " << (secs>0? total_kmers/secs : 0) << "\n";

    auto gt = read_ground_truth("FULL_counts_k31.csv");

    double sum_abs = 0.0, sum_sq = 0.0;
    int64_t max_abs = 0;
    size_t n = 0;

    for (const auto& e : gt) {
        double est = static_cast<double>(ts.estimate(e.kmer));
        double diff = est - e.count;
        double abs_err = std::fabs(diff);

        sum_abs += abs_err;
        sum_sq  += diff * diff;
        if (abs_err > max_abs) max_abs = static_cast<int64_t>(abs_err);
        ++n;
    }

    double mae  = sum_abs / n;
    double rmse = std::sqrt(sum_sq / n);

    std::cout << "[EVAL] N=" << n << "\n";
    std::cout << "[EVAL] MAE = " << mae << "\n";
    std::cout << "[EVAL] RMSE= " << rmse << "\n";
    std::cout << "[EVAL] MaxAbsError= " << max_abs << "\n";

    double phi = 1e-6;
    double thr = phi * total_kmers;
    auto gtHH = read_ground_truth("FULL_heavy_k31_phi1.000000e-06.csv");

    struct HH { std::string kmer; int64_t gt; double est; double est_freq; };
    std::vector<HH> hh;

    for (const auto& e : gt) {
        double est = static_cast<double>(ts.estimate(e.kmer));
        if (est >= thr){
            hh.push_back({ e.kmer, e.count, est, est / static_cast<double>(total_kmers) });
        }
    }
    std::unordered_map<std::string, int64_t> gtHHmap;
    for (const auto& e : gtHH) {
        gtHHmap[e.kmer] = e.count;
    }
    size_t TP = 0, FP = 0, FN = 0;
    double sumAbsErrHH = 0.0, sumRelErrHH = 0.0;
    size_t hhCount = 0;
    for (const auto& h : hh) {
        auto it = gtHHmap.find(h.kmer);
        if (it != gtHHmap.end()) {
            TP++;
            double absErr = std::fabs(h.est - it->second);
            double relErr = absErr / (it->second > 0 ? (double)it->second : 1.0);
            sumAbsErrHH += absErr;
            sumRelErrHH += relErr;
            hhCount++;
        } else {
            FP++;
        }
    }
    FN = gtHH.size() - TP;
    double precision = (TP + FP > 0) ? (double)TP / (TP + FP) : 0.0;
    double recall    = (TP + FN > 0) ? (double)TP / (TP + FN) : 0.0;
    double f1        = (precision + recall > 0) ? 2 * precision * recall / (precision + recall) : 0.0;

    double maeHH = (hhCount > 0) ? sumAbsErrHH / hhCount : 0.0;
    double mreHH = (hhCount > 0) ? sumRelErrHH / hhCount : 0.0;

    std::cout << "\n[HH METRICS]\n";
    std::cout << "TP=" << TP << " FP=" << FP << " FN=" << FN << "\n";
    std::cout << "Precision=" << precision << " Recall=" << recall << " F1=" << f1 << "\n";
    std::cout << "MAE(HH)=" << maeHH << "  MRE(HH)=" << mreHH << "\n";

    return 0;
}