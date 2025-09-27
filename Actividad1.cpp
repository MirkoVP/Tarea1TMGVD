// kmer_groundtruth_act1.cpp
// Adaptado para: k ∈ {21, 31} fijos y lectura automática desde carpeta Genomas/.
// Mantiene: --phi (umbral HH), --limit (N máx. de k-mers por k), --out-prefix (prefijo salidas).

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cctype>
#include <cstdint>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <filesystem>

#if defined(__linux__)
#include <sys/resource.h>
#include <unistd.h>
#endif

namespace fs = std::filesystem;

// -------- Utilidad memoria pico (Linux) --------
static size_t getPeakRSSBytes() {
#if defined(__linux__)
    std::ifstream in("/proc/self/status");
    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind("VmHWM:", 0) == 0) {
            std::istringstream iss(line.substr(6));
            size_t kb = 0; std::string unit;
            iss >> kb >> unit;
            return kb * 1024ULL;
        }
    }
    struct rusage ru{};
    if (getrusage(RUSAGE_SELF, &ru) == 0) return (size_t)ru.ru_maxrss * 1024ULL;
#endif
    return 0;
}

// -------- Mapping 2 bits/base --------
inline int encode_base(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return -1; // inválida (N u otra)
    }
}
inline char decode_base(int b) {
    static const char* ABC = "ACGT";
    return ABC[b & 3];
}
std::string decode_kmer(uint64_t code, int k) {
    std::string s(k, 'A');
    for (int i = k - 1; i >= 0; --i) {
        s[i] = decode_base(int(code & 3ULL));
        code >>= 2;
    }
    return s;
}

// -------- Rodante por k --------
struct KmerRolling {
    int k;
    uint64_t mask;
    uint64_t fwd=0, rev=0;
    int L=0;
    bool saturated=false;

    std::unordered_map<uint64_t, uint64_t> counts;
    uint64_t totalKmers = 0; // N

    explicit KmerRolling(int kk): k(kk) {
        mask = (k==32) ? ~0ULL : ((1ULL << (2*k)) - 1ULL);
    }

    bool push_base(int b, int cb, uint64_t &canon_out) {
        if (b < 0) { fwd=0; rev=0; L=0; return false; }
        fwd = ((fwd << 2) | (uint64_t)b) & mask;
        rev = (rev >> 2) | ( (uint64_t)cb << (2*(k-1)) );
        if (L < k) { ++L; if (L < k) return false; }
        canon_out = std::min(fwd, rev);
        return true;
    }
};

// -------- CLI simple (phi, limit, out-prefix) --------
struct Args {
    double phi = 1e-4;            // --phi
    uint64_t limitPerK = 0;       // --limit
    std::string outPrefix = "gt"; // --out-prefix
};
static void usage(const char* prog) {
    std::cerr <<
    "Uso: " << prog << " [--phi X] [--limit N] [--out-prefix NAME]\n"
    "  * k se fija en 21 y 31.\n"
    "  * Lee todos los archivos GCA_*.fna en carpeta Genomas/.\n";
}
static bool parse_args(int argc, char** argv, Args& a) {
    for (int i=1; i<argc; ++i) {
        std::string s = argv[i];
        if (s == "--phi" && i+1 < argc) {
            a.phi = std::stod(argv[++i]);
        } else if (s == "--limit" && i+1 < argc) {
            a.limitPerK = std::stoull(argv[++i]);
        } else if (s == "--out-prefix" && i+1 < argc) {
            a.outPrefix = argv[++i];
        } else if (!s.empty() && s[0]=='-') {
            std::cerr << "Opción no reconocida: " << s << "\n";
            return false;
        } else {
            std::cerr << "Argumento inesperado: " << s << "\n";
            return false;
        }
    }
    return true;
}

// -------- Escrituras --------
static void write_counts_csv(const std::string& path, const KmerRolling& kr) {
    std::ofstream out(path);
    out << "k,code,kmer,count,freq\n";
    for (const auto& [code, cnt] : kr.counts) {
        double freq = (kr.totalKmers==0) ? 0.0 : (double)cnt / (double)kr.totalKmers;
        out << kr.k << "," << code << "," << decode_kmer(code, kr.k) << ","
            << cnt << "," << std::setprecision(10) << freq << "\n";
    }
}
static void write_heavy_csv(const std::string& path, const KmerRolling& kr, double phi) {
    std::ofstream out(path);
    out << "k,code,kmer,count,freq\n";
    const uint64_t N = kr.totalKmers;
    const double thr = phi * (double)N;
    for (const auto& [code, cnt] : kr.counts) {
        if ((double)cnt >= thr) {
            double freq = (N==0) ? 0.0 : (double)cnt / (double)N;
            out << kr.k << "," << code << "," << decode_kmer(code, kr.k) << ","
                << cnt << "," << std::setprecision(10) << freq << "\n";
        }
    }
}

// -------- Recolectar archivos Genomas/GCA_*.fna --------
static std::vector<fs::path> collect_genomes() {
    std::vector<fs::path> files;
    fs::path root("Genomas");
    if (!fs::exists(root) || !fs::is_directory(root)) {
        std::cerr << "[ERROR] No existe carpeta 'Genomas/'.\n";
        return files;
    }
    for (auto &entry : fs::directory_iterator(root)) {
        if (!entry.is_regular_file()) continue;
        auto p = entry.path();
        // patrón: nombre que empiece con GCA_ y termine en .fna
        if (p.filename().string().rfind("GCA_", 0) == 0 && p.extension() == ".fna") {
            files.push_back(p);
        }
    }
    if (files.empty()) {
        std::cerr << "[ADVERTENCIA] No se encontraron archivos GCA_*.fna en Genomas/.\n";
    }
    std::sort(files.begin(), files.end());
    return files;
}

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    Args args;
    if (!parse_args(argc, argv, args)) { usage(argv[0]); return 1; }

    // k fijos
    const std::vector<int> KS = {21, 31};
    std::vector<KmerRolling> trackers;
    trackers.reserve(KS.size());
    for (int k : KS) trackers.emplace_back(k);

    // Recolectar entradas desde Genomas/
    auto inputs = collect_genomes();
    if (inputs.empty()) return 2;

    auto t0 = std::chrono::high_resolution_clock::now();

    // Procesar todos los .fna
    for (const auto& inpath : inputs) {
        std::ifstream in(inpath);
        if (!in) { std::cerr << "No se pudo abrir " << inpath << "\n"; return 3; }
        std::cerr << "[Leyendo] " << inpath << "\n";

        auto reset_all = [&]() {
            for (auto& kr : trackers) { kr.fwd=0; kr.rev=0; kr.L=0; }
        };

        std::string line;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') { reset_all(); continue; }
            for (char ch : line) {
                char c = std::toupper(static_cast<unsigned char>(ch));
                int b  = encode_base(c);
                int cb = (b >= 0) ? (b ^ 0b11) : -1;
                for (auto& kr : trackers) {
                    if (kr.saturated) continue;
                    uint64_t canon;
                    if (kr.push_base(b, cb, canon)) {
                        ++kr.totalKmers;
                        ++kr.counts[canon];
                        if (args.limitPerK > 0 && kr.totalKmers >= args.limitPerK) {
                            kr.saturated = true;
                        }
                    }
                }
            }
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    double secs = std::chrono::duration<double>(t1 - t0).count();
    size_t peak = getPeakRSSBytes();

    // Salidas
    std::cerr << std::fixed << std::setprecision(3);
    std::cerr << "[Tiempo] " << secs << " s\n";
    if (peak) std::cerr << "[Memoria pico] " << (peak / (1024.0*1024.0)) << " MB\n";
    std::cerr << std::scientific << std::setprecision(2);
    std::cerr << "[phi] " << args.phi << "\n";
    std::cerr << std::defaultfloat; 

    for (const auto& kr : trackers) {
        std::ostringstream cpath, hpath;
        cpath << args.outPrefix << "_counts_k" << kr.k << ".csv";
        hpath << args.outPrefix << "_heavy_k" << kr.k << "_phi" << std::scientific << args.phi << ".csv";
        write_counts_csv(cpath.str(), kr);
        write_heavy_csv(hpath.str(),  kr, args.phi);

        uint64_t hh = 0;
        double thr = args.phi * (double)kr.totalKmers;
        for (const auto& [_, cnt] : kr.counts) if ((double)cnt >= thr) ++hh;

        std::cerr << std::defaultfloat;
        std::cerr << "[k=" << kr.k << "] N=" << kr.totalKmers
                  << "  únicos=" << kr.counts.size()
                  << "  HH(phi)=" << hh << "  thr=" << thr << "\n";
    }
    return 0;
}
