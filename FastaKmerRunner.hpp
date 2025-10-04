#pragma once
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <functional>
#include <iostream>
#include <cctype>

namespace fs = std::filesystem;
inline bool has_fna_ext(const fs::path& p) {
    auto ext = p.extension().string();
    for (auto& ch : ext) ch = std::tolower(static_cast<unsigned char>(ch));
    return (ext == ".fna"); 
}

inline std::vector<fs::path> list_fna_files(const fs::path& root = "Genomas") {
    std::vector<fs::path> files;
    if (!fs::exists(root)) return files;
    for (auto const& entry : fs::directory_iterator(root)) {
        if (entry.is_regular_file() && has_fna_ext(entry.path())) {
            files.push_back(entry.path());
        }
    }
    return files;
}

inline void uppercase(std::string& s) {
    for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
}


inline void for_each_kmer_in_fasta(const fs::path& fna_path, size_t k,const std::function<void(const std::string&)>& on_kmer)
{
    std::ifstream in(fna_path);
    if (!in) {
        std::cerr << "[WARN] No se pudo abrir: " << fna_path << "\n";
        return;
    }

    std::string line;
    std::string buf;   

    auto flush_kmers = [&](const std::string& seq) {
        if (seq.size() < k) return;
        
        for (size_t i = 0; i + k <= seq.size(); ++i) {
            on_kmer(seq.substr(i, k));
        }
    };

    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '>') {
            if (!buf.empty()) {
                flush_kmers(buf);
                buf.clear();
            }
        } else {
            uppercase(line);
            size_t start = 0;
            for (size_t i = 0; i < line.size(); ++i) {
                char c = line[i];
                if (c!='A' && c!='C' && c!='G' && c!='T') {
                    // procesa tramo válido anterior
                    if (i > start) {
                        std::string chunk = line.substr(start, i - start);
                        buf.append(chunk);
                        flush_kmers(buf);
                        if (buf.size() >= k-1) buf = buf.substr(buf.size() - (k-1));
                        else buf.clear();
                    }
                    start = i + 1;
                }
            }
            
            if (start < line.size()) {
                buf.append(line.substr(start));
                flush_kmers(buf);
                if (buf.size() >= k-1) buf = buf.substr(buf.size() - (k-1));
                else buf.clear();
            }
        }
    }
    if (!buf.empty()) {
        //seguridad:
        if (buf.size() >= k) {
            for (size_t i = 0; i + k <= buf.size(); ++i) {
                on_kmer(buf.substr(i, k));
            }
        }
    }
}
