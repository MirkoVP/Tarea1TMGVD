#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <iostream>
#include <iomanip>
#include <climits>
#include "murmurhash32.hpp"

class CountMinCU {
public:
    CountMinCU(int d, int w, uint32_t seed = 0x12345678u)
    : d_(d), w_(w), seed_(seed), table_(d, std::vector<uint64_t>(w, 0ULL)) {
        if (d_ <= 0 || w_ <= 0) throw std::invalid_argument("d y w deben ser > 0");
    }

    void insert(uint64_t key, uint64_t weight = 1ULL) {
        if (weight == 0ULL) return;
        idx_.resize(d_);
        uint64_t minv = ULLONG_MAX;
        for (int i = 0; i < d_; ++i) {
            uint32_t h = bucket_hash(i, key) % static_cast<uint32_t>(w_);
            idx_[i] = h;
            uint64_t v = table_[i][h];
            if (v < minv) minv = v;
        }
        for (int i = 0; i < d_; ++i) {
            uint32_t h = idx_[i];
            if (table_[i][h] == minv) {
                table_[i][h] += weight;
            }
        }
    }
    uint64_t estimate_u64(uint64_t key) const {
        uint64_t est = ULLONG_MAX;
        for (int i = 0; i < d_; ++i) {
            uint32_t h = bucket_hash(i, key) % static_cast<uint32_t>(w_);
            uint64_t v = table_[i][h];
            if (v < est) est = v;
        }
        return est == ULLONG_MAX ? 0ULL : est;
    }

    static uint64_t encode_kmer_2bits(const std::string& kmer) {
        if (kmer.size() > 32)
            throw std::invalid_argument("k > 32 no cabe en 64 bits (2 bits por símbolo).");
        uint64_t code = 0ULL;
        for (char ch : kmer) {
            code <<= 2;
            switch (ch) {
                case 'A': case 'a': break;
                case 'C': case 'c': code |= 0b01ULL; break;
                case 'G': case 'g': code |= 0b10ULL; break;
                case 'T': case 't': code |= 0b11ULL; break;
                default: throw std::invalid_argument("Símbolo inválido en k-mer (usa A/C/G/T).");
            }
        }
        return code;
    }
    static uint64_t reverse_complement_code(uint64_t code, size_t k) {
        uint64_t rc = 0ULL;
        for (size_t i = 0; i < k; i++) {
            uint64_t bits = code & 0b11ULL;
            code >>= 2;
            bits ^= 0b11ULL;
            rc = (rc << 2) | bits;
        }
        return rc;
    }
    static uint64_t encode_kmer_canonical(const std::string& kmer) {
        uint64_t fwd = encode_kmer_2bits(kmer);
        uint64_t rev = reverse_complement_code(fwd, kmer.size());
        return std::min(fwd, rev);
    }

private:
    int d_, w_;
    uint32_t seed_;
    std::vector<std::vector<uint64_t>> table_;
    mutable std::vector<uint32_t> idx_;

    inline uint32_t bucket_hash(int i, uint64_t key) const {
        uint32_t s = seed_ ^ (0x9e3779b9u * (uint32_t)(i + 1));
        return murmurhash(&key, s);
    }
};


class TowerSketch_CMCU {
public:
    TowerSketch_CMCU(int d, int base_w, int L, uint32_t seed = 0xBADC0FFEu)
    : d_(d), base_w_(base_w), L_(L), seed_(seed) {
        if (d_ <= 0 || base_w_ <= 0 || L_ <= 0)
            throw std::invalid_argument("Parámetros inválidos: d, base_w, L deben ser > 0.");
        build_levels_();
    }
    void insert(uint64_t kmer_code, uint64_t weight = 1ULL) {
        for (int lvl = 0; lvl < L_; ++lvl) {
            if (promotes_to_level_(kmer_code, lvl)) {
                levels_[lvl].insert(kmer_code, weight);
            }
        }
    }

    uint64_t estimate_u64(uint64_t kmer_code) const {
        bool any = false;
        uint64_t best = ULLONG_MAX;
        for (int lvl = 0; lvl < L_; ++lvl) {
            if (promotes_to_level_(kmer_code, lvl)) {
                any = true;
                uint64_t e = levels_[lvl].estimate_u64(kmer_code);
                if (e < best) best = e;
            }
        }
        if (!any) return 0ULL;//seguridad
        return best == ULLONG_MAX ? 0ULL : best;
    }
    void insert(const std::string& kmer, uint64_t weight = 1ULL) {
        insert(CountMinCU::encode_kmer_canonical(kmer), weight);
    }

    uint64_t estimate(const std::string& kmer) const {
        return estimate_u64(CountMinCU::encode_kmer_canonical(kmer));
    }

private:
    int d_, base_w_, L_;
    uint32_t seed_;
    std::vector<CountMinCU> levels_;

    void build_levels_() {
        levels_.clear();
        levels_.reserve(L_);
        int w = base_w_;
        for (int lvl = 0; lvl < L_; ++lvl) {
            int w_lvl = std::max(1, base_w_ >> lvl);
            levels_.emplace_back(d_, w_lvl, seed_ ^ (0x9e3779b9u * (uint32_t)(lvl + 101)));
        }
    }
    inline bool promotes_to_level_(uint64_t key, int level) const {
        if (level == 0) return true;
        uint32_t s = (seed_ + 0x85ebca6bu) ^ (0xc2b2ae35u * (uint32_t)(level + 11));
        uint32_t h = murmurhash(&key, s);
        uint32_t mask = (level >= 32) ? 0xFFFFFFFFu : ((1u << level) - 1u);
        return (h & mask) == 0u;
    }
};