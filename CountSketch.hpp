#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include "murmurhash32.hpp"

class CountSketch {
public:
    CountSketch(int d, int w, uint32_t seed = 0x12345678u)
    : d_(d), w_(w), seed_(seed), table_(d, std::vector<int64_t>(w, 0)) {
        if (d_ <= 0 || w_ <= 0) throw std::invalid_argument("d y w deben ser > 0");
    }

    void insert(uint64_t kmer_code, int64_t weight = 1) {
        for (int i = 0; i < d_; ++i) {
            uint32_t h = bucket_hash(i, kmer_code) % static_cast<uint32_t>(w_);
            int s = sign_hash(i, kmer_code);
            table_[i][h] += s * weight;
        }
    }

    double estimate(uint64_t kmer_code) const {
        std::vector<int64_t> est(d_);
        for (int i = 0; i < d_; ++i) {
            uint32_t h = bucket_hash(i, kmer_code) % static_cast<uint32_t>(w_);
            int s = sign_hash(i, kmer_code);
            est[i] = s * table_[i][h];
        }
        return median(est);
    }
    void insert(const std::string& kmer, int64_t weight = 1) {
        insert(encode_kmer_canonical(kmer), weight);
    }

    double estimate(const std::string& kmer) const {
        return estimate(encode_kmer_canonical(kmer));
    }
    double size_MB() const {
        size_t bytes = 0;
        for (const auto& row : table_) bytes += row.capacity() * sizeof(int64_t);
        return static_cast<double>(bytes) / (1024.0 * 1024.0);
    }
    static uint64_t encode_kmer_2bits(const std::string& kmer) {
        if (kmer.size() > 32)
            throw std::invalid_argument("k > 32 no cabe en 64 bits (2 bits por símbolo).");

        uint64_t code = 0ULL;
        for (char ch : kmer) {
            code <<= 2;
            switch (ch) {
                case 'A': case 'a': code |= 0b00; break;
                case 'C': case 'c': code |= 0b01; break;
                case 'G': case 'g': code |= 0b10; break;
                case 'T': case 't': code |= 0b11; break;
                default:
                    throw std::invalid_argument("Símbolo inválido en k-mer (usa A/C/G/T).");
            }
        };
        return code;
    }

    static uint64_t reverse_complement_code(uint64_t code, size_t k) {
        uint64_t rc = 0ULL;
        for (size_t i = 0; i < k; i++) {
            uint64_t bits = code & 0b11ULL;
            code >>= 2;
            bits ^= 0b11ULL;
            rc <<= 2;
            rc |= bits;
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
    std::vector<std::vector<int64_t>> table_;
    inline uint32_t bucket_hash(int i, uint64_t key) const {
        uint32_t s = seed_ ^ (0x9e3779b9u * (uint32_t)(i + 1));
        return murmurhash(&key, s); 
    }
    inline int sign_hash(int i, uint64_t key) const {
        uint32_t s = (seed_ + 0x85ebca6bu) ^ (0xc2b2ae35u * (uint32_t)(i + 11));
        uint32_t h = murmurhash(&key, s);
        return (h & 1u) ? +1 : -1;
    }

    static double median(std::vector<int64_t> v) {
        const size_t n = v.size();
        if (n == 0) return 0.0;
        const size_t mid = n / 2;
        std::nth_element(v.begin(), v.begin() + mid, v.end());
        if (n & 1) return static_cast<double>(v[mid]);
        auto a = v[mid];
        std::nth_element(v.begin(), v.begin() + (mid - 1), v.begin() + mid);
        auto b = v[mid - 1];
        return 0.5 * (static_cast<double>(a) + static_cast<double>(b));
    }
};
