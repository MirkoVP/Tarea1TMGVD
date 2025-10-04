#ifndef MURMURHASH32_HPP
#define MURMURHASH32_HPP

#include <cstdint>
#include <cstddef>

// MurmurHash3 x86_32 — implementación estándar
// Fuente de referencia: Austin Appleby (Dominio público / MIT-like)
// Devuelve un hash de 32 bits para un buffer arbitrario (key,len) y seed.
static inline uint32_t murmurhash32_x86_32(const void* key, int len, uint32_t seed) {
    const uint8_t * data = (const uint8_t*)key;
    const int nblocks = len / 4;

    uint32_t h1 = seed;

    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;

    //----------
    // body

    const uint32_t * blocks = (const uint32_t *)(data + nblocks*4);

    for(int i = -nblocks; i; i++) {
        uint32_t k1 = blocks[i];

        k1 *= c1;
        k1 = (k1 << 15) | (k1 >> (32 - 15));
        k1 *= c2;

        h1 ^= k1;
        h1 = (h1 << 13) | (h1 >> (32 - 13));
        h1 = h1*5 + 0xe6546b64;
    }

    //----------
    // tail

    const uint8_t * tail = (const uint8_t*)(data + nblocks*4);

    uint32_t k1 = 0;

    switch(len & 3) {
      case 3: k1 ^= tail[2] << 16;
      case 2: k1 ^= tail[1] << 8;
      case 1: k1 ^= tail[0];
              k1 *= c1;
              k1 = (k1 << 15) | (k1 >> (32 - 15));
              k1 *= c2;
              h1 ^= k1;
    }

    //----------
    // finalization

    h1 ^= len;

    // fmix32
    h1 ^= h1 >> 16;
    h1 *= 0x85ebca6b;
    h1 ^= h1 >> 13;
    h1 *= 0xc2b2ae35;
    h1 ^= h1 >> 16;

    return h1;
}

// Wrapper 1: firma genérica que queremos usar desde el código
static inline uint32_t murmurhash(const void* key, int len, uint32_t seed) {
    return murmurhash32_x86_32(key, len, seed);
}

// Wrapper 2: compatibilidad con TU firma original (uint64_t*)
// Hashea los 8 bytes de la palabra de 64 bits.
static inline uint32_t murmurhash(const uint64_t* key, uint32_t seed) {
    return murmurhash32_x86_32((const void*)key, (int)sizeof(uint64_t), seed);
}

#endif // MURMURHASH32_HPP
