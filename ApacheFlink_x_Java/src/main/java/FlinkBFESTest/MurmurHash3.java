// //-----------------------------------------------------------------------------
// // MurmurHash3 was written by Austin Appleby, and is placed in the public
// // domain. The author hereby disclaims copyright to this source code.

// // Note - The x86 and x64 versions do _not_ produce the same results, as the
// // algorithms are optimized for their respective platforms. You can still
// // compile and run any of them on any platform, but your performance with the
// // non-native version will be less than optimal.
package FlinkBFESTest;

public class MurmurHash3 {

    // Platform-independent constants
    private static final int BIG_CONSTANT_32 = 0xff51afd7;
    private static final long BIG_CONSTANT_64_1 = 0xff51afd7ed558ccdL;
    private static final long BIG_CONSTANT_64_2 = 0xc4ceb9fe1a85ec53L;

    // Rotate left
    private static int rotl32(int x, int r) {
        return (x << r) | (x >>> (32 - r));
    }

    private static long rotl64(long x, int r) {
        return (x << r) | (x >>> (64 - r));
    }

    // Block of 32 and 64-bit reads
    private static int getBlock32(int[] key, int i) {
        return key[i];
    }

    private static long getBlock64(long[] key, int i) {
        return key[i];
    }

    // Finalization mix
    private static int fmix32(int h) {
        h ^= h >>> 16;
        h *= 0x85ebca6b;
        h ^= h >>> 13;
        h *= 0xc2b2ae35;
        h ^= h >>> 16;
        return h;
    }

    private static long fmix64(long k) {
        k ^= k >>> 33;
        k *= BIG_CONSTANT_64_1;
        k ^= k >>> 33;
        k *= BIG_CONSTANT_64_2;
        k ^= k >>> 33;
        return k;
    }

    public static long[] murmurHash3_x64_128(byte[] key, int seed) {
        int len = key.length;
        int nblocks = len / 16;

        long h1 = seed;
        long h2 = seed;

        final long c1 = 0x87c37b91114253d5L;
        final long c2 = 0x4cf5ad432745937fL;

        long[] blocks = new long[nblocks * 2];
        for (int i = 0; i < nblocks; i++) {
            blocks[i * 2] = ((long) (key[i * 16] & 0xFF)) |
                    ((long) (key[i * 16 + 1] & 0xFF) << 8) |
                    ((long) (key[i * 16 + 2] & 0xFF) << 16) |
                    ((long) (key[i * 16 + 3] & 0xFF) << 24) |
                    ((long) (key[i * 16 + 4] & 0xFF) << 32) |
                    ((long) (key[i * 16 + 5] & 0xFF) << 40) |
                    ((long) (key[i * 16 + 6] & 0xFF) << 48) |
                    ((long) (key[i * 16 + 7] & 0xFF) << 56);
            blocks[i * 2 + 1] = ((long) (key[i * 16 + 8] & 0xFF)) |
                    ((long) (key[i * 16 + 9] & 0xFF) << 8) |
                    ((long) (key[i * 16 + 10] & 0xFF) << 16) |
                    ((long) (key[i * 16 + 11] & 0xFF) << 24) |
                    ((long) (key[i * 16 + 12] & 0xFF) << 32) |
                    ((long) (key[i * 16 + 13] & 0xFF) << 40) |
                    ((long) (key[i * 16 + 14] & 0xFF) << 48) |
                    ((long) (key[i * 16 + 15] & 0xFF) << 56);
        }

        for (int i = 0; i < nblocks; i++) {
            long k1 = getBlock64(blocks, i * 2);
            long k2 = getBlock64(blocks, i * 2 + 1);

            k1 *= c1;
            k1 = rotl64(k1, 31);
            k1 *= c2;
            h1 ^= k1;

            h1 = rotl64(h1, 27);
            h1 += h2;
            h1 = h1 * 5 + 0x52dce729;

            k2 *= c2;
            k2 = rotl64(k2, 33);
            k2 *= c1;
            h2 ^= k2;

            h2 = rotl64(h2, 31);
            h2 += h1;
            h2 = h2 * 5 + 0x38495ab5;
        }

        // Tail
        long k1 = 0;
        long k2 = 0;
        int offset = nblocks * 16;
        switch (len & 15) {
            case 15:
                k2 ^= ((long) (key[offset + 14] & 0xFF)) << 48;
            case 14:
                k2 ^= ((long) (key[offset + 13] & 0xFF)) << 40;
            case 13:
                k2 ^= ((long) (key[offset + 12] & 0xFF)) << 32;
            case 12:
                k2 ^= ((long) (key[offset + 11] & 0xFF)) << 24;
            case 11:
                k2 ^= ((long) (key[offset + 10] & 0xFF)) << 16;
            case 10:
                k2 ^= ((long) (key[offset + 9] & 0xFF)) << 8;
            case 9:
                k2 ^= ((long) (key[offset + 8] & 0xFF));
                k2 *= c2;
                k2 = rotl64(k2, 33);
                k2 *= c1;
                h2 ^= k2;
            case 8:
                k1 ^= ((long) (key[offset + 7] & 0xFF)) << 56;
            case 7:
                k1 ^= ((long) (key[offset + 6] & 0xFF)) << 48;
            case 6:
                k1 ^= ((long) (key[offset + 5] & 0xFF)) << 40;
            case 5:
                k1 ^= ((long) (key[offset + 4] & 0xFF)) << 32;
            case 4:
                k1 ^= ((long) (key[offset + 3] & 0xFF)) << 24;
            case 3:
                k1 ^= ((long) (key[offset + 2] & 0xFF)) << 16;
            case 2:
                k1 ^= ((long) (key[offset + 1] & 0xFF)) << 8;
            case 1:
                k1 ^= ((long) (key[offset] & 0xFF));
                k1 *= c1;
                k1 = rotl64(k1, 31);
                k1 *= c2;
                h1 ^= k1;
        }

        // Finalization
        h1 ^= len;
        h2 ^= len;

        h1 += h2;
        h2 += h1;

        h1 = fmix64(h1);
        h2 = fmix64(h2);

        h1 += h2;
        h2 += h1;

        return new long[] { h1, h2 };
    }
}
