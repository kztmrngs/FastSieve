// FastSieve/src/sieve.cpp
// エラトステネスの篩の高速実装

// --- メモ群 ---
/*
### 素数候補 7系列 (j=7) の各式:
|-----------------------------|-------------------------------|
| パターン 0 (diff=0): 0 + 0m | パターン 4 (diff=12): 22 + 96m |
| パターン 1 (diff=4): 7 + 32m | パターン 5 (diff=16): 29 + 128m |
| パターン 2 (diff=6): 11 + 48m | パターン 6 (diff=22): 41 + 176m |
| パターン 3 (diff=10): 18 + 80m | パターン 7 (diff=24): 44 + 192m |

### 素数候補 11系列 (j=11) の各式:
|-----------------------------|-------------------------------|
| パターン 0 (diff=0): 0 + 0m | パターン 4 (diff=12): 35 + 96m |
| パターン 1 (diff=2): 6 + 16m | パターン 5 (diff=18): 53 + 144m |
| パターン 2 (diff=6): 17 + 48m | パターン 6 (diff=20): 58 + 160m |
| パターン 3 (diff=8): 23 + 64m | パターン 7 (diff=26): 76 + 208m |

### 素数候補 13系列 (j=13) の各式:
|-----------------------------|-------------------------------|
| パターン 0 (diff=0): 0 + 0m | パターン 4 (diff=16): 55 + 128m |
| パターン 1 (diff=4): 13 + 32m | パターン 5 (diff=18): 62 + 144m |
| パターン 2 (diff=6): 20 + 48m | パターン 6 (diff=24): 83 + 192m |
| パターン 3 (diff=10): 34 + 80m | パターン 7 (diff=28): 97 + 224m |

### 素数候補 17系列 (j=17) の各式:
|-----------------------------|-------------------------------|
| パターン 0 (diff=0): 0 + 0m | パターン 4 (diff=14): 63 + 112m |
| パターン 1 (diff=2): 9 + 16m | パターン 5 (diff=20): 90 + 160m |
| パターン 2 (diff=6): 27 + 48m | パターン 6 (diff=24): 108 + 192m |
| パターン 3 (diff=12): 54 + 96m | パターン 7 (diff=26): 117 + 208m |

### 素数候補 19系列 (j=19) の各式:
|-----------------------------|-------------------------------|
| パターン 0 (diff=0): 0 + 0m | パターン 4 (diff=18): 91 + 144m |
| パターン 1 (diff=4): 20 + 32m | パターン 5 (diff=22): 111 + 176m |
| パターン 2 (diff=10): 50 + 80m | パターン 6 (diff=24): 121 + 192m |
| パターン 3 (diff=12): 61 + 96m | パターン 7 (diff=28): 142 + 224m |

### 素数候補 23系列 (j=23) の各式:
|-----------------------------|-------------------------------|
| パターン 0 (diff=0): 0 + 0m | パターン 4 (diff=18): 110 + 144m |
| パターン 1 (diff=6): 36 + 48m | パターン 5 (diff=20): 122 + 160m |
| パターン 2 (diff=8): 49 + 64m | パターン 6 (diff=24): 147 + 192m |
| パターン 3 (diff=14): 85 + 112m | パターン 7 (diff=26): 159 + 208m |

### 素数候補 29系列 (j=29) の各式:
|-----------------------------|-------------------------------|
| パターン 0 (diff=0): 0 + 0m | パターン 4 (diff=14): 108 + 112m |
| パターン 1 (diff=2): 15 + 16m | パターン 5 (diff=18): 139 + 144m |
| パターン 2 (diff=8): 62 + 64m | パターン 6 (diff=20): 154 + 160m |
| パターン 3 (diff=12): 93 + 96m | パターン 7 (diff=24): 185 + 192m |

### 素数候補 31系列 (j=31) の各式:
|-----------------------------|-------------------------------|
| パターン 0 (diff=0): 0 + 0m | パターン 4 (diff=16): 132 + 128m |
| パターン 1 (diff=6): 49 + 48m | パターン 5 (diff=18): 149 + 144m |
| パターン 2 (diff=10): 82 + 80m | パターン 6 (diff=22): 182 + 176m |
| パターン 3 (diff=12): 99 + 96m | パターン 7 (diff=28): 231 + 224m |
*/


// --- インクルード群 ---
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>

// --- 型・定数定義群 ---
using ll = long long;
using ull = unsigned long long;
constexpr size_t CACHE_BYTES = 256 * 1024;
constexpr size_t SEGMENT_SIZE = CACHE_BYTES * 30;
static const char wheel30Offsets[8] = {1, 7, 11, 13, 17, 19, 23, 29};
static const char constWheel30Offsets[8] = {7, 11, 13, 17, 19, 23, 29, 31};
static const char wheel30Indexs[30] = {
    -1, 7, -1, -1, -1, -1, -1, 0,  // 0-7 (1->7: 31系列として扱う)
    -1, -1, -1, 1, -1, 2, -1, -1,  // 8-15
    -1, 3, -1, 4, -1, -1, -1, 5,   // 16-23
    -1, -1, -1, -1, -1, 6          // 24-29
};
static const char offsetMod8[8][8] = {
    {5, 4, 0, 7, 3, 2, 6, 1}, // 7系列
    {0, 6, 1, 7, 3, 5, 2, 4}, // 11系列
    {5, 2, 1, 7, 4, 3, 0, 6}, // 13系列
    {5, 6, 0, 3, 4, 7, 1, 2}, // 17系列
    {0, 4, 2, 5, 3, 7, 1, 6}, // 19系列
    {5, 1, 6, 2, 3, 7, 0, 4}, // 23系列
    {0, 7, 6, 5, 4, 3, 2, 1}, // 29系列
    {0, 1, 2, 3, 4, 5, 6, 7}  // 31系列
};
static const char offsetIdxConst[8][8] = {
    {0, 1, 2, 2, 3, 4, 5, 6}, // 7系列
    {0, 0, 2, 2, 4, 6, 7, 9}, // 11系列
    {0, 2, 3, 4, 7, 8, 11, 12}, // 13系列
    {0, 1, 4, 7, 8, 11, 14, 15}, // 17系列
    {0, 2, 6, 7, 11, 13, 15, 17}, // 19系列
    {0, 5, 6, 11, 14, 15, 19, 20}, // 23系列
    {0, 1, 7, 11, 13, 17, 19, 23}, // 29系列
    {0, 6, 10, 12, 16, 18, 22, 28} // 31系列
};
static const char offsetIdxCoeff[8][8] = {
    {0, 4, 6, 10, 12, 16, 22, 24}, // 7系列
    {0, 2, 6, 8, 12, 18, 20, 26}, // 11系列
    {0, 4, 6, 10, 16, 18, 24, 28}, // 13系列
    {0, 2, 6, 12, 14, 20, 24, 26}, // 17系列
    {0, 4, 10, 12, 18, 22, 24, 28}, // 19系列
    {0, 6, 8, 14, 18, 20, 24, 26}, // 23系列
    {0, 2, 8, 12, 14, 18, 20, 24}, // 29系列
    {0, 6, 10, 12, 16, 18, 22, 28} // 31系列
};

// --- ヘルパー関数群 ---
inline int getBit(const std::vector<unsigned char>& sieve, ull index) {
    ull byteIndex = index >> 3; // index / 8
    ull bitIndex = index & 7; // index % 8
    return (sieve[byteIndex] >> bitIndex) & 1;
}

// --- 主関数群 ---
// ベース篩
int baseSieve(ull hi, std::vector<ull>& basePrimes){
    // ベース篩用に配列を確保 8bitで作ってビットパッキングで扱う。 ホイール法のために8/30のサイズにする。
    size_t sieveSize = (hi / 30 * 8) / 8 + 1; // 30の倍数に切り捨て
    std::vector<unsigned char> sieve(sieveSize, 0xFF); // 全ビットを1に初期化 (この配列は[1, 7, 11, 13, 17, 19, 23, 29...]というふうに対応。)
    // 1は素数ではないので0にする
    sieve[0] &= ~(1 << 0); // 数字の1に対応するビットを0にする
    // メインの篩処理
    ull mul = 0, limitBit = sieveSize << 3, limitByte = sieveSize; // limitBit = sieveSize * 8
    ull prime, idx, prime8;
    // 開始位置は7の倍数なので、iを1からスタート
    for(ull i = 1; mul <= (ull)(sqrt(hi) / 30) + 1; i += 8, mul++){
        ull mul30 = 30 * mul;
        for(int j = 0; j < 8; j++){
            if(getBit(sieve, i + j)){
                // 素数の値を復元
                ull prime = mul30 + constWheel30Offsets[j];
                ull idx = (prime * prime) / 30;
                ull offsetIdx[8];
                char masks[8];
                // オフセットとマスクを取得
                for(int k = 0; k < 8; k++){
                    offsetIdx[k] = (ull)offsetIdxConst[j][k] + (ull)offsetIdxCoeff[j][k] * mul;
                    masks[k] = ~(1 << offsetMod8[j][k]);
                }
                // メインループ
                ull maxOffset = offsetIdx[7];
                ull safeLimit = (maxOffset < limitByte) ? (limitByte - maxOffset) : 0;
                for(; idx < safeLimit; idx += prime) {
                    sieve[idx + offsetIdx[0]] &= masks[0];
                    sieve[idx + offsetIdx[1]] &= masks[1];
                    sieve[idx + offsetIdx[2]] &= masks[2];
                    sieve[idx + offsetIdx[3]] &= masks[3];
                    sieve[idx + offsetIdx[4]] &= masks[4];
                    sieve[idx + offsetIdx[5]] &= masks[5];
                    sieve[idx + offsetIdx[6]] &= masks[6];
                    sieve[idx + offsetIdx[7]] &= masks[7];
                }
                //終端処理
                for(; idx < limitByte; idx += prime){
                    for(int j = 0; j < 8; j++){
                        ull currentIdx = idx + offsetIdx[j];
                        if(currentIdx < limitByte){
                            sieve[currentIdx] &= masks[j];
                        }
                    }
                }
            }
        }
    }
    //素数リストを作成
    size_t primeCount = 0;
    if(hi >= 2) basePrimes[primeCount++] = 2;
    if(hi >= 3) basePrimes[primeCount++] = 3;
    if(hi >= 5) basePrimes[primeCount++] = 5;
    int flag = 0;
    for (ull i = 0; i < sieveSize; i++) {
        if (sieve[i] == 0) continue; // 全ビットが0ならスキップ
        for (int bit = 0; bit < 8; bit++) {
            if (sieve[i] & (1 << bit)) {
                ull prime = 30 * i + (ull)wheel30Offsets[bit];
                if (prime <= hi) {
                    basePrimes[primeCount++] = prime;
                } else {
                    flag = 1;
                    break;
                }
            }
        }
        if(flag) break;
    }
    return primeCount;
}
// 区間篩
int segmentSieve(ull lo, ull hi, std::vector<ull>& segmentPrimes, ull* basePrimes, size_t basePrimesSize){ 
    // ベース篩と基本同じだがloによるオフセットがある。また最初のビット反転はしない。入力にベース篩からの素数リストが必要。
    ull segmentSize = (hi - lo + 29) / 30 + 1; // 30の倍数に切り捨て
    std::vector<unsigned char> sieve(segmentSize, 0xFF); // 全ビットを1に初期化 (この配列は[lo + 1, lo + 7, lo + 11, lo + 13, lo + 17, lo + 19, lo + 23, lo + 29...]というふうに対応。)
    ull loByte = lo / 30;
    ull limitByte = segmentSize;
    // basePrimesの4番目(7)からスタート
    for(size_t i = 3; i < basePrimesSize; i++){
        ull prime = basePrimes[i];
        ull prime2 = prime * prime;
        if (prime2 > hi) break;
        // 素数の系列を取得
        int k = wheel30Indexs[prime % 30];
        // 開始位置の計算
        ull startNum = prime2;
        if (prime2 < lo){
            ull step = prime * 30;
            startNum = prime2 + ((lo -prime2 + step - 1) / step) * step;
        }
        ll baseIdx = (ll)(startNum / 30) - (ll)(loByte);
        ull mul = prime / 30;
        if (k == 7) mul--;
        ull primeByte = prime;
        // オフセットとマスクの計算
        ll offsetIdx[8];
        char masks[8];
        for (int j = 0; j < 8; j++){
            offsetIdx[j] = (ll)offsetIdxConst[k][j] + (ll)offsetIdxCoeff[k][j] * mul;
            masks[j] = ~(1 << offsetMod8[k][j]);
        }
        // 先頭処理
        if (prime2 < lo) {
            ll prevBase = baseIdx - (ll)primeByte;
            for (int j = 0; j < 8; j++){
                ll pIdx = prevBase + offsetIdx[j];
                if (pIdx >= 0 && (ull)pIdx < limitByte){
                    sieve[pIdx] &= masks[j];
                }
            }
        }

        // メインループ
        ull maxOffset = offsetIdx[7];
        ull safeLimit = (maxOffset < limitByte) ? (limitByte - maxOffset) : 0;
        ll idx = baseIdx;
        for(; (ull)idx < safeLimit; idx += primeByte){
            sieve[idx + offsetIdx[0]] &= masks[0];
            sieve[idx + offsetIdx[1]] &= masks[1];
            sieve[idx + offsetIdx[2]] &= masks[2];
            sieve[idx + offsetIdx[3]] &= masks[3];
            sieve[idx + offsetIdx[4]] &= masks[4];
            sieve[idx + offsetIdx[5]] &= masks[5];
            sieve[idx + offsetIdx[6]] &= masks[6];
            sieve[idx + offsetIdx[7]] &= masks[7];
        }
        //終端処理
        for(; (ull)idx < limitByte; idx += primeByte){
            for(int j = 0; j < 8; j++){
                ll pIdx = idx + offsetIdx[j];
                if((ull)pIdx < limitByte){
                    sieve[pIdx] &= masks[j];
                }
            }
        }
    }
    // 素数リストを作成
    size_t primeCount = 0; // 素数の個数
    if(lo <= 2 && hi >= 2){
        segmentPrimes[primeCount++] = 2;
    }
    if(lo <= 3 && hi >= 3){
        segmentPrimes[primeCount++] = 3;
    }
    if(lo <= 5 && hi >= 5){
        segmentPrimes[primeCount++] = 5;
    }
    int flag = 0;
    for (ull i = 0; i < segmentSize; i++) {
        if (sieve[i] == 0) continue; // 全ビットが0ならスキップ
        for (int bit = 0; bit < 8; bit++) {
            if (sieve[i] & (1 << bit)) {
                ull prime = lo + (30 * i + (ull)wheel30Offsets[bit]);
                if (prime >= lo && prime <= hi) {
                    segmentPrimes[primeCount++] = prime;
                } else if (prime > hi) {
                    flag = 1;
                    break;
                }
            }
        }
        if(flag) break;
    }
    return primeCount;
}

// 篩の制御用
std::vector<ull> sieveCompute(ull lo, ull hi, int mode){
    // 配列を確保   
    ull limit = ((ull)sqrt((double)hi) / 30) * 30 + 30; // 30の倍数に切り捨て
    size_t basePrimesSize = (size_t)(limit / log(limit) * 1.25506 + 1); // 素数計数関数に基づく素数の上限個数
    std::vector<ull> basePrimes(basePrimesSize);

    ull segmentSize = SEGMENT_SIZE;// 30の倍数に切り捨て
    size_t segmentPrimesSize = (size_t)(segmentSize / log(segmentSize) * 1.25506 + 1); // 素数計数関数に基づく素数の上限個数
    std::vector<ull> segmentPrimes(segmentPrimesSize);
    std::vector<ull> foundPrimes;

    // モードに応じて篩を実行
    if(mode == 0 || lo < 3){ // 0: 2からhiまでの素数を求める。またloが3未満の場合変わりようがない。
        if (hi > 1) {
            foundPrimes.reserve((size_t)(hi / log(hi) * 1.25506 + 1));
        }
        if(961 < hi){ // hiが961より大きい場合は区間篩を使う。961は31^2。
            ull basePrimesCount = baseSieve(limit, basePrimes);
            // 30の倍数に揃える
            foundPrimes.insert(foundPrimes.end(), basePrimes.begin(), basePrimes.begin() + basePrimesCount);
            for (ull i = limit; i <= hi; i += segmentSize) {
                ull currentSegmentSize = std::min(i + segmentSize - 1, hi);
                ull segmentPrimesCount = segmentSieve(i, currentSegmentSize, segmentPrimes, basePrimes.data(), basePrimesCount);
                foundPrimes.insert(foundPrimes.end(), segmentPrimes.begin(), segmentPrimes.begin() + segmentPrimesCount);
            }
        } else {
            // hiが961以下の場合はベース篩のみで完結する
            if (hi > 1) {
                size_t newBasePrimesSize = (size_t)(hi / log(hi) * 1.25506 + 1);
                if (newBasePrimesSize > basePrimes.size()) {
                    basePrimes.resize(newBasePrimesSize);
                }
            }
            ull basePrimesCount = baseSieve(hi, basePrimes);
            foundPrimes.insert(foundPrimes.end(), basePrimes.begin(), basePrimes.begin() + basePrimesCount);
        }
    } else if(mode == 1){ // 1: loからhiまでの素数を求める
        if (hi > 1) {
            double lo_primes = (lo > 1) ? (lo / log(lo)) : 0;
            foundPrimes.reserve((size_t)((hi / log(hi) - lo_primes) * 1.25506 + 1));
        }
        lo = std::max(lo, 2ULL); // loが2未満の場合は2にする
        ull basePrimesCount = baseSieve(limit, basePrimes);
        ull lo30 = lo / 30 * 30;
        ull hi30 = hi / 30 * 30 + 30;
        // 30の倍数に揃える
        for (ull i = lo30; i <= hi30; i += segmentSize) {
            ull current_hi = std::min(i + segmentSize - 1, hi);
            ull segmentPrimesCount = segmentSieve(i, current_hi, segmentPrimes, basePrimes.data(), basePrimesCount);
            for(size_t j = 0; j < segmentPrimesCount; j++){
                if(segmentPrimes[j] >= lo && segmentPrimes[j] <= hi){
                    foundPrimes.push_back(segmentPrimes[j]);
                }
            }
            if (current_hi == hi) break;
        }
    } else {
        printf("sieveCompute: Invalid mode %d\n", mode);
        exit(1);
    }
    return foundPrimes;
}

// 実質io関数
int main(int argc, char *argv[]){
    // 配列を確保

    // 入力を受け取る。
    /*
    <実行ファイル> <上限値>
    例: ./FastSieve 1000000
    <実行ファイル> <下限値> <上限値>
    例: ./FastSieve 1000000 2000000
    */
    if(argc > 3){
        printf("Usage: %s <upper_limit> OR %s <lower_limit> <upper_limit>\n", argv[0], argv[0]);
        return 1;
    } else if(argc == 2){
        ull hi  = std::stoull(argv[1]);
        auto start = std::chrono::high_resolution_clock::now();
        auto primes = sieveCompute(2, hi, 0);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        ull primeCount = primes.size();
        printf("Primes up to %llu (total %llu):\n", hi, primeCount);
        for (const auto& p : primes) {
            printf("%llu ", p);
        }
        printf("\n");
        printf("Time: %f seconds\n", elapsed.count());
    } else if(argc == 3){
        ull lo = std::stoull(argv[1]);
        ull hi = std::stoull(argv[2]);
        auto start = std::chrono::high_resolution_clock::now();
        auto primes = sieveCompute(lo, hi, 1);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        ull primeCount = primes.size();
        printf("Primes from %llu to %llu (total %llu):\n", lo, hi, primeCount);
        for (const auto& p : primes) {
            printf("%llu ", p);
        }
        printf("\n");
        printf("Time: %f seconds\n", elapsed.count());
    } else {
        printf("Usage: %s <upper_limit> OR %s <lower_limit> <upper_limit>\n", argv[0], argv[0]);
        return 1;
    }
    return 0;
}
