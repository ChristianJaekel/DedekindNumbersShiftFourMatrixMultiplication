#include "IVLCalc.hpp"
#include "FDL.hpp"

#include <algorithm>
#include <cmath>

const std::array<uint64_t, 65> bitMasks = {
    0b0000000000000000000000000000000000000000000000000000000000000000,
    0b0000000000000000000000000000000000000000000000000000000000000001,
    0b0000000000000000000000000000000000000000000000000000000000000011,
    0b0000000000000000000000000000000000000000000000000000000000000111,
    0b0000000000000000000000000000000000000000000000000000000000001111,
    0b0000000000000000000000000000000000000000000000000000000000011111,
    0b0000000000000000000000000000000000000000000000000000000000111111,
    0b0000000000000000000000000000000000000000000000000000000001111111,
    0b0000000000000000000000000000000000000000000000000000000011111111,
    0b0000000000000000000000000000000000000000000000000000000111111111,
    0b0000000000000000000000000000000000000000000000000000001111111111,
    0b0000000000000000000000000000000000000000000000000000011111111111,
    0b0000000000000000000000000000000000000000000000000000111111111111,
    0b0000000000000000000000000000000000000000000000000001111111111111,
    0b0000000000000000000000000000000000000000000000000011111111111111,
    0b0000000000000000000000000000000000000000000000000111111111111111,
    0b0000000000000000000000000000000000000000000000001111111111111111,
    0b0000000000000000000000000000000000000000000000011111111111111111,
    0b0000000000000000000000000000000000000000000000111111111111111111,
    0b0000000000000000000000000000000000000000000001111111111111111111,
    0b0000000000000000000000000000000000000000000011111111111111111111,
    0b0000000000000000000000000000000000000000000111111111111111111111,
    0b0000000000000000000000000000000000000000001111111111111111111111,
    0b0000000000000000000000000000000000000000011111111111111111111111,
    0b0000000000000000000000000000000000000000111111111111111111111111,
    0b0000000000000000000000000000000000000001111111111111111111111111,
    0b0000000000000000000000000000000000000011111111111111111111111111,
    0b0000000000000000000000000000000000000111111111111111111111111111,
    0b0000000000000000000000000000000000001111111111111111111111111111,
    0b0000000000000000000000000000000000011111111111111111111111111111,
    0b0000000000000000000000000000000000111111111111111111111111111111,
    0b0000000000000000000000000000000001111111111111111111111111111111,
    0b0000000000000000000000000000000011111111111111111111111111111111,
    0b0000000000000000000000000000000111111111111111111111111111111111,
    0b0000000000000000000000000000001111111111111111111111111111111111,
    0b0000000000000000000000000000011111111111111111111111111111111111,
    0b0000000000000000000000000000111111111111111111111111111111111111,
    0b0000000000000000000000000001111111111111111111111111111111111111,
    0b0000000000000000000000000011111111111111111111111111111111111111,
    0b0000000000000000000000000111111111111111111111111111111111111111,
    0b0000000000000000000000001111111111111111111111111111111111111111,
    0b0000000000000000000000011111111111111111111111111111111111111111,
    0b0000000000000000000000111111111111111111111111111111111111111111,
    0b0000000000000000000001111111111111111111111111111111111111111111,
    0b0000000000000000000011111111111111111111111111111111111111111111,
    0b0000000000000000000111111111111111111111111111111111111111111111,
    0b0000000000000000001111111111111111111111111111111111111111111111,
    0b0000000000000000011111111111111111111111111111111111111111111111,
    0b0000000000000000111111111111111111111111111111111111111111111111,
    0b0000000000000001111111111111111111111111111111111111111111111111,
    0b0000000000000011111111111111111111111111111111111111111111111111,
    0b0000000000000111111111111111111111111111111111111111111111111111,
    0b0000000000001111111111111111111111111111111111111111111111111111,
    0b0000000000011111111111111111111111111111111111111111111111111111,
    0b0000000000111111111111111111111111111111111111111111111111111111,
    0b0000000001111111111111111111111111111111111111111111111111111111,
    0b0000000011111111111111111111111111111111111111111111111111111111,
    0b0000000111111111111111111111111111111111111111111111111111111111,
    0b0000001111111111111111111111111111111111111111111111111111111111,
    0b0000011111111111111111111111111111111111111111111111111111111111,
    0b0000111111111111111111111111111111111111111111111111111111111111,
    0b0001111111111111111111111111111111111111111111111111111111111111,
    0b0011111111111111111111111111111111111111111111111111111111111111,
    0b0111111111111111111111111111111111111111111111111111111111111111,
    0b1111111111111111111111111111111111111111111111111111111111111111};

IVLCalc::IVLCalc(unsigned n, unsigned difference) {
    mNrOfBits = static_cast<uint64_t>(std::pow(2, n));

    mDifference    = difference;
    mGeneratorsRed = n - difference;
    FreeDistLat fdl(n - difference);
    mFdlRed = fdl.getElements();

    initIndices();
    initIntervals();
}

void IVLCalc::initIndices() {
    mIndex.resize(mFdlRed.back() + 1);
    unsigned i = 0;
    for (const auto element : mFdlRed) {
        mIndex[element] = i;
        ++i;
    }
}

void IVLCalc::initIntervals() {
    mIntRed.resize(mFdlRed.size() * mFdlRed.size());
    for (auto a : mFdlRed)
        for (auto b : mFdlRed) {
            if (b == (a | b)) {
                uint64_t intervalSize = 0;
                for (unsigned int i = mIndex[a]; i <= mIndex[b]; ++i) {
                    if ((mFdlRed[i] == (mFdlRed[i] | a)) && (b == (mFdlRed[i] | b)))
                        intervalSize++;
                }

                mIntRed[mIndex[a] * mFdlRed.size() + mIndex[b]] = intervalSize;
            }
        }
}

std::vector<uint64_t> IVLCalc::generateIntervalMixed(uint64_t x, unsigned int ind) const {
    uint64_t y = mFdlRed[ind];

    std::vector<uint64_t> interval;

    for (unsigned i = mIndex[x]; i <= ind; ++i)
        if ((x == (x & mFdlRed[i])) && (y == (y | mFdlRed[i])))
            interval.push_back(mFdlRed[i]);

    return interval;
}

uint64_t IVLCalc::intervalLength1(uint64_t a, uint64_t b, uint64_t c, uint64_t d) const {
    uint64_t indD = mIndex[d];

    uint64_t ivl = 0;

    auto intervalAC = generateInterval(a, c);

    for (auto x : intervalAC)
        ivl += mIntRed[mIndex[b | x] * mFdlRed.size() + indD];

    return ivl;
}

uint64_t IVLCalc::intervalLength2(uint64_t a, uint64_t b, uint64_t c, uint64_t d, uint64_t e,
                                  uint64_t f, uint64_t g, uint64_t h) const {
    uint64_t indA = mIndex[a];
    uint64_t indH = mIndex[h];

    uint64_t ivl = 0;

    if ((b == c) && (f == g)) {
        auto intervalBF = generateInterval(b, f);

        for (unsigned alpha = 0; alpha < intervalBF.size(); ++alpha) {
            uint64_t x      = intervalBF[alpha];
            uint64_t xJoinE = x & e;
            uint64_t xMeetD = x | d;

            ivl += mIntRed[mIndex[xMeetD] * mFdlRed.size() + indH] *
                   mIntRed[indA * mFdlRed.size() + mIndex[xJoinE]];

            for (unsigned beta = alpha + 1; beta < intervalBF.size(); ++beta) {
                uint64_t y = intervalBF[beta];
                ivl += 2 * mIntRed[mIndex[xMeetD | y] * mFdlRed.size() + indH] *
                       mIntRed[indA * mFdlRed.size() + mIndex[xJoinE & y]];
            }
        }
    } else {
        auto intervalBF = generateInterval(b, f);
        auto intervalCG = generateInterval(c, g);

        for (auto x : intervalBF) {
            uint64_t xJoinE = x & e;
            uint64_t xMeetD = x | d;

            for (auto y : intervalCG) {
                ivl += mIntRed[mIndex[xMeetD | y] * mFdlRed.size() + indH] *
                       mIntRed[indA * mFdlRed.size() + mIndex[xJoinE & y]];
            }
        }
    }

    return ivl;
}

uint64_t IVLCalc::intervalLength3(uint64_t a, uint64_t b, uint64_t c, uint64_t d, uint64_t e,
                                  uint64_t f, uint64_t g, uint64_t h, uint64_t i, uint64_t j,
                                  uint64_t k, uint64_t l, uint64_t m, uint64_t n, uint64_t o,
                                  uint64_t p) const {
    uint64_t indA = mIndex[a] * mFdlRed.size();
    uint64_t indB = mIndex[b] * mFdlRed.size();
    uint64_t indC = mIndex[c] * mFdlRed.size();
    uint64_t indE = mIndex[e] * mFdlRed.size();
    uint64_t indL = mIndex[l];
    uint64_t indN = mIndex[n];
    uint64_t indO = mIndex[o];
    uint64_t indP = mIndex[p];

    uint64_t ivl = 0;

    // to detect symmetries
    bool FNeqGO = (f == g && indN == indO);
    bool FNeqDL = (f == d && indN == indL);
    bool GOeqDL = (g == d && indO == indL);
    bool JMKeq  = (j == m && j == k);
    bool BCEeq  = (b == c && b == e);

    auto intervalAI = generateInterval(a, i);

    for (auto s : intervalAI) {
        uint64_t indexBS = mIndex[b | s] * mFdlRed.size();
        uint64_t indexCS = mIndex[c | s] * mFdlRed.size();
        uint64_t indexES = mIndex[e | s] * mFdlRed.size();

        if (FNeqGO && FNeqDL && BCEeq && JMKeq) // all 3 intervals are equal
        {
            auto intervalFN = generateIntervalMixed(f | s, indN);

            for (unsigned alpha = 0; alpha < intervalFN.size(); ++alpha) {
                uint64_t x      = intervalFN[alpha];
                uint64_t xJoinJ = x & j;
                uint64_t xJoinK = x & k;
                uint64_t xJoinM = x & m;
                uint64_t xMeetH = x | h;

                uint64_t ivlTmp =
                    mIntRed[mIndex[xMeetH] * mFdlRed.size() + indP] // alpha = beta = gamma
                    * mIntRed[indexBS + mIndex[xJoinJ]] * mIntRed[indexCS + mIndex[xJoinK]];

                for (unsigned gamma = alpha + 1; gamma < intervalFN.size();
                     ++gamma) // alpha = beta != gamma
                {
                    uint64_t v = intervalFN[gamma];

                    ivlTmp += 3 * mIntRed[mIndex[xMeetH | v] * mFdlRed.size() + indP] *
                              mIntRed[indexBS + mIndex[xJoinJ & v]] *
                              mIntRed[indexCS + mIndex[xJoinK & v]];
                }

                ivl += ivlTmp * mIntRed[indexES + mIndex[xJoinM]];

                for (unsigned beta = alpha + 1; beta < intervalFN.size(); ++beta) {
                    uint64_t y           = intervalFN[beta];
                    uint64_t yJoinK      = y & k;
                    uint64_t xMeetYMeetH = xMeetH | y;

                    // alpha !=beta = gamma
                    ivl += 3 * mIntRed[mIndex[xMeetYMeetH] * mFdlRed.size() + indP] *
                           mIntRed[indexBS + mIndex[xJoinJ & y]] *
                           mIntRed[indexCS + mIndex[yJoinK]] *
                           mIntRed[indexES + mIndex[xJoinM & y]];

                    uint64_t ivlTmp = 0;
                    for (unsigned gamma = beta + 1; gamma < intervalFN.size();
                         ++gamma) // alpha !=beta != gamma
                    {
                        uint64_t v = intervalFN[gamma];

                        ivlTmp += 6 * mIntRed[mIndex[xMeetYMeetH | v] * mFdlRed.size() + indP] *
                                  mIntRed[indexBS + mIndex[xJoinJ & v]] *
                                  mIntRed[indexCS + mIndex[yJoinK & v]];
                    }

                    ivl += ivlTmp * mIntRed[indexES + mIndex[xJoinM & y]];
                }
            }
        } else if (FNeqGO && !FNeqDL && BCEeq && JMKeq) // only first two are equal
        {
            auto intervalFN = generateIntervalMixed(f | s, indN);
            auto intervalDL = generateIntervalMixed(d | s, indL);

            for (unsigned alpha = 0; alpha < intervalFN.size(); ++alpha) // x=y
            {
                uint64_t x      = intervalFN[alpha];
                uint64_t xJoinJ = x & j;
                uint64_t xJoinM = x & m;

                uint64_t xJoinK      = x & k; // alpha = beta
                uint64_t xMeetXMeetH = x | h;

                uint64_t ivlTmp = 0;

                for (auto v : intervalDL) {
                    ivlTmp += mIntRed[mIndex[xMeetXMeetH | v] * mFdlRed.size() + indP] *
                              mIntRed[indexBS + mIndex[xJoinJ & v]] *
                              mIntRed[indexCS + mIndex[xJoinK & v]];
                }

                ivl += ivlTmp * mIntRed[indexES + mIndex[xJoinM]];

                for (unsigned beta = alpha + 1; beta < intervalFN.size(); ++beta) // alpha != beta
                {
                    uint64_t y           = intervalFN[beta];
                    uint64_t yJoinK      = y & k;
                    uint64_t xMeetYMeetH = x | y | h;

                    uint64_t ivlTmp = 0;

                    uint64_t factor = 0;

                    for (auto v : intervalDL) {
                        ivlTmp += mIntRed[mIndex[xMeetYMeetH | v] * mFdlRed.size() + indP] *
                                  mIntRed[indexBS + mIndex[xJoinJ & v]] *
                                  mIntRed[indexCS + mIndex[yJoinK & v]];
                    }

                    ivl += 2 * ivlTmp * mIntRed[indexES + mIndex[xJoinM & y]];
                }
            }
        } else if (!FNeqGO && GOeqDL && BCEeq && JMKeq) // only last two are equal
        {
            auto intervalFN = generateIntervalMixed(f | s, indN);
            auto intervalGO = generateIntervalMixed(g | s, indO);

            for (auto x : intervalFN) {
                uint64_t xJoinJ = x & j;
                uint64_t xJoinM = x & m;

                for (unsigned int alpha = 0; alpha < intervalGO.size(); ++alpha) {
                    uint64_t y           = intervalGO[alpha];
                    uint64_t yJoinK      = y & k;
                    uint64_t xMeetYMeetH = x | y | h;

                    // case alpha = beta
                    ivl += mIntRed[mIndex[xMeetYMeetH | y] * mFdlRed.size() + indP] *
                           mIntRed[indexBS + mIndex[xJoinJ & y]] *
                           mIntRed[indexCS + mIndex[yJoinK]] *
                           mIntRed[indexES + mIndex[xJoinM & y]];

                    uint64_t ivlTmp = 0;

                    for (unsigned int beta = alpha + 1; beta < intervalGO.size(); ++beta) {
                        uint64_t v = intervalGO[beta];
                        ivlTmp += mIntRed[mIndex[xMeetYMeetH | v] * mFdlRed.size() + indP] *
                                  mIntRed[indexBS + mIndex[xJoinJ & v]] *
                                  mIntRed[indexCS + mIndex[yJoinK & v]];
                    }

                    ivl += 2 * ivlTmp * mIntRed[indexES + mIndex[xJoinM & y]];
                }
            }

        } else if (!FNeqGO && FNeqDL && BCEeq && JMKeq) // only first and last are equal
        {
            auto intervalFN = generateIntervalMixed(f | s, indN);
            auto intervalGO = generateIntervalMixed(g | s, indO);

            for (unsigned alpha = 0; alpha < intervalFN.size(); ++alpha) {
                uint64_t x                = intervalFN[alpha];
                uint64_t xJoinJ           = x & j;
                uint64_t xJoinM           = x & m;
                uint64_t intervalBSXJoinJ = mIntRed[indexBS + mIndex[xJoinJ]];

                for (auto y : intervalGO) {
                    uint64_t yJoinK      = y & k;
                    uint64_t xMeetYMeetH = x | y | h;

                    uint64_t ivlTmp =
                        mIntRed[mIndex[xMeetYMeetH] * mFdlRed.size() + indP] // alpha = beta
                        * intervalBSXJoinJ * mIntRed[indexCS + mIndex[yJoinK & x]];

                    uint64_t ivlTmp2 = 0;
                    for (unsigned beta = alpha + 1; beta < intervalFN.size();
                         ++beta) // alpha !=beta
                    {
                        uint64_t v = intervalFN[beta];

                        ivlTmp2 += mIntRed[mIndex[xMeetYMeetH | v] * mFdlRed.size() + indP] *
                                   mIntRed[indexBS + mIndex[xJoinJ & v]] *
                                   mIntRed[indexCS + mIndex[yJoinK & v]];
                    }

                    ivl += (2 * ivlTmp2 + ivlTmp) * mIntRed[indexES + mIndex[xJoinM & y]];
                }
            }

        } else {
            auto intervalFN = generateIntervalMixed(f | s, indN);
            auto intervalGO = generateIntervalMixed(g | s, indO);
            auto intervalDL = generateIntervalMixed(d | s, indL);

            for (auto x : intervalFN) {
                uint64_t xJoinJ = x & j;
                uint64_t xJoinM = x & m;

                for (auto y : intervalGO) {
                    uint64_t yJoinK      = y & k;
                    uint64_t xMeetYMeetH = x | y | h;

                    uint64_t ivlTmp = 0;

                    for (auto v : intervalDL) {
                        ivlTmp += mIntRed[mIndex[xMeetYMeetH | v] * mFdlRed.size() + indP] *
                                  mIntRed[indexBS + mIndex[xJoinJ & v]] *
                                  mIntRed[indexCS + mIndex[yJoinK & v]];
                    }

                    ivl += ivlTmp * mIntRed[indexES + mIndex[xJoinM & y]];
                }
            }
        }
    }

    return ivl;
}

// internally a dual order is used!
uint64_t IVLCalc::intervalLength(const std::vector<uint64_t>& v) const {
    if (std::adjacent_find(v.begin(), v.end(), std::not_equal_to<>()) == v.end())
        return 1;

    switch (mDifference) {
        case 3: {
            switch (v.size()) {
                case 4: {
                    if (v[2] != (v[0] | v[2]) || v[3] != (v[1] | v[3])) // assure that y<=x
                        return 0;

                    auto     split4X = splitIn4(v[0]);
                    uint64_t a       = split4X[0];
                    uint64_t b       = split4X[1];
                    uint64_t c       = split4X[2];
                    uint64_t d       = split4X[3];

                    auto     split4Y = splitIn4(v[1]);
                    uint64_t e       = split4Y[0];
                    uint64_t f       = split4Y[1];
                    uint64_t g       = split4Y[2];
                    uint64_t h       = split4Y[3];

                    auto     split4Z = splitIn4(v[2]);
                    uint64_t i       = split4Z[0];
                    uint64_t j       = split4Z[1];
                    uint64_t k       = split4Z[2];
                    uint64_t l       = split4Z[3];

                    auto     split4W = splitIn4(v[3]);
                    uint64_t m       = split4W[0];
                    uint64_t n       = split4W[1];
                    uint64_t o       = split4W[2];
                    uint64_t p       = split4W[3];

                    return intervalLength3(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p);

                    break;
                }
                default: // 2
                {
                    if (v[1] != (v[0] | v[1])) // assure that y<=x
                        return 0;

                    auto     splitX = splitIn8(v[0]);
                    uint64_t a      = splitX[0];
                    uint64_t b      = splitX[1];
                    uint64_t c      = splitX[2];
                    uint64_t d      = splitX[3];
                    uint64_t e      = splitX[4];
                    uint64_t f      = splitX[5];
                    uint64_t g      = splitX[6];
                    uint64_t h      = splitX[7];

                    auto     splitY = splitIn8(v[1]);
                    uint64_t i      = splitY[0];
                    uint64_t j      = splitY[1];
                    uint64_t k      = splitY[2];
                    uint64_t l      = splitY[3];
                    uint64_t m      = splitY[4];
                    uint64_t n      = splitY[5];
                    uint64_t o      = splitY[6];
                    uint64_t p      = splitY[7];

                    return intervalLength3(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p);
                    break;
                }
            }
        }
        case 2: {
            switch (v.size()) {
                case 4: {
                    if (v[2] != (v[0] | v[2]) || v[3] != (v[1] | v[3])) // assure that y<=x
                        return 0;

                    auto     splitW = splitIn2(v[0]);
                    uint64_t a      = splitW[0];
                    uint64_t b      = splitW[1];

                    auto     splitX = splitIn2(v[1]);
                    uint64_t c      = splitX[0];
                    uint64_t d      = splitX[1];

                    auto     splitY = splitIn2(v[2]);
                    uint64_t e      = splitY[0];
                    uint64_t f      = splitY[1];

                    auto     splitZ = splitIn2(v[3]);
                    uint64_t g      = splitZ[0];
                    uint64_t h      = splitZ[1];

                    return intervalLength2(a, b, c, d, e, f, g, h);

                    break;
                }
                default: // 2
                {
                    if (v[1] != (v[0] | v[1])) // assure that y<=x
                        return 0;

                    auto     split4X = splitIn4(v[0]);
                    uint64_t a       = split4X[0];
                    uint64_t b       = split4X[1];
                    uint64_t c       = split4X[2];
                    uint64_t d       = split4X[3];

                    auto     split4Y = splitIn4(v[1]);
                    uint64_t e       = split4Y[0];
                    uint64_t f       = split4Y[1];
                    uint64_t g       = split4Y[2];
                    uint64_t h       = split4Y[3];

                    return intervalLength2(a, b, c, d, e, f, g, h);

                    break;
                }
            }
        }
        default: // 1
        {
            switch (v.size()) {
                case 4: {
                    if (v[2] != (v[0] | v[2]) || v[3] != (v[1] | v[3])) // assure that y<=x
                        return 0;

                    return intervalLength1(v[0], v[1], v[2], v[3]);

                    break;
                }
                default: // 2
                {
                    if (v[1] != (v[0] | v[1])) // assure that y<=x
                        return 0;

                    auto     splitX = splitIn2(v[0]);
                    uint64_t a      = splitX[0];
                    uint64_t b      = splitX[1];

                    auto     splitY = splitIn2(v[1]);
                    uint64_t c      = splitY[0];
                    uint64_t d      = splitY[1];

                    return intervalLength1(a, b, c, d);

                    break;
                }
            }
        }
    }
}

std::array<uint64_t, 2> IVLCalc::splitIn2(uint64_t x) const {
    uint64_t a = x >> mNrOfBits / 2;
    uint64_t b = x & bitMasks[mNrOfBits / 2];

    return {a, b};
}

std::array<uint64_t, 4> IVLCalc::splitIn4(uint64_t x) const {
    auto     split2X = splitIn2(x);
    uint64_t a       = split2X[0] >> mNrOfBits / 4;
    uint64_t b       = split2X[0] & bitMasks[mNrOfBits / 4];

    uint64_t c = split2X[1] >> mNrOfBits / 4;
    uint64_t d = split2X[1] & bitMasks[mNrOfBits / 4];

    return {a, b, c, d};
}

std::array<uint64_t, 8> IVLCalc::splitIn8(uint64_t x) const {
    auto split4X = splitIn4(x);

    uint64_t a = split4X[0] >> mNrOfBits / 8;
    uint64_t b = split4X[0] & bitMasks[mNrOfBits / 8];

    uint64_t c = split4X[1] >> mNrOfBits / 8;
    uint64_t d = split4X[1] & bitMasks[mNrOfBits / 8];

    uint64_t e = split4X[2] >> mNrOfBits / 8;
    uint64_t f = split4X[2] & bitMasks[mNrOfBits / 8];

    uint64_t g = split4X[3] >> mNrOfBits / 8;
    uint64_t h = split4X[3] & bitMasks[mNrOfBits / 8];

    return {a, b, c, d, e, f, g, h};
}

std::vector<uint64_t> IVLCalc::generateInterval(uint64_t x, uint64_t y) const {
    std::vector<uint64_t> interval;

    for (unsigned i = mIndex[x]; i <= mIndex[y]; ++i)
        if ((x == (x & mFdlRed[i])) && (y == (y | mFdlRed[i])))
            interval.push_back(mFdlRed[i]);

    return interval;
}
