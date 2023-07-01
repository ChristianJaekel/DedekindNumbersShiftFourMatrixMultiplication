#include "Device.h"
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/transform_reduce.h>

static void checkCublasErrors(cublasStatus_t err) {
    switch (err) {
        case CUBLAS_STATUS_SUCCESS: {
            break;
        }
        case CUBLAS_STATUS_NOT_INITIALIZED: {
            std::cout << "CUBLAS_STATUS_NOT_INITIALIZED" << std::endl;
            break;
        }
        case CUBLAS_STATUS_ALLOC_FAILED: {
            std::cout << "CUBLAS_STATUS_ALLOC_FAILED" << std::endl;
            break;
        }
        case CUBLAS_STATUS_INVALID_VALUE: {
            std::cout << "CUBLAS_STATUS_INVALID_VALUE" << std::endl;
            break;
        }
        case CUBLAS_STATUS_ARCH_MISMATCH: {
            std::cout << "CUBLAS_STATUS_ARCH_MISMATCH" << std::endl;
            break;
        }
        case CUBLAS_STATUS_MAPPING_ERROR: {
            std::cout << "CUBLAS_STATUS_MAPPING_ERROR" << std::endl;
            break;
        }
        case CUBLAS_STATUS_EXECUTION_FAILED: {
            std::cout << "CUBLAS_STATUS_EXECUTION_FAILED" << std::endl;
            break;
        }
        case CUBLAS_STATUS_INTERNAL_ERROR:
            std::cout << "CUBLAS_STATUS_INTERNAL_ERROR" << std::endl;
    }
}

DeviceEnumerator::DeviceEnumerator(const std::vector<uint8_t>& idx) {
    // only one device is used
    mNrOfDevices = 1;

    // std::cout << "Device Memory: " << getAvailableGPUMemory() << std::endl;
    // std::cout << std::endl;

    mAs.resize(mNrOfDevices);
    mBs.resize(mNrOfDevices);
    mCs.resize(mNrOfDevices);

    mIdx.resize(mNrOfDevices);
    mIdxCalc.resize(mNrOfDevices);

    for (unsigned i = 0; i < mNrOfDevices; ++i) {
        cudaError_t err{cudaSetDevice(i)};
        if (err != cudaSuccess) {
            std::cout << cudaGetErrorString(err) << std::endl;
        }

        mIdx[i]     = toDevice(idx, i);
        mIdxCalc[i] = IndexCalc(mIdx[i]);
    }

    mDistTop.resize(mNrOfDevices);
    mDistBottom.resize(mNrOfDevices);
    mElements.resize(mNrOfDevices);
    mResults.resize(mNrOfDevices);

    mHandles = std::vector<cublasHandle_t>(mNrOfDevices);
    for (unsigned i = 0; i < mNrOfDevices; ++i) {
        cublasCreate(&mHandles[i]);
    }
}

void DeviceEnumerator::setTop(const std::vector<uint64_t>& distTopValues) {
    for (unsigned i = 0; i < mNrOfDevices; ++i) {
        mDistTop[i] = toDevice(distTopValues, i);
    }
}

void DeviceEnumerator::setBottom(const std::vector<uint64_t>& distBottomValues) {
    for (unsigned i = 0; i < mNrOfDevices; ++i) {
        mDistBottom[i] = toDevice(distBottomValues, i);
    }
}

void DeviceEnumerator::setElements(const std::vector<uint16_t>& elements) {
    for (unsigned i = 0; i < mNrOfDevices; ++i) {
        mElements[i] = toDevice(elements, i);
    }
}

void DeviceEnumerator::freeElements() {
    for (unsigned i = 0; i < mNrOfDevices; ++i) {
        cudaFree(mElements[i]);
    }
}

void DeviceEnumerator::freeTop() {
    for (unsigned i = 0; i < mNrOfDevices; ++i) {
        cudaFree(mDistTop[i]);
    }
}

void DeviceEnumerator::freeBottom() {
    for (unsigned i = 0; i < mNrOfDevices; ++i) {
        cudaFree(mDistBottom[i]);
    }
}

DeviceEnumerator::~DeviceEnumerator() {
    for (unsigned i = 0; i < mNrOfDevices; ++i) {
        cudaFree(mIdx[i]);
        cublasDestroy(mHandles[i]);
    }
}

size_t DeviceEnumerator::getAvailableGPUMemory() {
    size_t mf, ma;
    cudaMemGetInfo(&mf, &ma);
    return ma;
}

double* DeviceEnumerator::initOnDevice(size_t s, unsigned int deviceId) {
    double* d;
    cudaSetDevice(deviceId);
    cudaMalloc((void**)&d, s * sizeof(double));

    return d;
}

unsigned int DeviceEnumerator::getNumberOfDevices() {
    int i;
    cudaGetDeviceCount(&i);
    return static_cast<unsigned int>(i);
}

void DeviceEnumerator::initMatrices(size_t nrOfElements, size_t matrixBatchSize) {
    const size_t s =
        (((nrOfElements + 1) * nrOfElements) / 2 + numberOfThreads - 1) / numberOfThreads;
    for (unsigned i = 0; i < mNrOfDevices; ++i) {
        cudaSetDevice(i);
        mAs[i] = initOnDevice(matrixBatchSize, i);
        mBs[i] = initOnDevice(matrixBatchSize, i);
        mCs[i] = initOnDevice(matrixBatchSize, i);

        uint128T* d;
        cudaMalloc((void**)&d, s * sizeof(uint128T));
        mResults[i] = thrust::device_pointer_cast(d);
    }

    cudaError_t err{cudaGetLastError()};
    if (err != cudaSuccess) {
        std::cout << cudaGetErrorString(err) << std::endl;
    }
}

void DeviceEnumerator::freeMatrices() {
    for (unsigned i = 0; i < mNrOfDevices; ++i) {
        cudaSetDevice(i);
        cudaFree(mAs[i]);
        cudaFree(mBs[i]);
        cudaFree(mCs[i]);
        cudaFree(mResults[i].get());
    }
}

void DeviceEnumerator::doMatMulStridedBatched(size_t s, int batchCount, double* C, double* A,
                                              double* B, cublasHandle_t& handle) {
    const double alpha = 1.0;
    const double beta  = 0.0;

    cublasStatus_t err =
        cublasDgemmStridedBatched(handle, CUBLAS_OP_N, CUBLAS_OP_N, s, s, s, &alpha, A, s, s * s, B,
                                  s, s * s, &beta, C, s, s * s, batchCount);

    checkCublasErrors(err);
}

// - computes the trace of C^2
// - lineare indicees are mapped to the upper triangular matrix
__global__ void kernelTrace(size_t shift, size_t s, double* C, uint128T* result) {
    __shared__ uint128T tmpVals[numberOfThreads];

    if (threadIdx.x == 0) {
        for (unsigned int i = 0; i < numberOfThreads; ++i)
            tmpVals[i] = 0;
    }
    __syncthreads();

    const uint64_t index = static_cast<uint64_t>(blockIdx.x) * static_cast<uint64_t>(blockDim.x) +
                           static_cast<uint64_t>(threadIdx.x);
    uint64_t indAC = static_cast<uint64_t>((-1 + sqrt(static_cast<double>(8 * index + 1))) / 2);
    uint64_t indBD = index - indAC * (indAC + 1) / 2;

    if (indAC < s && indBD < s) {
        const uint128T symWeight = (static_cast<uint64_t>(indAC != indBD) + 1);

        const uint128T sum1 = C[shift + indAC * s + indBD];
        const uint128T sum2 = C[shift + indBD * s + indAC];

        tmpVals[threadIdx.x] = sum1 * sum2 * symWeight;
    }

    __syncthreads();

    if (threadIdx.x == 0) {
        uint128T sum = 0;
        for (unsigned i = 0; i < numberOfThreads; ++i)
            sum += tmpVals[i];

        result[blockIdx.x] = sum;
    }
}

static uint128T computeTrace(double* C, size_t start, size_t end,
                             const std::vector<std::vector<uint64_t>>& abValues, size_t s,
                             unsigned int deviceId, uint128T* results) {
    uint128T sum{};

    uint64_t numOfBlocks = (((s + 1) * s) / 2 + numberOfThreads - 1) / numberOfThreads;

    dim3 dimGrid(numOfBlocks);
    dim3 dimBlock(numberOfThreads);

    cudaSetDevice(deviceId);

    for (size_t batch = 0; batch < (end - start); ++batch) {
        kernelTrace<<<dimGrid, dimBlock>>>(batch * s * s, s, C, results);

        const auto trace = thrust::reduce(thrust::device, results, results + numOfBlocks,
                                          static_cast<uint128T>(0), thrust::plus<uint128T>());

        sum += trace * static_cast<uint128T>(abValues[start + batch][2]);

        cudaError_t err{cudaGetLastError()};
        if (err != cudaSuccess) {
            std::cout << cudaGetErrorString(err) << std::endl;
        }
    }

    return sum;
}

// - fills matrices A and B with values
// - lineare indicees are mapped to the upper triangular matrix
// - variable names don't match the ones in the paper
__global__ void fillMatrices(uint16_t a, uint16_t b, uint16_t* elements, uint64_t* distTop,
                             uint64_t* distBottom, IndexCalc idxCalc, size_t s, double* A,
                             double* B, size_t shift) {
    const uint64_t index = static_cast<uint64_t>(blockIdx.x) * static_cast<uint64_t>(blockDim.x) +
                           static_cast<uint64_t>(threadIdx.x);
    const uint64_t ind1 =
        static_cast<uint64_t>((-1 + sqrt(static_cast<double>(8 * index + 1))) / 2);
    const uint64_t ind2 = index - ind1 * (ind1 + 1) / 2;

    if (ind1 < s && ind2 < s) {
        const auto c   = elements[ind1];
        const auto aMc = c & a;
        const auto aJc = c | a;
        const auto bMc = c & b;
        const auto bJc = c | b;

        const auto d     = elements[ind2];
        const auto aMcMd = d & aMc;
        const auto bJcJd = d | bJc;
        const auto bMcMd = d & bMc;
        const auto aJcJd = d | aJc;

        const auto alpha = distBottom[idxCalc.index(aMcMd)] * distTop[idxCalc.index(bJcJd)];
        const auto beta  = distBottom[idxCalc.index(bMcMd)] * distTop[idxCalc.index(aJcJd)];

        A[shift + ind1 * s + ind2] = alpha;
        B[shift + ind1 * s + ind2] = beta;
        A[shift + ind2 * s + ind1] = alpha;
        B[shift + ind2 * s + ind1] = beta;
    }
}

static void fillMatrices(size_t start, size_t end,
                         const std::vector<std::vector<uint64_t>>& abValues, size_t s,
                         unsigned int deviceId, double* A, double* B, uint16_t* elements,
                         uint64_t* distTop, uint64_t* distBottom, IndexCalc idxCalc) {
    uint64_t numOfBlocks = (((s + 1) * s) / 2 + numberOfThreads - 1) / numberOfThreads;

    dim3 dimGrid(numOfBlocks);
    dim3 dimBlock(numberOfThreads);

    for (size_t batch = 0; batch < (end - start); ++batch) {
        cudaSetDevice(deviceId);
        fillMatrices<<<dimGrid, dimBlock>>>(static_cast<uint16_t>(abValues[start + batch][0]),
                                            static_cast<uint16_t>(abValues[start + batch][1]),
                                            elements, distTop, distBottom, idxCalc, s, A, B,
                                            static_cast<size_t>(s * s * batch));
    }
};

uint128T DeviceEnumerator::doEnumerationGPU(size_t start, size_t end,
                                            const std::vector<std::vector<uint64_t>>& abValues,
                                            size_t s, unsigned int deviceId) {
    const size_t batchCount = end - start;

    fillMatrices(start, end, abValues, s, deviceId, mAs[deviceId], mBs[deviceId],
                 mElements[deviceId], mDistTop[deviceId], mDistBottom[deviceId],
                 mIdxCalc[deviceId]);

    doMatMulStridedBatched(s, batchCount, mCs[deviceId], mAs[deviceId], mBs[deviceId],
                           mHandles[deviceId]);

    cudaError_t err{cudaGetLastError()};
    if (err != cudaSuccess) {
        std::cout << cudaGetErrorString(err) << std::endl;
    }

    return computeTrace(mCs[deviceId], start, end, abValues, s, deviceId, mResults[deviceId].get());
}
