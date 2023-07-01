#pragma once
#include <Uint128.hpp>
#include <thrust/device_ptr.h>
#include <vector>

// CUDA runtime
#include <cublas_v2.h>
#include <cuda_runtime.h>

const int numberOfThreads = 128;

/// @brief Index calculation on device. Here it is just a lookup,
/// but other means of calculation are possible.
struct IndexCalc {
    IndexCalc() = default;

    IndexCalc(uint8_t* idx) { mIdx = idx; };

    __device__ unsigned int index(uint64_t x) { return mIdx[x]; };

    uint8_t* mIdx;
};

/// @brief Class that facilitates the enumeration of the fee distributive lattice on GPU
class DeviceEnumerator {
public:
    /// @brief C'tor that initializes the device enumerator.
    /// @param idx index map for the elements
    DeviceEnumerator(const std::vector<uint8_t>& idx);

    /// @brief D'tor that frees all allocated memory.
    ~DeviceEnumerator();

    /// @brief Sets values for the top operator.
    void setTop(const std::vector<uint64_t>& distTopValues);

    /// @brief Sets values for the bottom operator.
    void setBottom(const std::vector<uint64_t>& distBottomValues);

    /// @brief Sets the elements.
    void setElements(const std::vector<uint16_t>& elements);

    /// @brief Frees the elements.
    void freeElements();

    /// @brief Frees the top operator.
    void freeTop();

    /// @brief Frees the bottom operator.
    void freeBottom();

    /// @brief Initializes the matrices.
    void initMatrices(size_t nrOfElements, size_t matrixBatchSize);

    /// @brief Frees the matrices.
    void freeMatrices();

    /// @brief Performs the enumeration on GPU. Matrices A,B are filled, multiplied C=A*B and the
    /// trace of C^2 is computed.
    /// @param start start index of ab-values
    /// @param end end index of ab-values
    /// @param abValues ab-values vector
    /// @param s interval size
    /// @param deviceID device to use
    uint128T doEnumerationGPU(size_t start, size_t end,
                              const std::vector<std::vector<uint64_t>>& abValues, size_t s,
                              unsigned int deviceID = 0);

    /// @brief Returns the number of devices.
    unsigned int getNumberOfDevices();

    /// @brief Returns the available GPU memory in bytes.
    size_t getAvailableGPUMemory();

private:
    /// @brief Copies a vector to device and returns a pointer to the device vector.
    template <typename T> T* toDevice(const std::vector<T>& v, unsigned int deviceID = 0) const {
        T* d;
        cudaSetDevice(deviceID);
        cudaMalloc((void**)&d, (v.size()) * sizeof(T));
        cudaMemcpy(d, v.data(), (v.size()) * sizeof(T), cudaMemcpyHostToDevice);

        return d;
    }

    /// @brief Multiplies a batch of matrices on GPU.
    void doMatMulStridedBatched(size_t s, int batchCount, double* C, double* A, double* B,
                                cublasHandle_t& handle);

    /// @brief Initializes device momory s*sizeof(double).
    /// @param s size of memory to allocate
    /// @param devideID device to use
    double* initOnDevice(size_t s, unsigned int deviceID);

    // matrix batches  for the computation C=A*B
    std::vector<double*> mAs{};
    std::vector<double*> mBs{};
    std::vector<double*> mCs{};

    // memory to compute the trace of C^2
    std::vector<thrust::device_ptr<uint128T>> mResults{};

    // CUDA handles
    std::vector<cublasHandle_t> mHandles{};

    // number of available devices
    unsigned int mNrOfDevices{};

    // index map for the elements
    std::vector<uint8_t*> mIdx{};

    // index calculation
    std::vector<IndexCalc> mIdxCalc;

    // bottom operator
    std::vector<uint64_t*> mDistBottom{};
    // top operator
    std::vector<uint64_t*> mDistTop{};
    // elements of FDL
    std::vector<uint16_t*> mElements{};
};
