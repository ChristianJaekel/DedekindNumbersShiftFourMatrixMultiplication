This repository contains two implementations of the algorithm described in the preprint [A computation of the ninth Dedekind Number](https://arxiv.org/abs/2304.00895) by Christian Jäkel, 2023. Matrix multiplication is utilized to compute the ninth Dedekind Number. A CPU implementation based on Eigen and a GPU implementation based on CUDA is provided. Both have Intel TBB dependency.

CPU:
- navigate into the CPU folder
- Eigen and Intel TBB have to be installed
- execute:
    - mkdir build
    - cd build
    - cmake -DCMAKE_BUILD_TYPE=Release ..
    - make
- benchmarks: navigate into the CPU/build/benchmark folder and execute "dedekind_benchmarks"
- tests: navigate into the CPU/build/test folder and execute "dedekind_tests"  

GPU:
- navigate into the GPU folder
- Intel TBB has to be installed
- a recent version of CUDA is required
- run make  

For CPU and GPU, the executable is called "dedekind". It computes the 6th, 7th and 8th [dedekind number](https://en.wikipedia.org/wiki/Dedekind_number).

Timings to compute the 8th Dedekind number with different GPUs.

Nvidia Quadro M2200: 8.78s  
Nvidia RTX A6000: 2.95s  
Nvidia A10: 2.76s  
Nvidia H100: 2.75s  

Timings to compute the 8th Dedekind number with different CPUs.  

Intel Core i5-2430M: 12.8s  
Intel Core i7-7920HQ: 2.1s  
Intel Core i7-10850H: 1.92s  

Timings to compute the 8th Dedekind number with different CPUs with 1 thread.  

Intel Core i5-2430M: 23.9s  
Intel Core i7-7920HQ: 7.8s  (@4.1GHz)  
Intel Core i7-10850H: 11.46s (@2.7GHz)  

The benefits of the GPU implementation can only be experienced for the computation of the 9th Dedekind number, since the involved matrices are much larger and numerous. For this computation, approximately 20 GB of precomputed data is required.

On some systems, compilation issues with clang++ occurred. Installing "libstdc++-12-dev" solved the problem.
