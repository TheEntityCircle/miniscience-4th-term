#include <iostream>
#include <math.h>

// CUDA Kernel function to add the elements of two arrays on the GPU
__global__
void add(int n, float *x, float *y)
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (int i = index; i < n; i += stride)
    y[i] = x[i] + y[i];
}

int main(void)
{
  int N = 1e8;

  // Allocate host mem
  float *x = new float[N];
  float *y = new float[N];

  // initialize x and y arrays on the host
  for (int i = 0; i < N; i++) {
    x[i] = 1.0f;
    y[i] = 2.0f;
  }

  // Prepare GPU mem
  float *dev_a, *dev_b;
  cudaMalloc((void**)&dev_a, N*sizeof(float));
  cudaMalloc((void**)&dev_b, N*sizeof(float));
  cudaMemcpy(dev_a, x, N*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_b, y, N*sizeof(float), cudaMemcpyHostToDevice);

  // Run kernel on 1M elements on the GPU
  int blockSize = 256;
  int numBlocks = (N + blockSize - 1) / blockSize;
  add<<<numBlocks, blockSize>>>(N, dev_a, dev_b);

  // Wait for GPU to finish before accessing on host
  cudaDeviceSynchronize();

  // Get the results back
  cudaMemcpy(y, dev_b, N*sizeof(float), cudaMemcpyDeviceToHost);

  // Free GPU mem
  cudaFree(dev_a);
  cudaFree(dev_b);

  // Check for errors (all values should be 3.0f)
//  float maxError = 0.0f;
//  for (int i = 0; i < N; i++)
//    maxError = fmax(maxError, fabs(y[i]-3.0f));
//  std::cout << "Max error: " << maxError << std::endl;

  // Free memory
  delete [] x;
  delete [] y;

  return 0;
}
