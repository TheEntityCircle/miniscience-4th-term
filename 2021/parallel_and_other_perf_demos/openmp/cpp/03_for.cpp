#include <iostream>
#include <stdio.h>
#include <omp.h>

int main()
{
    long n = 10;
    long long sum = 0;
    #pragma omp parallel
    {
        long long local_sum = 0;
        #pragma omp for 
        for (long i = 0; i < n; ++i) {
            local_sum += i;
            for (long j = 0; j < 3; ++j) {
                printf("Thread %d values %d and %d\n", omp_get_thread_num(), i, j);
            }
        }
        #pragma omp atomic
        sum += local_sum;
    }

    std::cout << sum << std::endl;
    return 0;
}

