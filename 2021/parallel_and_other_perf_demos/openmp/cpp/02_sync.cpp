#include <iostream>

int main()
{
    #pragma omp parallel
    {
        #pragma omp critical
        {
            std::cout << "Hello World" << std::endl;
        }
    }
    return 0;
}
