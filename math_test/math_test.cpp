#include <iostream>
#include "matrix.h"

// testing math library
int run_all_matrix_tests();
int run_all_vector_tests();

using namespace base_math;

int main() {
    run_all_matrix_tests();
    run_all_vector_tests();
    
    return 0;
}
