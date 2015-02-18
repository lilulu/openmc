/* A test case for Lulu's C++ CMFD solver. */

#include <cstdio>

/* A simple test program that */
static void test_cxx_code(void) __attribute__((constructor));

void test_cxx_code(void)
{
    std::fprintf(stderr, "Called a C++ function!\n");
}
