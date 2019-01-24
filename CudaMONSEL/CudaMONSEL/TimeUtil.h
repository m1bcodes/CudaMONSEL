#include <ctime>
#include <iostream>

#define TimeFunc(func) { std::time_t start = time(0); func; std::time_t end = time(0); double timePassed = difftime(end, start) * 1000.0; std::cout << timePassed << std::endl; }
