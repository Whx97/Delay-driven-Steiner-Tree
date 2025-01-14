#include <iostream>
#include <fstream>
#define COUT_RED_START std::cout << "\033[1;31m";
#define COUT_GREEN_START std::cout << "\033[1;32m";
#define COUT_YELLOW_START std::cout << "\033[1;33m";
#define COUT_BLUE_START std::cout << "\033[1;34m";
#define COUT_PURPLE_START std::cout << "\033[1;35m";
#define COUT_CYAN_START std::cout << "\033[1;36m";

#define COUT_COLOR_END std::cout << "\033[0m" << std::endl;

template <typename... T>
void printGreenLog(T... t) {
    COUT_GREEN_START;
    (std::cout << ... << t);
    COUT_COLOR_END;
}