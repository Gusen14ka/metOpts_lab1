#include "LPproblem.hpp"
#include "Logger.hpp"
#include "simplex.hpp"
#include <Windows.h>
#include <iostream>

#define LOG Logger::instance()

int main(){
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    setlocale(LC_ALL, "ru_RU.UTF-8");
    LOG.set_log_to_console(true);
    LPproblem problem("data.txt");
    std::cout << "--- 1. Ввод задачи --- \n" << std::endl;
    //problem.print();

    std::cout << "\n --- 2. Приведение к остальным формам ---\n" << std::endl;
    std::cout << "Симметрическая" << std::endl;
    auto symprbl = problem.get_symetric();
    //symprbl.print();

    std::cout << "\nКаноническая" << std::endl;
    auto canonprbl = problem.get_canon();
    //canonprbl.print();

    std::cout << "\n --- 3. Двойственные задачи ---\n" << std::endl;
    std::cout << "Двойственная общая" << std::endl;
    auto dualprbl = problem.get_dual();
    //dualprbl.print();

    std::cout << "\nДвойственная симметрическая" << std::endl;
    auto dualsymprbl = symprbl.get_dual();
    //dualsymprbl.print();

    std::cout << "\nДвойственная каноническая" << std::endl;
    auto dualcanonprbl = canonprbl.get_dual();
    //dualcanonprbl.print();

    runSimplex(canonprbl);
}