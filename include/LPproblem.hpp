#pragma once

#include "problemForm.hpp"
#include <cstddef>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <relations.hpp>

constexpr size_t NofVar = 5;
constexpr size_t NofCon = 4; 


struct Constraint{
    std::vector<double> coefs = std::vector<double>(NofVar, 0);
    Relation rel;
    double rhs;

    void multipl(double x);
    Constraint(size_t vars) : coefs(std::vector<double>(vars, 0)){}
};

struct TargetFunc{
    std::vector<double> coefs = std::vector<double>(NofVar, 0);
    bool minimaze;

    TargetFunc(size_t vars) : coefs(std::vector<double>(vars, 0)) {}
};

class LPproblem{
public:
    LPproblem(std::string const & filename);
    LPproblem(size_t vars, size_t cons);
    void print() const;
    ProblemForm get_type() const;
    LPproblem get_symetric() const;
    LPproblem get_canon() const;
    LPproblem get_dual() const;
private:
    std::vector<size_t> nonnegs = std::vector<size_t>(NofVar, 0);
    std::vector<Constraint> constraints = std::vector<Constraint>(NofCon, Constraint(NofVar));
    TargetFunc targetFunc;
    std::unordered_map<size_t, std::pair<size_t, size_t>> switches; // Если пустой значит задача в той же форме, что и была по условию
    std::vector<size_t> inequals; // Массив индексов ограничений, которые являются неравенствами (или были ими до приведения к канон форме)

    void parser(std::ifstream & file);

    // Меняет <= на >= если задача на min или наоборот если на max
    static void fixCons(LPproblem& problem);

    friend void runSimplex(LPproblem problem);
};