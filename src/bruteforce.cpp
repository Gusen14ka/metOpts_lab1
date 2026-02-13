#include "bruteforce.hpp"
#include "problemForm.hpp"
#include <cstddef>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

const double EPS = 1e-9;

// Функция для решения СЛАУ методом Гаусса
// B - матрица коэффициентов (m x m), b - вектор правых частей
static bool solveSLAU(int m, vector<vector<double>>& B, vector<double>& b, vector<double>& res) {
    for (int i = 0; i < m; ++i) {
        int pivot = i;
        for (int j = i + 1; j < m; ++j) {
            if (abs(B[j][i]) > abs(B[pivot][i])) pivot = j;
        }
        if (abs(B[pivot][i]) < EPS) return false;
        swap(B[i], B[pivot]);
        swap(b[i], b[pivot]);

        double div = B[i][i];
        for (int j = i; j < m; ++j) B[i][j] /= div;
        b[i] /= div;

        for (int j = 0; j < m; ++j) {
            if (i != j) {
                double mult = B[j][i];
                for (int k = i; k < m; ++k) B[j][k] -= mult * B[i][k];
                b[j] -= mult * b[i];
            }
        }
    }
    res = b;
    return true;
}

void runBruteForce(LPproblem problem) {
    if (problem.get_type() != ProblemForm::Canon) throw std::runtime_error("runBruteForce: ожидается Canon");

    // max(f) = min(-f)
    if(!problem.targetFunc.minimaze){
        for(auto& el: problem.targetFunc.coefs){
            el *= -1;
        }
    }

    // Делаем неотрицательной правую часть ограничений - но вообще не обязательно
    for(auto& el: problem.constraints){
        if(el.rhs < 0) {
            el.multipl(-1);
        }
    }

    // Для удобства создадим общую матрицу
    size_t nVars = problem.constraints[0].coefs.size(); // число исходных переменных
    size_t nCons = problem.constraints.size();           // число ограничений (m)
    size_t rows = nCons + 1;                             // m строк + строка целевой функции
    size_t cols = nVars + 1;                     // n + RHS

    std::vector<std::vector<double>> mat(
        rows, std::vector<double>(cols, 0.0)
    );

    // Запонлняем матрицу
    // 0-строка - строка целевой функции
    // 0 
    for (size_t i = 0; i < nCons + 1; ++i) {
        if(i == nCons){
            // Целевая функция
            for (size_t j = 0; j < nVars; ++j)
                mat[i][j] = problem.targetFunc.coefs[j];
        }
        else{
            // Коэфы в ограничениях
            for (size_t j = 0; j < nVars; ++j)
                mat[i][j] = problem.constraints[i].coefs[j];

            // RHS
            mat[i][cols - 1] = problem.constraints[i].rhs;
        }
    }

    size_t n = nVars;
    size_t m = nCons;



    // Индексы столбцов переменных: 0, 1, ..., n-1
    vector<int> colIndices(n);
    for (int i = 0; i < n; ++i) colIndices[i] = i;

    // Массив для генерации сочетаний (m из n)
    vector<bool> v(n);
    fill(v.end() - m, v.end(), true); 

    vector<double> bestSolution(n, 0);
    double minF = problem.targetFunc.minimaze ? std::numeric_limits<double>::infinity() : -std::numeric_limits<double>::infinity();
    bool found = false;

    cout << "--- Перебор вершин (RHS в последнем столбце, Цель в последней строке) ---" << endl;

    do {
        vector<int> currentBasis;
        for (int i = 0; i < n; ++i) {
            if (v[i]) currentBasis.push_back(colIndices[i]);
        }

        // 1. Формируем матрицу B (m x m) и вектор b (m)
        vector<vector<double>> B(m, vector<double>(m));
        vector<double> b_vec(m);
        for (int i = 0; i < m; ++i) {
            b_vec[i] = mat[i][n]; // RHS из последнего столбца n
            for (int j = 0; j < m; ++j) {
                B[i][j] = mat[i][currentBasis[j]];
            }
        }

        // 2. Решаем СЛАУ для базисных переменных
        vector<double> x_basis(m);
        if (solveSLAU(m, B, b_vec, x_basis)) {
            // 3. Проверка допустимости (x >= 0) [163, Определение 4.3]
            bool feasible = true;
            vector<double> currentX(n, 0);
            for (int i = 0; i < m; ++i) {
                if (x_basis[i] < -EPS) feasible = false;
                currentX[currentBasis[i]] = x_basis[i];
            }

            if (feasible) {
                found = true;
                // 4. Расчет f = sum(c_j * x_j), где c_j в строке m
                double currentF = 0;
                for (int j = 0; j < n; ++j) {
                    currentF += mat[m][j] * currentX[j];
                }

                currentF *= !problem.targetFunc.minimaze ? -1 : 1;
                cout << "Крайняя точка: f = " << setw(10) << currentF << " | x = (";

                // Делаем обратные замены если они есть
                std::vector<double> real_ans(NofVar, 0.0);
                for(size_t i = 0; i < NofVar; ++i){
                    real_ans[i] = currentX[i];
                }
                if(!problem.switches.empty()){
                    for(auto& el: problem.switches){
                        real_ans[el.first] = currentX[el.second.first] - currentX[el.second.second];
                    }
                }

                for (int j = 0; j < real_ans.size(); ++j) cout << (abs(real_ans[j]) < EPS ? 0 : real_ans[j]) << (j == real_ans.size() - 1 ? "" : ", ");
                cout << ")" << endl;

                if ((problem.targetFunc.minimaze && currentF < minF) || (!problem.targetFunc.minimaze && currentF > minF)) {
                    minF = currentF;
                    for (int j = 0; j < n; ++j) bestSolution[j] = currentX[j];
                }
            }
        }
    } while (next_permutation(v.begin(), v.end()));

    if (found) {
        cout << "\nИТОГ: f = " << minF << " при x* = (";

        // Делаем обратные замены если они есть
        std::vector<double> real_ans(NofVar, 0.0);
        for(size_t i = 0; i < NofVar; ++i){
            real_ans[i] = bestSolution[i];
        }
        if(!problem.switches.empty()){
            for(auto& el: problem.switches){
                real_ans[el.first] = bestSolution[el.second.first] - bestSolution[el.second.second];
            }
        }

        for (int i = 0; i < real_ans.size(); ++i) cout << real_ans[i] << (i == real_ans.size() - 1 ? "" : ", ");
        cout << ")" << endl;
    } else {
        cout << "Допустимых опорных векторов нет." << endl;
    }
}
