#include "simplex.hpp"
#include "LPproblem.hpp"
#include "problemForm.hpp"
#include <cstddef>
#include <vector>
#include <iostream>
#include <iomanip>

void runSimplex(LPproblem problem){
    const double EPS = 1e-12;

    // Работаем только с канонической формой
    if(problem.get_type() != ProblemForm::Canon) throw std::runtime_error("runSimplex: ожидается Canon");

    // Делаем неотрицательной правую часть ограничений
    for(auto& el: problem.constraints){
        if(el.rhs < 0) {
            el.multipl(-1);
        }
    }

    // max(f) = min(-f)
    if(!problem.targetFunc.minimaze){
        for(auto& el: problem.targetFunc.coefs){
            el *= -1;
        }
    }

    // Для удобства создадим симлекс матрицу
    // | A | I (искусственные) | b |
    size_t nVars = problem.constraints[0].coefs.size(); // число исходных переменных
    size_t nCons = problem.constraints.size();           // число ограничений (m)
    size_t rows = nCons + 1;                             // m строк + строка оценок
    size_t cols = nVars + nCons + 1;                     // n + m столбцов переменных + RHS

    std::vector<std::vector<double>> mat(
        rows, std::vector<double>(cols, 0.0)
    );

    // Ограничения: заполняем A, единичную матрицу искусственных и RHS
    for (size_t i = 0; i < nCons; ++i) {
        for (size_t j = 0; j < nVars; ++j)
            mat[i][j] = problem.constraints[i].coefs[j];

        // искусственный столбец: позиция nVars + i
        mat[i][nVars + i] = 1.0;

        // RHS
        mat[i][cols - 1] = problem.constraints[i].rhs;
    }

    // Инициализируем базис: изначально в базисе все искусственные переменные
    std::vector<size_t> basis(nCons);
    for (size_t i = 0; i < nCons; ++i) basis[i] = nVars + i;

    // Считаем первоначальные оценки (Phase I): Δ_j = - sum_i a_ij
    for(size_t j = 0; j < nVars; ++j){
        double sum = 0.0;
        for(size_t i = 0; i < nCons; ++i){
            sum += mat[i][j];
        }
        mat[rows-1][j] = -sum;
    }

    // Опционально: можно печатать начальную таблицу для отладки

    auto printTable = [&](int iter){
        std::cout << "=== Phase I iter " << iter << " ===\n";
        std::cout << std::setw(6) << "bas";
        for (size_t j = 0; j < cols - 1; ++j) std::cout << std::setw(8) << ("x" + std::to_string(j+1));
        std::cout << std::setw(10) << "RHS" << std::setw(12) << "resid\n";

        // Печатаем таблицу (включая RHS)
        for (size_t i = 0; i < nCons; ++i) {
            std::cout << std::setw(6) << ("x" + std::to_string(basis[i]+1));
            for (size_t j = 0; j < cols; ++j) std::cout << std::setw(8) << std::setprecision(6) << mat[i][j];
            // временно оставим перенос строки, невязку посчитаем ниже
            std::cout << "\n";
        }
        std::cout << std::setw(6) << "Δ";
        for (size_t j = 0; j < cols; ++j) std::cout << std::setw(8) << std::setprecision(6) << mat[rows-1][j];
        std::cout << "\n";

        // --- Вычислим текущий базисный вектор x (для исходных nVars переменных) ---
        std::vector<double> cur_x(nVars, 0.0);
        for (size_t i = 0; i < nCons; ++i) {
            if (basis[i] < nVars) {
                cur_x[basis[i]] = mat[i][cols - 1];
            }
        }

        // --- Вычислим невязки для каждого ограничения: residual = A*x - b ---
        double max_abs_resid = 0.0;
        std::cout << "Residuals (A*x - b) per constraint:\n";
        for (size_t i = 0; i < nCons; ++i) {
            double lhs = 0.0;
            for (size_t j = 0; j < nVars; ++j) lhs += problem.constraints[i].coefs[j] * cur_x[j];
            double resid = lhs - problem.constraints[i].rhs;
            max_abs_resid = std::max(max_abs_resid, std::fabs(resid));
            std::cout << " row " << (i+1) << ": resid = " << std::setprecision(8) << resid
                      << "   (LHS=" << lhs << " RHS=" << problem.constraints[i].rhs << ")\n";
        }
        std::cout << " max |resid| = " << std::setprecision(8) << max_abs_resid << "\n\n";
    };
    printTable(0);

    // Основной цикл Phase I
    const int MAX_ITERS = 10000;
    int iter = 0;
    while(true){
        ++iter;
        if(iter > MAX_ITERS) throw std::runtime_error("runSimplex: превышен лимит итераций Phase I");

        // 1) Выбор входящего столбца (самый отрицательный Δ). Если нет — оптимум Phase I.
        double bestDelta = 0.0;
        size_t enterIdx = cols; // признак "не найден"
        for(size_t j = 0; j < cols - 1; ++j){
            if(mat[rows-1][j] < bestDelta - EPS){
                bestDelta = mat[rows-1][j];
                enterIdx = j;
            }
        }
        if(enterIdx == cols){
            // все Δ >= 0 => оптимум Phase I
            printTable(iter);
            break;
        }
        // 2) Ratio test — выбираем строку с минимальным положительным отношением b_i / a_i,enter
        double minRel = std::numeric_limits<double>::infinity();
        size_t leaveRow = rows; // признак "не найден"
        for(size_t i = 0; i < nCons; ++i){
            double a = mat[i][enterIdx];
            if(a > EPS){ // Только позитивные коэффициенты допустимы
                double rel = mat[i][cols - 1] / a;
                if (rel < minRel - EPS) {
                    minRel = rel;
                    leaveRow = i;
                }
            }
        }
        if (leaveRow == rows) {
            throw std::runtime_error("Ошибка: целевая функция не ограничена снизу");
        }

        // 3) Пересчет симплекс-таблицы (Метод Гаусса-Жордана)
        basis[leaveRow] = enterIdx;
        double pivotVal = mat[leaveRow][enterIdx];

        // Нормализация ведущей строки
        for (size_t j = 0; j < cols; ++j) mat[leaveRow][j] /= pivotVal;

        // Исключение переменной из остальных строк (включая строку оценок)
        for (size_t i = 0; i < rows; ++i) {
            if (i != leaveRow) {
                double factor = mat[i][enterIdx];
                for (size_t j = 0; j < cols; ++j) {
                    mat[i][j] -= factor * mat[leaveRow][j];
                }
            }
        }

        // Печатаем таблицу и невязки после pivot'а
        printTable(iter);
    }
  // Конец цикла Phase I

    // --- Анализ результатов Phase I ---
    // Вспомогательная функция минимизировала сумму искусственных переменных. 
    // Согласно источникам, если минимум > 0, исходная задача не имеет допустимых решений [1].
    // Текущее значение функции (сумма y_i) находится в mat[rows-1][cols-1] (с инверсией знака).
    const double W_EPS = 1e-8; // можно увеличить до 1e-7 при численной нестабильности
    double W = 0.0;
    for (size_t i = 0; i < nCons; ++i) {
        size_t col = basis[i];
        // искусственные столбцы имеют индексы от nVars до nVars + nCons - 1
        if (col >= nVars && col < nVars + nCons) {
            // значение переменной = RHS этой строки
            W += mat[i][cols - 1];
        }
    }
    if (W > W_EPS) {
        std::cout << "Задача не имеет допустимых решений (Phase I: W = " << W << " > " << W_EPS << ").\n";
        // диагностическая печать полезна при отладке:
        std::cerr << "DEBUG: basis and RHS after Phase I:\n";
        for (size_t i = 0; i < nCons; ++i) {
            std::cerr << " row " << i << " basis x" << basis[i]+1 << " = " << mat[i][cols - 1] << "\n";
        }
        // при желании — вывести всю таблицу
        return;
    }

    // --- Переход к Phase II ---
    
    // СМЕНА БАЗИСА (Вытеснение искусственных переменных) [1, 2]
    // Если Фаза I завершилась успехом (L=0), но в базисе остались искусственные переменные 
    // (их значения в mat[i][cols-1] равны 0), их нужно заменить на реальные переменные.
    for (size_t i = 0; i < nCons; ++i) {
        if (basis[i] >= nVars) { // Индекс искусственной переменной
            bool replaced = false;
            for (size_t j = 0; j < nVars; ++j) {
                // Ищем в этой строке ненулевой элемент для реальной переменной [2]
                if (std::abs(mat[i][j]) > EPS) {
                    // Выполняем Pivot (шаг Жордана-Гаусса) для замены базиса [5]
                    double pivotVal = mat[i][j];
                    for (size_t k = 0; k < cols; ++k) mat[i][k] /= pivotVal;
                    for (size_t k = 0; k < rows; ++k) {
                        if (k != i) {
                            double factor = mat[k][j];
                            for (size_t l = 0; l < cols; ++l) mat[k][l] -= factor * mat[i][l];
                        }
                    }
                    basis[i] = j;
                    replaced = true;
                    break;
                }
            }
            // Если replaced == false, значит строка i полностью нулевая для реальных переменных.
            // Согласно [2], это означает, что данное ограничение является лишним (линейно зависимым).
        }
    }

    // 2. Формируем новую строку оценок для исходной функции цели f(x) = c*x
    // Теперь все basis[i] гарантированно < nVars (если задача не имела лишних ограничений)
    for (size_t j = 0; j < nVars; ++j) {
        double z_j = 0.0;
        for (size_t i = 0; i < nCons; ++i) {
            z_j += problem.targetFunc.coefs[basis[i]] * mat[i][j];
        }
        mat[rows-1][j] = problem.targetFunc.coefs[j] - z_j;
    }

    // Считаем текущее значение f(x)
    double current_f = 0.0;
    for (size_t i = 0; i < nCons; ++i) {
        current_f += problem.targetFunc.coefs[basis[i]] * mat[i][cols-1];
    }
    mat[rows-1][cols-1] = current_f;

    // --- Итерации Phase II ---
    int iter2 = 0;
    while (true) {
        ++iter2;
        // 1) Выбор входящего столбца
        size_t enterIdx = cols;
        double minDelta = -EPS;
        for (size_t j = 0; j < nVars; ++j) {
            if (mat[rows - 1][j] < minDelta) {
                minDelta = mat[rows - 1][j];
                enterIdx = j;
            }
        }

        if (enterIdx == cols) break; // Оптимум найден [6]

        // 2) Ratio test
        double minRatio = std::numeric_limits<double>::infinity();
        size_t leaveRow = rows;

        for (size_t i = 0; i < nCons; ++i) {
            double a_ij = mat[i][enterIdx];
            if (a_ij > EPS) {
                double ratio = mat[i][cols - 1] / a_ij;
                // Если mat[i][cols-1] около 0, ratio будет 0.
                // Это инициирует "смену базиса" (pivot с нулевым шагом) [3, 4].
                if (ratio < minRatio - EPS) {
                    minRatio = ratio;
                    leaveRow = i;
                }
            }
        }

        if (leaveRow == rows) {
            throw std::runtime_error("Функция не ограничена снизу [7]");
        }

        // 3) Pivot (пересчет всей таблицы)
        basis[leaveRow] = enterIdx;
        double pivot = mat[leaveRow][enterIdx];
        for (size_t j = 0; j < cols; ++j) mat[leaveRow][j] /= pivot;
        for (size_t i = 0; i < rows; ++i) {
            if (i != leaveRow) {
                double factor = mat[i][enterIdx];
                for (size_t j = 0; j < cols; ++j) mat[i][j] -= factor * mat[leaveRow][j];
            }
        }
    }
    

    // Вывод результата
    /*
    std::cout << "Оптимальное значение: " << mat[rows - 1][cols - 1] << std::endl;
    for (size_t i = 0; i < nCons; ++i) {
        std::cout << "x" << basis[i] + 1 << " = " << mat[i][cols - 1] << std::endl;
    }
    */

    // Получение решения только для исходных переменных
    std::vector<double> x(nVars, 0.0);
    for (size_t i = 0; i < nCons; ++i) {
        if (basis[i] < nVars) {
            x[basis[i]] = mat[i][cols - 1]; // значение базисной исходной переменной
        }
    }

    /*
    // Вывод исходных переменных
    std::cout << "Solution (original variables):\n";
    for (size_t j = 0; j < nVars; ++j) {
        std::cout << "x" << (j+1) << " = " << x[j] << "\n";
    }

    // Вывод значений оставшихся базисных (вспомогательных) переменных
    std::cout << "Base auxiliary/slack/artificial variables:\n";
    for (size_t i = 0; i < nCons; ++i) {
        size_t col = basis[i];
        if (col >= nVars) {
            std::cout << "var x" << (col+1) << " (aux) = " << mat[i][cols - 1] << "\n";
        }
    }
    */
    /*
    // Проверка ограничений: A*x ?= b (дополнительная печать невязок финального решения)
    std::cout << "Final constraints residuals (post-solution):\n";
    for (size_t i = 0; i < nCons; ++i) {
        double lhs = 0.0;
        for (size_t j = 0; j < nVars; ++j) lhs += problem.constraints[i].coefs[j] * x[j];
        double rhs = problem.constraints[i].rhs;
        double resid = lhs - rhs;
        std::cout << "Row " << i+1 << ": LHS=" << lhs << " RHS=" << rhs
                  << " residual=" << std::setprecision(12) << resid << "\n";
    }
    */

    // Делаем обратные замены если они есть
    std::vector<double> real_ans(NofVar, 0.0);
    for(size_t i = 0; i < NofVar; ++i){
        real_ans[i] = x[i];
    }
    if(!problem.switches.empty()){
        for(auto& el: problem.switches){
            real_ans[el.first] = x[el.second.first] - x[el.second.second];
        }
    }

    // Вывод оптимальных значений
    for(size_t i = 0; i < real_ans.size(); ++i){
        std::cout << "x" << (i+1) << " = " << real_ans[i] << "\n";
    }

    // Считаем целевую функцию по оригинальным коэффициентам
    double obj = 0.0;
    for (size_t j = 0; j < nVars; ++j) obj += problem.targetFunc.coefs[j] * x[j];
    if (!problem.targetFunc.minimaze) {
        obj *= -1;
    }
    std::cout << "Значание функции цели: " << obj << "\n";
}