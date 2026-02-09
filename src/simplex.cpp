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
        std::cout << std::setw(10) << "RHS\n";
        for (size_t i = 0; i < nCons; ++i) {
            std::cout << std::setw(6) << ("x" + std::to_string(basis[i]+1));
            for (size_t j = 0; j < cols; ++j) std::cout << std::setw(8) << std::setprecision(6) << mat[i][j];
            std::cout << "\n";
        }
        std::cout << std::setw(6) << "Δ";
        for (size_t j = 0; j < cols; ++j) std::cout << std::setw(8) << std::setprecision(6) << mat[rows-1][j];
        std::cout << "\n\n";
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
            //printTable(iter);
            break;
        }

        // 2) Ratio test — выбираем строку с минимальным положительным отношением b_i / a_i,enter
        double minRel = std::numeric_limits<double>::infinity();
        size_t leaveRow = rows; // признак "не найден"
        for(size_t i = 0; i < nCons; ++i){
            double a = mat[i][enterIdx];
            if(a > EPS){ // Только позитивные коэффициенты допустимы
                double rel = mat[i][cols - 1] / a;
                // Вырождение: если rel == minRel — выбираем меньший индекс строки (Bland-like)
                if(rel + EPS < minRel || (std::abs(rel - minRel) <= EPS && i < leaveRow)){
                    minRel = rel;
                    leaveRow = i;
                }
            }
        }

        if(leaveRow == rows){
            // Нет положительных элементов в выбранном столбце => направление убывания целевой неограничено.
            // В Phase I это означает, что W можно уменьшать без ограничения — на практике это сигнал о проблеме построения,
            // но корректно — бросить понятное исключение.
            throw std::runtime_error("runSimplex: no positive entries in entering column -> unbounded direction (Phase I)");
        }

        // 3) Pivot: нормализация ведущей строки и обнуление столбца
        double pivot = mat[leaveRow][enterIdx];
        if(std::abs(pivot) < EPS) throw std::runtime_error("runSimplex: pivot ~= 0 (numerical issue)");

        // нормализуем ведущую строку
        for(size_t j = 0; j < cols; ++j){
            mat[leaveRow][j] /= pivot;
        }

        // обнуляем столбец enterIdx во всех остальных строках (включая строку Δ)
        for(size_t i = 0; i < rows; ++i){
            if(i == leaveRow) continue;
            double factor = mat[i][enterIdx];
            if(std::abs(factor) <= EPS) continue;
            for(size_t j = 0; j < cols; ++j){
                mat[i][j] -= factor * mat[leaveRow][j];
                // числовая чистка
                if(std::abs(mat[i][j]) < EPS) mat[i][j] = 0.0;
            }
        }

        // 4) Обновление базиса
        basis[leaveRow] = enterIdx;

        // Опционально: вывод после pivot
        //printTable(iter);
    }

    // Phase I завершён — проверяем W (сумма искусственных). Надёжнее считать сумму значений искусственных переменных:
    double W = 0.0;
    for(size_t i = 0; i < nCons; ++i){
        size_t col = basis[i];
        if(col >= nVars && col < nVars + nCons){ // если базисная переменная — искусственная
            W += mat[i][cols - 1]; // её значение = b_i
        }
    }
    // также можно проверить через последнюю ячейку строки Δ: mat[rows-1][cols-1] = -sum(rows) по инициализации,
    // но при pivot-ах явное суммирование базисных искусственных безопаснее.

    if(std::abs(W) > 1e-8){
        // нет допустимого решения
        throw std::runtime_error("runSimplex: Phase I failed, problem infeasible (W > 0).");
    }

    std::cout << "runSimplex: Phase I completed successfully (W == 0). Base variables:\n";
    for(size_t i = 0; i < nCons; ++i){
        std::cout << " row " << i << " -> var x" << (basis[i] + 1) << " = " << mat[i][cols - 1] << "\n";
    }

    // Phase II

    // ---------- НОВЫЙ БЛОК: выталкиваем искусственные из базиса перед удалением столбцов ----------
    // Для каждой строки, где в базисе искусственная переменная (basis[i] >= nVars),
    // пытаемся найти ненулевой столбец в зоне реальных переменных (0..nVars-1) и сделать pivot.
    for (size_t i = 0; i < nCons; ++i) {
        size_t bj = basis[i];
        if (bj >= nVars && bj < nVars + nCons) {
            // искусственная в базисе — пробуем заменить
            bool replaced = false;
            for (size_t j = 0; j < nVars; ++j) {
                if (std::abs(mat[i][j]) > EPS) {
                    // Сделаем pivot в (i, j) — аналогично циклу pivot'а
                    double pv = mat[i][j];
                    // Нормализуем строку i
                    for (size_t col = 0; col < cols; ++col) mat[i][col] /= pv;
                    // Обнуляем столбец j во всех остальных строках
                    for (size_t r = 0; r < rows; ++r) {
                        if (r == i) continue;
                        double factor = mat[r][j];
                        if (std::abs(factor) <= EPS) continue;
                        for (size_t col = 0; col < cols; ++col) {
                            mat[r][col] -= factor * mat[i][col];
                            if (std::abs(mat[r][col]) < EPS) mat[r][col] = 0.0;
                        }
                    }
                    // обновляем базис
                    basis[i] = j;
                    replaced = true;
                    break;
                }
            }
            if (!replaced) {
                // Не нашли ненулевой столбец в реальной части.
                // Если RHS ~= 0, строка линейно зависима — допустимо (оставим строку как есть).
                // Если RHS != 0 — это означает, что Phase I дал корявый результат (редко),
                // и корректнее сообщить об ошибке.
                if (std::abs(mat[i][cols - 1]) > 1e-8) {
                    throw std::runtime_error("runSimplex: cannot remove artificial from basis, inconsistent row after Phase I");
                } else {
                    // строка нулевая, искусственная остается, но её значение 0 (можно удалить столбец позже)
                    // Для отладки — можно напечатать предупреждение:
                    std::cerr << "runSimplex: warning — artificial remains in basis for row " << i
                              << " but RHS == 0; row is dependent and will be ignored.\n";
                }
            }
        }
    }
    // ---------- конец блока вытолкования искусственных ----------

    // Удаляем столбцы искусственных переменных
    size_t newCols = nVars + 1; // только x и RHS
    std::vector<std::vector<double>> mat2(rows, std::vector<double>(newCols, 0.0));

    for(size_t i = 0; i < rows; ++i){
        for(size_t j = 0; j < nVars; ++j){
            mat2[i][j] = mat[i][j];
        }
        mat2[i][newCols - 1] = mat[i][cols - 1]; // RHS
    }

    mat.swap(mat2);
    cols = newCols;

    // Обнуляем строку оценок
    for(size_t j = 0; j < cols; ++j)
        mat[rows-1][j] = 0.0;

    // Записываем c_j
    for(size_t j = 0; j < nVars; ++j){
        mat[rows-1][j] = problem.targetFunc.coefs[j];
    }

    // Корректируем по базису: Δ = c - c_B * A
    for(size_t i = 0; i < nCons; ++i){
        size_t bj = basis[i];
        double cB = problem.targetFunc.coefs[bj];

        for(size_t j = 0; j < cols; ++j){
            mat[rows-1][j] -= cB * mat[i][j];
        }
    }


    while(true){
        // Входящий столбец
        double best = 0.0;
        size_t enter = cols;

        for(size_t j = 0; j < cols - 1; ++j){
            if(problem.targetFunc.minimaze){
                if(mat[rows-1][j] < best - EPS){  // для min: ищем < 0
                    best = mat[rows-1][j];
                    enter = j;
                }
            } else {
                if(mat[rows-1][j] > best + EPS){  // для max: ищем > 0
                    best = mat[rows-1][j];
                    enter = j;
                }
            }
        }

        if(enter == cols) break; // оптимум

        // Симплекс-отношение
        double minRel = std::numeric_limits<double>::infinity();
        size_t leave = rows;

        for(size_t i = 0; i < nCons; ++i){
            if(mat[i][enter] > 0){
                double r = mat[i][cols-1] / mat[i][enter];
                if(r < minRel){
                    minRel = r;
                    leave = i;
                }
            }
        }

        if(leave == rows)
            throw std::runtime_error("Unbounded solution");

        // Pivot
        double pivot = mat[leave][enter];
        for(size_t j = 0; j < cols; ++j)
            mat[leave][j] /= pivot;

        for(size_t i = 0; i < rows; ++i){
            if(i == leave) continue;
            double f = mat[i][enter];
            for(size_t j = 0; j < cols; ++j)
                mat[i][j] -= f * mat[leave][j];
        }

        basis[leave] = enter;
    }


    
    std::vector<double> x(nVars, 0.0);
    for(size_t i = 0; i < nCons; ++i){
        if(basis[i] < nVars)
            x[basis[i]] = mat[i][cols-1];
    }

    double value = mat[rows-1][cols-1];
    // Снова учитываем тип задачи
    if(problem.targetFunc.minimaze)
        value = -value;

    std::cout << "Результат Симплекс метода:" << std::endl;
    for(size_t i = 0; i < x.size(); ++i){
        std::cout << "x" << i+1 << " = " << x[i] << std::endl;
    }
    std::string text = problem.targetFunc.minimaze ? "Минимальное " : "Максимальное ";
    std::cout << text + "значение целевой функции = " << value << std::endl;
}
