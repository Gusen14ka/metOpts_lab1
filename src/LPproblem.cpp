#include "LPproblem.hpp"
#include "Logger.hpp"
#include "relations.hpp"
#include <algorithm>
#include <cstddef>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <vector>

#define LOG Logger::instance()

LPproblem::LPproblem(std::string const & filename) : targetFunc(TargetFunc(NofVar)){
    std::ifstream file(filename);
    if(!file || !file.is_open()){
        LOG.error("Bad file or way", "LPproblem");
        throw;
    }
    try{
        parser(file);
    }
    catch(std::exception& e){
        LOG.error(e.what(), "LPproblem");
        throw;
    }
}

void LPproblem::parser(std::ifstream & file){
    std::string tok;
    for(int i = 0; i < NofCon; ++i){
        for(int j = 0; j < NofVar; ++j){
            file >> constraints[i].coefs[j];
        }
        file >> tok;
        constraints[i].rel = rel_from_str(tok);
        file >> constraints[i].rhs;
    }

    for(int i = 0; i < NofVar; ++i){
        int x;
        file >> x;
        nonnegs[i] = x != 0 ? 1 : 0;
    }

    file >> tok;
    if(tok != "min" && tok != "max"){
        throw std::invalid_argument("Unknown type of target function");
    }
    targetFunc.minimaze = tok == "min";

    for(size_t i = 0; i < NofVar; ++i){
        file >> targetFunc.coefs[i];
    }
}

void LPproblem::print() const {
    using std::cout;
    using std::endl;

    cout << "Целевая фунция:" << endl;
    cout << (targetFunc.minimaze ? "min " : "max ");

    for (size_t i = 0; i < targetFunc.coefs.size(); ++i) {
        double c = targetFunc.coefs[i];
        if (i > 0 && c >= 0) cout << "+ ";
        cout << c << "*x" << i + 1 << " ";
    }
    cout << endl << endl;

    cout << "Ограничения:" << endl;

    for (size_t i = 0; i < constraints.size(); ++i) {
        const Constraint &c = constraints[i];

        for (size_t j = 0; j < c.coefs.size(); ++j) {
            double a = c.coefs[j];
            if (j > 0 && a >= 0) cout << "+ ";
            cout << a << "*x" << j + 1 << " ";
        }

        cout << rel_to_str(c.rel) << " " << c.rhs << endl;
    }

    cout << endl;

    cout << "Ограничения на знак:" << endl;
    for (size_t i = 0; i < nonnegs.size(); ++i) {
        if (nonnegs[i] != 0) {
            cout << "x" << i + 1 << " >= 0" << endl;
        } else {
            cout << "x" << i + 1 << " free" << endl;
        }
    }

    cout << "\nФорма задачи: " << problemForm_to_str(get_type()) << endl;
}

ProblemForm LPproblem::get_type() const{
    if(!std::all_of(nonnegs.begin(), nonnegs.end(), [](int a){if(a != 1) return false; else return true;})){
        return ProblemForm::Common;
    }
    if(std::all_of(constraints.begin(), constraints.end(), [](Constraint const & c){
            if(c.rel == Relation::EQ) return true; else return false;})){
        return ProblemForm::Canon;
    }
    if(std::all_of(constraints.begin(), constraints.end(), [](Constraint const & c){
            if(c.rel == Relation::GE) return true; else return false;}) || 
            std::all_of(constraints.begin(), constraints.end(), [](Constraint const & c){
            if(c.rel == Relation::LE) return true; else return false;})){
        return ProblemForm::Symetric;
    }
    else return ProblemForm::Common;
}

void Constraint::multipl(double x){
    if(x < 0 && rel != Relation::EQ){
        rel = rel == Relation::GE ? Relation::LE : Relation::GE;
    }
    for(auto & el: coefs){
        el *= x;
    }
    rhs *= x;
}

LPproblem::LPproblem(size_t vars, size_t cons) : targetFunc(TargetFunc(vars)){
    nonnegs = std::vector<size_t>(vars, 0);
    constraints = std::vector<Constraint>(cons, Constraint(vars));
}

void LPproblem::fixCons(LPproblem& problem){
    std::vector<size_t> toSwitch;
    toSwitch.reserve(5);
    for(size_t i = 0; i < problem.constraints.size(); ++i){
        if((problem.targetFunc.minimaze && problem.constraints[i].rel == Relation::LE) || 
                (!problem.targetFunc.minimaze && problem.constraints[i].rel == Relation::GE)){
            problem.constraints[i].multipl(-1);
        }
        if(problem.constraints[i].rel == Relation::EQ){
            problem.constraints[i].rel = problem.targetFunc.minimaze ? Relation::GE : Relation::LE;
            toSwitch.push_back(i);
        }
    }
    for(auto& el: toSwitch){
        problem.constraints.emplace_back(problem.constraints[el]);
        problem.constraints[problem.constraints.size() - 1].multipl(-1);
        problem.constraints[problem.constraints.size() - 1].rel = problem.targetFunc.minimaze ? Relation::GE : Relation::LE;
    }
}

LPproblem LPproblem::get_symetric() const{
    // Копируем задачу
    LPproblem problem = *this;
    problem.switches.clear();
    if(get_type() == ProblemForm::Symetric) return problem;
    // Меняем <= на >= если задача на min или наоборот если на max
    fixCons(problem);

    // Подготовка индексов свободных переменных
    std::vector<size_t> toSwitch;
    toSwitch.reserve(5);
    for(size_t i = 0; i < nonnegs.size(); ++i){
        if(nonnegs[i] != 1) toSwitch.push_back(i);
    }
    // Записываем замену в словарь для обратной замены
    for(size_t i = 0; i < toSwitch.size(); ++i){
        problem.switches[toSwitch[i]] = {toSwitch[i], nonnegs.size() - 1 + toSwitch[i]};
    }
    // Делаем замену в ограничениях
    for(auto& el: problem.constraints){
        for(auto& idx: toSwitch){
            el.coefs.push_back(el.coefs[idx] * -1);
        }
    }
    // Делаем замену в целевой функции
    for(auto& idx: toSwitch){
        problem.targetFunc.coefs.push_back(problem.targetFunc.coefs[idx] * -1);
    }

    // Меняем сами ограничения на знак

    for(auto& idx: toSwitch){
        problem.nonnegs[idx] = 1;
        problem.nonnegs.push_back(1);
    }
    return problem;
}

LPproblem LPproblem::get_canon() const{
    // Копируем задачу
    LPproblem problem = *this;
    problem.switches.clear();
    if(get_type() == ProblemForm::Canon) return problem;    

    // Подготовка индексов свободных переменных
    std::vector<size_t> toSwitch;
    toSwitch.reserve(5);
    for(size_t i = 0; i < nonnegs.size(); ++i){
        if(nonnegs[i] != 1) toSwitch.push_back(i);
    }
    // Записываем замену в словарь для обратной замены
    for(size_t i = 0; i < toSwitch.size(); ++i){
        problem.switches[toSwitch[i]] = {toSwitch[i], nonnegs.size() - 1 + toSwitch[i]};
    }
    // Делаем замену в ограничениях
    for(auto& el: problem.constraints){
        for(auto& idx: toSwitch){
            el.coefs.push_back(el.coefs[idx] * -1);
        }
    }
    // Делаем замену в целевой функции
    for(auto& idx: toSwitch){
        problem.targetFunc.coefs.push_back(problem.targetFunc.coefs[idx] * -1);
    }

    // Меняем сами ограничения на знак
    for(auto& idx: toSwitch){
        problem.nonnegs[idx] = 1;
        problem.nonnegs.push_back(1);
    }
    toSwitch.clear();
    /*
    Находим неравенства для их замены на равенства путем добавления избыточных и слабительных переменных
    Далее будем добавлять их коэффиценты в концы массивов
    */
    std::vector<size_t> toSwitchLE, toSwitchGE;
    toSwitchGE.reserve(5); toSwitchLE.reserve(5);
    for(size_t i = 0; i < constraints.size(); ++i){
        if(constraints[i].rel == Relation::GE) toSwitchGE.push_back(i);
        else if(constraints[i].rel == Relation::LE) toSwitchLE.push_back(i);
    }
    // Учитываем добавочные переменные в порядке возрастания индексов их "родных" ограничений
    problem.inequals.reserve(5);
    auto itGE = toSwitchGE.begin();
    auto itLE = toSwitchLE.begin();
    while(itGE != toSwitchGE.end() || itLE != toSwitchLE.end()){
        if((itGE == toSwitchGE.end() &&  itLE != toSwitchLE.end()) || *itGE > *itLE){
            problem.inequals.push_back(*itLE);
            for(size_t i = 0; i < problem.constraints.size(); ++i){
                if(*itLE == i){
                    problem.constraints[i].coefs.push_back(1);
                    problem.constraints[i].rel = Relation::EQ;
                }
                else{
                    problem.constraints[i].coefs.push_back(0);
                }
            }
            ++itLE;
        }
        else if (itGE != toSwitchGE.end()){
            problem.inequals.push_back(*itGE);
            for(size_t i = 0; i < problem.constraints.size(); ++i){
                if(*itGE == i){
                    problem.constraints[i].coefs.push_back(-1);
                    problem.constraints[i].rel = Relation::EQ;
                }
                else{
                    problem.constraints[i].coefs.push_back(0);
                }
            }
            ++itGE;
        }
    }
    for(size_t i = 0; i < toSwitchGE.size() + toSwitchLE.size(); ++i) {
        problem.nonnegs.push_back(1);
        problem.targetFunc.coefs.push_back(0);
    }
    return problem;
}

LPproblem LPproblem::get_dual() const {
    int fact = targetFunc.minimaze ? 1 : -1;
    auto fixedprbl = *this;
    fixCons(fixedprbl);

    LPproblem problem = LPproblem(fixedprbl.constraints.size(), fixedprbl.targetFunc.coefs.size());
    // Задача становится обратной
    problem.targetFunc.minimaze = !fixedprbl.targetFunc.minimaze;
    // Правые части ограничений становятся коэфициентами функции цели
    for(size_t i = 0; i < fixedprbl.constraints.size(); ++i){
        problem.targetFunc.coefs[i] = fixedprbl.constraints[i].rhs;
    }
    // Коэфы функции цели становятся правыми частями ограничений
    for(size_t i = 0; i < fixedprbl.targetFunc.coefs.size(); ++i){
        problem.constraints[i].rhs = fixedprbl.targetFunc.coefs[i];
    }

    // Строки ограничений становятся столбцами ограничений
    for(size_t i = 0; i < fixedprbl.constraints.size(); ++i){
        for(size_t j = 0; j < fixedprbl.constraints[i].coefs.size(); ++j){
            problem.constraints[j].coefs[i] = fixedprbl.constraints[i].coefs[j];
        }
        // Знак неравенства переходит в ограничение по знаку для переменной
        if((targetFunc.minimaze && fixedprbl.constraints[i].rel == Relation::GE) || 
                (!targetFunc.minimaze && fixedprbl.constraints[i].rel == Relation::LE)){
            problem.nonnegs[i] = 1;
        }
    }

    // Ограничение по знаку для переменной переходит в знак неравенства
    for(size_t i = 0; i < fixedprbl.nonnegs.size(); ++i){
        if(nonnegs[i] == 1) problem.constraints[i].rel = fixedprbl.targetFunc.minimaze ? Relation::LE : Relation::GE;
        else problem.constraints[i].rel = Relation::EQ;
    }

    return problem;
}