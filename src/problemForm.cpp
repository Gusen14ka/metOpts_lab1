#include "problemForm.hpp"

std::string problemForm_to_str(ProblemForm const & r){
    switch (r) {
        case ProblemForm::Common: return "Общая";
        case ProblemForm::Canon: return "Канноническая";
        case ProblemForm::Symetric: return "Симметрическая";
    }
    return "?";
}