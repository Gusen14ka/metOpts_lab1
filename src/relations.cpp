#include "relations.hpp"

Relation rel_from_str (const std::string &s){
    if(s == "<=") return Relation::LE;
    if(s == ">=") return Relation::GE;
    if(s == "=") return Relation::EQ;
    else{
        throw;
    }
}

std::string rel_to_str(Relation const & r){
    switch (r) {
        case Relation::LE: return "<=";
        case Relation::GE: return ">=";
        case Relation::EQ: return "=";
    }
    return "?";
}