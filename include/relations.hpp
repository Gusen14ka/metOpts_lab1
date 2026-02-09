#pragma once
#include <string>


enum class Relation {LE, GE, EQ}; // <=, =>, =

Relation rel_from_str(std::string const &s);

std::string rel_to_str(Relation const & r);