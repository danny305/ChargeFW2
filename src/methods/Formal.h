//
// Created by krab1k on 12/11/18.
//

#pragma once


#include <boost/config.hpp>

#include "../Method.h"

class Formal : public Method {

public:
    explicit Formal() : Method({}, {}, {}) {};

    std::vector<double> calculate_charges(const Molecule &molecule) override;
};

extern "C" BOOST_SYMBOL_EXPORT Formal method;
Formal method;