//
// Created by krab1k on 24/10/18.
//

#pragma once

#include <string>

#include "Reader.h"

class SDF: public Reader {

public:
    MoleculeSet read_file(const std::string &filename) override;
};