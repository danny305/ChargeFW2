//
// Created by krab1k on 12/11/18.
//

#include "Formal.h"

std::vector<double> Formal::calculate_charges(const Molecule &molecule) {
    std::vector<double> res;
    res.reserve(molecule.atoms().size());
    for (const auto &atom: molecule.atoms()) {
        res.push_back(atom.formal_charge());
    }
    return res;
}
