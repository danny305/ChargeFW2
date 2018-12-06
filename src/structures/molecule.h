//
// Created by krab1k on 23/10/18.
//

#pragma once

#include <iostream>
#include <utility>
#include <string>
#include <vector>
#include <tuple>
#include <map>
#include <memory>

#include "atom.h"
#include "bond.h"

class Molecule {
    std::string name_;
    std::unique_ptr<std::vector<Atom> > atoms_;
    std::unique_ptr<std::vector<Bond> > bonds_;
    std::vector<char> bond_info_{};

public:
    const std::vector<Atom> &atoms() const { return *atoms_; }

    const std::vector<Bond> &bonds() const { return *bonds_; }

    const std::string &name() const { return name_; }

    bool bonded(const Atom &atom1, const Atom &atom2) const;

    Molecule() = default;

    Molecule(std::string name, std::unique_ptr<std::vector<Atom> > atoms, std::unique_ptr<std::vector<Bond> > bonds,
             const std::map<int, int> &charges);

    friend std::ostream &operator<<(std::ostream &str, const Molecule &molecule);

    friend class MoleculeSet;
};