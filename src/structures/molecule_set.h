//
// Created by krab1k on 24/10/18.
//

#pragma once

#include <vector>
#include <memory>
#include <utility>
#include <tuple>
#include <string>

#include "molecule.h"
#include "../method.h"
#include "../classifier.h"
#include "../parameters.h"

class AtomClassifier;

class BondClassifier;

class MoleculeSet {
    std::vector<std::tuple<std::string, std::string, std::string>> atom_types_{};
    std::vector<std::tuple<std::string, std::string, std::string, std::string>> bond_types_{};
    std::unique_ptr<std::vector<Molecule> > molecules_{nullptr};

    size_t classify_atoms_from_parameters(const Parameters &parameters, bool remove_unclassified=true);

    size_t classify_bonds_from_parameters(const Parameters &parameters, bool remove_unclassified=true);

public:
    explicit MoleculeSet(std::unique_ptr<std::vector<Molecule> > molecules);

    void info() const;

    const std::vector<Molecule> &molecules() const { return *molecules_; }

    void classify_atoms(const AtomClassifier &cls);

    void classify_bonds(const BondClassifier &cls);

    size_t classify_set_from_parameters(const Parameters &parameters, bool remove_unclassified=true);

    std::vector<std::tuple<std::string, std::string, std::string>> atom_types() const { return atom_types_; }

    std::vector<std::tuple<std::string, std::string, std::string, std::string>>
    bond_types() const { return bond_types_; }

    void fulfill_requirements(const std::vector<RequiredFeatures> &features);
};
