//
// Created by krab1k on 28.1.19.
//

#include <string>
#include <map>
#include <sstream>
#include <fstream>
#include <boost/algorithm/string.hpp>

#include "mmcif.h"
#include "config.h"
#include "../periodic_table.h"


MoleculeSet mmCIF::read_file(const std::string &filename) {
    std::ifstream file(filename);
    if (!file) {
        fmt::print(stderr, "Cannot open file: {}\n", filename);
        exit(EXIT_FILE_ERROR);
    }

    auto molecules = std::make_unique<std::vector<Molecule> >();

    std::string line;
    std::string name;
    try {
        auto atoms = std::make_unique<std::vector<Atom>>();
        do {
            std::getline(file, line);
            if (boost::starts_with(line, "_entry.id")) {
                std::stringstream ss(line);
                ss >> name >> name;
            }
        } while (not boost::starts_with(line, "_atom_site."));

        std::map<std::string, size_t> record_positions;
        size_t category_idx = 0;
        do {
            auto it = line.find('.');
            auto category = line.substr(it + 1);
            boost::trim(category);

            record_positions[category] = category_idx;
            category_idx++;
            std::getline(file, line);
        } while (boost::starts_with(line, "_atom_site."));

        size_t idx = 0;
        do {
            std::stringstream ss(line);
            std::vector<std::string> records{std::istream_iterator<std::string>(ss),
                                             std::istream_iterator<std::string>{}};

            auto it = record_positions.find("Cartn_x");
            auto x = std::stod(records[it->second]);

            it = record_positions.find("Cartn_y");
            auto y = std::stod(records[it->second]);

            it = record_positions.find("Cartn_z");
            auto z = std::stod(records[it->second]);

            it = record_positions.find("type_symbol");
            auto element_symbol = records[it->second];
            auto element = PeriodicTable::pte().getElement(element_symbol);

            it = record_positions.find("label_atom_id");
            auto atom_name = records[it->second];

            it = record_positions.find("label_seq_id");
            auto residue_id = std::stoul(records[it->second]);

            it = record_positions.find("label_comp_id");
            auto residue = records[it->second];

            atoms->emplace_back(idx, element, x, y, z, atom_name, residue_id, residue);

            idx++;
            std::getline(file, line);
        } while (boost::starts_with(line, "  ATOM"));

        auto bonds = std::make_unique<std::vector<Bond>>();
        std::map<size_t, int> charges;
        molecules->emplace_back(name, std::move(atoms), std::move(bonds), charges);

    } catch (const std::exception &) {
        fmt::print(stderr, "Invalid PDB/PQR file\n");
        exit(EXIT_FILE_ERROR);
    }

    return MoleculeSet(std::move(molecules));
}
