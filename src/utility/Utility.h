//
// Created by krab1k on 6.11.18.
//

#pragma once

#include <iostream>
#include <vector>

template<typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec) {
    for (const auto &val: vec) {
        os << val << " ";
    }
    return os;
}
