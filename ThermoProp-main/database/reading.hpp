#pragma once

#include "../includes.hpp"

auto split_string (const std::string &s) -> std::vector<std::string>;

auto get_param(std::ifstream &file) -> real;

auto read_database(std::vector<real> &Tc, std::vector<real> &Pc, 
  std::vector<real> &omega, std::string databasePath, std::string components) -> void;