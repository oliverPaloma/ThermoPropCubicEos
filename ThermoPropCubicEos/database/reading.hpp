#pragma once

#include "../includes.hpp"

auto split_string (const std::string &s) -> std::vector<std::string>;

auto get_param(std::ifstream &file) -> double;

auto read_database(std::vector<double> &Tc, std::vector<double> &Pc, 
  std::vector<double> &omega, std::string databasePath, std::string components) -> void;