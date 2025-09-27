#include "reading.hpp"

auto split_string (const std::string &s) -> std::vector<std::string>
{
  std::vector<std::string> result;
  std::vector<char> line (s.begin(),s.end());
  int n = -1;
  for (unsigned i = 0; i < line.size(); i++)
  {
    if (std::isgraph(line[i]))
    {
      if (i==0) 
      {
        result.push_back(std::string{line[i]}); 
        n++;
      }
      else if ((i!=0) && ((line[i-1]==' ') || (line[i-1]=='\t')))
      {
        result.push_back(std::string{line[i]});
        n++;
      }
      else
      {
        result[n].append(std::string{line[i]});
      } 
    } 
  }
  return result;
}

auto get_param(std::ifstream &file) -> double
{
  std::string line;
  std::getline (file,line);
  auto parameters = split_string(line);
  return std::stod(parameters[1]);
}

auto read_database(std::vector<double> &Tcr, std::vector<double> &Pcr, 
  std::vector<double> &omega, std::string databasePath, std::string components) -> void
{
  // Tc, Pc, and omega
  auto n_param = 3; 

  // Getting component names
  auto species = split_string(components);
  auto ncomp = species.size();  

  Tcr.resize(ncomp); Pcr.resize(ncomp); omega.resize(ncomp); 
  
  std::ifstream input_file (databasePath);
  if (!input_file.is_open()) 
  {
    std::cout << "Error opening database file: `" << databasePath << "`." << std::endl; 
    exit(0);
  }

  auto i = 0;
  bool test;
  std::string line;
  while(getline (input_file,line))
  { 
    // getting species name
    auto SpeciesName = split_string(line)[2];

    test = false;
    for (auto j = 0; j < ncomp; j++)
    {
      if (!SpeciesName.compare(species[j]))
      { 
        Tcr[j] = get_param(input_file); 
        Pcr[j] = get_param(input_file);
        omega[j] = get_param(input_file);
        i++; 
        test = true;
        break;
      } 
    }

    if (i==ncomp) break;
    if (!test) {for (int k=0;k<n_param;k++) input_file.ignore(100,'\n');}
  }
  input_file.close();

  // Testing if all input species are present in the database
  for (auto j = 0; j < ncomp; j++) 
  {
    if (Tcr[j] == 0. || Pcr[j] == 0. ||omega[j] == 0.)
      std::cout << "\nSpecies index " << j << " is not present in the database" << std::endl;
  }
}
