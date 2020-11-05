#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <vector>
#include <map>

#include "typedef.h"



// Count Matrix of gRNAs
// The is the gRNA lib ID
typedef std::string gRNAlibraryID;
typedef std::map<gRNAlibraryID, std::vector<std::string> > Model;



typedef std::vector<std::string> Tags;// for csv
const char* model_tags[] = {"libID", "geneName"};

Model count_matrix;
Tags csv_tags;


// Parses gRNA information
void read_model(std::string filename)
{

  std::ifstream fp(filename.c_str());
  std::string line;
  if (!fp.good())
  {
    std::cerr << filename <<" :file not found" << std::endl;
  }

  // discard header of input file
  std::getline(fp, line);
  
  // Insert tags as mageck expects them
  csv_tags.insert(csv_tags.begin(), model_tags, model_tags + 2);
  
  // Read lines 
  while (std::getline(fp, line))
  {
    std::vector<std::string> fields;
    split(fields, line, is_any_of(","));

    auto clibID = fields[2];
    if (fields.size() != 4)
    {
      std::cerr <<"Something is really wrong with that line:" << line << std::endl;
      continue;
    }

    // fields.pop_back();

    
    count_matrix[clibID] = std::vector<std::string>();
    // fields[0] gRNA
    count_matrix[clibID].push_back(fields[2]);
    count_matrix[clibID].push_back(fields[1]);
  }
  fp.close();
}


void read_csv(std::string filename)
{

  std::ifstream fp(filename.c_str());
  std::string line;
  
  if (!fp.good())
  {
    std::cerr << filename <<" :file not found" << std::endl;
  }
  
  if (!std::getline(fp, line))
  {
      std::cerr << filename <<" :empty file" << std::endl;
  }
 
 
  std::vector<std::string> header;
  split(header, line, is_any_of(","));
  auto sample_name = header[3];
  std::cout << "Reading " << sample_name << " from " << filename << std::endl;
  std::cout << "Header Format" << line << std::endl;

  csv_tags.push_back(sample_name);
  
  while (std::getline(fp, line))
  {
    std::vector<std::string> fields;
    split(fields, line, is_any_of(","));
    if (fields.size() != 4)
    {
      std::cerr << "Something is wrong with that line:" << line << std::endl;
      continue;
    }
    
    count_matrix[fields[2]].push_back(fields[3]);
  }
  fp.close();
}

void export_csv(std::string path, const Tags& columns, const Model& data)
{
  std::ofstream out(path.c_str());

  out << boost::algorithm::join(columns, ",") << std::endl;
  // assert(!data.empty() && data.front().size() == columns.size() - 1);

  for (auto const& it : data)
  {
    auto values = join(it.second, ",");
    out << values << std::endl;
  }
  
  out.close();
}


int main(int argc, char* argv[])
{

  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0] << " <outfile-csv> <infile1-csv> .. <infile2-csv>" << std::endl;
    std::cerr << "Merges csv multiple files in the format guide,geneID,libID,sample_name" << std::endl;
    return -1;
  }
  
  // Iterate through all the files 
  for ( int i = 2; i < argc; i++){
    if (count_matrix.empty())
    {
      read_model(argv[2]);
    }
    read_csv(argv[i]);
  }
  export_csv(argv[1], csv_tags, count_matrix);
  
  return 0;

}
