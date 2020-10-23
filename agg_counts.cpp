#include <iostream>
#include <fstream>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>

#include <string>
#include <cassert>
#include <vector>
#include <map>

using boost::algorithm::join;
using boost::algorithm::split;
using boost::is_any_of;

typedef std::map<std::string, std::vector<std::string> > Model;

typedef std::vector<std::string> Tags;// for csv
const char* model_tags[] = {"libID", "gRNA", "geneID"};

Model count_matrix;
Tags csv_tags;


// Reads the first file and peeks the gRNAs
void read_model(std::string filename)
{

  std::ifstream fp(filename.c_str());
  std::string line;
  if (!fp.good())
  {
    std::cerr << filename <<" :file not found" << std::endl;
  }
  std::getline(fp, line);

  while (std::getline(fp, line))
  {
    std::vector<std::string> fields;
    split(fields, line, is_any_of(","));
    if (fields.size() != 4)
    {
      std::cerr <<"Something is really wrong with that line:" << line << std::endl;
      continue;
    }

    // fields.pop_back();
 
    count_matrix[fields[2]] = std::vector<std::string>();
    count_matrix[fields[2]].push_back(fields[0]);
    count_matrix[fields[2]].push_back(fields[1]);
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
  
  if(!std::getline(fp, line))
  {
      std::cerr << filename <<" :empty file" << std::endl;
  }
  else
  {
      std::vector<std::string> header;
      split(header, line, is_any_of(","));
      std::cout << "Reading " << header[3] << std::endl;
      csv_tags.push_back(header[3]);
  }
  
  while (std::getline(fp, line))
  {
    std::vector<std::string> fields;
    split(fields, line, is_any_of(","));
    if (fields.size() != 4)
    {
      std::cerr << "Something is really wrong with that line:" << line << std::endl;
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
    out << it.first << ","<<values << std::endl;
  }
  
  out.close();
}


int main(int argc, char* argv[])
{

  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0] << "<outfile-csv> <infile1-csv> .. <infile2-csv>" << std::endl;
    std::cerr << "Merges csv multiple files in the format guide,geneID,libID,sample_name" << std::endl;
    return -1;
  }
  
  csv_tags.insert(csv_tags.begin(), model_tags, model_tags + 3 );
  read_model(argv[2]);
  for ( int i = 2; i < argc; i++){
          read_csv(argv[i]);
  }
  export_csv(argv[1], csv_tags, count_matrix);
  
  return 0;

}
