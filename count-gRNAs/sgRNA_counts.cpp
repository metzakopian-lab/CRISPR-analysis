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

typedef std::vector<std::string> Model;
typedef std::map <std::string, unsigned int> CountTable;
typedef std::vector<std::string> Tags;// for csv
const char* model_tags[] = {"guide","geneID","libID"};


std::vector<Model> library_index;

CountTable counts;
unsigned short guide_length = 0;

void read_library(std::string filename)
{

  std::ifstream fp(filename.c_str());
  std::string line;
  if (!fp.good())
  {
    std::cerr << filename <<" :file not found" << std::endl;
  }
  while (std::getline(fp, line))
  {
    std::vector<std::string> fields;
    split(fields, line, is_any_of(","));
    if (fields.size() != 3)
    {
      std::cerr <<"something is really wrong with that line:" << line << std::endl;
      continue;
    }
    
    library_index.push_back(fields);
    // std::cout << "gRNA" << fields[0] << std::endl;
    // std::cout << join(fields, ",") <<std::endl;
    counts[fields[0]] = 0;
  }
  guide_length = counts.begin()->first.size();
  std::cout <<" Sucessfully read index with " << library_index.size() << " gRNA of " << guide_length  << std::endl;

}

int map_sequence(CountTable& guide_index, std::string sequence, const unsigned short window_length)
{
  int subseq = sequence.size() - window_length;
  assert(subseq >= 0);
  // std::cout << std::endl << "Testing ";
  for (int i = 0; i <= subseq; i++)
  { 
    auto search = guide_index.find(sequence.substr(i, window_length));
    if (search != guide_index.end())
    {
      // increase the count table and exit
      search->second++;
      return 0;
    }
  }
  // std::cout<< "read" << sequence << std::endl;
  return 1;
}


void read_fq(const std::string filename)
{
  int total_reads = 0, unmapped_reads = 0;
  
  std::ifstream fp(filename.c_str());
  fp.ignore(65536, '\n'); // skip first id
  std::string seq;
  while (std::getline(fp, seq) && fp.good() && !fp.eof())
  {
    
    for(int i = 0; i < 3 && fp.good() && !fp.eof(); i++)
    {
      fp.ignore(65536,'\n');
    }
    
    if(map_sequence(counts, seq, guide_length))
    {
      unmapped_reads++;
    }
        
    if(++total_reads % 500000 == 0)
    {
      std::cout << "Processed reads " << total_reads << 
        ", unmapped percentage: " << (100.0 * unmapped_reads) / total_reads << "\r\n";
    }
  }


  std::cout << "Total processed reads " << total_reads << 
        ", unmapped percentage: " << (100.0 * unmapped_reads) / total_reads << std::endl;

}


void export_csv(std::string path, const Tags& columns, const std::vector<Model>& data, const CountTable& guide_index)
{
  std::ofstream out(path.c_str());

  out << boost::algorithm::join(columns, ",") << std::endl;
  assert(!data.empty() && data.front().size() == columns.size() - 1);

  for (auto& it : data)
  {
    auto gRNA_count  = guide_index.find(it[0])->second;
    auto values = join(it, ",");




    out << values <<"," << gRNA_count << std::endl;
    // out << values << "," <<  << std::endl;
    // out << gRNA_count <<std::endl;
  }
  
  out.close();
}





int main(int argc, char* argv[])
{

  if (argc < 5)
  {
    std::cerr << "Usage: " << argv[0] << "<gRNA-csv-file> <case-name> <fq-file> <csv-outfile>" << std::endl;
    return -1;
  }
  
  Tags csv_tags(model_tags, model_tags + 3 );
  csv_tags.push_back(argv[2]);
  read_library(argv[1]);
  read_fq(argv[3]);
  export_csv(argv[4], csv_tags, library_index, counts);
  
  return 0;

}
