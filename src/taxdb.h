/* 
 * Original work Copyright 2013 David Ainsworth
 * Modified work copyright 2017 Florian Breitwieser 
 *
 * The original file is part of SLAM
 * The modified file is part of KrakenUniq
 *
 * SLAM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SLAM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.

 * You should have received a copy of the GNU Affero General Public License
 * along with SLAM.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef TAXD_DB_H_
#define TAXD_DB_H_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <unordered_set>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include "report-cols.h"
//#include "readcounts.hpp"


using namespace std;
//using kraken::ReadCounts;

namespace patch
{
  template < typename T > std::string to_string( const T& n )
  {
    std::ostringstream stm ;
    stm << n ;
    return stm.str() ;
  }
}


void log_msg (const std::string& s);

template<typename T> uint64_t string_to_T(std::string str);

std::vector<std::string> in_betweens(const std::string &s, const char start_char, const char end_char, size_t start_at = 0);

std::vector<std::string> tokenise(const std::string &s, const std::string& delimiter, size_t max_fields = 0, size_t end_chars = 0);

std::vector<std::string> get_fields(const std::string &s, const std::string& delimiter, std::vector<size_t> fields); 

// TODO: Consider using TaxRank instead of string in TaxonomyEntry
//       However, then it would not be possible to define custom ranks..
struct TaxRank {
  // All ranks that appear in the NCBI taxonomy database,
  //  plus 'sequence', 'assembly', and 'root'
  //static constexpr vector<string> rank_strings = {
  // "no rank", "sequence", "assembly",
  // "subspecies", "species", "subgenus", "genus", "tribe", "subfamily",
  //"family", "superfamily", "parvorder", "infraorder", "suborder",
  //"order", "superorder", "parvclass", "infraclass", "subclass",
  //"class", "superclass", "subphylum", "phylum", "kingdom",
  //"superkingdom", "root"};

  enum RANK { unknown, no_rank, sequence, assembly,
    subspecies, species, species_subgroup, species_group, subgenus, genus, tribe, subfamily,
    family, superfamily, parvorder, infraorder, suborder,
    order, superorder, parvclass, infraclass, subclass,
    class_, superclass, subphylum, phylum, kingdom,
    superkingdom, root
  };

  static const unordered_map<string, RANK> string_to_rank;

  static RANK toRank(const string& rank) {
    const auto& it = string_to_rank.find(rank);
    if (it == string_to_rank.end()) {
      cerr << "ERROR: Could not find rank " << rank << endl;
      exit(1);
    }
    return it->second;
  }

  static const char* toString(const TaxRank::RANK& rank) {
    switch(rank) {
      case RANK::unknown:           return "unknown";
      case RANK::no_rank:          return "no rank";
      case RANK::sequence:         return "sequence";
      case RANK::assembly:         return "assembly";
      case RANK::subspecies:       return "subspecies";
      case RANK::species:          return "species";
      case RANK::species_subgroup: return "species subgroup";
      case RANK::species_group:    return "species group";
      case RANK::subgenus:         return "subgenus";
      case RANK::genus:            return "genus";
      case RANK::tribe:            return "tribe";
      case RANK::subfamily:        return "subfamily";
      case RANK::family:           return "family";
      case RANK::superfamily:      return "superfamily";
      case RANK::parvorder:        return "parvorder";
      case RANK::infraorder:       return "infraorder";
      case RANK::suborder:         return "suborder";
      case RANK::order:            return "order";
      case RANK::superorder:       return "superorder";
      case RANK::parvclass:        return "parvclass";
      case RANK::infraclass:       return "infraclass";
      case RANK::subclass:         return "subclass";
      case RANK::class_:            return "class";
      case RANK::superclass:       return "superclass";
      case RANK::subphylum:        return "subphylum";
      case RANK::phylum:           return "phylum";
      case RANK::kingdom:          return "kingdom";
      case RANK::superkingdom:     return "superkingdom";
      case RANK::root:             return "root";
      default:
           log_msg("Invalid rank!\n");
    }
    return "NA";
  }

};

const unordered_map<string, TaxRank::RANK> TaxRank::string_to_rank = {
  {"unknown", TaxRank::unknown},
  {"no rank", TaxRank::no_rank}, 
  {"sequence", TaxRank::sequence},
  {"assembly", TaxRank::assembly},
  {"subspecies", TaxRank::subspecies},
  {"species", TaxRank::species},
  {"species subgroup", TaxRank::species_subgroup},
  {"species group", TaxRank::species_group},
  {"subgenus", TaxRank::subgenus},
  {"genus", TaxRank::genus},
  {"tribe", TaxRank::tribe},
  {"subfamily", TaxRank::subfamily},
  {"family", TaxRank::family},
  {"superfamily", TaxRank::superfamily},
  {"parvorder", TaxRank::parvorder},
  {"infraorder", TaxRank::infraorder},
  {"suborder", TaxRank::suborder},
  {"order", TaxRank::order},
  {"superorder", TaxRank::superorder},
  {"parvclass", TaxRank::parvclass},
  {"infraclass", TaxRank::infraclass},
  {"subclass", TaxRank::subclass},
  {"class", TaxRank::class_},
  {"superclass", TaxRank::superclass},
  {"subphylum", TaxRank::subphylum},
  {"phylum", TaxRank::phylum},
  {"kingdom", TaxRank::kingdom},
  {"superkingdom", TaxRank::superkingdom},
  {"root", TaxRank::root}
};


template<typename TAXID>
class TaxonomyEntry {
  public:
    TAXID taxonomyID;
    TaxonomyEntry<TAXID>* parent;
    std::vector<TaxonomyEntry*> children;

    string rank;
    std::string scientificName;
    uint64_t genomeSize;
    uint64_t genomeSizeOfChildren;

    TaxonomyEntry() : taxonomyID(0), parent(NULL), genomeSize(0), genomeSizeOfChildren(0) {}

    TaxonomyEntry(TAXID taxonomyID_, TaxonomyEntry<TAXID>* parent_, std::string rank_, std::string scientificName_, uint64_t genomeSize_ = 0, uint64_t genomeSizeOfChildren_ = 0) :
      taxonomyID(taxonomyID_), parent(parent_), rank(rank_), scientificName(scientificName_),
      genomeSize(genomeSize_), genomeSizeOfChildren(genomeSizeOfChildren_) {

  if (parent_ != NULL) {
    parent->children.push_back(this);
  }

      }

    inline bool operator==(const TaxonomyEntry& other) const; 

    friend std::ostream &operator<<(std::ostream &os, const TaxonomyEntry<TAXID> &m) { 
      TAXID parentTaxonomyID = (m.parent == NULL)? m.taxonomyID : m.parent->taxonomyID;
      os << '[' << m.taxonomyID << ";parent="<< parentTaxonomyID << ";name=" << m.scientificName << ";rank=" << m.rank << ']';
      return os;
    }

};




//template<>
//TaxonomyEntry<uint32_t, uint64_t>::TaxonomyEntry () {
//  readCounts = 0;
//  readCountsOfChildren = 0;
//}

/*
   template<typename TAXID>
   struct TaxonomyEntryPtr_comp {
   bool operator() ( const TaxonomyEntry<TAXID>* a, const TaxonomyEntry<TAXID>* b) const;
   };
   */


template<typename TAXID>
class TaxonomyDB {
  public:
    TaxonomyDB(const std::string namesDumpFileName, const std::string nodesDumpFileName);
    TaxonomyDB(const std::string inFileName, bool hasGenomeSizes = false);
    TaxonomyDB();

    TaxonomyDB(TaxonomyDB&& rhs) : entries(std::move(rhs.entries)) {
    }

    TaxonomyDB& operator=(TaxonomyDB&& rhs) {
      entries = std::move(rhs.entries);
      return *this;
    }


    void writeTaxonomyIndex(std::ostream & outs) const;
    void readTaxonomyIndex(const std::string inFileName, bool hasGenomeSizes);

    TAXID getTaxIDAtRank(const TAXID taxID, const std::string& rank) const;
    std::string getScientificName(const TAXID taxID) const;
    std::string getRank(const TAXID taxID) const;
    TAXID getLowestCommonAncestor(const std::vector<TAXID>& taxIDs) const;
    pair<TAXID,int> getLowestCommonAncestor(TAXID a, TAXID b) const;
    string getNextProperRank(TAXID a) const;
    TAXID getTaxIDAtNextProperRank(TAXID a) const;

    TAXID getParentTaxID(const TAXID taxID) const;
    std::unordered_map<TAXID, TAXID> getParentMap() const;
    TAXID getByScientificName(string name) const;
    std::unordered_map<std::string, TAXID> getScientificNameMap() const;
    std::string getLineage(TAXID taxonomyID) const;
    std::string getMetaPhlAnLineage(TAXID taxonomyID) const;
    TaxonomyEntry<TAXID> getEntry(TAXID taxID) const;

    bool insert(TAXID taxonomyID_, TAXID parentTaxonomyID_, std::string rank_, std::string scientificName_);
    bool hasTaxon(TAXID taxonomyID_);

    size_t distance(TAXID taxID1, TAXID taxID2) const;

    bool isSubSpecies(TAXID taxonomyID) const;
    int isBelowInTree(TAXID upper, TAXID lower) const;

    void setGenomeSizes(const std::unordered_map<TAXID, uint64_t> & genomeSizes);
    void readGenomeSizes(string file);
    void setGenomeSize(const TAXID taxid, const uint64_t genomeSize);

    void printReport();

    std::unordered_map<TAXID, TaxonomyEntry<TAXID> > entries;
    bool genomeSizes_are_set;

  private:

    std::unordered_map<TAXID, TaxonomyEntry<TAXID> >
      readTaxonomyIndex_(const std::string inFileName, bool hasGenomeSizes);
};


template<typename TAXID, typename READCOUNTS>
class TaxReport {
  private:
    std::ostream& _reportOfb;
    const TaxonomyDB<TAXID> & _taxdb;
    const std::unordered_map<TAXID, READCOUNTS>& _taxCounts; // set in constructor, from classification
    //std::unordered_map<const TaxonomyEntry<TAXID>*, vector<const TaxonomyEntry<TAXID>*> > _children;
    std::unordered_map<const TaxonomyEntry<TAXID>*, vector<const READCOUNTS*> > _children;
    std::unordered_map<const TaxonomyEntry<TAXID>*, READCOUNTS> _cladeCounts; // consider accessing by TaxEntry*
    uint64_t _total_n_reads = 0;
    bool _show_zeros;
    void printLine(const TaxonomyEntry<TAXID>& tax, const READCOUNTS& rc, unsigned depth);
    READCOUNTS setCladeCounts(const TaxonomyEntry<TAXID>* tax, unordered_map<const TaxonomyEntry<TAXID>*, unordered_set<const TaxonomyEntry<TAXID>*> >& _children);

  public:
    TaxReport(std::ostream& _reportOfb, const TaxonomyDB<TAXID> & taxdb, const std::unordered_map<TAXID, READCOUNTS>&, bool _show_zeros);
    void printReport(const std::string & format);
    void printReport(const TaxonomyEntry<TAXID>& tax, unsigned depth);
    void setReportCols(const std::vector<std::string>& names);

    std::vector<std::string> _report_col_names;
    std::vector<REPORTCOLS> _report_cols;
};

template<typename K,typename V>
inline
V find_or_use_default(const std::unordered_map<K, V>& my_map, const K& query, const V default_value);

//////////////////////////// DEFINITIONS
void log_msg (const std::string& s) {
  std::cerr << s;
}

template<typename T>
uint64_t string_to_T(string str) {
  stringstream stream(str);
  T result;
  stream >> result;
  return result;
}

std::vector<std::string> in_betweens(const std::string &s, const char start_char, const char end_char, size_t start_at) {
  std::vector<std::string> tokens;
  size_t i = 0;
  size_t next_end = start_at-1;

  for (size_t next_start = s.find(start_char, next_end + 1); \
      next_start != string::npos;
      next_start = s.find(start_char, next_end + 1), ++i) {

    next_end = s.find(end_char, next_start + 1);
    if (next_end == string::npos) {
      cerr << "unmatched start and end!";
      exit(1);
    }

    tokens.push_back(s.substr(next_start+1, next_end-1));
  }

  return tokens;
}



std::vector<std::string> tokenise(const std::string &s, const std::string& delimiter, size_t max_fields, size_t end_chars) {
  std::vector<std::string> tokens(max_fields);
  size_t delim_length = delimiter.length();
  size_t last = 0;
  size_t i = 0;

  for (size_t next = s.find(delimiter, last);
      (max_fields > 0 && i < max_fields) && next != string::npos;
      next = s.find(delimiter, last), ++i) {
    tokens[i] = s.substr(last, next-last);
    last = next + delim_length;
  }
  if (max_fields > 0 && i < max_fields) {
    tokens[max_fields-1] = s.substr(last, s.length()-last-end_chars);
  }

  return tokens;
}

std::vector<std::string> get_fields(const std::string &s, const std::string& delimiter, vector<size_t> fields) {
  std::vector<std::string> tokens;
  tokens.reserve(fields.size());
  size_t delim_length = delimiter.length();
  size_t last = 0;
  size_t i = 0;
  size_t current_field = 0;

  for (size_t next = s.find(delimiter, last);
      tokens.size() < fields.size() && next != string::npos;
      next = s.find(delimiter, last), ++i) {
    if (i == fields[current_field]) {
      tokens.push_back(s.substr(last, next-last));
      ++current_field;
    }
    last = next + delim_length;
  }

  return tokens;
}


template<typename TAXID>
unordered_map<TAXID, TAXID> TaxonomyDB<TAXID>::getParentMap() const {
  unordered_map<TAXID, TAXID> Parent_map;
  //for (const auto & tax : entries) {
  for (auto tax_it = entries.begin(); tax_it != entries.end(); ++tax_it) {
    if (tax_it->first == 0) 
      continue;
    if (tax_it->second.parent == NULL) {
      //cerr << "Parent for " << tax.first << " is 0\n";
      Parent_map[tax_it->first] = 0; // for kraken::lca
    } else {
      Parent_map[tax_it->first] = tax_it->second.parent->taxonomyID;
    }
  }
  return Parent_map;
}

template<typename TAXID>
TaxonomyEntry<TAXID> TaxonomyDB<TAXID>::getEntry(TAXID taxID) const {
  auto it = entries.find(taxID);
  if (it == entries.end()) {
    TaxonomyEntry<TAXID> ti { 0, 0, "NA"};
    return ti;
  } else {
    return it->second;
  }
}

template<typename TAXID>
void createPointers(
    std::unordered_map<TAXID, TaxonomyEntry<TAXID> >& entries, 
    const std::unordered_map<TAXID, TAXID>& parentMap) {
  for (auto entry_it = entries.begin(); entry_it != entries.end(); ++entry_it) {
    TAXID taxonomyID = entry_it->first;
    auto parent_it = parentMap.find(taxonomyID);
    if (parent_it == parentMap.end()) {
      //cerr << "Cannot find parent for " << taxonomyID << endl;
    } else {
      TAXID parentTaxonomyID = parent_it->second;
      if (taxonomyID != parentTaxonomyID) {
        auto parent_ptr = entries.find(parentTaxonomyID);
        if (parent_ptr != entries.end()) {
          entry_it->second.parent = &parent_ptr->second;
          parent_ptr->second.children.push_back(&entry_it->second);
       } else {
         cerr << "Could not find parent with taxonomy ID " << parentTaxonomyID << " for taxonomy ID " << taxonomyID << endl;
       }
      }
    }
  }
}

template<typename TAXID>
TaxonomyDB<TAXID>::TaxonomyDB() : genomeSizes_are_set(false) { }

template<typename TAXID>
TaxonomyDB<TAXID>::TaxonomyDB(const std::string inFileName, bool hasGenomeSizes) :
  entries( readTaxonomyIndex_(inFileName, hasGenomeSizes) ), genomeSizes_are_set(hasGenomeSizes)
{ }

template<typename TAXID>
unordered_map<TAXID, TaxonomyEntry<TAXID>> readDumps(const std::string namesDumpFileName, const std::string nodesDumpFileName) {
  std::unordered_map<TAXID, TaxonomyEntry<TAXID> > entries;
  log_msg("Building taxonomy index from " + nodesDumpFileName + " and " + namesDumpFileName);
  unordered_map<TAXID, TAXID> parentMap = parseNodesDump(nodesDumpFileName, entries);
  createPointers(entries, parentMap);
  parseNamesDump(namesDumpFileName, entries);
  log_msg(". Done, got " + patch::to_string(entries.size()) + " taxa\n");
  return(entries);
}

template<typename TAXID>
TaxonomyDB<TAXID>::TaxonomyDB(const std::string namesDumpFileName, const std::string nodesDumpFileName) : 
  entries(readDumps<TAXID>(namesDumpFileName, nodesDumpFileName)) {
  }

template<typename TAXID>
std::unordered_map<TAXID,TAXID> parseNodesDump(const std::string nodesDumpFileName, std::unordered_map<TAXID, TaxonomyEntry<TAXID> >& entries) {
  std::ifstream nodesDumpFile(nodesDumpFileName);
  if (!nodesDumpFile.is_open())
    throw std::runtime_error("unable to open nodes file");

  std::string line;

  TAXID taxonomyID;
  TAXID parentTaxonomyID;
  std::string rank;
  char delim;
  std::unordered_map<TAXID,TAXID> parentMap;

  while (nodesDumpFile >> taxonomyID >> delim >> parentTaxonomyID) {
    nodesDumpFile.ignore(3);
    getline(nodesDumpFile, rank, '\t');
    auto entryIt = entries.find(taxonomyID);
    if (entryIt == entries.end()) {
      entries[taxonomyID] = TaxonomyEntry<TAXID>(taxonomyID, NULL, rank, "");
      parentMap[taxonomyID] = parentTaxonomyID;
    } else {
      parentMap[taxonomyID] = parentTaxonomyID;
      entryIt->second.rank = rank;
    }

    nodesDumpFile.ignore(2560, '\n');
  }
  return parentMap;
}

template<typename TAXID>
void parseNamesDump(const std::string namesDumpFileName, std::unordered_map<TAXID, TaxonomyEntry<TAXID> >& entries) {
  std::ifstream namesDumpFile(namesDumpFileName);
  if (!namesDumpFile.is_open())
    throw std::runtime_error("unable to open names file");
  std::string line;

  TAXID taxonomyID;
  std::string scientificName, type;
  while (namesDumpFile.good()) {
    namesDumpFile >> taxonomyID;
    namesDumpFile.ignore(3);
    getline(namesDumpFile, scientificName, '\t');
    namesDumpFile.ignore(3);
    namesDumpFile.ignore(256, '|');
    namesDumpFile.ignore(1);
    getline(namesDumpFile, type, '\t');

    if (type == "scientific name") {
      auto entryIt = entries.find(taxonomyID);
      if (entryIt == entries.end()) {
  cerr << "Entry for " << taxonomyID << " does not exist - it should!" << '\n';
  //entries[taxonomyID] = TaxonomyEntry<TAXID>(taxonomyID, NULL, "", scientificName);
      } else {
  entryIt->second.scientificName = scientificName;
      }
    }
    namesDumpFile.ignore(2560, '\n');
  }
}

template<typename KeyType, typename ValueType>
std::vector<KeyType> getSortedKeys(const std::unordered_map<KeyType, ValueType>& my_unordered_map) {
  std::vector<KeyType> keys;
  keys.reserve (my_unordered_map.size());
  for (auto it = my_unordered_map.begin(); it != my_unordered_map.end(); ++it) {
    keys.push_back(it->first);
  }
  std::sort (keys.begin(), keys.end());
  return keys;
}

template<typename TAXID>
void TaxonomyDB<TAXID>::writeTaxonomyIndex(std::ostream & outs) const {
  std::vector<TAXID> sorted_keys = getSortedKeys(entries);
  for (size_t i = 0; i < sorted_keys.size(); ++i) {
    TAXID taxonomyID = sorted_keys[i];
    const auto& entry = entries.at(taxonomyID);
    TAXID parentTaxonomyID = (entry.parent==NULL? taxonomyID : entry.parent->taxonomyID);
    outs << taxonomyID << '\t' << parentTaxonomyID << '\t'
      << entry.scientificName << '\t' << entry.rank;
    if (genomeSizes_are_set) {
      outs << '\t' << entry.genomeSize << '\t' << entry.genomeSizeOfChildren;
    }
    outs << '\n';
  }
  outs.flush();
}

template<typename TAXID>
void TaxonomyDB<TAXID>::setGenomeSizes(const std::unordered_map<TAXID, uint64_t> & genomeSizes) {
  for (auto it = genomeSizes.begin(); it != genomeSizes.end(); ++it) {
    setGenomeSize(it->first, it->second);
  }
  genomeSizes_are_set = true;
}

template<typename TAXID>
void TaxonomyDB<TAXID>::readTaxonomyIndex(const std::string inFileName, bool hasGenomeSizes) {
  entries = readTaxonomyIndex_(inFileName, hasGenomeSizes);
  genomeSizes_are_set = hasGenomeSizes;
}

template<typename TAXID>
std::unordered_map<TAXID, TaxonomyEntry<TAXID> > 
TaxonomyDB<TAXID>::readTaxonomyIndex_(const std::string inFileName, bool hasGenomeSizes) {
  log_msg("Reading taxonomy index from " + inFileName);
  std::ifstream inFile(inFileName);
  if (!inFile.is_open())
    throw std::runtime_error("unable to open taxonomy index file " + inFileName);

  std::unordered_map<TAXID, TaxonomyEntry<TAXID> > entries;
  std::unordered_map<TAXID, TAXID> parentMap;
  TAXID taxonomyID, parentTaxonomyID;
  std::string scientificName, rank;
  uint64_t genomeSize = 0;
  uint64_t genomeSizeOfChildren = 0;

  std::string line;
  while (!inFile.eof()) {
    inFile >> taxonomyID >> parentTaxonomyID;
    if (taxonomyID > 1 && taxonomyID == parentTaxonomyID) {
      cerr << "ERROR: the parent of " << taxonomyID << " is itself. Should not happend for taxa other than the root.\n";
      exit(1);
    }
    inFile.get(); // read tab
    std::getline(inFile, scientificName, '\t');
    if (hasGenomeSizes) {
      std::getline(inFile, rank, '\t');
      inFile >> genomeSize >> genomeSizeOfChildren;
    } else {
      std::getline(inFile, rank, '\n');
    }
    TaxonomyEntry<TAXID> newEntry(taxonomyID, NULL, rank, scientificName, genomeSize, genomeSizeOfChildren);

    //auto insert_res = entries.insert({ taxonomyID, newEntry });
    entries.insert({ taxonomyID, newEntry });
    parentMap[taxonomyID] = parentTaxonomyID;
  }
  entries.insert({0, {0, NULL, "no rank", "unclassified" }});
  //entries.insert({-1, {-1, 0, "no rank", "uncategorized" }});
  createPointers(entries, parentMap);
  log_msg(". Done.\n");
  //log_msg(". Done, read " + patch::to_string(entries.size()) + " taxa.\n");
  return(entries);
}

template<typename TAXID>
string TaxonomyDB<TAXID>::getNextProperRank(TAXID a) const {
  if (a == 0) {
    return "NA";
  }
  while (getRank(a) == "no rank" && a != getParentTaxID(a)) {
    a = getParentTaxID(a);
  }
  if ( a == 1 ) {
    return "root";
  }
  return getRank(a);
}

template<typename TAXID>
TAXID TaxonomyDB<TAXID>::getTaxIDAtNextProperRank(TAXID a) const {
  if (a == 0 || a == 1) {
    return 0;
  }
  while (getRank(a) == "no rank" && a != getParentTaxID(a)) {
    a = getParentTaxID(a);
  }
  return a;
}

template<typename TAXID>
pair<TAXID,int> TaxonomyDB<TAXID>::getLowestCommonAncestor(TAXID a, TAXID b) const {
  if (a == 0 || b == 0) {
    return a ? pair<TAXID,int>(a,-1) : pair<TAXID,int>(b,-1); 
  }

  // create a path from a to the root
  std::unordered_set<uint32_t> a_path;
  int distA = 0;
  while (a > 0 && a != getParentTaxID(a)) {
    if (a == b)
      return pair<TAXID,int>{a, distA};
    a_path.insert(a);
    a = getParentTaxID(a);
    ++distA;
  }

  int distB = 0;
  // search for b in the path from a to the root
  while (b > 0 && b != getParentTaxID(b)) {
    auto it = a_path.find(b);
    if (it != a_path.end()) {
      return pair<TAXID,int>(b, distB + std::distance(a_path.begin(), it));
    }
    b = getParentTaxID(b);
    ++distB;
  }
  return pair<TAXID,int>(1, distA+distB);
}

/*

   template<typename TAXID>
   TAXID TaxonomyDB<TAXID>::getLowestCommonAncestor(
   const std::vector<TAXID>& taxIDs) const {
   if (taxIDs.size() == 0) {
   return 0;
   }
   std::vector<std::vector<READCOUNTS> > paths;
   for (auto& taxID : taxIDs) {
   bool good = true;
   std::vector<READCOUNTS> path;
   TAXID tempTaxID = taxID;
   while (tempTaxID != 0) {
   path.push_back(tempTaxID);
   tempTaxID = getParentTaxID(tempTaxID);
   }
   if (good) paths.push_back(path);
   }
   if (paths.size() == 0) {
   return 0;
   }
   for (auto& path : paths)
   std::reverse(path.begin(), path.end());
   std::sort(paths.begin(), paths.end(),
   [](std::vector<READCOUNTS> i, std::vector<READCOUNTS> j) {
   return i.size() < j.size();
   });
   TAXID consensus = 0;
// assumes equal paths lengths??
for (unsigned i = 0; i < paths[0].size(); i++) {
TAXID temp = 0;
for (auto& path : paths) {
if (temp == 0)
temp = path[i];
else if (temp != path[i]) {
return consensus;
}
}
consensus = temp;
}
return consensus;
}
*/

template<typename TAXID>
bool TaxonomyDB<TAXID>::hasTaxon(TAXID taxonomyID_) {
  return entries.find(taxonomyID_) != entries.end();
}

template<typename TAXID>
bool TaxonomyDB<TAXID>::insert(TAXID taxonomyID_, TAXID parentTaxonomyID_, 
    std::string rank_, std::string scientificName_) {
  if (parentTaxonomyID_ == taxonomyID_) {
    return false;
  }

  auto parentIt = entries.find(parentTaxonomyID_);
  if (parentIt == entries.end()) {
    cerr << "ERROR with taxon [" << taxonomyID_  <<";"<<rank_<<";"<<scientificName_<<"] - parent taxon " << parentTaxonomyID_ << " not in database!" << endl;
    return false;
  }

  TaxonomyEntry<TAXID> newEntry(taxonomyID_, &parentIt->second, rank_, scientificName_, 0, 0);

  newEntry.parent = &(parentIt->second);
  auto insert_res = entries.insert({taxonomyID_, newEntry});
  if (insert_res.second) {
    parentIt->second.children.push_back(&insert_res.first->second);
  }
  return insert_res.second;

}

template<typename TAXID>
TAXID TaxonomyDB<TAXID>::getParentTaxID(const TAXID taxID) const {
  auto entry = entries.find(taxID);
  if (entry != entries.end() && entry->second.parent != NULL)
    return entry->second.parent->taxonomyID;
  else
    return 0;
}

template<typename TAXID>
std::string TaxonomyDB<TAXID>::getScientificName(const TAXID taxID) const {
  auto entry = entries.find(taxID);
  if (entry != entries.end()) {
    return entry->second.scientificName;
  } else
    return std::string();
}

template<typename TAXID>
std::string TaxonomyDB<TAXID>::getRank(const TAXID taxID) const {
  auto entry = entries.find(taxID);
  if (entry != entries.end()) {
    return entry->second.rank;
  } else
    return std::string();
}

template<typename TAXID>
std::string TaxonomyDB<TAXID>::getLineage(TAXID taxonomyID) const {
  std::string lineage;
  while (true) {
    // 131567 = Cellular organisms
    if (taxonomyID != 131567) {
      if (lineage.size()) lineage.insert(0, "; ");
      lineage.insert(0, getScientificName(taxonomyID));
      if (getRank(taxonomyID) == "species") lineage.clear();
    }
    taxonomyID = getParentTaxID(taxonomyID);
    if (taxonomyID == 0) {
      if (lineage.size()) lineage.append(".");
      break;
    }
  }
  return lineage;
}

template<typename TAXID>
std::string TaxonomyDB<TAXID>::getMetaPhlAnLineage(TAXID taxonomyID) const {
  std::string rank = getRank(taxonomyID);
  if (rank == "superphylum") return std::string();
  std::string lineage;
  while (true) {
    // 131567 = Cellular organisms
    //if (taxonomyID != 131567) {
      std::string rank = getRank(taxonomyID);
      if (rank == "species") {
  lineage.insert(0, "|s__");
  lineage.insert(4, getScientificName(taxonomyID));
      } else if (rank == "genus") {
  lineage.insert(0, "|g__");
  lineage.insert(4, getScientificName(taxonomyID));
      } else if (rank == "family") {
  lineage.insert(0, "|f__");
  lineage.insert(4, getScientificName(taxonomyID));
      } else if (rank == "order") {
  lineage.insert(0, "|o__");
  lineage.insert(4, getScientificName(taxonomyID));
      } else if (rank == "class") {
  lineage.insert(0, "|c__");
  lineage.insert(4, getScientificName(taxonomyID));
      } else if (rank == "phylum") {
  lineage.insert(0, "|p__");
  lineage.insert(4, getScientificName(taxonomyID));
      } else if (rank == "superkingdom") {
  lineage.insert(0, "|k__");
  lineage.insert(4, getScientificName(taxonomyID));
      } else {
  lineage.insert(0, "|-__");
  lineage.insert(4, getScientificName(taxonomyID));

	 // }
    }
    taxonomyID = getParentTaxID(taxonomyID);
    if (taxonomyID == 0) {
      break;
    }
  }
  std::replace(lineage.begin(), lineage.end(), ' ', '_');
  return lineage;
}

template<typename TAXID>
TAXID TaxonomyDB<TAXID>::getTaxIDAtRank(const TAXID taxID,
    const std::string& rank) const {
  if (taxID == 0 || taxID == 1)
    return 0;
  auto entry_it = entries.find(taxID);
  // cerr << "getTaxIDAtRank(" << taxID << "," << rank << ")" << endl;
  if (entry_it != entries.end()) {
    const TaxonomyEntry<TAXID>* entry_ptr = &entry_it->second;
    while (entry_ptr != NULL
  && entry_ptr->parent != NULL) {
      // cerr << "Checking rank of " << entry->second.taxonomyID << ": " << entry->second.rank << endl;
      if (entry_ptr->rank == rank) {
  return entry_ptr->taxonomyID;
      } else {
  entry_ptr = entry_ptr->parent;
      }
    }
  }
  return 0;
}


template<typename TAXID>
void TaxonomyDB<TAXID>::setGenomeSize(const TAXID taxid, const uint64_t genomeSize) {
  auto it = entries.find(taxid);
  if (it == entries.end()) {
    cerr << "No taxonomy entry for " << taxid << "!!" << endl;
    return;
  }
  TaxonomyEntry<TAXID>* tax = &it->second;
  tax->genomeSize += genomeSize;

  while (tax->parent != NULL) {
    tax = tax->parent;
    //std::cerr << "setting genomeSizeOfChildren of parent" << std::endl;
    tax->genomeSizeOfChildren += genomeSize;
  }
}

template<typename TAXID>
void TaxonomyDB<TAXID>::readGenomeSizes(string file) {
  //for (auto entry_it = entries.begin(); entry_it != entries.end(); ++entry_it) {
  //  entry_it->second.genomeSize = 0;
  //  entry_it->second.genomeSizeOfChildren = 0;
  //}
  cerr << "Reading genome sizes from " << file << " ...";
  std::ifstream inFile(file);
  if (!inFile.is_open())
    throw std::runtime_error("unable to open file " + file);
  TAXID taxonomyID;
  uint64_t size;
  while (!inFile.eof()) {
    inFile >> taxonomyID >> size;
    setGenomeSize(taxonomyID, size);
  }

  cerr << " done" << endl;
}

/*
   template<typename TAXID>
   void TaxonomyDB<TAXID>::setReadCounts(const unordered_map<TAXID>& readCounts) {
   for (auto& elem : readCounts) {
   addReadCount(elem.first, elem.second);
   }

   for (auto& tax : entries) {
   std::sort(tax.second.children.begin(), tax.second.children.end(),TaxonomyEntryPtr_comp<TAXID>());
   }
   }
   */

/*
template<typename TAXID, typename READCOUNTS>
READCOUNTS TaxReport<TAXID,READCOUNTS>::setCladeCounts(const TaxonomyEntry<TAXID>* tax, unordered_map<const TaxonomyEntry<TAXID>*, unordered_set<const TaxonomyEntry<TAXID>*> >& _children) {
  auto itt = _taxCounts.find(tax->taxonomyID);
  auto itc = _children.find(tax);
  if (itc == _children.end()) {
    //this is a leaf node, return its taxCounts
    if (itt == _taxCounts.end()) {
      // cerr << "This leaf node [taxid "<< tax->taxonomyID <<"] has no taxon count" << endl;
      return(READCOUNTS());
    }
    _cladeCounts[tax] = itt->second;
  } else {
    auto c = itc->second.begin();

      READCOUNTS rc = setCladeCounts(*c, _children);

    for (++c; c != itc->second.end(); ++c) {
      rc += setCladeCounts(*c, _children);
    }
    if (itt != _taxCounts.end())
      rc += itt->second;

    _cladeCounts[tax] = rc;
  }
  return _cladeCounts.at(tax);
}*/

template<typename TAXID, typename READCOUNTS>
TaxReport<TAXID,READCOUNTS>::TaxReport(std::ostream& reportOfb, const TaxonomyDB<TAXID>& taxdb, 
    const std::unordered_map<TAXID, READCOUNTS>& readCounts,
    bool show_zeros) : _reportOfb(reportOfb), _taxdb(taxdb), _taxCounts(readCounts), _show_zeros(show_zeros) {

  cerr << "Setting values in the taxonomy tree ...";
  //unordered_map<const TaxonomyEntry<TAXID>*, unordered_set<const TaxonomyEntry<TAXID>*> > _children1;
  for (auto it = _taxCounts.begin(); it != _taxCounts.end(); ++it) {
    auto tax_it = taxdb.entries.find(it->first);
    if (tax_it == taxdb.entries.end()) {
      cerr << "No entry for " << it->first << " in database!" << endl;
    } else {
      const TaxonomyEntry<TAXID>* tax = &(tax_it->second);
      while (tax != NULL) {
        _children[tax].push_back(&(it->second));
        /*
        auto res = (_children1[tax->parent]).insert(tax);
        if (!res.second)
          break;
        */
        tax = tax->parent;
      }
    }
  }

  //cerr << " Nr children: " << _children.size() << endl;


#ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic, 50)
#endif
  for (size_t i = 0; i < _children.size(); ++i) {

    auto cit = _children.begin();
    advance(cit, i);
    READCOUNTS rc = *(cit->second.front());
    for (size_t j = 1; j < cit->second.size(); ++j)
       rc += *(cit->second[j]);

#ifdef _OPENMP
    #pragma omp critical(update_clade_counts)
#endif
    {
      _cladeCounts.insert( std::make_pair( cit->first, std::move(rc) ) );
    }
  }
  
  cerr << " done" << endl;


  _report_cols = {REPORTCOLS::PERCENTAGE, REPORTCOLS::NUM_READS_CLADE, REPORTCOLS::NUM_READS, 
    REPORTCOLS::NUM_KMERS_CLADE, REPORTCOLS::NUM_UNIQUE_KMERS_CLADE, 
    REPORTCOLS::NUM_KMERS_IN_DATABASE_CLADE, REPORTCOLS::TAX_RANK, REPORTCOLS::TAX_ID, 
    REPORTCOLS::SPACED_NAME};
}


template<typename TAXID, typename READCOUNTS>
void TaxReport<TAXID,READCOUNTS>::setReportCols(const std::vector<std::string>& names) {
  _report_cols.clear();
  for (size_t i = 0; i< names.size(); ++i) {
    auto& s = names[i];
    auto it = report_col_name_map.find(s);
    if (it == report_col_name_map.end()) {
      throw std::runtime_error(s + " is not a valid report column name");
    }
    _report_cols.push_back(it->second);
  }
  _report_col_names = names;
}

template<typename TAXID, typename READCOUNTS>
void TaxReport<TAXID,READCOUNTS>::printReport(const std::string& format) {

  cerr << "Printing classification report ... ";
  _total_n_reads = 0;
  for (int i : { 0, 1, -1 } ) {
     auto it = _taxdb.entries.find(i);
     if (it != _taxdb.entries.end()) {
      const auto it2 = _cladeCounts.find(&(it->second));
      if (it2 != _cladeCounts.end()) 
        //_total_n_reads += reads(it2->second);
        _total_n_reads += it2->second.readCount();
     }
   }

  if (_total_n_reads == 0) {
    std::cerr << "total number of reads is zero - not creating a report!" << endl;
    return;
  }
  if (_report_cols.size() == _report_col_names.size()) {
    // print header
    _reportOfb << _report_col_names[0];
    for (size_t i=1; i < _report_col_names.size(); ++i) {
      _reportOfb << '\t' << _report_col_names[i];
    }
    _reportOfb << '\n';
  }

  if (format == "kraken") {
    // A: print number of unidentified reads
    // B: print classified results
    // C: Print Unclassified stuff
    for (int i : { 0, 1, -1 } ) {
      auto it = _taxdb.entries.find(i);
      if (it != _taxdb.entries.end()) {
        printReport(it->second,0u);
      }
    }
  } else {
    // print stuff at a certain level ..
    //_uid_abundance;
    //_taxinfo

  }

  cerr << " done" << endl;
}

template<typename TAXID, typename READCOUNTS>
void TaxReport<TAXID,READCOUNTS>::printReport(const TaxonomyEntry<TAXID>& tax, unsigned depth) {

    const auto taxit_ptr = _cladeCounts.find(&tax);
    if (taxit_ptr == _cladeCounts.end())
      return;
    const auto & cladecounts = taxit_ptr->second;
    if (!_show_zeros && cladecounts.readCount() == 0)
      return;

    printLine(tax, cladecounts, depth);


    // Sort children
    vector<size_t> pos;
    unordered_map<size_t, READCOUNTS*> rc;
    for (size_t i =0; i < tax.children.size(); ++i) {
      auto it = _cladeCounts.find(tax.children[i]);
      if (it != _cladeCounts.end()) {
        pos.push_back(i);
        rc[i] = &(it->second);
      }
    }
    std::sort(pos.begin(), pos.end(), [&](size_t a, size_t b) { return *(rc.at(b)) < *(rc.at(a)) ;} );

    for (size_t i=0; i < rc.size(); ++i) {
      auto child_it = tax.children[ pos[i] ];
      printReport(*child_it, depth+1);
    }
}

template<typename TAXID, typename READCOUNTS>
void TaxReport<TAXID,READCOUNTS>::printLine(const TaxonomyEntry<TAXID>& tax, const READCOUNTS& rc, unsigned depth) {
  const auto r_it = _taxCounts.find(tax.taxonomyID);
  const bool has_tax_data = r_it != _taxCounts.end();

  const uint64_t unique_kmers_for_clade = rc.uniqueKmerCount();
  double genome_size = double(tax.genomeSize+tax.genomeSizeOfChildren);

  for (size_t i = 0; i< _report_cols.size(); ++i) {
    auto& col = _report_cols[i];
    switch (col) {
      case REPORTCOLS::NAME:              _reportOfb << tax.scientificName ; break;
      case REPORTCOLS::SPACED_NAME:       _reportOfb << string(2*depth, ' ') + tax.scientificName; break;
      case REPORTCOLS::TAX_ID:            _reportOfb << (tax.taxonomyID == (uint32_t)-1? -1 : (int32_t) tax.taxonomyID); break;
      case REPORTCOLS::DEPTH:             _reportOfb << depth; break;
      case REPORTCOLS::PERCENTAGE:       _reportOfb << setprecision(4) << 100.0*rc.readCount()/_total_n_reads; break;
           //case REPORTCOLS::ABUNDANCE:      _reportOfb << 100*counts.abundance[0]; break;
           //case REPORTCOLS::ABUNDANCE_LEN:  _reportOfb << 100*counts.abundance[1]; break;
      case REPORTCOLS::NUM_READS:        _reportOfb << (has_tax_data? r_it->second.readCount() : 0); break;
      case REPORTCOLS::NUM_READS_CLADE:  _reportOfb << rc.readCount(); break;
      case REPORTCOLS::NUM_UNIQUE_KMERS: _reportOfb << (has_tax_data? r_it->second.kmerCount() : 0); break;
      case REPORTCOLS::NUM_UNIQUE_KMERS_CLADE:  _reportOfb << unique_kmers_for_clade; break;
      case REPORTCOLS::NUM_KMERS:        _reportOfb << (has_tax_data? r_it->second.kmerCount() : 0); break;
      case REPORTCOLS::NUM_KMERS_CLADE:  _reportOfb << rc.kmerCount(); break;
      case REPORTCOLS::NUM_KMERS_IN_DATABASE: _reportOfb << tax.genomeSize; break;
      case REPORTCOLS::CLADE_KMER_COVERAGE: 
                if (genome_size == 0) { 
            _reportOfb << "NA"; 
                } else {
            _reportOfb << setprecision(4) << (unique_kmers_for_clade  / genome_size); 
                }; break;
      case REPORTCOLS::CLADE_KMER_DUPLICITY: _reportOfb << setprecision(3) << ( double(rc.kmerCount()) / unique_kmers_for_clade ); break;
      case REPORTCOLS::NUM_KMERS_IN_DATABASE_CLADE: _reportOfb << tax.genomeSize + tax.genomeSizeOfChildren; break;
                //case REPORTCOLS::GENOME_SIZE: ; break;
                //case REPORTCOLS::NUM_WEIGHTED_READS: ; break;
                //case REPORTCOLS::SUM_SCORE: ; break;
      case REPORTCOLS::TAX_RANK: _reportOfb << tax.rank; break;
      default: _reportOfb << "NA";
    }
    if (&col == &_report_cols.back()) {
      _reportOfb << '\n';
    } else {
      _reportOfb << '\t';
    }
  }
}


template<typename K,typename V>
inline
V find_or_use_default(const std::unordered_map<K, V>& my_map, const K& query, const V default_value) {
  auto itr = my_map.find(query);

  if (itr == my_map.end()) {
    return default_value;
  }

  return itr->second;
}



#endif /* TAXD_DB_H_ */

