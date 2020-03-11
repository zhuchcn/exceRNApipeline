#include "TaxaCounter.h"
#include "utilities.hpp"
#include <algorithm>
#include <utility>
#include <fstream>
#include <iostream>

void TaxaCounter::loadNCBITaxDump(const std::string nodes_dmp, const std::string names_dmp, bool reduce){
    _db = NCBITaxonomy::loadNCBITaxDump(nodes_dmp, names_dmp);
    if(reduce) _db->reduceTaxLevels();
}

void TaxaCounter::readReads(const std::string& path){
    if(endswith(path, ".gz")) {
        igzstream io (path.c_str());
        std::istream_iterator<Line> end;
        std::istream_iterator<Line> it (io);
        _readReads(it, end);
    } else {
        std::ifstream io (path.c_str());
        std::istream_iterator<Line> end;
        std::istream_iterator<Line> it (io);
        _readReads(it, end);
    }
}

void TaxaCounter::_readReads(std::istream_iterator<Line> it, std::istream_iterator<Line> end) {
    std::cout << "Start reading reads." << std::endl;
    while(it != end) {
        std::vector<std::string> line = split(*(it++), "\t");
        std::string& read_id = line[0];
        std::string& tax_name = line[1];
        for (int i = 0; i < int(tax_name.length()); ++i) {
            if(tax_name[i] == '_') tax_name[i] = ' ';
        }
        std::transform(tax_name.begin(), tax_name.end(), tax_name.begin(), ::tolower);
        if (!(_db->hasTax(tax_name))) {
            std::vector<std::string> tax_name_split = split(tax_name, " ");
            if(tax_name_split.size() >= 2) {
                tax_name = tax_name_split[0] + tax_name_split[1];
                if(!_db->hasTax(tax_name)) {
                    tax_name = tax_name_split[0];
                    if(!_db->hasTax(tax_name)) {
                        continue;
                    }
                }
            } else {
                continue;
            }
        }
        int tax_id = _db->getTaxId(tax_name)[0];
        if(_reads.find(read_id) == _reads.end()) {
            NodePair nodePair = std::make_pair(tax_id, _db->getNode(tax_id));
            _reads.insert({read_id, nodePair});
            // _reads[read_id] = std::make_pair(tax_id, _db->getNode(tax_id));
        } else {
            _reads[read_id] = _db->getCommonAncestor(_reads[read_id].first, tax_id);
        }
        
        while(std::find(_taxonomy_levels.begin(), _taxonomy_levels.end(), 
                _reads[read_id].second->taxonomy_level) == _taxonomy_levels.end()) {
            _reads[read_id] = _db->getParentNode(_reads[read_id].second->parent_id);
        }
    }
    std::cout << "Reading reads done." << std::endl;
}

void TaxaCounter::countReads() {
    std::cout << "Start counting reads." << std::endl;
    for(auto it = _reads.begin(); it != _reads.end(); ++it) {
        NodePair& nodePair = it->second;
        while (_counts.find(nodePair.second->taxonomy_level) != _counts.end()) {
            const std::string& level = nodePair.second->taxonomy_level;
            const std::string& name_txt = nodePair.second->name_txt;
            if(level == "root") {
                if(_counts["root"].size() == 0) {
                    _counts["root"]["root"] = 1;
                } else {
                    _counts["root"]["root"] ++;
                }
                break;
            }
            if(_counts[level].find(name_txt) == _counts[level].end()) {
                _counts[level][name_txt] = 1;
            } else {
                _counts[level][name_txt]++;
            }

            nodePair = _db->getParentNode(nodePair.first);
            while(std::find(_taxonomy_levels.begin(), _taxonomy_levels.end(),
                    nodePair.second->taxonomy_level) == _taxonomy_levels.end()) {
                nodePair = _db->getParentNode(nodePair.first);
            }
        }
    }
    std::cout << "Counting reads done." << std::endl;
}

void TaxaCounter::writeCounts(const std::string& prefix) {
    std::cout << "Start writing counts." << std::endl;
    for(auto it = _counts.begin(); it != _counts.end(); ++it) {
        std::string output_path = prefix + it->first + ".txt";
        std::ofstream os(output_path.c_str());
        for(auto cit = it->second.begin(); cit != it->second.end(); ++cit) {
            os << cit->first << "\t" << cit->second << std::endl;
        }
        os.close();
    }
    std::cout << "Writing counts done." << std::endl;
}

void taxaCounter(const std::string& input_path, const std::string& output_prefix,
const std::string& names_dmp, const std::string& nodes_dmp) {
    TaxaCounter tc;
    tc.loadNCBITaxDump(nodes_dmp, names_dmp, true);
    tc.readReads(input_path);
    tc.countReads();
    tc.writeCounts(output_prefix);
}