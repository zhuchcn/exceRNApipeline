#include "NCBITaxonomy.h"
#include "NCBITaxonomyRanks.hpp"
#include "utilities.hpp"
#include <utility>
#include <algorithm>

const std::unordered_map<std::string, int>* NCBITaxonomy::_taxonomy_ranks = &_NCBI_TAXONOMY_RANKS_;

const NCBITaxonomy::Node* NCBITaxonomy::getNode(int tax_id) const {
    if(_nodes.find(tax_id) == _nodes.end()) {
        throw std::runtime_error("The tax_id doesn't exist.");
    }
    return &_nodes.at(tax_id);
}

void NCBITaxonomy::insertName(std::string name, int tax_id) {
    auto it = _names.find(name);
    if(it != _names.end()) {
        std::vector<int> & second = it->second;
        if(std::find(second.begin(), second.end(), tax_id) != second.end()){
            std::cerr << "The taxonomy and tax_id pair (" << name << ", " << tax_id << ") already exists." << std::endl;
            return;
        }
        second.push_back(tax_id);
    } else {
        _names[name] = {tax_id};
    }
}

void NCBITaxonomy::insertNode(int tax_id, std::string name_txt, int parent_id, std::string rank) {
    if(_nodes.find(tax_id) != _nodes.end()) {
        std::cerr << "The tax_id " << tax_id << "already exists" << std::endl;
        return;
    }
    _nodes.insert({tax_id, {name_txt, parent_id, rank}});
}

std::vector<int> NCBITaxonomy::getTaxId(std::string name_txt) const{
    auto it = _names.find(name_txt);
    if(it == _names.end()) {
        std::cerr << "The tax name " << name_txt << " wasn't found." << std::endl;
        throw std::runtime_error("NCBITaxonomy::getTaxId(): key " + name_txt + " was not found");
    }
    return it->second;
}

bool NCBITaxonomy::hasTax(std::string name_txt) const {
    return _names.find(name_txt) != _names.end();
}

std::pair<int, NCBITaxonomy::Node*> NCBITaxonomy::getParentNode(int tax_id) {
    Node* node = &_nodes.at(tax_id);
    if(_nodes.find(node->parent_id) == _nodes.end()) {
        std::string msg = "The parent (id = " + std::to_string(node->parent_id) + 
            ") of entry " + std::to_string(tax_id) + " does not exist.";
        throw std::runtime_error(msg);
    }
    return std::make_pair(node->parent_id, &_nodes.at(node->parent_id));
}

std::pair<int, NCBITaxonomy::Node*> NCBITaxonomy::getCommonAncestor(int tax_id1, int tax_id2){
    
    struct TaxNode {
        std::string* name_txt;
        std::string* taxonomy_level;
        int* parent_id;
        int rank;
        TaxNode(std::string* name, std::string* level, int* parent, int rnk):
            name_txt(name), taxonomy_level(level), parent_id(parent), rank(rnk) {};
    };

    auto getTaxNode = [&](int& tax_id) {
        Node* node = &_nodes.at(tax_id);
        std::string* name_txt = &node->name_txt;
        std::string* level = &node->taxonomy_level;
        int* parent = &node->parent_id;
        return (TaxNode (name_txt, level, parent, _taxonomy_ranks->at(*level)));
    };

    auto node1 = getTaxNode(tax_id1);
    auto node2 = getTaxNode(tax_id2);

    while(*node1.name_txt != *node2.name_txt){
        while (node1.rank != node2.rank) {
            if(node1.rank > node2.rank) {
                tax_id1 = *node1.parent_id;
                node1 = getTaxNode(tax_id1);
            } else {
                tax_id2 = *node2.parent_id;
                node2 = getTaxNode(tax_id2);
            }
        }
        tax_id1 = *node1.parent_id;
        tax_id2 = *node2.parent_id;
        node1 = getTaxNode(*node1.parent_id);
        node2 = getTaxNode(*node2.parent_id);
    }
    return std::make_pair(tax_id1, &_nodes.at(tax_id1));
};

void NCBITaxonomy::reduceTaxLevels() {
    for(auto it = _nodes.begin(); it != _nodes.end(); ++it) {
        int tax_id = it->first;
        Node* node = &it->second;
        
        if(std::find(_taxonomy_levels.begin(), _taxonomy_levels.end(), node->taxonomy_level) 
                == _taxonomy_levels.end()) 
            continue;

        auto parent = getParentNode(tax_id);
        tax_id = parent.first;
        node = parent.second;

        while(std::find(_taxonomy_levels.begin(), _taxonomy_levels.end(), node->taxonomy_level) 
                == _taxonomy_levels.end()) {
            auto parent = getParentNode(tax_id);
            tax_id = parent.first;
            node = parent.second;
        }
        it->second.parent_id = tax_id;
    }

    for(auto it = _nodes.begin(); it != _nodes.end();) {
        if(std::find(_taxonomy_levels.begin(), _taxonomy_levels.end(), it->second.taxonomy_level) 
                == _taxonomy_levels.end()) {
            it = _nodes.erase(it);
            continue;
        }
        ++it;
    }

    for(auto it = _names.begin(); it != _names.end();) {
        for(auto vit = it->second.begin(); vit < it->second.end();) {
            if(_nodes.find(*vit) == _nodes.end()) {
                vit = it->second.erase(vit);
                continue;
            }
            ++vit;
        }
        if(it->second.size() == 0) {
            it = _names.erase(it);
            continue;
        }
        ++it;
    }
}

NCBITaxonomy* NCBITaxonomy::loadNCBITaxDump(std::string nodes_dmp, std::string names_dmp) {
    NCBITaxonomy* taxa = new NCBITaxonomy;

    // Read names.dmp
    std::cout << "loading " << names_dmp << std::endl;
    std::fstream names_io (names_dmp);
    std::istream_iterator<Line> end;
    std::istream_iterator<Line> it1 (names_io);

    std::unordered_map<int, std::string> idNameDict;

    while(it1 != end) {
        std::string l (*(it1++));
        std::regex e ("\t\\|$");
        l = std::regex_replace(l, e, "");

        std::vector<std::string> row = split(l, "\t|\t");
        if(row[row.size()-1] != "scientific name") continue;
        
        e = "[^A-Za-z0-9]";
        row[1] = std::regex_replace(row[1], e, " ");
        e = " +";
        row[1] = std::regex_replace(row[1], e, " ");
        std::transform(row[1].begin(), row[1].end(), row[1].begin(), ::tolower);
        idNameDict[std::stoi(row[0])] = row[1];
    }

    // Read nodes.dmp
    std::cout << "loading " << nodes_dmp << std::endl;
    std::fstream nodes_io (nodes_dmp);
    std::istream_iterator<Line> it2 (nodes_io);

    while(it2 != end) {
        std::string l (*(it2++));
        std::vector<std::string> row = split(l, "\t|\t");
        if(row[0] == "1" && row[2] == "no rank") row[2] = "root";
        int id = std::stoi(row[0]);
        int parentId = std::stoi(row[1]);
        taxa->insertNode(id, idNameDict[id], parentId, row[2]);
    }

    std::cout << "saving names" << std::endl;
    for(auto it = idNameDict.begin(); it != idNameDict.end(); ++it) {
        // std::cout << "inserting: (" << it->second << ", " << it->first << ")" << std::endl;
        taxa->insertName(it->second, it->first);
    }

    return taxa;
}