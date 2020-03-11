#ifndef _NCBI_TAXONOMY_H_
#define _NCBI_TAXONOMY_H_

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <exception>
#include <fstream>
#include <iterator>
#include <regex>

const std::vector<std::string> _DEFAULT_TAXONOMY_LEVELS ({
    "root", "kingdom", "phylum", "class", "order", "family", "genus", "species" });

// The main NCBITaxonomy class
class NCBITaxonomy {
public:
    // Contains the name_txt, parent_id, and the taxonomy rank of the current
    // node.
    class Node {
    public:
        std::string name_txt;
        int parent_id;
        std::string taxonomy_level;

        Node() {};
        Node(std::string name, int parent, std::string level):
            name_txt(name), parent_id(parent), taxonomy_level(level) {};
    };

private:
    std::unordered_map<int, Node> _nodes;
    std::unordered_map<std::string, std::vector<int>> _names;
    std::vector<std::string> _taxonomy_levels;
    static const std::unordered_map<std::string, int>* _taxonomy_ranks;
    
public:
    NCBITaxonomy(): _nodes(), _names(), _taxonomy_levels(_DEFAULT_TAXONOMY_LEVELS) {};
    virtual ~NCBITaxonomy(){};

    const Node* getNode(int) const;
    
    void insertNode(int tax_id, std::string name_txt, int parent_id, std::string rank);
    void insertName(std::string name, int tax_id);

    std::vector<int> getTaxId(std::string name_txt) const;
    bool hasTax(std::string name_txt) const;
    std::pair<int, Node*> getParentNode(int tax_id);

    // Given any two tax_ids, this method returns the tax label that is above 
    // both of them, aka, the common acestor.
    std::pair<int, Node*> getCommonAncestor(int tax_id1, int tax_id2);

    // The NCBI's taxdump is very detailed thus it contains also some taxonomy 
    // levels that are outside of the common levels (kingdom, phylum, class,
    // order, family, genus, and species). This method reduces the tax levels
    // by removing anything that are not in the taxonomy_levels, in order to 
    // provide a better performance of getCommonAncestor.
    void reduceTaxLevels();

    // Create a NCBITaxonomy from the ncbi tax-dump files.
    static  NCBITaxonomy* loadNCBITaxDump(std::string nodes_dmp, std::string names_dmp);
};

#endif