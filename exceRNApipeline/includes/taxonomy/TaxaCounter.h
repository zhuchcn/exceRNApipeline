#ifndef _TAXA_COUNTER_H_
#define _TAXA_COUNTER_H_

#include "NCBITaxonomy.h"
#include "utilities.hpp"
#include "../gzstream/gzstream.h"
#include <unordered_map>
#include <vector>
#include <string>


class TaxaCounter {
public:
    using NodePair = std::pair<int, const NCBITaxonomy::Node*>;

private:
    std::unordered_map<std::string, NodePair> _reads;
    std::unordered_map<std::string, std::unordered_map<std::string, int>> _counts;
    std::vector<std::string> _taxonomy_levels;
    NCBITaxonomy* _db;

    void _readReads(std::istream_iterator<Line> it, std::istream_iterator<Line> end);

public:
    TaxaCounter(): _reads(), _counts(), _taxonomy_levels(_DEFAULT_TAXONOMY_LEVELS), _db() {
        for (auto it = _taxonomy_levels.begin(); it < _taxonomy_levels.end(); ++it) {
            _counts.insert({*it, {}});
        }
    };
    virtual ~TaxaCounter() {};

    void loadNCBITaxDump(std::string nodes_dmp, std::string names_dmp, bool reduce);

    void readReads(const std::string& path);

    void countReads();

    void writeCounts(const std::string& prefix);
};


void taxaCounter(const std::string& input_path, const std::string& output_prefix,
const std::string& names_dmp, const std::string& nodes_dmp);


#endif