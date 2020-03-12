%module exceRNApipeline_taxaCounter
%{
#include "TaxaCounter.h"
%}

%include "std_vector.i"
%include "std_string.i"
%include "std_pair.i"
%include "std_unordered_map.i"

void taxaCounter(const std::string& input_path, const std::string& output_prefix,
const std::string& names_dmp, const std::string& nodes_dmp);