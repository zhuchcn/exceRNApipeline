#ifndef CATCH_CONFIG_MAIN
#define CATCH_CONFIG_MAIN

#include "../NCBITaxonomy.h"
#include "../TaxaCounter.h"
#include "../../catch/catch.hpp"
#include <iostream>
#include <string>

TEST_CASE("Test NCBITaxonomy class", "[NCBITaxonomy]") {
    SECTION("Testing names") {
        NCBITaxonomy tax;
        tax.insertName("root", 1);
        tax.insertName("Azorhizobium caulinodans", 7);
        tax.insertName("Azorhizobium caulinodans", 7);
        REQUIRE(tax.getTaxId("Azorhizobium caulinodans")[0] == 7);
        REQUIRE(tax.getTaxId("Azorhizobium caulinodans").size() == 1);
        REQUIRE(tax.hasTax("Azorhizobium caulinodans"));
        REQUIRE(!tax.hasTax("Azorhizobium caulin"));
    }

    SECTION("Testing nodes") {
        NCBITaxonomy tax;
        tax.insertNode(10, "Cellvibrio", 1706371, "genus");
    }

    SECTION("Testing load from dump") {
        NCBITaxonomy* taxa = NCBITaxonomy::loadNCBITaxDump("tests/nodes_test.dmp", "tests/names_test.dmp");
        REQUIRE(taxa->hasTax("methylophilus methylotrophus"));
        std::vector<int> ids = taxa->getTaxId("methylophilus methylotrophus");
        REQUIRE(std::find(ids.begin(), ids.end(), 17) != ids.end());
    }

    SECTION("Testing common ancestor") {
        NCBITaxonomy* taxa = NCBITaxonomy::loadNCBITaxDump("tests/nodes_test.dmp", "tests/names_test.dmp");
        auto res = taxa->getCommonAncestor(23, 24);
        REQUIRE(res.first == 22);
        REQUIRE(res.second->name_txt == "shewanella");
    }

    SECTION("Testing reduce tax levels") {
        NCBITaxonomy* taxa = NCBITaxonomy::loadNCBITaxDump("tests/nodes_test2.dmp", "tests/names_test2.dmp");
        taxa->reduceTaxLevels();
        // delta/epsilon subdivisions is a subphylum which should be removed.
        REQUIRE(!taxa->hasTax("delta/epsilon subdivisions"));
        // Buchnera is a genus so should be kept.
        REQUIRE(taxa->hasTax("buchnera"));
    }
}

TEST_CASE("Test taxa counter") {
    SECTION("Testing taxa counter") {
        taxaCounter("tests/sample_readCounts.txt", "tests/sample_readCounts_", "tests/names_test2.dmp", "tests/nodes_test2.dmp");
    }
}

#endif