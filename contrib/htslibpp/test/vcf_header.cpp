#include <gmock/gmock.h>
#include "../htslibpp.h"
#include "../htslibpp_variant.h"

#include <string>
#include <fstream>
#include <streambuf>
#include <algorithm>
#include <cstdio>

using namespace YiCppLib::HTSLibpp;

class VcfHeader: public testing::Test {
    public:
        YiCppLib::HTSLibpp::htsFile htsFileHandler = htsOpen("datasets/brca2.platnium-trio.vcf", "r");
};


// Tests for READ

TEST_F(VcfHeader, CanCreateBcfHeaderObject) {
    auto header = htsHeader<bcfHeader>::read(htsFileHandler);
    ASSERT_NE(header.get(), nullptr);
}

TEST_F(VcfHeader, CanIterateOverRecordsSequencially) {
    auto header = htsHeader<bcfHeader>::read(htsFileHandler);
    int32_t header_line_ct = 0;
    std::for_each(std::cbegin(header), std::cend(header), [&header_line_ct](const auto& p) { header_line_ct++; });
    ASSERT_EQ(header_line_ct, 59);
}

TEST_F(VcfHeader, CanIterateOverIDDictionary) {
    auto header = htsHeader<bcfHeader>::read(htsFileHandler);
    auto dict = htsHeader<bcfHeader>::DictType::ID;

    int32_t header_filter_ct = 0;
    int32_t header_info_ct = 0;
    int32_t header_format_ct = 0;

    std::for_each(htsHeader<bcfHeader>::dictBegin(header, dict), htsHeader<bcfHeader>::dictEnd(header, dict), [&](const auto& p){

            auto proxy = htsProxy(p);

            if(proxy.hasValueForLineType(htsHeader<bcfHeader>::LineType::FILTER)) header_filter_ct++;
            if(proxy.hasValueForLineType(htsHeader<bcfHeader>::LineType::INFO)) header_info_ct++;
            if(proxy.hasValueForLineType(htsHeader<bcfHeader>::LineType::FORMAT)) header_format_ct++;

        });

    ASSERT_EQ(header_filter_ct, 1);
    ASSERT_EQ(header_info_ct, 43);
    ASSERT_EQ(header_format_ct, 8);
}

TEST_F(VcfHeader, CanIterateOverContigDictionary) {
    auto header = htsHeader<bcfHeader>::read(htsFileHandler);
    auto dict = htsHeader<bcfHeader>::DictType::CONTIG;

    int32_t contig_ct = 0;

    std::for_each(htsHeader<bcfHeader>::dictBegin(header, dict), htsHeader<bcfHeader>::dictEnd(header, dict), [&](const auto& p) {
            auto proxy = HTSProxyIDPairContig(p);
            if(proxy.hasValueForLineType(htsHeader<bcfHeader>::LineType::CONTIG)) contig_ct++;
        });

    ASSERT_EQ(contig_ct, 1);
}

TEST_F(VcfHeader, CanIterateOverSampleDictionary) {
    auto header = htsHeader<bcfHeader>::read(htsFileHandler);
    auto dict = htsHeader<bcfHeader>::DictType::SAMPLE;

    int32_t sample_ct = 0;

    std::for_each(htsHeader<bcfHeader>::dictBegin(header, dict), htsHeader<bcfHeader>::dictEnd(header, dict), [&](const auto& p) {
            auto proxy = htsProxy(p);
            ASSERT_EQ(strncmp(proxy.key(), "NA128", 5), 0);  // The sample names should all start with NA128
            sample_ct++;
        });

    ASSERT_EQ(sample_ct, 3);
}

// Tests for MODIFY

/* operations right now supported in HTSLib:
 * merge(dst, src): ret_val <- dst + src
 * add_sample(hdr, sample): hdr <- hdr + sample
 * set(hdr, fname): hdr <- hdr_read(fname)
 * remove(hdr, type, key): hdr <- hdr - (key of type id)
 *
 * HREC funcs:
 * hrec_add_key(hrec, str, len): add key of 'str' 
 */



// Tests for CREATE
// Tests for WRITE

TEST_F(VcfHeader, CanWriteHeaderUsingStreamOperator) {
    auto header = htsHeader<bcfHeader>::read(htsFileHandler);

    // forced scope so outHandle closes
    {
        auto outHandle = htsOpen("test-can-write-header-using-stream-operator.vcf", "w");
        outHandle<<header;
    }

    std::ifstream t("test-can-write-header-using-stream-operator.vcf");
    std::string str((std::istreambuf_iterator<char>(t)), std::istreambuf_iterator<char>());
    ASSERT_NE(str.length(), 0);

    auto n = std::count(str.cbegin(), str.cend(), '\n');
    ASSERT_EQ(n, 60);

    remove("test-can-write-header-using-stream-operator.vcf");
}
