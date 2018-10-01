#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm>
#include <tuple>
#include <math.h>
#include "log.h"
//#include "config.h" TODO: Ask YQ
#include "../contrib/htslibpp/htslibpp.h"
#include "../contrib/htslibpp/htslibpp_variant.h"
#include "/usr/local/Cellar/gsl/2.5/include/gsl/gsl_randist.h"
#include "/usr/local/Cellar/gsl/2.5/include/gsl/gsl_randist.h"
#include "/usr/local/Cellar/gsl/2.5/include/gsl/gsl_cdf.h"
#include "/usr/local/Cellar/gsl/2.5/include/gsl/gsl_sf_gamma.h"

using namespace YiCppLib::HTSLibpp;

int main(int argc, const char *argv[]) {

    // Ensure we have at least one file coming in
    int numFiles = argc - 1;
    if (numFiles < 2) {
        std::cout << "VtCombiner requires at least two files to combine. Please try again." << std::endl;
        exit(1);
    }

    // Open write stream to stdout
    auto htsOutHandle = htsOpen("-", "w");

    // Declare containers for holding necessary info for each file
    std::vector<htsIterator> currRecords;               // One iterator per record for current position
    std::vector<htsIterator> nextRecords;               // One iterator per record for next position

    std::vector<int> currCoords(numFiles, 0);           // Coordinates for current position
    std::vector<int> nextCoords(numFiles, 0);           // Coordinates for next position

    std::vector<std::string> currAlts(numFiles, "");    // Variants at current position
    std::vector<std::string> nextAlts(numFiles, "");    // Variants at next position

    std::vector<int*> currCounts(numFiles, nullptr);    // Counts for current position
    std::vector<int*> nextCounts(numFiles, nullptr);    // Counts for next position

    // ***HELPER FXNS
    htsIterator advanceAndFill(htsIterator record, std::vector<int> &coords, std::vector<std::string> &alts, std::vector<int*> &counts, int arrIndex, int advanceNum);
    void sortAndCombine(std::vector<htsIterator> &currRecords, std::vector<htsIterator> &nextRecords, std::vector<int> &currCoords, std::vector<int> &nextCoords,
            std::vector<std::string> &currAlts, std::vector<std::string> &nextAlts, std::vector<int*> &currCounts, std::vector<int*> &nextCounts);

    // ***GET ITERATORS FOR EACH FILE AND POPULATE ABOVE ARRAYS TO PRIME MAIN LOOP
    for (int i = 0; i < numFiles; i++) {
        auto file = htsOpen(argv[i + 1], "r");
        auto header = htsHeader<bcfHeader>::read(file);
        auto currRec = begin(file, header);             // Get first line in file

        auto nextRec = begin(file, header);             // Get next line in file
        nextRec++;

        // Get current fields and put in respective arrays
        currCoords[i] = currRec.pos + 1;
        // TODO: get alt string(s) and put in array
        // TODO: get counts from INFO field and put in array

        // Advance iterator and get next fields
        nextCoords[i] = nextRec.pos + 1;
        // TODO: get alt string(s) and put in array
        // TODO: get counts from INFO field and put in array

        // Move iterator back to first position and assign control of iterator to vector
        currRecords.push_back(std::move(currRec));      // **Why doesn't this work?
        nextRecords.push_back(std::move(nextRec));      // **Why doesn't this work?
    }

    // ***MAIN LOOP
    int finishCount = 0;
    while (finishCount < numFiles) {
        for (int i = 0; i < numFiles; i++) {

            // Stable sort the arrays and combine any duplicates
            sortAndCombine(currRecords, nextRecords, currCoords, nextCoords, currAlts, nextAlts, currCounts,
                           nextCounts);

            // Get first record and write out
            auto rec = std::move(currRecords[i]);
            htsOutHandle << rec;                        // **Why doesn't this work?

            // Advance to the next-next record and fill the coordinate, alt, and count arrays
            rec = std::move(advanceAndFill(std::move(rec), currCoords, currAlts, currCounts, i, 2));

            // Find next smallest coordinate
            if (nextCoords[i] < currCoords[i + 1]) {
                rec = std::move(nextRecords[i]);
                htsOutHandle << rec;
                rec = std::move(advanceAndFill(std::move(rec), nextCoords, nextAlts, nextCounts, i, 2));
            } else if (currCoords[i + 1] < nextCoords[i]) {
                rec = std::move(currRecords[i + 1]);
                htsOutHandle << rec;
                rec = std::move(advanceAndFill(std::move(rec), currCoords, currAlts, currCounts, i + 1, 2));
            }
            // If we get here, we have a tie for INF, meaning we're at the end of both files
            else {
                finishCount += 2;
            }
        }
    }
    return 0;
}

/* Advances given iterator by given value. Overwrites the information in the given arrays as appropriate at the given index.
 * Returns control of iterator back to caller. */
htsIterator advanceAndFill(htsIterator record, std::vector<int> &coords, std::vector<std::string> &alts, std::vector<int*> &counts, int arrIndex, int advanceNum) {
    // TODO: implement
    // Check for null/end of file here - if so, put +inf in corresponding coord position
    return record;
}

/* Takes in current and next arrays. Stable sorts all arrays in same relative order. Then searches for duplicates and combines as necessary. */
void sortAndCombine() {
    // TODO: implement - this will be the pain in the butt to get correct
}

/* Takes in record, coordinate, alternate, and count arrays, along with indices to swap. Performs swap. */
void swap() {
    // TODO: implement
}






// *** OTHER APPROACHES / MUSINGS...
//    MAY NEED TO USE NON-WRAPPED VERSION
//    std::vector<::htsFile*> files(numFiles, nullptr);
//    std::vector<bcf_hdr_t*> headers(numFiles, nullptr);
//    std::vector<bcf1_t*> lines(numFiles, nullptr);
//
//    for (int i = 0; i < numFiles; i++) {
//        files[i] = HTSLIB_VCF_H::vcf_open(argv[i+1], "r");
//        headers[i] = HTSLIB_VCF_H::bcf_hdr_read(files[i]);
//
//        // Read next line
//        HTSLIB_VCF_H::vcf_read(files[i], headers[i], lines[i]);
//
//        // TEST: Output line
//        htsOutHandle << lines[i];
//    }

//  TODO: can't seem to get Yi's lib to work - ask him if I'm reinventing the wheel
//    std::vector<YiCppLib::HTSLibpp::htsFile> files;
//    std::vector<YiCppLib::HTSLibpp::htsHeader<bcfHeader>> headers;
//    std::vector<YiCppLib::HTSLibpp::htsReader<bcfRecord>> lines;

//    for (int i = 0; i < numFiles; i++) {
//        auto currFile = YiCppLib::HTSLibpp::htsOpen(argv[i+1], "r");
//        //auto currHdr = YiCppLib::HTSLibpp::htsHeader<bcfHeader>::read(currFile);
//        auto currHdr = YiCppLib::HTSLibpp::htsHeader<bcfHeader>::read(currFile);
//        headers[i] = std::move(currHdr);
//    }

//    // Pipe any single header out
//    auto header = htsHeader<bcfHeader>::read(files[0]);
//    if (header.get() == nullptr) {
//        std::cerr << "Unable to read header from first file. Please try again." << std::endl;
//        exit(1);
//    }
//    htsOutHandle << header;

// TODO: implement multi-loop iterations
//    std::vector<bcf1_t*> linePtrs(numFiles, nullptr);
//    for (int i = 0; i < numFiles; i++) {
//
//        // Get pointers to first line for each file
//        header = htsHeader<bcfHeader>::read(htsInHandles[i]);
//        auto rec = htsReader<bcfRecord>::read(htsInHandles[i], header);
//
//
//    }