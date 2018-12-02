// SJG & YQ Dec2018

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

    std::vector<bcfRecord> bcf2vec(htsReader<bcfRecord>::iterator begin, htsReader<bcfRecord>::iterator end);


    // ***LOAD RECORDS AND HEADERS INTO MEMORY
    std::vector<bcfHeader> headers;                     // Vector of headers, one per file - stay in static positions
    std::vector<std::vector<bcfRecord>> recordArr;      // Vector of vectors of records - each record is a smart ptr (bcfRecord)
    for (int i = 0; i < numFiles; i++) {
        auto file = htsOpen(argv[i + 1], "r");

        auto header = htsHeader<bcfHeader>::read(file);
        bcfHeader dupHdr{bcf_hdr_dup(header.get())};
        headers.push_back(std::move(dupHdr));

        auto records = bcf2vec(begin(file, header), end(file, header));
        recordArr.push_back(records);
    }

    // ***COMBINE HEADER INFO AND OUTPUT
    // TODO: need all phase info here to coordinate with chips, for now just outputting first header
    htsOutHandle << headers[0];


    // ***COMBINE ALL RECORD ARRAYS INTO SINGLE ARRAY
    std::vector<int> lineIndex(numFiles, 0);            // Array of line index for which record looking at in each file
    std::vector<bcfRecord*> currRecs(numFiles);         // Array of the current record being compared for each file
    std::vector<bcfRecord> combinedRecs;                // Giant array of ordered, combined lines

    // Initialize index array and put first line of each file in current record array
    for (int i = 0; i < numFiles; i++) {
        lineIndex.push_back(0);
        currRecs.push_back(&(recordArr[i][0]));
    }

    // Combine lines into single array
    std::vector<bool> fileFinished(numFiles, false);
    bool allFilesFinished = false;
    int lowIndex = 0;
    int lowCoord = INFINITY;
    while (!allFilesFinished) {
        for (int i = 0; i < numFiles; i++) {
            if (!fileFinished[i] && (*currRecs[i])->pos+1 <= lowCoord) {
                // TODO: if we have variant at same position but diff allele, order by allele
                lowCoord = (*currRecs[i])->pos+1;
                lowIndex = i;
            }
        }

        // Add lowest record & fill in spot
        combinedRecs.push_back(std::move(*currRecs[lowIndex]));             // Transfer control of pointer to combined array
        lineIndex[lowIndex] = lineIndex[lowIndex]++;                        // Advance current line in file we just pulled from
        if (lineIndex[lowIndex] < recordArr[lowIndex].size()) {
            currRecs[lowIndex] = &recordArr[lowIndex][lineIndex[lowIndex]]; // Fill in next line from file
        } else {
            fileFinished[lowIndex] = true;
        }

        // Reset coord for next run through
        lowCoord = INFINITY;

        // Check if we've reached the end of all files
        allFilesFinished = true;
        for (int i = 0; i < numFiles; i++) {
            allFilesFinished &= fileFinished[i];
        }
    }


    // ***REMOVE DUPLICATES FROM COMBINED, ORDERED ARRAY
    int currLine = 0;
    int nextLine = 1;
    while (nextLine < combinedRecs.size()) {
        if (combinedRecs[currLine]->pos+1 == combinedRecs[nextLine]->pos+1) {
            if (combinedRecs[currLine]->alt == combinedRecs[nextLine]->alt) {   // TODO: figure out how to get allele info here
            }
            else {
                htsOutHandle << bcfHdrRecPair(headers[0], combinedRecs[currLine]);  // TODO: use combined header here
                currLine = nextLine;
            }
            nextLine++;
        }
    }
}









// *** OTHER APPROACHES / MUSINGS...


//    std::vector<int> currMap(numFiles, 0);              // Lookup for recordArr - static position matches coord array, dynamic value matches record position
//    std::vector<int> nextMap(numFiles, 0);
//
//    std::vector<int> currCoords(numFiles, 0);           // Coordinates for current position
//    std::vector<int> nextCoords(numFiles, 0);           // Coordinates for next position
//
//    std::vector<std::string> currAlts(numFiles, "");    // Variants at current position
//    std::vector<std::string> nextAlts(numFiles, "");    // Variants at next position
//
//    std::vector<std::vector<int>> currCounts;           // Counts for current position
//    std::vector<std::vector<int>> nextCounts;           // Counts for next position


// ***HELPER FXNS
//htsIterator advanceAndFill(htsIterator record, std::vector<int> &coords, std::vector<std::string> &alts, std::vector<int*> &counts, int arrIndex, int advanceNum);

//    void sortAndCombine(std::vector<htsIterator> &currRecords, std::vector<htsIterator> &nextRecords, std::vector<int> &currCoords, std::vector<int> &nextCoords,
//            std::vector<std::string> &currAlts, std::vector<std::string> &nextAlts, std::vector<int*> &currCounts, std::vector<int*> &nextCounts);



//    // ***INITIALIZE ARRAYS
//    for (int i = 0; i < numFiles; i++) {
//
//    }
//
//
//        // TODO: might make below into a callable loop w/ line # to retrieve records by (instead of hard coded 0 and 1)
//        // Get current fields and put in respective arrays
//        currCoords[i] = records[0]->pos + 1;
//
//        // Get alt strings and put in array
//        // TODO: have to call vcf_unpack here
//
//        // Get counts and put in array
//        int numAlts = (records[0]->n_allele) - 1;
//        int countArrLength = (numAlts * 2) + 2;
//        std::vector<int> counts((countArrLength * 2), 0);
//        int status = bcf_get_info_int32(header.get(), &(records[0]), "ENRICH_COUNTS", &counts[0], (countArrLength * 2));     // TODO: why doesn't this work?
//        if (status < 0) {
//            logger(LOGLV_WARN) << "Could not obtain enrichment counts at " << bcf_hdr_id2name(header.get(), records[0]->rid)
//                << ":"
//                << records[0]->pos + 1 << std::endl;
//        }
//        currCounts[0] = counts;
//
//        // Advance iterator and get next fields
//        nextCoords[i] = records[1]->pos + 1;
//        // TODO: get alt string(s) and put in array - copy above when complete
//        // TODO: get counts from INFO field and put in array - copy above when complete
//
//
//    // ***MAIN LOOP
//    int finishCount = 0;
//    while (finishCount < numFiles) {
//        for (int i = 0; i < numFiles; i++) {
//
//            // Stable sort the arrays and combine any duplicates
//            sortAndCombine(currRecords, nextRecords, currCoords, nextCoords, currAlts, nextAlts, currCounts,
//                           nextCounts);
//
//            // Get first record and write out
//            auto rec = std::move(currRecords[i]);
//            htsOutHandle << rec;                        // **Why doesn't this work?
//
//            // Advance to the next-next record and fill the coordinate, alt, and count arrays
//            rec = std::move(advanceAndFill(std::move(rec), currCoords, currAlts, currCounts, i, 2));
//
//            // Find next smallest coordinate
//            if (nextCoords[i] < currCoords[i + 1]) {
//                rec = std::move(nextRecords[i]);
//                htsOutHandle << rec;
//                rec = std::move(advanceAndFill(std::move(rec), nextCoords, nextAlts, nextCounts, i, 2));
//            } else if (currCoords[i + 1] < nextCoords[i]) {
//                rec = std::move(currRecords[i + 1]);
//                htsOutHandle << rec;
//                rec = std::move(advanceAndFill(std::move(rec), currCoords, currAlts, currCounts, i + 1, 2));
//            }
//            // If we get here, we have a tie for INF, meaning we're at the end of both files
//            else {
//                finishCount += 2;
//            }
//        }
//    }
//    return 0;
//}
//
///* Retrieves control of unique pointer to which lineNum corresponds. Obtains relative information and overwrites current info in given arrays as appropriate at the given index.
// * Returns control of unique record pointer to array.
// *
// * NOTE: ASSUMES ONLY VALID INDICES ARE PROVIDED
// */
//void advanceAndFill(int index, int lineNum, std::vector<std::vector<bcfRecord>> recordArr, std::vector<bcfHeader> headers,
//        std::vector<int> &map, std::vector<int> &coords, std::vector<std::string> &alts, std::vector<int*> &counts) {
//
//    // Get record that map points to
//    int fileIndex = map[index];
//    auto rec = std::move((recordArr[fileIndex])[lineNum]);
//
//    // Get header that map points to
//    auto hdr = std::move(headers[fileIndex]);
//
//
//}
//
///* Takes in current and next arrays. Stable sorts all arrays in same relative order. Then searches for duplicates and combines as necessary. */
//void sortAndCombine() {
//    // TODO: implement - this will be the pain in the butt to get correct
//}
//
///* Takes in record, coordinate, alternate, and count arrays, along with indices to swap. Performs swap. */
//void swap() {
//    // TODO: implement
//}
//
///* Takes in iterator to single file, copies all lines as smart pointers to vector, returns vector and control of smart pointers. */
//std::vector<bcfRecord> bcf2vec(htsReader<bcfRecord>::iterator begin, htsReader<bcfRecord>::iterator end) {
//    std::vector<bcfRecord> records;
//    std::for_each(std::move(begin), std::move(end), [&records](const auto &rec) {
//       bcfRecord dupRec{bcf_dup(rec.get())};
//       records.push_back(std::move(dupRec));        // Record pointer only exists in array - pass reference around
//    });
//    return std::move(records);
//}






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