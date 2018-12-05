// SJG & YQ Dec2018
// Combines variants from multiple files
// Puts multi-allelic sites into a single line
// Combines the following fields:
// ALT, INFO-> GT, AF,

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


    // ***PRIME LOOP
    // Load first record from each file into our working array
    std::vector<bcfRecord*> currRecs(numFiles);                 // Vector of individual lines - one per file and in relative order to files/headers
    for (int i = 0; i < numFiles; i++) {
        currRecs.push_back(&(recordArr[i][0]));
    }
    bool advanceMin = false;        // True if we've seen all lines at or below the current coordMin
    int coordMin = INFINITY;        // The minimum coordinate in currRecs per iteration
    bcfHeader* outHdr = nullptr;       // The combined header for the output line
    bcfRecord* outLine = nullptr;      // The next line to be written to the output stream
    bool outOfLines = false;        // True if we have any remaining lines to be processed in any files


    // ***PROCESS LINES
    // Find minimum coordinate
    while (!outOfLines) {
        for (int i = 0; i < numFiles; i++) {
            int currCoord = (*currRecs[i])->pos + 1;
            if (currCoord < coordMin) {
                coordMin = currCoord;
            }
        }
        // Retrieve lines with the minimum coordinate
        while(!advanceMin) {
            bool noLines = true;
            bool noLineAtMin = true;
            for (int i = 0; i < numFiles; i++) {
                int currCoord = (*currRecs[i])->pos + 1;

                // Found coordinate match
                if (currCoord == coordMin) {
                    noLineAtMin = false;
                    // If first line at current coordinate
                    if (outLine == NULL) {
                        outHdr = &headers[i];
                        outLine = currRecs[i];
                        // TODO: removeInfo(outHdr, outLine); - remove all INFO fields besides genotypes for simplicity & browser memory issues? [optional]
                    // Combine lines
                    } else {
                        // Get genotypes out of outLine
                        int status, *out_gt_arr = nullptr, out_ngt_arr = 0;
                        status = bcf_get_genotypes((*outHdr).get(), (*outLine).get(), &out_gt_arr, &out_ngt_arr);
                        if (status <= 0) {
                            logger(LOGLV_WARN) << "skipping record at " << bcf_hdr_id2name((*outHdr).get(), (*outLine)->rid)
                                               << ":" << (*outLine)->pos + 1 << " - could not obtain genotype" << std::endl;
                            continue;
                        }
                        std::unique_ptr<int> uniq_out_gt_arr(out_gt_arr);

                        // Get genotypes out of current line
                        int *curr_gt_arr = nullptr, curr_ngt_arr = 0;
                        status = bcf_get_genotypes(headers[i].get(), (*currRecs[i]).get(), &curr_gt_arr, &curr_ngt_arr);
                        if (status <= 0) {
                            logger(LOGLV_WARN) << "skipping record at " << bcf_hdr_id2name(headers[i].get(), (*currRecs[i])->rid)
                                               << ":" << (*currRecs[i])->pos + 1 << " - could not obtain genotype" << std::endl;
                            continue;
                        }
                        std::unique_ptr<int> uniq_curr_gt_arr(curr_gt_arr);

                        // Combine genotypes and header info
                        // TODO
                    }
                    // Put next line in spot we just pulled out
                    currCoord = getNextLine(currRecs[i]);
                    // TODO: noLines &= endOfFile

                // Sanity check for testing
                } else if (currCoord < coordMin) {
                    logger(LOGLV_WARN) << "ERROR: minimum is not minimum";
                }
            }
            advanceMin = noLineAtMin;
            outOfLines = noLines;
        }

        // Output line - guaranteed to have all files processed at min here
        htsOutHandle << bcfHdrRecPair(*outHdr, *outLine);

        // Reset loop parameters
        outLine = nullptr;
        outHdr = nullptr;
        coordMin = INFINITY;
    }




    // <editor-fold desc="Combine All Records">
    // ***COMBINE ALL RECORD ARRAYS INTO SINGLE ARRAY
    /*std::vector<int> lineIndex(numFiles, 0);                // Array of line index for which record looking at in each file
    std::vector<bcfRecord *> currRecs(numFiles);            // Array of the current record being compared for each file
    std::vector<bcfRecord> combinedRecs;                    // Giant array of ordered, combined lines

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
            if (!fileFinished[i] && (*currRecs[i])->pos + 1 <= lowCoord) {
                lowCoord = (*currRecs[i])->pos + 1;
                lowIndex = i;
            }
        }

        // Add lowest record & fill in spot
        combinedRecs.push_back(
                std::move(*currRecs[lowIndex]));                    // Transfer control of pointer to combined array
        lineIndex[lowIndex] = lineIndex[lowIndex]++;                // Advance current line in file we just pulled from
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
    }*/
    //</editor-fold>


    // ***COMBINE ALL CALLS FOR A SINGLE SITE
    int currLine = 0;
    int nextLine = 1;
    bool matchFlag = false;
    // declare start line array of enrichment counts here

    while (nextLine < combinedRecs.size()) {
        if (combinedRecs[currLine]->pos + 1 == combinedRecs[nextLine]->pos + 1) {
            matchFlag = true;

            // Pull out gt array for start line like below - store in outside prop



            // Pull out gt array for next line
            int status, *next_gt_arr = nullptr, next_ngt_arr = 0;

            // TODO: have to have header lookup here
            status = bcf_get_genotypes(header.get(), combinedRecs[nextLine].get(), &next_gt_arr, &next_ngt_arr);
            if (status <= 0) {
                logger(LOGLV_WARN) << "skipping record at " << bcf_hdr_id2name(header.get(), rec->rid)
                                   << ":" << rec->pos + 1 << " - could not obtain genotype" << std::endl;
                continue;
            }
            std::unique_ptr<int> uniq_gt_arr(next_gt_arr);


            // for each alt in next_gt_arr
                // if (start_gt_arr.contains(alt)
                    // find alt index of this alt in start
                    // combine counts from next at pulled index - helper fxn
                // else
                    // append alt to start_gt_arr
                    // combine enrichment counts for alt into start enrichment counts - helper fxn

        } else {
            if (matchFlag) {
                // Output combined line at start             htsOutHandle << bcfHdrRecPair(headers[0], combinedRecs[currLine]);  // TODO: use combined header here
                matchFlag = false;
                currLine = nextLine;
            } else {
                // Just output line at start? TODO: flesh this out
            }
        }
        nextLine++;
    }
}


// SUBLOOP
// foreach line in currLines
    // if (line.pos == globalMin)
        // if (outputLine empty) outputLine = line & clear out INFO fields other than gt;
        // else combine genotypes from line into genotypes from outputLine;
        // advance line
// get nextMin of lines
    // if nextMin > currMin
        // calculate enrichment for line and append to INFO field
        // write outputLine
        // outputLine - NULL
    // else
        // goto foreach
