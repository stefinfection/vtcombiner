#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm>
#include <tuple>
#include <math.h>
#include "log.h"
//#include "config.h" TODO: ask Yi about this
#include "../contrib/htslibpp/htslibpp.h"
#include "../contrib/htslibpp/htslibpp_variant.h"
#include "/usr/local/Cellar/gsl/2.5/include/gsl/gsl_randist.h"
#include "/usr/local/Cellar/gsl/2.5/include/gsl/gsl_cdf.h"
#include "/usr/local/Cellar/gsl/2.5/include/gsl/gsl_sf_gamma.h"

using namespace YiCppLib::HTSLibpp;

int main(int argc, const char *argv[]) {

    // Helper fxns
    double getPValFishers(int a, int b, int c, int d, int e, int f);
    double getPValCochArm(int caseHomRef, int caseHet, int caseHomAlt, int ctrlHomRef, int ctrlHet, int ctrlHomAlt);

    // Ensure we have at least one file coming in
    if (argc < 2) {
        std::cout << "No files found to process." << std::endl;
        exit(1);
    }

    // Open write stream to stdout
    auto htsOutHandle = htsOpen("-", "w");

    // Single file case
    if (argc == 2) {
        std::cout << "Only one file provided" << std::endl;

        // Write out updated header
        auto htsInHandle = htsOpen(argv[1], "r");
        auto header = htsHeader<bcfHeader>::read(htsInHandle);
        bcf_hdr_append(header.get(),
                "##INFO=<ID=ENRICH_PVAL,Number=A,Type=Float,Description=\"P-value calculated from Fisher Exact Test or Cochrane Armitage Trend Test as appropriate\">");
        htsOutHandle << header;

        std::for_each(begin(htsInHandle, header), end(htsInHandle, header), [&](auto &rec) {
            int numAlts = rec->n_allele - 1;
            int offset = (2 * numAlts) + 2;

            // Get enrich counts from record
            int *countArr = nullptr, nCountArr = 0;
            int status = bcf_get_info_int32(header.get(), rec.get(), "ENRICH_COUNTS", &countArr, &nCountArr);
            if (status < 0) {
                std::cout << "Skipping record at " << bcf_hdr_id2name(header.get(), rec->rid)
                    << ":" << rec->pos + 1 << " - could not obtain enrichment counts" << std::endl;
                return;
            }

            // Create smart pointer for counts
            std::unique_ptr<int> uniqCountArr(countArr);

            // Assign p-value function according to sample size
            int numSamples = 0;
            for (int i = 0; i < offset; i++) {
                numSamples += countArr[i];
            }
            double (*pFunc)(int a, int b, int c, int d, int e, int f) = getPValCochArm;
            if (numSamples < 1000) {
                pFunc = getPValFishers;
            }

            // Get p-values for each alt
            std::vector<double> pVals(numAlts, 0);
            for (int i = 1; i < (numAlts * 2) - 1; i = i + 2) {
                int proHomRef = countArr[0];
                int subHomRef = countArr[0 + offset];
                int nonSubHomRef = proHomRef - subHomRef;

                int proHet = countArr[i];
                int subHet = countArr[i + offset];
                int nonSubHet = proHet - subHet;

                int proHomAlt = countArr[i + 1];
                int subHomAlt = countArr[i + 1 + offset];
                int nonSubHomAlt = proHomAlt - subHomAlt;

                pVals[i - 1] = pFunc(subHomRef, subHet, subHomAlt, nonSubHomRef, nonSubHet, nonSubHomAlt);
            }

            // Insert p-values into line and output
            bcf_update_info_int32(header.get(), rec.get(), "ENRICH_PVAL", &pVals[0], numAlts);
            bcf_subset(header.get(), rec.get(), 0, nullptr);
            htsOutHandle << bcfHdrRecPair(header, rec);
        });
        exit(0);
    }

    // Open read streams for inputs
    std::vector<YiCppLib::HTSLibpp::htsFile> htsInHandles;
    htsInHandles.reserve(argc-1);
    for (int i = 0; i < argc-1; i++) {
        // Get handler for file
        htsInHandles[i] = htsOpen(argv[i+1], "r");
    }

    // Pipe any single header out
    auto header = htsHeader<bcfHeader>::read(htsInHandles[0]);
    if (header.get() == nullptr) {
        std::cerr << "unable to read header from input stream" << std::endl;
        exit(1);
    }
    htsOutHandle << header;


    //








    return 0;
}

/* Takes in values A, B, C, D representing a 2x2 contingency table:
 *
 *     case    control
 *      --      --
 *      a    |   b      |   ref. allele count
 *      c    |   d      |   alt. allele count
 *
 * Performs a Fisher's Exact Test and returns the corresponding p-value.
 *
 * NOTE: utilizes GSL Gamma library
 * NOTE 2: e & f parameters inserted for compatibility with fxn pointer
 */
double getPValFishers(int caseHomRef, int caseHet, int caseHomAlt, int ctrlHomRef, int ctrlHet, int ctrlHomAlt) {
    const size_t a = caseHomRef * 2 + caseHet;    // caseRefCount
    int c = caseHomAlt * 2 + caseHet;    // caseAltCount
    const size_t b = ctrlHomRef * 2 + ctrlHet;    // ctrlRefCount
    int d = ctrlHomAlt * 2 + ctrlHet;    // ctrlAltCount
    int n = a + b + c + d;

    double num = gsl_sf_fact(a + b) * gsl_sf_fact(c + d) * gsl_sf_fact(a + c) * gsl_sf_fact(b + d);
    double denom = gsl_sf_fact(a) + gsl_sf_fact(b) + gsl_sf_fact(c) + gsl_sf_fact(d) + gsl_sf_fact(n);
    return num / denom;
}

/* Performs a Cochran Armitage Trend test for the given values. Takes in the following parameters:
*  homozygous reference counts (AA), heterozygous alternate counts (Aa), homozygous alternate counts (aa) from an
*  experimental group (aka subset cohort), and a control group (aka non-subset cohort), respectively.
*  Assumes a codominant model for weighting of (0,1,2).
*  Returns a p-value based on a single degree of freedom, or one independent variable.
*
*  For more info: http://statgen.org/wp-content/uploads/2012/08/armitage.pdf
*
* NOTE: utilizes GSL CDF and Random Numbers libraries
*/
double getPValCochArm(int caseHomRef, int caseHet, int caseHomAlt, int ctrlHomRef, int ctrlHet, int ctrlHomAlt) {
    // Sum matrix totals
    int ctrlTotal = ctrlHomRef + ctrlHet + ctrlHomAlt;  // S
    int caseTotal = caseHomRef + caseHet + caseHomAlt;  // R
    int sampleTotal = ctrlTotal + caseTotal;            // N
    int hetTotal = ctrlHet + caseHet;                   // n_1
    int homAltTotal = ctrlHomAlt + caseHomAlt;          // n_2

    // Calculate z-value using Saseini notation and weighted vector of (0,1,2)
    // TODO: double check we won't get overflow here
    int testStat = (sampleTotal * (caseHet + (2 * caseHomAlt))) - (caseTotal * (hetTotal + (2 * homAltTotal)));     // U
    int testVarA = ctrlTotal * caseTotal;                                                                           // SR
    int testVarB = sampleTotal * (hetTotal + (4 * homAltTotal));                                                    // N(n_1 + 4n_2)
    int testVarC = pow((hetTotal + (2 * homAltTotal)), 2.0);                                                        // (n_1 + 2n_2)^2
    int num = sampleTotal * (pow(testStat, 2.0));                                                                   // N * U^2
    int denom = testVarA * (testVarB - testVarC);                                                                   // A(B-C)
    double z = sqrt((num / denom));

    // Calculate p-value from z-value using one degree of freedom
    return (gsl_cdf_chisq_Q(z, 1) - gsl_cdf_chisq_P(z, 1)) * 2;
}