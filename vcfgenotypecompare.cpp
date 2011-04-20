#include "Variant.h"
#include "split.h"
#include <string>
#include <iostream>
#include <sstream>

using namespace std;
using namespace vcf;

map<string, int> decomposeGenotype(string& genotype) {
    string splitter = "/";
    if (genotype.find("|") != string::npos) {
        splitter = "|";
    }
    vector<string> haps = split(genotype, splitter);
    map<string, int> decomposed;
    for (vector<string>::iterator h = haps.begin(); h != haps.end(); ++h) {
        ++decomposed[*h];
    }
    return decomposed;
}

bool isHet(map<string, int>& genotype) {
    return genotype.size() > 1;
}

bool isHom(map<string, int>& genotype) {
    return genotype.size() == 1;
}

bool hasNonRef(map<string, int>& genotype) {
    return genotype.find("1") != genotype.end();
}

bool isHomRef(map<string, int>& genotype) {
    return isHom(genotype) && !hasNonRef(genotype);
}

bool isHomNonRef(map<string, int>& genotype) {
    return isHom(genotype) && hasNonRef(genotype);
}

bool isNull(map<string, int>& genotype) {
    return genotype.find(".") != genotype.end();
}

string genotypeSpec(map<string, int>& genotype) {
    string gtspec;
    if (isNull(genotype)) {
        gtspec = "NN";
    } else if (isHom(genotype)) {
        if (hasNonRef(genotype)) {
            gtspec = "AA";
        } else {
            gtspec = "RR";
        }
    } else {
        gtspec = "AR";
    }
    return gtspec;
}

int main(int argc, char** argv) {

    if (argc != 3) {
        cerr << "usage: " << argv[0] << " <other-genotype-tag> <vcf file>" << endl
             << "adds statistics to the INFO field of the vcf file describing the" << endl
             << "amount of discrepancy between the genotypes (GT) in the vcf file and the" << endl
             << "genotypes reported in the <other-genotype-tag>.  use this after" << endl
             << "vcfannotategenotypes to get correspondence statistics for two vcfs." << endl;
        return 1;
    }

    string otherGenoTag = argv[1];
    string filename = argv[2];

    VariantCallFile variantFile;
    if (filename == "-") {
        variantFile.open(std::cin);
    } else {
        variantFile.open(filename);
    }

    if (!variantFile.is_open()) {
        return 1;
    }

    vector<string> specs;
    specs.push_back("AA_AA");
    specs.push_back("AA_AR");
    specs.push_back("AA_RR");
    specs.push_back("AR_AR");
    specs.push_back("AR_RR");
    specs.push_back("RR_RR");
    specs.push_back("AA_NN");
    specs.push_back("AR_NN");
    specs.push_back("RR_NN");
    specs.push_back("NN_NN");

    for (vector<string>::iterator spec = specs.begin(); spec != specs.end(); ++spec) {
        string line = "##INFO=<ID=" + otherGenoTag + ".genotypes." + *spec
            + ",Number=1,Type=Integer,Description=\"Number of genotypes with "
            + *spec + " relationship with " + otherGenoTag + "\">";
        variantFile.addHeaderLine(line);
    }

    string line;

    line = "##INFO=<ID=" + otherGenoTag + ".genotypes.count,Number=1,Type=Integer,Description=\"Count of genotypes under comparison.\">";
    variantFile.addHeaderLine(line);

    line = "##INFO=<ID=" + otherGenoTag + ".genotypes.alternate_count,Number=1,Type=Integer,Description=\"Count of alternate genotypes in the first file.\">";
    variantFile.addHeaderLine(line);

    line = "##INFO=<ID=" + otherGenoTag
        + ".site.alternate_positive_discrepancy,Number=1,Type=Integer,Description=\"Estimated positive discrepancy rate of "
        + otherGenoTag + " genotypes, where positive discrepancies are all cases where an alternate allele is called GT "
        + " but none is represented in " + otherGenoTag + " or " + otherGenoTag + " is null/no-call\">";
    variantFile.addHeaderLine(line);

    line = "##INFO=<ID=" + otherGenoTag
        + ".site.alternate_negative_discrepancy,Number=1,Type=Integer,Description=\"Estimated negative discrepancy rate of "
        + otherGenoTag + " genotypes, where negative discrepancies are all cases where no alternate allele is called in "
        + " GT but an alternate is represented in " + otherGenoTag + ", including no-calls or partly null genotypes\">";
    variantFile.addHeaderLine(line);

    line = "##INFO=<ID=" + otherGenoTag
        + ".site.alternate_null_discrepancy,Number=1,Type=Integer,Description=\"Estimated null discrepancy rate of "
        + otherGenoTag + " genotypes, where null discrepancies are all cases where GT is specified and contains an alternate but "
        + otherGenoTag + " is null.  Cases where GT is null or partly null are excluded.\">";
    variantFile.addHeaderLine(line);

    line = "##INFO=<ID=" + otherGenoTag
        + ".site.call_discrepancy,Number=1,Type=Integer,Description=\"Estimated call discrepancy rate of "
        + otherGenoTag + " genotypes (het->hom, hom->het) between " + otherGenoTag + " and GT\">";
    variantFile.addHeaderLine(line);

    line = "##INFO=<ID=" + otherGenoTag
        + ".site.call_concordance,Number=1,Type=Integer,Description=\"Estimated call concorndance rate of "
        + otherGenoTag + " genotypes between " + otherGenoTag + " and GT\">";
    variantFile.addHeaderLine(line);

    line = "##INFO=<ID=" + otherGenoTag
        + ".site.non_reference_discrepancy,Number=1,Type=Float,Description=\"Estimated non-reference discrepancy relative to "
        + otherGenoTag + " genotypes,\">";
    variantFile.addHeaderLine(line);

    line = "##INFO=<ID=" + otherGenoTag
        + ".site.non_reference_discrepancy.count,Number=1,Type=Int,Description=\"non-reference discrepancy normalizer relative to "
        + otherGenoTag + " genotypes,\">";
    variantFile.addHeaderLine(line);

    line = "##INFO=<ID=" + otherGenoTag
        + ".site.non_reference_discrepancy.normalizer,Number=1,Type=Int,Description=\"non-reference discrepancy count relative to "
        + otherGenoTag + " genotypes,\">";
    variantFile.addHeaderLine(line);

    line = "##INFO=<ID=" + otherGenoTag
        + ".site.non_reference_sensitivity,Number=1,Type=Float,Description=\"Estimated non-reference sensitivity relative to "
        + otherGenoTag + " genotypes,\">";
    variantFile.addHeaderLine(line);

    line = "##INFO=<ID=" + otherGenoTag
        + ".site.non_reference_sensitivity.count,Number=1,Type=Int,Description=\"non-reference sensitivity normalizer relative to "
        + otherGenoTag + " genotypes,\">";
    variantFile.addHeaderLine(line);

    line = "##INFO=<ID=" + otherGenoTag
        + ".site.non_reference_sensitivity.normalizer,Number=1,Type=Int,Description=\"non-reference sensitivity count relative to "
        + otherGenoTag + " genotypes,\">";
    variantFile.addHeaderLine(line);

    cout << variantFile.header << endl;

    Variant var(variantFile);

    while (variantFile.getNextVariant(var)) {
        // for each sample, check GT against <other-genotype-tag>
        // tally stats, and append to info
        map<string, map<string, string> >::iterator s     = var.samples.begin();
        map<string, map<string, string> >::iterator sEnd  = var.samples.end();

        map<string, int> genotypeComparisonCounts;
        int gtCount = var.samples.size();
        int gtAltCount = 0; // number of alternate-containing genotypes in the first file
        int pdCount = 0; // positive discrepancy count
        int ndCount = 0; // negative discrepancy count
        int nnCount = 0; // null discrepancy count
        int cdCount = 0; // call discrepancy count
        int ccCount = 0; // call concordance count
        int nrdCount = 0; // non-reference discrepancy count
        int nrdNormalizer = 0; // divisor for nrd rate
        int nrsCount = 0; // non-reference sensitivity count
        int nrsNormalizer = 0; // divisor for nrs rate

        for (; s != sEnd; ++s) {
            map<string, string>& sample = s->second;
            const string& name = s->first;
            // decompose genotypes into counts of strings
            // to facilitate comparison
            string& gtA = sample["GT"];
            string& gtB = sample[otherGenoTag];
            map<string, int> genotypeA = decomposeGenotype(gtA);
            map<string, int> genotypeB = decomposeGenotype(gtB);

            string gtspecA = genotypeSpec(genotypeA);
            string gtspecB = genotypeSpec(genotypeB);
            //cout << gtA << " " << gtB << endl;
            //cout << gtspecA << " " << gtspecB << endl;
            ++genotypeComparisonCounts[gtspecA + "_" + gtspecB];

            if (hasNonRef(genotypeA)) {
                ++gtAltCount;
            }

            if (genotypeA != genotypeB) {
                if (isNull(genotypeA)) {
                    // TODO handle this somehow, maybe via a different flag?
                    if (!isNull(genotypeB)) {
                        ++nnCount;  // null discrepancy, the second set makes a call, this one does not
                    }
                } else if (hasNonRef(genotypeA)) {
                    if (hasNonRef(genotypeB)) { // they cannot be the same, but they both represent an alternate
                        ++cdCount;  // the calls are discrepant
                    } else { // the other call does not have an alternate
                        ++pdCount;
                        // it is also null
                        if (isNull(genotypeB)) {
                            ++nnCount;
                        }
                    }
                } else { // the current genotype has no non-ref alternate
                    if (hasNonRef(genotypeB)) {
                        ++ndCount;
                        if (isNull(genotypeB)) {
                            ++nnCount;
                        }
                    }
                }
            } else {
                if (!isNull(genotypeA)) {
                    ++ccCount;
                }
            }


            if (!(isNull(genotypeA) || isNull(genotypeB))
                    && !(isHomRef(genotypeA) && isHomRef(genotypeB))) {
                ++nrdNormalizer;
                if (genotypeA != genotypeB) {
                    ++nrdCount;
                }
            }

            if (!(isNull(genotypeB) || isHomRef(genotypeB))) {
                ++nrsNormalizer;
                if (!(isNull(genotypeA) || isHomRef(genotypeA))) {
                    ++nrsCount;
                }
            }

        }

        for (map<string, int>::iterator g = genotypeComparisonCounts.begin();
                g != genotypeComparisonCounts.end(); ++g) {
            stringstream c;
            c << g->second;
            var.info[otherGenoTag + ".genotypes." + g->first] = c.str();
        }

        stringstream gtc;
        gtc << gtCount;
        var.info[otherGenoTag + ".genotypes.count"] = gtc.str();

        stringstream gtac;
        gtac << gtAltCount;
        var.info[otherGenoTag + ".genotypes.alternate_count"] = gtac.str();

        stringstream pd;
        pd << pdCount;
        var.info[otherGenoTag + ".site.alternate_positive_discrepancy"] = pd.str();

        stringstream nd;
        nd << ndCount;
        var.info[otherGenoTag + ".site.alternate_negative_discrepancy"] = nd.str();

        stringstream nn;
        nn << nnCount;
        var.info[otherGenoTag + ".site.alternate_null_discrepancy"] = nn.str();

        stringstream cd;
        cd << cdCount;
        var.info[otherGenoTag + ".site.call_discrepancy"] = cd.str();

        stringstream cc;
        cc << ccCount;
        var.info[otherGenoTag + ".site.call_concordance"] = cc.str();

        stringstream nrdc;
        nrdc << nrdCount;
        var.info[otherGenoTag + ".site.non_reference_discrepancy.count"] = nrdc.str();

        stringstream nrdn;
        nrdn << nrdNormalizer;
        var.info[otherGenoTag + ".site.non_reference_discrepancy.normalizer"] = nrdn.str();

        if (nrdNormalizer > 0) {
            stringstream nrd;
            nrd << (double) nrdCount / (double) nrdNormalizer;
            var.info[otherGenoTag + ".site.non_reference_discrepancy"] = nrd.str();
        }

        stringstream nrsc;
        nrsc << nrsCount;
        var.info[otherGenoTag + ".site.non_reference_sensitivity.count"] = nrsc.str();

        stringstream nrsn;
        nrsn << nrsNormalizer;
        var.info[otherGenoTag + ".site.non_reference_sensitivity.normalizer"] = nrsn.str();

        if (nrsNormalizer > 0) {
            stringstream nrs;
            nrs << (double) nrsCount / (double) nrsNormalizer;
            var.info[otherGenoTag + ".site.non_reference_sensitivity"] = nrs.str();
        }

        cout << var << endl;

    }

    return 0;

}

