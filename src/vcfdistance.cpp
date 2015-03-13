#include "Variant.h"
#include "split.h"
#include <string>
#include <iostream>

using namespace std;
using namespace vcflib;

int main(int argc, char** argv) {

    if (argc > 1) {
        cerr << "usage: " << argv[0] << " <[vcf file]" << endl
             << "adds a tag (BasesToClosestVariant) to each variant record which indicates" << endl
             << "the distance to the nearest variant" << endl;
        return 1;
    }

    VariantCallFile variantFile;
    variantFile.open(std::cin);

    if (!variantFile.is_open()) {
        return 1;
    }

    Variant varA(variantFile);
    Variant varB(variantFile);
    Variant varC(variantFile);

    vector<Variant*> vars;
    vars.push_back(&varA);
    vars.push_back(&varB);
    vars.push_back(&varC);
    
    for (vector<Variant*>::iterator v = vars.begin(); v != vars.end(); ++v) {
        variantFile.getNextVariant(**v);
    }

    string tag = "BasesToClosestVariant";
    string line = "##INFO=<ID=" + tag + ",Number=1,Type=Integer,Description=\"" \
        + "Number of bases to the closest variant in the file.\">";
    variantFile.addHeaderLine(line);

    cout << variantFile.header << endl;

    if (!vars.at(0)->sequenceName.empty()) {
        if (!vars.at(1)->sequenceName.empty()) {
            // at least two variants, so calculate the first distance
            if (vars.at(0)->sequenceName == vars.at(1)->sequenceName) {
                vars.at(0)->info[tag].push_back(convert(vars.at(1)->position - vars.at(0)->position));
            }
            cout << *vars.at(0) << endl;

            if (!vars.at(2)->sequenceName.empty()) {
                // at least three variants, so starting with the first three,
                // calculate the middle variant's closest distance, and then
                // slide the window forward one.
                do {
                    if (vars.at(1)->sequenceName == vars.at(0)->sequenceName &&
                        vars.at(1)->sequenceName == vars.at(2)->sequenceName) {
                        vars.at(1)->info[tag].push_back(convert(min(vars.at(1)->position - vars.at(0)->position,
                                                                    vars.at(2)->position - vars.at(1)->position)));
                    } else if (vars.at(1)->sequenceName == vars.at(0)->sequenceName) {
                        vars.at(1)->info[tag].push_back(convert(vars.at(1)->position - vars.at(0)->position));
                    } else if (vars.at(2)->sequenceName == vars.at(1)->sequenceName) {
                        vars.at(1)->info[tag].push_back(convert(vars.at(2)->position - vars.at(1)->position));
                    } else {
                        // don't add the tag
                    }
                    cout << *vars.at(1) << endl;
                    // rotate
                    Variant* v = vars.at(0);
                    vars.at(0) = vars.at(1);
                    vars.at(1) = vars.at(2);
                    vars.at(2) = v;
                } while (variantFile.getNextVariant(*vars.back()));
            }

            // assign the last distance and output the last variant
            if (vars.at(0)->sequenceName == vars.at(1)->sequenceName) {
                vars.at(1)->info[tag].push_back(convert(vars.at(1)->position - vars.at(0)->position));
            }
            cout << *vars.at(1) << endl;
        } else {
            // output the lone variant line untouched
            cout << *vars.at(0) << endl;
        }
    }

    return 0;

}

