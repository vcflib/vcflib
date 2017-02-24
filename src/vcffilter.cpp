#include "Variant.h"
#include "split.h"
#include <getopt.h>

using namespace std;
using namespace vcflib;

void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] <vcf file>" << endl
         << endl
         << "options:" << endl
         << "    -f, --info-filter     specifies a filter to apply to the info fields of records," << endl
         << "                          removes alleles which do not pass the filter" << endl
         << "    -g, --genotype-filter specifies a filter to apply to the genotype fields of records" << endl
         << "    -k, --keep-info       used in conjunction with '-g', keeps variant info, but removes genotype" << endl
         << "    -s, --filter-sites    filter entire records, not just alleles" << endl
         << "    -t, --tag-pass        tag vcf records as positively filtered with this tag, print all records" << endl
         << "    -F, --tag-fail        tag vcf records as negatively filtered with this tag, print all records" << endl
         << "    -A, --append-filter   append the existing filter tag, don't just replace it" << endl
         << "    -a, --allele-tag      apply -t on a per-allele basis.  adds or sets the corresponding INFO field tag" << endl
         << "    -v, --invert          inverts the filter, e.g. grep -v" << endl
         << "    -o, --or              use logical OR instead of AND to combine filters" << endl
         << "    -r, --region          specify a region on which to target the filtering, requires a BGZF" << endl
         << "                          compressed file which has been indexed with tabix.  any number of" << endl
         << "                          regions may be specified." << endl
         << endl
         << "Filter the specified vcf file using the set of filters." << endl
         << "Filters are specified in the form \"<ID> <operator> <value>:" << endl
         << " -f \"DP > 10\"  # for info fields" << endl
         << " -g \"GT = 1|1\" # for genotype fields" << endl
         << " -f \"CpG\"  # for 'flag' fields" << endl
         << endl
         << "Operators can be any of: =, !, <, >, |, &" << endl
         << endl
         << "Any number of filters may be specified.  They are combined via logical AND" << endl
         << "unless --or is specified on the command line.  Obtain logical negation through" << endl
         << "the use of parentheses, e.g. \"! ( DP = 10 )\"" << endl
         << endl
         << "For convenience, you can specify \"QUAL\" to refer to the quality of the site, even" << endl
         << "though it does not appear in the INFO fields." << endl
         << endl;
    exit(0);
}

bool passesFilters(Variant& var,
		   vector<VariantFilter>& filters,
		   bool logicalOr, string alt = "") {
    for (vector<VariantFilter>::iterator f = filters.begin();
	 f != filters.end(); ++f) {
      string s = "";
      if (logicalOr) {
	if (alt.empty()) {
	  if (f->passes(var, s)) return true;
	} else {
	  if (f->passes(var, s, alt)) return true;
	}
      } else {
	if (alt.empty()) {
	  if (!f->passes(var, s)) return false;
	} else {
	  if (!f->passes(var, s, alt)) return false;
	}
      }
    }
    if (logicalOr)
      return false;
    else
      return true;
}


int main(int argc, char** argv) {

    int c;
    bool invert = false;
    bool logicalOr = false;
    bool filterSites = false;
    bool keepInfo = false;
    vector<string> infofilterStrs;
    vector<VariantFilter> infofilters;
    vector<string> genofilterStrs;
    vector<VariantFilter> genofilters;
    string tagPass = "";
    string tagFail = "";
    string filterSpec;
    string alleleTag;
    vector<string> regions;
    bool replaceFilter = true;

    if (argc == 1)
        printSummary(argv);

    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"help", no_argument, 0, 'h'},
                {"filter-sites", no_argument, 0, 's'},
                {"info-filter",  required_argument, 0, 'f'},
                {"genotype-filter",  required_argument, 0, 'g'},
                {"tag-pass", required_argument, 0, 't'},
                {"tag-pass", required_argument, 0, 'F'},
                {"append-filter", no_argument, 0, 'A'},
                {"allele-tag", required_argument, 0, 'a'},
                {"invert", no_argument, 0, 'v'},
                {"or", no_argument, 0, 'o'},
                {"region", required_argument, 0, 'r'},
                {"keep-info", no_argument, 0, 'k'},
                //{"length",  no_argument, &printLength, true},
                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hvAsof:g:kt:F:r:a:",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
                break;
            printf ("option %s", long_options[option_index].name);
            if (optarg)
                printf (" with arg %s", optarg);
            printf ("\n");
            break;

        case 'f':
            filterSpec += " " + string(optarg);
            infofilterStrs.push_back(string(optarg));
            break;

        case 's':
            filterSites = true;
            break;

        case 'a':
            alleleTag = optarg;
            break;

        case 'g':
            filterSpec += " genotypes filtered with: " + string(optarg);
            genofilterStrs.push_back(string(optarg));
            break;

        case 't':
            tagPass = optarg;
            break;

        case 'F':
            tagFail = optarg;
            break;

        case 'A':
            replaceFilter = false;
            break;

        case 'h':
            printSummary(argv);
            exit(0);
            break;

        case 'v':
            invert = true;
            break;

        case 'o':
            logicalOr = true;
            break;

        case 'r':
            regions.push_back(optarg);
            break;

        case 'k':
        	keepInfo = true;
        	break;

        case '?':
            /* getopt_long already printed an error message. */
            printSummary(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    if(filterSpec.empty()){
      cout << "Argument -g or -f is required" << endl;
      exit(1);
    }

    if (genofilterStrs.size() == 0 && keepInfo) {
      cout << "argument '-k' (--keep-info) requires a Genotype filter: ('-g')" << endl
	   << "i.e.: -g \"GT = 1|1\" -k" << endl;
      exit(1);
    }

    filterSpec = filterSpec.substr(1); // strip leading " "

    VariantCallFile variantFile;
    string inputFilename;
    if (optind == argc - 1) {
      inputFilename = argv[optind];
      variantFile.open(inputFilename);
    } else {
      variantFile.open(std::cin);
    }

    if (!variantFile.is_open()) {
        return 1;
    }

    for (vector<string>::iterator f = infofilterStrs.begin();
	 f != infofilterStrs.end(); ++f) {
      infofilters.push_back(VariantFilter(*f, VariantFilter::RECORD,
					  variantFile.infoTypes));
    }

    for (vector<string>::iterator f = genofilterStrs.begin();
	 f != genofilterStrs.end(); ++f) {
        genofilters.push_back(VariantFilter(*f, VariantFilter::SAMPLE,
					    variantFile.formatTypes));
    }

    vector<string> headerlines = split(variantFile.header, "\n");
    variantFile.header.clear();
    for (vector<string>::iterator l = headerlines.begin();
	 l != headerlines.end(); ++l) {
      if (!filterSpec.empty()
	  && (l->find("INFO")
	      != string::npos || l + 1 == headerlines.end())) {
	variantFile.header += "##filter=\"" + filterSpec + "\"\n";
	filterSpec.clear();
      }
      variantFile.header += *l + ((l + 1 == headerlines.end()) ? "" : "\n");
    }

    if (!tagPass.empty()) {
      variantFile.addHeaderLine("##FILTER=<ID="+ tagPass +",Description=\"Record passes the filters: " + filterSpec + ".\">");
    }

    if (!tagFail.empty()) {
      variantFile.addHeaderLine("##FILTER=<ID="+ tagFail +",Description=\"Record fails the filters: " + filterSpec + ".\">");
    }

    if (!alleleTag.empty()) {
        if (tagFail.empty()) {
            tagFail = "";
        }
        if (tagPass.empty()) {
	  tagPass = "PASS";
        }
        variantFile.addHeaderLine("##INFO=<ID="+ alleleTag +",Number=A,Type=String,Description=\"" + tagPass + " if this allele passes the filters, " + tagFail + " if not, filters are: " + filterSpec + ".\">");
    }

    cout << variantFile.header << endl;

    /*
      if (genofilters.empty()) {
      variantFile.parseSamples = false;
      }
    */

    if (filterSites) {
      variantFile.parseSamples = false;
    }

    Variant var(variantFile);

    vector<string>::iterator regionItr = regions.begin();

    do {
      if (!inputFilename.empty() && !regions.empty()) {
	string regionStr = *regionItr++;
	variantFile.setRegion(regionStr);
      }

      while (variantFile.getNextVariant(var)) {

	if (!genofilters.empty()) {
	  for (vector<VariantFilter>::iterator f
		 = genofilters.begin(); f != genofilters.end(); ++f) {
	    f->removeFilteredGenotypes(var, keepInfo);
	  }
	}
	if (!infofilters.empty()) {
	  if (filterSites) {
	    bool passes = passesFilters(var, infofilters, logicalOr);
	    if (invert) {
	      passes = !passes;
	    }
	    if (passes) {
	      if (!tagPass.empty()) {
		if (alleleTag.empty()) {
		  if (replaceFilter) {
		    var.filter.clear();
		    var.addFilter(tagPass);
		  } else {
		    var.addFilter(tagPass);
		  }
		} else {
		  var.info[alleleTag].clear();
		  for (vector<string>::iterator
			 a = var.alt.begin(); a != var.alt.end(); ++a) {
		    var.info[alleleTag].push_back(tagPass);
		  }
		}
	      } else {
		if (!var.originalLine.empty()) {
		  cout << var.originalLine << endl;
		} else {
		  cout << var << endl;
		}
	      }
	    } else {
	      if (!tagFail.empty()) {
		if (alleleTag.empty()) {
		  if (replaceFilter) {
		    var.filter.clear();
		    var.addFilter(tagFail);
		  } else {
		    var.addFilter(tagFail);
		  }
		} else {
		  var.info[alleleTag].clear();
		  for (vector<string>::iterator
			 a = var.alt.begin(); a != var.alt.end(); ++a) {
		    var.info[alleleTag].push_back(tagFail);
		  }
		}
	      }
	    }
	    if (passes && !tagPass.empty()) {
	      cout << var << endl;
	    } else if (!tagFail.empty()) {
	      cout << var << endl;
	    }
	  } else { // filter out alleles which pass
	    // removes the failing alleles
	    vector<string> failingAlts;
	    vector<string> passingAlts;
	    vector<bool> passes;
	    for (vector<string>::iterator a = var.alt.begin();
		 a != var.alt.end(); ++a) {
	      if (invert xor passesFilters(var, infofilters, logicalOr, *a)) {
		passingAlts.push_back(*a);
		passes.push_back(true);
	      } else {
		failingAlts.push_back(*a);
		passes.push_back(false);
	      }
	    }
	    if (tagPass.empty()) { // if there is no specified tag, just remove the failing alts
	      if (failingAlts.size() < var.alt.size()) {
		for (vector<string>::iterator a = failingAlts.begin(); a != failingAlts.end(); ++a) {
		  var.removeAlt(*a);
		}
		cout << var << endl;
	      }
	    } else { // otherwise, apply the tag
	      if (alleleTag.empty()) {
		if (!passingAlts.empty()) {
		  if (replaceFilter) {
		    var.filter.clear();
		    var.addFilter(tagPass);
		  } else {
		    var.addFilter(tagPass);
		  }
		} else {
		  if (replaceFilter) {
		    var.filter.clear();
		    if (!tagFail.empty()) {
		      var.addFilter(tagFail);
		    }
		  } else {
		    if (!tagFail.empty()) {
		      var.addFilter(tagFail);
		    }
		  }
		}
	      } else {
		var.info[alleleTag].clear();
		for (vector<bool>::iterator p = passes.begin();
		     p != passes.end(); ++p) {
		  if (*p) {
		    var.info[alleleTag].push_back(tagPass);
		  } else {
		    var.info[alleleTag].push_back(tagFail);
		  }
		}
	      }
	      // TODO
	      // here, if we don't use genotype filters, we shouldn't re-print the samples
	      // we haven't done anything to this part of the input.
	      cout << var << endl;
	    }
	  }
	} else {
	  if (genofilters.empty()) {
	    cout << variantFile.line << endl;
	  } else {
	    cout << var << endl;
	  }
	}
      }

    } while (regionItr != regions.end());

    return 0;

}

