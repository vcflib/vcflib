#include "Variant.h"
#include "split.h"
#include <getopt.h>

using namespace std;
using namespace vcflib;

void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] <vcf file>" << endl
         << endl
         << "options:" << endl 
         << "    -n, --fix-null-genotypes   only apply to null and partly-null genotypes" << endl
         << endl
         << "Set genotypes using the maximum genotype likelihood for each sample." << endl
         << endl;
    exit(0);
}


int main(int argc, char** argv) {

    bool fixNull = false;
    int c;

    while (true) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"help", no_argument, 0, 'h'},
            {"fix-null-genotypes",  no_argument, 0, 'n'},
            //{"length",  no_argument, &printLength, true},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hn",
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

          case 'n':
	      fixNull = true;
	      break;
 
          case 'h':
            printSummary(argv);
            exit(0);
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

    cout << variantFile.header << endl;

    map<pair<int, int>, list<list<int> > > glOrderCache;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        if (find(var.format.begin(), var.format.end(), "GL") == var.format.end()) {
            cout << var << endl;
            continue;
        }
        if (find(var.format.begin(), var.format.end(), "GT") == var.format.end()) {
            var.format.push_back("GT");
            reverse(var.format.begin(), var.format.end());
        }
        for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin();
             s != var.samples.end(); ++s) {
            map<string, vector<string> >& sample = s->second;
            map<string, vector<string> >::iterator g = sample.find("GT");
            map<string, vector<string> >::iterator l = sample.find("GL");
            if (l != sample.end()) {
                if (g == sample.end()) {
                    sample["GT"].push_back("./.");
                    g = sample.find("GT");
                }

                string& gt = g->second.front();
                // if we are fixing null but the genotype is fully specified, continue
                if (fixNull && gt.find(".") == string::npos) continue;
                string splitter = "/";
                if (gt.find("|") != string::npos) splitter = "|";
                int samplePloidy = split(gt, splitter).size();
                int numAlleles = var.alt.size() + 1; // including reference

                // get the gt GL ordering
                pair<int, int> pa = make_pair(samplePloidy, numAlleles);
                map<pair<int, int>, list<list<int> > >::iterator order = glOrderCache.find(pa);
                if (order == glOrderCache.end()) {
                    glOrderCache[pa] = glorder(samplePloidy, numAlleles);
                }
                list<list<int> >& glOrdering = glOrderCache[pa];

                // find the gl max
                vector<string>& gls = l->second;
                vector<string>::iterator p = gls.begin();
                double maxGl;
                convert(*p, maxGl); ++p;
                int i = 1, maxindex = 0;
                for (; p != gls.end(); ++p, ++i) {
                    double cgl;
                    convert(*p, cgl);
                    if (cgl > maxGl) {
                        maxGl = cgl;
                        maxindex = i; // prefers == gls in order of listing
                    }
                }
		
                // determine which genotype it represents
                // modify, if the GT is part-null
                vector<string>& gtv = g->second;
                list<list<int> >::iterator b = glOrdering.begin();
                advance(b, maxindex);
                /*
                  cout << "changing sample " << s->first << " gt from " << gt << " to " << join(*b, "/")
                  << " gls are ";
                  int q = 0;
                  for (list<list<int> >::iterator i = glOrdering.begin(); i != glOrdering.end(); ++i, ++q) {
                  cout << join(*i, "/") << ":" << sample["GL"].at(q) << ", ";
                  }
                  cout << endl;
                */

                gtv.clear();
                gtv.push_back(join(*b, "/"));
            }
        }
        cout << var << endl;
    }

    return 0;

}

