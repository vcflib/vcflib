#include "Variant.h"
#include "split.h"
#include "convert.h"
#include <string>
#include <iostream>
#include <set>
#include <sys/time.h>
#include "fsom/fsom.h"
#include <getopt.h>
#include <cmath>

using namespace std;
using namespace vcflib;

double mean(const vector<double>& data) {
    double total = 0;
    for (vector<double>::const_iterator i = data.begin(); i != data.end(); ++i) {
        total += *i;
    }
    return total/data.size();
}

double median(vector <double>& data) {
    double median;
    size_t size = data.size();
    // ascending order
    sort(data.begin(), data.end());
    // get middle value
    if (size % 2 == 0) {
        median = (data[size/2-1] + data[size/2]) / 2;
    } else {
        median = data[size/2];
    }
    return median;
}

double variance(const vector <double>& data, const double mean) {
    double total = 0;
    for (vector <double>::const_iterator i = data.begin(); i != data.end(); ++i) {
        total += (*i - mean)*(*i - mean);
    }
    return total / (data.size());
}

double standard_deviation(const vector <double>& data, const double mean) {
    return sqrt(variance(data, mean));
}

struct Stats {
    double mean;
    double stdev;
    Stats(void) : mean(0), stdev(1) { }
};

bool load_som_metadata(string& som_metadata_file, int& x, int& y, vector<string>& fields, map<string, Stats>& stats) {
    ifstream in(som_metadata_file.c_str());
    if (!in.is_open()) {
        return false;
    }
    string linebuf;
    getline(in, linebuf);
    vector<string> xy = split(linebuf, "\t ");
    convert(xy.front(), x);
    convert(xy.back(), y);
    while (getline(in, linebuf)) {
        // format is: field_name, mean, stdev
        vector<string> m = split(linebuf, "\t ");
        fields.push_back(m[0]);
        Stats& s = stats[m[0]];
        convert(m[1], s.mean);
        convert(m[2], s.stdev);
    }
    in.close();
    return true;
}

bool save_som_metadata(string& som_metadata_file, int x, int y, vector<string>& fields, map<string, Stats>& stats) {
    ofstream out(som_metadata_file.c_str());
    if (!out.is_open()) {
        return false;
    }
    out << x << "\t" << y << endl;
    for (vector<string>::iterator f = fields.begin(); f != fields.end(); ++f) {
        Stats& s = stats[*f];
        out << *f << "\t" << s.mean << "\t" << s.stdev << endl;
    }
    out.close();
    return true;
}

void normalize_inputs(vector<double>& record, vector<string>& fields, map<string, Stats>& stats) {
    vector<double>::iterator r = record.begin();
    for (vector<string>::iterator f = fields.begin(); f != fields.end(); ++f, ++r) {
        Stats& s = stats[*f];
        *r = (*r - s.mean) / s.stdev;
    }
}

void read_fields(Variant& var, int ai, vector<string>& fields, vector<double>& record) {
    double td;
    vector<string>::iterator j = fields.begin();
    for (; j != fields.end(); ++j) {
        if (*j == "QUAL") { // special handling...
            td = var.quality;
        } else {
            if (var.info.find(*j) == var.info.end()) {
                td = 0;
            } else {
                if (var.vcf->infoCounts[*j] == 1) { // for non Allele-variant fields
                    convert(var.info[*j][0], td);
                } else {
                    convert(var.info[*j][ai], td);
                }
            }
        }
        record.push_back(td);
    }
}

struct SomPaint {
    int true_count;
    int false_count;
    double prob_true;
    SomPaint(void) : true_count(0), false_count(0), prob_true(0) { }
};

static unsigned long prev_uticks = 0;

static unsigned long get_uticks(){
    struct timeval ts;
    gettimeofday(&ts,0);
    return ((ts.tv_sec * 1000000) + ts.tv_usec);
}

static void start_timer(){
    prev_uticks = get_uticks();
}

static void print_timing( const char *msg ){
#define MS_DELTA (1000.0)
#define SS_DELTA (MS_DELTA * 1000.0)
#define MM_DELTA (SS_DELTA * 60.0)
#define HH_DELTA (MM_DELTA * 60.0)

    double ticks = get_uticks() - prev_uticks;

    if( ticks < MS_DELTA ){
        fprintf(stderr, "%s\t : %lf us\n", msg, ticks );
    }
    else if( ticks < SS_DELTA ){
        fprintf(stderr, "%s\t : %lf ms\n", msg, ticks / MS_DELTA );
    }
    else if( ticks < MM_DELTA ){
        fprintf(stderr, "%s\t : %lf s\n", msg, ticks / SS_DELTA );
    }
    else if( ticks < HH_DELTA ){
        fprintf(stderr, "%s\t : %lf m\n", msg, ticks / MM_DELTA );
    }
    else{
        fprintf(stderr, "%s\t : %lf h\n", msg, ticks / HH_DELTA );
    }

    start_timer();
}


void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] [vcf file]" << endl
         << endl
         << "training: " << endl
         << "    " << argv[0] << " -s output.som -x 20 -y 20 -f \"AF DP ABP\" training.vcf" << endl
         << endl
         << "application: " << endl
         << "    " << argv[0] << " -a output.som test.vcf >results.vcf" << endl
         << endl
         << argv[0] << "trains and/or applies a self-organizing map to the input VCF data" << endl
         << "on stdin, adding two columns for the x and y coordinates of the winning" << endl
         << "neuron in the network and an optional euclidean distance from a given" << endl
         << "node (--center)." << endl
         << endl
         << "If a map is provided via --apply, it will be applied to input without" << endl
         << "training.  A .meta file describing network parameters and input parameter" << endl
         << "distributions is used to automatically setup the network." << endl
         << endl
         << "options:" << endl
         << endl
         << "    -h, --help             this dialog" << endl
         << endl
         << "training:" << endl
         << endl
         << "    -f, --fields \"FIELD ...\"  INFO fields to provide to the SOM" << endl
         << "    -a, --apply FILE       apply the saved map to input data to FILE" << endl
         << "    -s, --save  FILE       train on input data and save the map to FILE" << endl
         << "    -p, --print-training-results" << endl
         << "                           print results of SOM on training input" << endl
         << "                           (you can also just use --apply on the same input)" << endl
         << "    -x, --width X          width in columns of the output array" << endl
         << "    -y, --height Y         height in columns of the output array" << endl
         << "    -i, --iterations N     number of training iterations or epochs" << endl
         << "    -d, --debug            print timing information" << endl
         << endl
         << "recalibration:" << endl
         << endl
         << "    -c, --center X,Y       annotate with euclidean distance from center" << endl
         << "    -T, --paint-true VCF   use VCF file to annotate true variants (multiple)" << endl
         << "    -F, --paint-false VCF  use VCF file to annotate false variants (multiple)" << endl
         << "    -R, --paint-tag TAG    provide estimated FDR% in TAG in variant INFO" << endl
         << "    -N, --false-negative   replace FDR% (false detection) with FNR% (false negative)" << endl;

}


int main(int argc, char** argv) {

    int width = 100;
    int height = 100;
    int num_dimensions = 2;
    int iterations = 1000;
    string som_file;
    string som_metadata_file;
    bool apply = false;
    bool train = false;
    bool apply_to_training_data = false; // print results against training data
    bool debug = false;
    vector<string> fields;
    vector<string> centerv;
    int centerx;
    int centery;
    string trueVCF;
    string falseVCF;
    bool normalize = true;

    int c;

    if (argc == 1) {
        printSummary(argv);
        exit(1);
    }

    while (true) {
        static struct option long_options[] =
        {  
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"help", no_argument, 0, 'h'},
            {"iterations", required_argument, 0, 'i'},
            {"width", required_argument, 0, 'x'},
            {"height", required_argument, 0, 'y'},
            {"apply", required_argument, 0, 'a'},
            {"save", required_argument, 0, 's'},
            {"fields", required_argument, 0, 'f'},
            {"print-training-results", no_argument, 0, 'p'},
            {"center", required_argument, 0, 'c'},
            {"paint-true", required_argument, 0, 'T'},
            {"paint-false", required_argument, 0, 'F'},
            {"debug", no_argument, 0, 'd'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hpdi:x:y:a:s:f:c:T:F:",
                         long_options, &option_index);

        if (c == -1)
            break;

        string field;

        switch (c)
        {

            case 'x':
                if (!convert(optarg, width)) {
                    cerr << "could not parse --width, -x" << endl;
                    exit(1);
                }
                break;

            case 'y':
                if (!convert(optarg, height)) {
                    cerr << "could not parse --height, -y" << endl;
                    exit(1);
                }
                break;

            case 'i':
                if (!convert(optarg, iterations)) {
                    cerr << "could not parse --iterations, -i" << endl;
                    exit(1);
                }
                break;

            case 'p':
                apply_to_training_data = true;
                break;

            case 'T':
                trueVCF = optarg;
                break;

            case 'F':
                falseVCF = optarg;
                break;

            case 'd':
                debug = true;
                break;

            case 'a':
                som_file = optarg;
                apply = true;
                break;
                
            case 's':
                som_file = optarg;
                train = true;
                break;

            case 'f':
                fields = split(string(optarg), ' ');
                break;

            case 'c':
                centerv = split(string(optarg), ',');
                convert(centerv.at(0), centerx);
                convert(centerv.at(1), centery);
                break;

            case 'h':
                printSummary(argv);
                exit(0);
                break;

            default:
                break;
        }
    }

    size_t i, j;
    som_network_t *net = NULL;
    vector<string> inputs;
    vector<vector<double> > data;
    map<string, Stats> stats;

    string line;
    stringstream ss;

    VariantCallFile variantFile;
    bool usingstdin = false;
    string inputFilename;
    if (optind == argc - 1) {
        inputFilename = argv[optind];
        variantFile.open(inputFilename);
    } else {
        variantFile.open(std::cin);
        usingstdin = true;
    }

    if (!variantFile.is_open()) {
        cerr << "could not open VCF file" << endl;
        return 1;
    }

    som_metadata_file = som_file + ".meta";

    Variant var(variantFile);

    variantFile.addHeaderLine("##INFO=<ID=SOMX,Number=A,Type=Integer,Description=\"X position of best neuron for variant in self-ordering map defined in " + som_file + "\">");
    variantFile.addHeaderLine("##INFO=<ID=SOMY,Number=A,Type=Integer,Description=\"Y position of best neuron for variant in self-ordering map defined in " + som_file + "\">");
    if (!centerv.empty()) {
        variantFile.addHeaderLine("##INFO=<ID=SOMD,Number=A,Type=Float,Description=\"Euclidean distance from "
                                  + convert(centerx) + "," + convert(centery) + " as defined by " + som_file + "\">");
    }
    if (!trueVCF.empty() && !falseVCF.empty()) {
        variantFile.addHeaderLine("##INFO=<ID=SOMP,Number=A,Type=Float,Description=\"Estimated probability the variant is true using som "
                                  + som_file + ", true variants from " + trueVCF + ", and false variants from " + falseVCF + "\">");
    }

    if (debug) start_timer();
    
    vector<Variant> variants;
    if (train) {
        map<string, pair<double, double> > normalizationLimits;
        while (variantFile.getNextVariant(var)) {
            variants.push_back(var);
            int ai = 0;
            vector<string>::iterator a = var.alt.begin();
            for ( ; a != var.alt.end(); ++a, ++ai) {
                vector<double> record;
                double td;
                vector<string>::iterator j = fields.begin();
                for (; j != fields.end(); ++j) {
                    if (*j == "QUAL") { // special handling...
                        td = var.quality;
                    } else {
                        if (var.info.find(*j) == var.info.end()) {
                            td = 0;
                        } else {
                            if (variantFile.infoCounts[*j] == 1) { // for non Allele-variant fields
                                convert(var.info[*j][0], td);
                            } else {
                                convert(var.info[*j][ai], td);
                            }
                        }
                    }
                    if (normalize) {
                        pair<double, double>& limits = normalizationLimits[*j];
                        if (td < limits.first) limits.first = td;
                        if (td > limits.second) limits.second = td;
                    }
                    record.push_back(td);
                }
                data.push_back(record);
            }
        }
        // normalize inputs
        if (normalize) {
            // get normalization vector
            // goal is normalization at 0, sd=1
            int i = 0;
            for (vector<string>::iterator f = fields.begin(); f != fields.end(); ++f, ++i) {
                vector<double> fv;
                for (vector<vector<double> >::iterator d = data.begin(); d != data.end(); ++d) {
                    fv.push_back(d->at(i));
                }
                Stats& s = stats[*f];
                // get normalization constants
                s.mean = mean(fv);
                s.stdev = standard_deviation(fv, s.mean);
                // normalize
                for (vector<vector<double> >::iterator d = data.begin(); d != data.end(); ++d) {
                    double v = d->at(i);
                    d->at(i) = (v - s.mean) / s.stdev;
                }
            }
        }
    }

    vector<double*> dataptrs (data.size());
    for (unsigned i=0, e=dataptrs.size(); i<e; ++i) {
        dataptrs[i] = &(data[i][0]); // assuming !thing[i].empty()
    }

    if (debug) print_timing( "Input Processing" );

    if (apply) {
        if (! (net = som_deserialize(som_file.c_str()))) {
            cerr << "could not load SOM from " << som_file << endl;
            return 1;
        }
        if (!fields.empty()) {
            cerr << "fields specified, but a SOM is to be applied, and metadata should be stored at " << som_metadata_file << endl;
            return 1;
        }
        if (!load_som_metadata(som_metadata_file, width, height, fields, stats)) {
            cerr << "could not load SOM metadata from " << som_metadata_file << endl;
            return 1;
        }
    } else {

        net = som_network_new(data[0].size(), height, width);
	
        if ( !net )	{
            printf( "ERROR: som_network_new failed.\n" );
            return 1;
        }
    }

    if (debug) print_timing( "Network Creation" );

    if (train) {
        if (debug) cerr << "Training using " << data.size() << " input vectors" << endl;
        som_init_weights ( net, &dataptrs[0], data.size() );
        som_train ( net, &dataptrs[0], data.size(), iterations );
    }

    if (debug) print_timing( "Network Training" );

    // open and calibrate using the true and false datasets

    if (train && apply_to_training_data) {
        // currently disabled
        /*
        cout << variantFile.header << endl;
        vector<Variant>::iterator v = variants.begin(); int di = 0;
        for ( ; v != variants.end() && di < data.size(); ++v) {
            var.info["SOMX"].clear();
            var.info["SOMY"].clear();
            var.info["SOMP"].clear();
            var.info["SOMD"].clear();
            for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a, ++di) {
                som_set_inputs ( net, dataptrs[di] );
                size_t x=0, y=0;
                som_get_best_neuron_coordinates ( net, &x, &y );
                v->info["SOMX"].push_back(convert(x));
                v->info["SOMY"].push_back(convert(y));
                if (!centerv.empty()) {
                    float distance = sqrt(pow(abs((float)centerx - (float)x), 2)
                                          + pow(abs((float)centery - (float)y), 2));
                    var.info["SOMD"].clear();
                    var.info["SOMD"].push_back(convert(distance));
                }
            }
            cout << *v << endl;
        }
        */
    } else if (apply) {

        // if we have true and false sets, use them to "paint" the map
        vector<vector<SomPaint> > paintedSOM;
        paintedSOM.resize(width);
        for (vector<vector<SomPaint> >::iterator t = paintedSOM.begin();
             t != paintedSOM.end(); ++t) {
            t->resize(height);
        }

        // handle trues
        if (!trueVCF.empty()) {
            VariantCallFile trueVariantFile;
            trueVariantFile.open(trueVCF);
            Variant v(trueVariantFile);
            while (trueVariantFile.getNextVariant(v)) {
                int ai = 0;
                vector<string>::iterator a = v.alt.begin();
                for ( ; a != v.alt.end(); ++a, ++ai) {
                    vector<double> record;
                    read_fields(v, ai, fields, record);
                    if (normalize) {
                        normalize_inputs(record, fields, stats);
                    }
                    som_set_inputs ( net, &record[0] );
                    size_t x=0, y=0;
                    som_get_best_neuron_coordinates ( net, &x, &y );
                    paintedSOM[x][y].true_count += 1;
                }
            }
        }

        // get falses
        if (!falseVCF.empty()) {
            VariantCallFile falseVariantFile;
            falseVariantFile.open(falseVCF);
            Variant v(falseVariantFile);
            while (falseVariantFile.getNextVariant(v)) {
                int ai = 0;
                vector<string>::iterator a = v.alt.begin();
                for ( ; a != v.alt.end(); ++a, ++ai) {
                    vector<double> record;
                    read_fields(v, ai, fields, record);
                    if (normalize) {
                        normalize_inputs(record, fields, stats);
                    }
                    som_set_inputs ( net, &record[0] );
                    size_t x=0, y=0;
                    som_get_best_neuron_coordinates ( net, &x, &y );
                    paintedSOM[x][y].false_count += 1;
                }
            }
        }

        // estimate probability of each node using true and false set
        for (vector<vector<SomPaint> >::iterator t = paintedSOM.begin();
             t != paintedSOM.end(); ++t) {
            for (vector<SomPaint>::iterator p = t->begin(); p != t->end(); ++p) {
                //cout << "count at node " << t - paintedSOM.begin() << "," << p - t->begin()
                //     << " is " << p->true_count << " true, " << p->false_count << " false" << endl;
                if (p->true_count + p->false_count > 0) {
                    p->prob_true = (double) p->true_count / (double) (p->true_count + p->false_count);
                } else {
                    // for nodes without training data, could we estimate from surrounding nodes?
                    // yes, TODO, but for now we can be conservative and say "0"
                    p->prob_true = 0;
                }
            }
        }

        cout << variantFile.header << endl;
        while (variantFile.getNextVariant(var)) {
            var.info["SOMX"].clear();
            var.info["SOMY"].clear();
            var.info["SOMP"].clear();
            var.info["SOMD"].clear();
            int ai = 0;
            vector<string>::iterator a = var.alt.begin();
            for ( ; a != var.alt.end(); ++a, ++ai) {
                vector<double> record;
                read_fields(var, ai, fields, record);
                if (normalize) {
                    normalize_inputs(record, fields, stats);
                }
                som_set_inputs ( net, &record[0] );
                size_t x=0, y=0;
                som_get_best_neuron_coordinates ( net, &x, &y );
                if (!trueVCF.empty() && !falseVCF.empty()) {
                    SomPaint& paint = paintedSOM[x][y];
                    var.info["SOMP"].push_back(convert(paint.prob_true));
                }
                var.info["SOMX"].push_back(convert(x));
                var.info["SOMY"].push_back(convert(y));
                if (!centerv.empty()) {
                    float distance = sqrt(pow(abs((float)centerx - (float)x), 2)
                                          + pow(abs((float)centery - (float)y), 2));
                    var.info["SOMD"].push_back(convert(distance));
                }
            }
            cout << var << endl;
        }
    }

    if (debug) print_timing( "Input Recognition" );

    if (train) {
        if (!save_som_metadata(som_metadata_file, width, height, fields, stats)) {
            cerr << "could not save metadata to " << som_metadata_file << endl;
        }
        som_serialize(net, som_file.c_str());
    }

    som_network_destroy ( net );

    if (debug) print_timing( "Network Destruction" );

    return 0;

}
