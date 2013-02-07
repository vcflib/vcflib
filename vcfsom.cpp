#include "Variant.h"
#include "split.h"
#include "convert.h"
#include <string>
#include <iostream>
#include <set>
#include <sys/time.h>
#include "fsom/fsom.h"
#include <getopt.h>

using namespace std;
using namespace vcf;


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
         << "Trains and applies a self-organizing map to the input tab-separated data on" << endl
	 << "stdin, adding two columns for the x and y coordinates of the winning neuron" << endl
	 << "in the network." << endl
	 << endl
	 << "If a map is provided, training will be skipped and the map will be applied to" << endl
	 << "the input.  Maps may be saved for later use." << endl
	 << endl
         << "options:" << endl
         << "    -h, --help             this dialog" << endl
	 << "    -f, --f \"FIELD ...\"  INFO fields to provide to the SOM" << endl
	 << "    -a, --apply FILE       apply the saved map to input data to FILE" << endl
	 << "    -s, --save  FILE       train on input data and save the map to FILE" << endl
         << "    -x, --width X          width in columns of the output array" << endl
         << "    -y, --height Y         height in columns of the output array" << endl
         << "    -i, --iterations N     number of training iterations or epochs" << endl;
}


int main(int argc, char** argv) {

    int width = 100;
    int height = 100;
    int num_dimensions = 2;
    int iterations = 1000;
    string som_file;
    bool apply = false;
    bool train = false;
    vector<string> fields;

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
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hi:x:y:a:s:f:",
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

    Variant var(variantFile);

    variantFile.addHeaderLine("##INFO=<ID=SOMX,Number=A,Type=Integer,Description=\"X position of best neuron for variant in self-ordering map defined in " + som_file + "\">");
    variantFile.addHeaderLine("##INFO=<ID=SOMY,Number=A,Type=Integer,Description=\"Y position of best neuron for variant in self-ordering map defined in " + som_file + "\">");

    start_timer();
    
    vector<Variant> variants;
    while (variantFile.getNextVariant(var)) {
	variants.push_back(var);
	int ai = 0;
	vector<string>::iterator a = var.alt.begin();
	 for ( ; a != var.alt.end(); ++a, ++ai) {
	     vector<double> record;
	     double td;
	     vector<string>::iterator j = fields.begin();
	     for (; j != fields.end(); ++j) {
		 convert(var.info[*j][ai], td);
		record.push_back(td);
	    }
	    data.push_back(record);
	}
    }

    vector<double*> dataptrs (data.size());
    for (unsigned i=0, e=dataptrs.size(); i<e; ++i) {
	dataptrs[i] = &(data[i][0]); // assuming !thing[i].empty()
    }

    print_timing( "Input Processing" );

    if (apply) {
	cerr << "Loading ... "  << endl;
	if (! (net = som_deserialize(som_file.c_str()))) {
	    cerr << "could not load SOM from " << som_file << endl;
	    return 1;
	}
    } else {

	net = som_network_new(data[0].size(), height, width);
	
	if ( !net )	{
	    printf( "ERROR: som_network_new failed.\n" );
	    return 1;
	}
    }

    print_timing( "Network Creation" );

    if (train) {
	cerr << "Training using " << data.size() << " input vectors" << endl;
	som_init_weights ( net, &dataptrs[0], data.size() );
	som_train ( net, &dataptrs[0], data.size(), iterations );
    }

    print_timing( "Network Training" );

    cout << variantFile.header << endl;

    vector<Variant>::iterator v = variants.begin(); int di = 0;
    for ( ; v != variants.end() && di < data.size(); ++v) {
	for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a, ++di) {
	    som_set_inputs ( net, dataptrs[di] );
	    size_t x=0, y=0;
	    som_get_best_neuron_coordinates ( net, &x, &y );
	    var.info["SOMX"].push_back(convert(x));
	    var.info["SOMY"].push_back(convert(y));
	}
	cout << var << endl;
    }

    print_timing( "Input Recognition" );

    if (train) {
	som_serialize(net, som_file.c_str());
    }

    som_network_destroy ( net );

    print_timing( "Network Destruction" );

    return 0;

}
