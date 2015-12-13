void advance_state ( int k );
bool antithetic_get ( );
void antithetic_memory ( int i, bool &value );
void antithetic_set ( bool value );
void cg_get ( int g, int &cg1, int &cg2 );
void cg_memory ( int i, int g, int &cg1, int &cg2 );
void cg_set ( int g, int cg1, int cg2 );
int cgn_get ( );
void cgn_memory ( int i, int &g );
void cgn_set ( int g );
void get_state ( int &cg1, int &cg2 );
int i4_uni ( );
void ig_get ( int g, int &ig1, int &ig2 );
void ig_memory ( int i, int g, int &ig1, int &ig2 );
void ig_set ( int g, int ig1, int ig2 );
void init_generator ( int t );
void initialize ( );
bool initialized_get ( );
void initialized_memory ( int i, bool &initialized );
void initialized_set ( );
void lg_get ( int g, int &lg1, int &lg2 );
void lg_memory ( int i, int g, int &lg1, int &lg2 );
void lg_set ( int g, int lg1, int lg2 );
int multmod ( int a, int s, int m );
float r4_uni_01 ( );
double r8_uni_01 ( );
void set_initial_seed ( int ig1, int ig2 );
void set_seed (  int cg1, int cg2 );
#if !defined(TIMESTAMP)
#define TIMESTAMP
void timestamp ( );
#endif
