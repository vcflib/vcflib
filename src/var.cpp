
#include "var.hpp"

genotype::~genotype(){}

zvar::~zvar(){}

void zvar::setPopName(string  popName){
  name = popName;
}

pooled::~pooled(){}

double pooled::bound(double v){
  if(v <= 0.00001){
    return  0.00001;
  }
  if(v >= 0.99999){
    return 0.99999;
  }
  return v;
}

gl::~gl(){}

gl::gl(void){
  nalt  = 0;
  nref  = 0;
  af    = 0;
  nhomr = 0;
  nhoma = 0;
  nhet  = 0;
  ngeno = 0;
  fis   = 0;
  hfrq  = 0;

  alpha = 0.01;
  beta  = 0.01;
}

pl::~pl(){}

pl::pl(void){
  nalt  = 0;
  nref  = 0;
  af    = 0;
  nhomr = 0;
  nhoma = 0;
  nhet  = 0;
  ngeno = 0;
  fis   = 0;
  hfrq  = 0;

  alpha = 0.01;
  beta  = 0.01;

}

gp::~gp(){}

gp::gp(void){
  nalt  = 0;
  nref  = 0;
  af    = 0;
  nhomr = 0;
  nhoma = 0;
  nhet  = 0;
  ngeno = 0;
  fis   = 0;
  hfrq  = 0;

  alpha = 0.01;
  beta  = 0.01;

}

gt::~gt(){}

gt::gt(void){
  nalt  = 0;
  nref  = 0;
  af    = 0;
  nhomr = 0;
  nhoma = 0;
  nhet  = 0;
  ngeno = 0;
  fis   = 0;
  hfrq  = 0;

  alpha = 0.01;
  beta  = 0.01;

}

pooled::pooled(void){

  npop   = 0;
  afsum  = 0;
  nalt   = 0;
  nref   = 0;
  af     = 0;
  npop   = 0;
  ntot   = 0;

  alpha  = 0.01;
  beta   = 0.01;

}

// polymorphism for GL and PL

double gt::unphred(map< string, vector<string> > & geno, int index){
  return -1;
}

double gl::unphred(map< string, vector<string> > & geno, int index){
 
  double unphreded = atof(geno["GL"][index].c_str());
  return unphreded;
}

double gp::unphred(map< string, vector<string> > & geno, int index){
 
  double unphreded = atof(geno["GP"][index].c_str());
  return log(unphreded) ;
}

double pl::unphred(map< string, vector<string> > & geno, int index){

   double unphreded = atof(geno["PL"][index].c_str());

   unphreded = log(pow(10, (-unphreded/10)));
     
   return unphreded;
}

void pooled::loadPop(vector< map< string, vector<string> > > & group, string seqid, long int position){
  vector< map< string, vector<string> > >::iterator targ_it = group.begin();


  for(; targ_it != group.end(); targ_it++){

    string genotype = (*targ_it)["GT"].front();

    if(genotype == "./."){
      continue;
    }

    string allelecounts = (*targ_it)["AD"].front();

    vector<string> ac   = (*targ_it)["AD"];
    
    npop += 1;
    

    double af = atof(ac[1].c_str()) / ( atof(ac[0].c_str()) + atof(ac[1].c_str()) );

    if(atof(ac[1].c_str()) == 0){
      af = 0;
    }
    if(atof(ac[0].c_str()) == 0){
      af = 1;
    }
        
    afsum += af;
				 
    afs.push_back(af);

    nrefs.push_back(atof(ac[0].c_str()));
    nalts.push_back(atof(ac[1].c_str()));

    nref += atof(ac[0].c_str());
    ntot += atof(ac[0].c_str());
    nalt += atof(ac[1].c_str());
    ntot += atof(ac[1].c_str());
  }
  if(npop < 1){
    af = -1;
  }
  else{
    af = afsum / npop;
  }
}

void genotype::estimatePosterior(void){

  int ng = genoIndex.size();

  for(int i = 0 ; i < ng; i++){

    if(genoIndex[i] == -1){
      continue;
    }

    double aa = genoLikelihoods[i][0] ;
    double ab = genoLikelihoods[i][1] ;
    double bb = genoLikelihoods[i][2] ;

    alpha += exp(ab);
    beta  += exp(ab);

    alpha += 2 * exp(aa);
    beta  += 2 * exp(bb);
  }
  
  
}

void pooled::estimatePosterior(void){
  
  if(npop < 2){
    cerr << "FATAL: not enough pooled populations in the target or background\n";
    exit(1);
  }
  
  double xbar = af;
  double ss = 0;
  
  for(int i = 0 ; i < npop; i++){
    ss += pow(( afs[i] - xbar),2);
  }
  
  double var = (1/(npop-1))*ss;
  
  xbar = bound(xbar);
  if(var < 0.01){
    var = 0.01;
  }
  if(var < xbar*(1-xbar)){
    alpha = xbar * (((xbar*(1-xbar))/var) -1);
    beta  = (1 - xbar) * (((xbar*(1-xbar))/var) -1);
  }
  else{
    alpha = -1;
    beta  = -1;
  }
}

void genotype::loadPop( vector< map< string, vector<string> > >& group, string seqid, long int position){

  seqid = seqid;
  pos   = position  ;

  vector< map< string, vector<string> > >::iterator targ_it = group.begin();

  for(; targ_it != group.end(); targ_it++){
        
    string genotype = (*targ_it)["GT"].front();

    if(genotype == "./0" || genotype == "./1" || genotype == "."){
      genotype = "./.";
    }

    gts.push_back(genotype);

    vector<double> phreds;
    vector<double> phredsCDF;

    double sum  = 0;
    if(genotype != "./."){

      double pa  =  unphred( (*targ_it), 0)  ;
      double pab =  unphred( (*targ_it), 1)  ;
      double pbb =  unphred( (*targ_it), 2)  ;

      double norm = log(exp(pa) + exp(pab) + exp(pbb))  ;

      phreds.push_back(pa  - norm);
      phreds.push_back(pab - norm);
      phreds.push_back(pbb - norm);

    
      sum += exp(pa - norm);
      //      cerr << sum << endl;
      phredsCDF.push_back(sum);
      sum += exp(pab - norm);
      //      cerr << sum << endl;
      phredsCDF.push_back(sum);
      sum += exp(pbb - norm);
      //      cerr << sum << endl;
      //      cerr << endl;
      phredsCDF.push_back(sum);
    }
    else{
      phreds.push_back(log(1/3));
      phreds.push_back(log(1/3));
      phreds.push_back(log(1/3));
      phredsCDF.push_back(1/3);
      phredsCDF.push_back(2/3);
      phredsCDF.push_back(1);
    }

    genoLikelihoods.push_back(phreds);
    genoLikelihoodsCDF.push_back(phredsCDF);

    
    while(1){
      if(genotype == "./."){
        genoIndex.push_back(-1);
        break;
      }
      if(genotype == "0/0"){
        ngeno += 1;
        nhomr += 1;
        nref  += 2;
        genoIndex.push_back(0);
        break;
      }
      if(genotype == "0/1"){
        ngeno += 1;
        nhet  += 1;
        nref  += 1;
        nalt  += 1;
        genoIndex.push_back(1);
        break;
      }
      if(genotype == "1/0"){
        ngeno += 1;
        nhet  += 1;
	nref  += 1;
        nalt  += 1;
        genoIndex.push_back(1);
        break;
      }
      if(genotype == "1/1"){
	ngeno += 1;
        nhoma += 1;
        nalt  += 2;
        genoIndex.push_back(2);
        break;
      }
      if(genotype == "0|0"){
        ngeno += 1;
        nhomr += 1;
        nref  += 2;
        genoIndex.push_back(0);
        break;
      }
      if(genotype == "0|1"){
        ngeno += 1;
        nhet  += 1;
        nref  += 1;
        nalt  += 1;
        genoIndex.push_back(1);
        break;
      }
      if(genotype == "1|0"){
        ngeno += 1;
        nhet  += 1;
        nref  += 1;
        nalt  += 1;
        genoIndex.push_back(1);
        break;
      }
      if(genotype == "1|1"){
        ngeno += 1;
        nhoma += 1;
        nalt  += 2;
        genoIndex.push_back(2);
        break;
      }
      cerr << "FATAL: unknown genotype: " << genotype << endl;
      exit(1);
    }
  }
  if(nalt == 0 && nref == 0){
    af = -1;
  }
  else{
    af  = (nalt / (nref + nalt));

    if(nhet > 0){
      fis = ( 1 - ((nhet/ngeno) / (2*af*(1 - af))));
    }
    else{
      fis = 1;
    }
    if(fis < 0){
      fis = 0.00001;
    }
  }
  hfrq = nhet / ngeno;
  npop = ngeno;
}
