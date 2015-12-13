# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "pdflib.hpp"
# include "rnglib.hpp"

//****************************************************************************80

double i4_binomial_pdf ( int n, double p, int k )

//****************************************************************************80
//
//  Purpose:
//
//    I4_BINOMIAL_PDF evaluates the binomial PDF.
//
//  Discussion:
//
//    pdf(n,p,k) = C(n,k) p^k (1-p)^(n-k)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 June 2013
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, int N, the number of binomial trials.
//    0 < N.
//
//    Input, double P, the probability of a success in one trial.
//
//    Input, int K, the number of successes.
//
//    Output, double I4_BINOMIAL_PDF, the probability of K successes
//    in N trials with a per-trial success probability of P.
//
{
  double value;

  if ( k < 0 )
  {
    value = 0.0;
  }
  else if ( k <= n )
  {
    value = r8_choose ( n, k ) * pow ( p, k ) * pow ( 1.0 - p, k );
  }
  else
  {
    value = 0.0;
  }
  return value;
}
//****************************************************************************80

int i4_binomial_sample ( int n, double pp )

//****************************************************************************80
//
//  Purpose:
//
//    I4_BINOMIAL_SAMPLE generates a binomial random deviate.
//
//  Discussion:
//
//    This procedure generates a single random deviate from a binomial
//    distribution whose number of trials is N and whose
//    probability of an event in each trial is P.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN77 version by Barry Brown, James Lovato.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Voratas Kachitvichyanukul, Bruce Schmeiser,
//    Binomial Random Variate Generation,
//    Communications of the ACM,
//    Volume 31, Number 2, February 1988, pages 216-222.
//
//  Parameters:
//
//    Input, int N, the number of binomial trials, from which a
//    random deviate will be generated.
//    0 < N.
//
//    Input, double PP, the probability of an event in each trial of
//    the binomial distribution from which a random deviate is to be generated.
//    0.0 < PP < 1.0.
//
//    Output, int I4_BINOMIAL_SAMPLE, a random deviate from the
//    distribution.
//
{
  double al;
  double alv;
  double amaxp;
  double c;
  double f;
  double f1;
  double f2;
  double ffm;
  double fm;
  double g;
  int i;
  int ix;
  int ix1;
  int k;
  int m;
  int mp;
  double p;
  double p1;
  double p2;
  double p3;
  double p4;
  double q;
  double qn;
  double r;
  double t;
  double u;
  double v;
  int value;
  double w;
  double w2;
  double x;
  double x1;
  double x2;
  double xl;
  double xll;
  double xlr;
  double xm;
  double xnp;
  double xnpq;
  double xr;
  double ynorm;
  double z;
  double z2;

  if ( pp <= 0.0 || 1.0 <= pp )
  {
    cerr << "\n";
    cerr << "I4_BINOMIAL_SAMPLE - Fatal error!\n";
    cerr << "  PP is out of range.\n";
    exit ( 1 );
  }

  p = r8_min ( pp, 1.0 - pp );
  q = 1.0 - p;
  xnp = ( double ) ( n ) * p;

  if ( xnp < 30.0 )
  {
    qn = pow ( q, n );
    r = p / q;
    g = r * ( double ) ( n + 1 );

    for ( ; ; )
    {
      ix = 0;
      f = qn;
      u = r8_uniform_01_sample ( );

      for ( ; ; )
      {
        if ( u < f )
        {
          if ( 0.5 < pp )
          {
            ix = n - ix;
          }
          value = ix;
          return value;
        }

        if ( 110 < ix )
        {
          break;
        }
        u = u - f;
        ix = ix + 1;
        f = f * ( g / ( double ) ( ix ) - r );
      }
    }
  }
  ffm = xnp + p;
  m = ffm;
  fm = m;
  xnpq = xnp * q;
  p1 = ( int ) ( 2.195 * sqrt ( xnpq ) - 4.6 * q ) + 0.5;
  xm = fm + 0.5;
  xl = xm - p1;
  xr = xm + p1;
  c = 0.134 + 20.5 / ( 15.3 + fm );
  al = ( ffm - xl ) / ( ffm - xl * p );
  xll = al * ( 1.0 + 0.5 * al );
  al = ( xr - ffm ) / ( xr * q );
  xlr = al * ( 1.0 + 0.5 * al );
  p2 = p1 * ( 1.0 + c + c );
  p3 = p2 + c / xll;
  p4 = p3 + c / xlr;
//
//  Generate a variate.
//
  for ( ; ; )
  {
    u = r8_uniform_01_sample ( ) * p4;
    v = r8_uniform_01_sample ( );
//
//  Triangle
//
    if ( u < p1 )
    {
      ix = xm - p1 * v + u;
      if ( 0.5 < pp ) 
      {
        ix = n - ix;
      }
      value = ix;
      return value;
    }
//
//  Parallelogram
//
    if ( u <= p2 )
    {
      x = xl + ( u - p1 ) / c;
      v = v * c + 1.0 - fabs ( xm - x ) / p1;

      if ( v <= 0.0 || 1.0 < v )
      {
        continue;
      }
      ix = x;
    }
    else if ( u <= p3 )
    {
      ix = xl + log ( v ) / xll;
      if ( ix < 0 )
      {
        continue;
      }
      v = v * ( u - p2 ) * xll;
    }
    else
    {
      ix = xr - log ( v ) / xlr;
      if ( n < ix )
      {
        continue;
      }
      v = v * ( u - p3 ) * xlr;
    }
    k = abs ( ix - m );

    if ( k <= 20 || xnpq / 2.0 - 1.0 <= k )
    {
      f = 1.0;
      r = p / q;
      g = ( n + 1 ) * r;

      if ( m < ix )
      {
        mp = m + 1;
        for ( i = mp; i <= ix; i++ )
        {
          f = f * ( g / i - r );
        }
      }
      else if ( ix < m )
      {
        ix1 = ix + 1;
        for ( i = ix1; i <= m; i++ )
        {
          f = f / ( g / i - r );
        }
      }

      if ( v <= f )
      {
        if ( 0.5 < pp )
        {
          ix = n - ix;
        }
        value = ix;
        return value;
      }
    }
    else
    {
      amaxp = ( k / xnpq ) * ( ( k * ( k / 3.0 
        + 0.625 ) + 0.1666666666666 ) / xnpq + 0.5 );
      ynorm = - ( double ) ( k * k ) / ( 2.0 * xnpq );
      alv = log ( v );

      if ( alv < ynorm - amaxp )
      {
        if ( 0.5 < pp )
        {
          ix = n - ix;
        }
        value = ix;
        return value;
      }

      if ( ynorm + amaxp < alv )
      {
        continue;
      }

      x1 = ( double ) ( ix + 1 );
      f1 = fm + 1.0;
      z = ( double ) ( n + 1 ) - fm;
      w = ( double ) ( n - ix + 1 );
      z2 = z * z;
      x2 = x1 * x1;
      f2 = f1 * f1;
      w2 = w * w;

      t = xm * log ( f1 / x1 ) + ( n - m + 0.5 ) * log ( z / w ) 
        + ( double ) ( ix - m ) * log ( w * p / ( x1 * q )) 
        + ( 13860.0 - ( 462.0 - ( 132.0 - ( 99.0 - 140.0 
        / f2 ) / f2 ) / f2 ) / f2 ) / f1 / 166320.0 
        + ( 13860.0 - ( 462.0 - ( 132.0 - ( 99.0 - 140.0 
        / z2 ) / z2 ) / z2 ) / z2 ) / z / 166320.0 
        + ( 13860.0 - ( 462.0 - ( 132.0 - ( 99.0 - 140.0 
        / x2 ) / x2 ) / x2 ) / x2 ) / x1 / 166320.0 
        + ( 13860.0 - ( 462.0 - ( 132.0 - ( 99.0 - 140.0 
        / w2 ) / w2 ) / w2 ) / w2 ) / w / 166320.0;

      if ( alv <= t )
      {
        if ( 0.5 < pp )
        {
          ix = n - ix;
        }
        value = ix;
        return value;
      }
    }
  }
  return value;
}
//****************************************************************************80

double i4vec_multinomial_pdf ( int n, double p[], int m, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MULTINOMIAL_PDF evaluates the multinomial PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 June 2013
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, int N, the number of trials.
//
//    Input, double P[M], the probability of each outcome
//    on any single trial.
//
//    Input, int M, the number of possible outcomes
//    of a single trial.
//
//    Input, int X[M], the results of N trials,
//    with X(I) the number of times outcome I occurred.
//
//    Output, double I4VEC_MULTINOMIAL_PDF, the probability
//    density function evaluated at X.
//
{
  int bot;
  int c;
  int i;
  int j;
  double pdf;
  int top;
//
//  The combinatorial coefficient is an integer.
//
  c = 1;
  top = n;
  for ( i = 0; i < m; i++ )
  {
    bot = 1;
    for ( j = 0; j < x[i]; j++ )
    {
      c = ( c * top ) / bot;
      top = top - 1;
      bot = bot + 1;
    }
  }

  pdf = ( double ) ( c );
  for ( i = 0; i < m; i++ )
  {
    pdf = pdf * pow ( p[i], x[i] );
  }

  return pdf;
}
//****************************************************************************80

int *i4vec_multinomial_sample ( int n, double p[], int ncat )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MULTINOMIAL_SAMPLE generates a multinomial random deviate.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN77 version by Barry Brown, James Lovato.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Luc Devroye,
//    Non-Uniform Random Variate Generation,
//    Springer, 1986,
//    ISBN: 0387963057,
//    LC: QA274.D48.
//
//  Parameters:
//
//    Input, int N, the number of events, which will be
//    classified into one of the NCAT categories.
//
//    Input, double P[NCAT-1].  P(I) is the probability that an event
//    will be classified into category I.  Thus, each P(I) must be between 
//    0.0 and 1.0.  Only the first NCAT-1 values of P must be defined since 
//    P(NCAT) would be 1.0 minus the sum of the first NCAT-1 P's.
//
//    Input, int NCAT, the number of categories.
//
//    Output, int I4VEC_MULTINOMIAL_SAMPLE[NCAT], a random observation from 
//    the multinomial distribution.  All IX(i) will be nonnegative and their 
//    sum will be N.
//
{
  int i;
  int icat;
  int *ix;
  int ntot;
  double prob;
  double ptot;

  if ( n < 0 )
  {
    cerr << "\n";
    cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
    cerr << "  N < 0\n";
    exit ( 1 );
  }

  if ( ncat <= 1 )
  {
    cerr << "\n";
    cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
    cerr << "  NCAT <= 1\n";
    exit ( 1 );
  }

  for ( i = 0; i < ncat - 1; i++ )
  {
    if ( p[i] < 0.0 )
    {
      cerr << "\n";
      cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
      cerr << "  Some P(i) < 0.\n";
      exit ( 1 );
    }

    if ( 1.0 < p[i] )
    {
      cerr << "\n";
      cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
      cerr << "  Some 1 < P(i).\n";
      exit ( 1 );
    }
  }

  ptot = 0.0;
  for ( i = 0; i < ncat - 1; i++ )
  {
    ptot = ptot + p[i];
  }

  if ( 0.99999 < ptot ) 
  {
    cerr << "\n";
    cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
    cerr << "  1.0 < Sum of P().\n";
    exit ( 1 );
  }
//
//  Initialize variables.
//
  ntot = n;
  ptot = 1.0;

  ix = new int[ncat];
  for ( i = 0; i < ncat; i++ )
  {
    ix[i] = 0;
  }
//
//  Generate the observation.
//
  for ( icat = 0; icat < ncat - 1; icat++ )
  {
    prob = p[icat] / ptot;
    ix[icat] = i4_binomial_sample ( ntot, prob );
    ntot = ntot - ix[icat];
    if ( ntot <= 0 )
    {
      return ix;
    }
    ptot = ptot - p[icat];
  }

  ix[ncat-1] = ntot;

  return ix;
}
//****************************************************************************80

double r8_beta_pdf ( double alpha, double beta, double rval )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BETA_PDF evaluates the PDF of a beta distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double ALPHA, BETA, shape parameters.
//    0.0 < ALPHA, BETA.
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_BETA_PDF, the value of the PDF at RVAL.
//
{
  double temp;
  double value;

  if ( alpha <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_BETA_PDF - Fatal error!\n";
    cerr << "  Parameter ALPHA is not positive.\n";
    exit ( 1 );
  }

  if ( beta <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_BETA_PDF- Fatal error!\n";
    cerr << "  Parameter BETA is not positive.\n";
    exit ( 1 );
  }

  if ( rval < 0.0 || 1.0 < rval )
  {
    value = 0.0;
  }
  else
  {
    temp = r8_gamma_log ( alpha + beta ) - r8_gamma_log ( alpha ) 
      - r8_gamma_log ( beta );

    value = exp ( temp ) * pow ( rval, alpha - 1.0 )
      * pow ( 1.0 - rval, beta - 1.0 );
  }
  return value;
}
//****************************************************************************80

double r8_beta_sample ( double aa, double bb )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BETA_SAMPLE generates a beta random deviate.
//
//  Discussion:
//
//    This procedure returns a single random deviate from the beta distribution
//    with parameters A and B.  The density is
//
//      x^(a-1) * (1-x)^(b-1) / Beta(a,b) for 0 < x < 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN77 version by Barry Brown, James Lovato.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Russell Cheng,
//    Generating Beta Variates with Nonintegral Shape Parameters,
//    Communications of the ACM,
//    Volume 21, Number 4, April 1978, pages 317-322.
//
//  Parameters:
//
//    Input, double AA, the first parameter of the beta distribution.
//    0.0 < AA.
//
//    Input, double BB, the second parameter of the beta distribution.
//    0.0 < BB.
//
//    Output, double R8_BETA_SAMPLE, a beta random variate.
//
{
  double a;
  double alpha;
  double b;
  double beta;
  double delta;
  double gamma;
  double k1;
  double k2;
  const double log4 = 1.3862943611198906188;
  const double log5 = 1.6094379124341003746;
  double r;
  double s;
  double t;
  double u1;
  double u2;
  double v;
  double value;
  double w;
  double y;
  double z;

  if ( aa <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_BETA_SAMPLE - Fatal error!\n";
    cerr << "  AA <= 0.0\n";
    exit ( 1 );
  }

  if ( bb <= 0.0 ) 
  {
    cerr << "\n";
    cerr << "R8_BETA_SAMPLE - Fatal error!\n";
    cerr << "  BB <= 0.0\n";
    exit ( 1 );
  }
//
//  Algorithm BB
//
  if ( 1.0 < aa && 1.0 < bb )
  {
    if ( aa < bb )
    {
      a = aa;
      b = bb;
    }
    else
    {
      a = bb;
      b = aa;
    }
    alpha = a + b;
    beta = sqrt ( ( alpha - 2.0 ) / ( 2.0 * a * b - alpha ) );
    gamma = a + 1.0 / beta;

    for ( ; ; )
    {
      u1 = r8_uniform_01_sample ( );
      u2 = r8_uniform_01_sample ( );
      v = beta * log ( u1 / ( 1.0 - u1 ) );
      w = a * exp ( v );

      z = u1 * u1 * u2;
      r = gamma * v - log4;
      s = a + r - w;

      if ( 5.0 * z <= s + 1.0 + log5 )
      {
        break;
      }

      t = log ( z );
      if ( t <= s )
      {
        break;
      }

      if ( t <= ( r + alpha * log ( alpha / ( b + w ) ) ) )
      {
        break;
      }
    }
  }
//
//  Algorithm BC
//
  else
  {
    if ( aa < bb )
    {
      a = bb;
      b = aa;
    }
    else
    {
      a = aa;
      b = bb;
    }
    alpha = a + b;
    beta = 1.0 / b;
    delta = 1.0 + a - b;
    k1 = delta * ( 1.0 / 72.0 + b / 24.0 ) 
      / ( a / b - 7.0 / 9.0 );
    k2 = 0.25 + ( 0.5 + 0.25 / delta ) * b;

    for ( ; ; )
    {
      u1 = r8_uniform_01_sample ( );
      u2 = r8_uniform_01_sample ( );

      if ( u1 < 0.5 )
      {
        y = u1 * u2;
        z = u1 * y;

        if ( k1 <= 0.25 * u2 + z - y )
        {
          continue;
        }
      }
      else
      {
        z = u1 * u1 * u2;

        if ( z <= 0.25 )
        {
          v = beta * log ( u1 / ( 1.0 - u1 ) );
          w = a * exp ( v );

          if ( aa == a )
          {
            value = w / ( b + w );
          }
          else
          {
            value = b / ( b + w );
          }
          return value;
        }

        if ( k2 < z )
        {
          continue;
        }
      }

      v = beta * log ( u1 / ( 1.0 - u1 ) );
      w = a * exp ( v );

      if ( log ( z ) <= alpha * ( log ( alpha / ( b + w ) ) + v ) - log4 )
      {
        break;
      }
    }
  }

  if ( aa == a )
  {
    value = w / ( b + w );
  }
  else
  {
    value = b / ( b + w );
  }
  return value;
}
//****************************************************************************80

double r8_chi_pdf ( double df, double rval )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CHI_PDF evaluates the PDF of a chi-squared distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C  by John Burkardt.
//
//  Parameters:
//
//    Input, double DF, the degrees of freedom.
//    0.0 < DF.
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_CHI_PDF, the value of the PDF at RVAL.
//
{
  double temp1;
  double temp2;
  double value;

  if ( df <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_CHI_PDF - Fatal error!\n";
    cerr << "  Degrees of freedom must be positive.\n";
    exit ( 1 );
  }
      
  if ( rval <= 0.0 )
  {
    value = 0.0;
  }
  else
  {
    temp2 = df * 0.5;

    temp1 = ( temp2 - 1.0 ) * log ( rval ) - 0.5 * rval 
      - temp2 * log ( 2.0 ) - r8_gamma_log ( temp2 );

    value = exp ( temp1 );
  }
  return value;
}
//****************************************************************************80

double r8_chi_sample ( double df )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CHI_SAMPLE generates a Chi-Square random deviate.
//
//  Discussion:
//
//    This procedure generates a random deviate from the chi square distribution
//    with DF degrees of freedom random variable.
//
//    The algorithm exploits the relation between chisquare and gamma.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN77 version by Barry Brown, James Lovato.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double DF, the degrees of freedom.
//    0.0 < DF.
//
//    Output, double R8_CHI_SAMPLE, a random deviate from the distribution.
//
{
  double arg1;
  double arg2;
  double value;

  if ( df <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_CHI_SAMPLE - Fatal error!\n";
    cerr << "  DF <= 0.\n";
    cerr << "  Value of DF: " << df << "\n";
    exit ( 1 );
  }

  arg1 = 1.0;
  arg2 = df / 2.0;

  value = 2.0 * r8_gamma_sample ( arg1, arg2 );

  return value;
}
//****************************************************************************80

double r8_choose ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
//
//  Discussion:
//
//    The value is calculated in such a way as to avoid overflow and
//    roundoff.  The calculation is done in R8 arithmetic.
//
//    The formula used is:
//
//      C(N,K) = N! / ( K! * (N-K)! )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 June 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    ML Wolfson, HV Wright,
//    Algorithm 160:
//    Combinatorial of M Things Taken N at a Time,
//    Communications of the ACM,
//    Volume 6, Number 4, April 1963, page 161.
//
//  Parameters:
//
//    Input, int N, K, the values of N and K.
//
//    Output, double R8_CHOOSE, the number of combinations of N
//    things taken K at a time.
//
{
  int i;
  int mn;
  int mx;
  double value;

  if ( k < n - k )
  {
    mn = k;
    mx = n - k;
  }
  else
  {
    mn = n - k;
    mx = k;
  }

  if ( mn < 0 )
  {
    value = 0.0;
  }
  else if ( mn == 0 )
  {
    value = 1.0;
  }
  else
  {
    value = ( double ) ( mx + 1 );

    for ( i = 2; i <= mn; i++ )
    {
      value = ( value * ( double ) ( mx + i ) ) / ( double ) i;
    }
  }

  return value;
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  const double value = 2.220446049250313E-016;

  return value;
}
//****************************************************************************80

double r8_exponential_pdf ( double beta, double rval )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EXPONENTIAL_PDF evaluates the PDF of an exponential distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double BETA, the scale value.
//    0.0 < BETA.
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_EXPONENTIAL_PDF, the value of the PDF at RVAL.
//
{
  double value;

  if ( beta <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_EXPONENTIAL_PDF - Fatal error!\n";
    cerr << "  BETA parameter must be positive.\n";
    exit ( 1 );
  }

  if ( rval < 0.0 )
  {
    value = 0.0;
  }
  else
  {
    value = exp ( - rval / beta ) / beta;
  }

  return value;
}
//****************************************************************************80

double r8_exponential_sample ( double lambda )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EXPONENTIAL_SAMPLE samples the exponential PDF.
//
//  Discussion:
//
//    Note that the parameter LAMBDA is a multiplier.  In some formulations,
//    it is used as a divisor instead.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double LAMBDA, the parameter of the PDF.
//
//    Output, double R8_EXPONENTIAL_SAMPLE, a sample of the PDF.
//
{
  double r;
  double value;

  r = r8_uniform_01_sample ( );

  value = - log ( r ) * lambda;

  return value;
}
//****************************************************************************80

double r8_exponential_01_pdf ( double rval )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EXPONENTIAL_01_PDF: PDF of the standard exponential distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 June 2013
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_EXPONENTIAL_01_PDF, the value of the PDF at RVAL.
//
{
  double value;

  if ( rval < 0.0 )
  {
    value = 0.0;
  }
  else
  {
    value = exp ( - rval );
  }

  return value;
}
//****************************************************************************80

double r8_exponential_01_sample ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EXPONENTIAL_01_SAMPLE samples the standard exponential PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EXPONENTIAL_01_SAMPLE, a sample of the PDF.
//
{
  double r;
  double value;

  r = r8_uniform_01_sample ( );

  value = - log ( r );

  return value;
}
//****************************************************************************80

double r8_gamma_log ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA_LOG evaluates the logarithm of the gamma function.
//
//  Discussion:
//
//    This routine calculates the LOG(GAMMA) function for a positive real
//    argument X.  Computation is based on an algorithm outlined in
//    references 1 and 2.  The program uses rational functions that
//    theoretically approximate LOG(GAMMA) to at least 18 significant
//    decimal digits.  The approximation for X > 12 is from reference
//    3, while approximations for X < 12.0 are similar to those in
//    reference 1, but are unpublished.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 April 2013
//
//  Author:
//
//    Original FORTRAN77 version by William Cody, Laura Stoltz.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody, Kenneth Hillstrom,
//    Chebyshev Approximations for the Natural Logarithm of the
//    Gamma Function,
//    Mathematics of Computation,
//    Volume 21, Number 98, April 1967, pages 198-203.
//
//    Kenneth Hillstrom,
//    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
//    May 1969.
//
//    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
//    Charles Mesztenyi, John Rice, Henry Thatcher,
//    Christoph Witzgall,
//    Computer Approximations,
//    Wiley, 1968,
//    LC: QA297.C64.
//
//  Parameters:
//
//    Input, double X, the argument of the function.
//
//    Output, double R8_GAMMA_LOG, the value of the function.
//
{
  double c[7] = {
    -1.910444077728E-03, 
     8.4171387781295E-04, 
    -5.952379913043012E-04,
     7.93650793500350248E-04, 
    -2.777777777777681622553E-03, 
     8.333333333333333331554247E-02, 
     5.7083835261E-03 };
  double corr;
  const double d1 = -5.772156649015328605195174E-01;
  const double d2 = 4.227843350984671393993777E-01;
  const double d4 = 1.791759469228055000094023;
  const double frtbig = 2.25E+76;
  int i;
  double p1[8] = {
    4.945235359296727046734888, 
    2.018112620856775083915565E+02, 
    2.290838373831346393026739E+03, 
    1.131967205903380828685045E+04, 
    2.855724635671635335736389E+04, 
    3.848496228443793359990269E+04, 
    2.637748787624195437963534E+04, 
    7.225813979700288197698961E+03 };
  double p2[8] = { 
    4.974607845568932035012064, 
    5.424138599891070494101986E+02, 
    1.550693864978364947665077E+04, 
    1.847932904445632425417223E+05, 
    1.088204769468828767498470E+06, 
    3.338152967987029735917223E+06, 
    5.106661678927352456275255E+06, 
    3.074109054850539556250927E+06 };
  double p4[8] = {
    1.474502166059939948905062E+04, 
    2.426813369486704502836312E+06, 
    1.214755574045093227939592E+08, 
    2.663432449630976949898078E+09, 
    2.940378956634553899906876E+10, 
    1.702665737765398868392998E+11, 
    4.926125793377430887588120E+11, 
    5.606251856223951465078242E+11 };
  double q1[8] = { 
    6.748212550303777196073036E+01, 
    1.113332393857199323513008E+03, 
    7.738757056935398733233834E+03, 
    2.763987074403340708898585E+04, 
    5.499310206226157329794414E+04, 
    6.161122180066002127833352E+04, 
    3.635127591501940507276287E+04, 
    8.785536302431013170870835E+03 };
  double q2[8] = { 
    1.830328399370592604055942E+02, 
    7.765049321445005871323047E+03, 
    1.331903827966074194402448E+05, 
    1.136705821321969608938755E+06, 
    5.267964117437946917577538E+06, 
    1.346701454311101692290052E+07, 
    1.782736530353274213975932E+07, 
    9.533095591844353613395747E+06 };
  double q4[8] = { 
    2.690530175870899333379843E+03, 
    6.393885654300092398984238E+05, 
    4.135599930241388052042842E+07, 
    1.120872109616147941376570E+09, 
    1.488613728678813811542398E+10, 
    1.016803586272438228077304E+11, 
    3.417476345507377132798597E+11, 
    4.463158187419713286462081E+11 };
  double res;
  const double sqrtpi = 0.9189385332046727417803297;
  const double xbig = 2.55E+305;
  double xden;
  const double xinf = 1.79E+308;
  double xm1;
  double xm2;
  double xm4;
  double xnum;
  double y;
  double ysq;

  y = x;

  if ( 0.0 < y && y <= xbig )
  {
    if ( y <= r8_epsilon ( ) )
    {
      res = - log ( y );
    }
//
//  EPS < X <= 1.5.
//
    else if ( y <= 1.5 )
    {
      if ( y < 0.6796875 )
      {
        corr = -log ( y );
        xm1 = y;
      }
      else
      {
        corr = 0.0;
        xm1 = ( y - 0.5 ) - 0.5;
      }

      if ( y <= 0.5 || 0.6796875 <= y )
      {
        xden = 1.0;
        xnum = 0.0;
        for ( i = 0; i < 8; i++ )
        {
          xnum = xnum * xm1 + p1[i];
          xden = xden * xm1 + q1[i];
        }
        res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) );
      }
      else
      {
        xm2 = ( y - 0.5 ) - 0.5;
        xden = 1.0;
        xnum = 0.0;
        for ( i = 0; i < 8; i++ )
        {
          xnum = xnum * xm2 + p2[i];
          xden = xden * xm2 + q2[i];
        }
        res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) );
      }
    }
//
//  1.5 < X <= 4.0.
//
    else if ( y <= 4.0 )
    {
      xm2 = y - 2.0;
      xden = 1.0;
      xnum = 0.0;
      for ( i = 0; i < 8; i++ )
      {
        xnum = xnum * xm2 + p2[i];
        xden = xden * xm2 + q2[i];
      }
      res = xm2 * ( d2 + xm2 * ( xnum / xden ) );
    }
//
//  4.0 < X <= 12.0.
//
    else if ( y <= 12.0 )
    {
      xm4 = y - 4.0;
      xden = -1.0;
      xnum = 0.0;
      for ( i = 0; i < 8; i++ )
      {
        xnum = xnum * xm4 + p4[i];
        xden = xden * xm4 + q4[i];
      }
      res = d4 + xm4 * ( xnum / xden );
    }
//
//  Evaluate for 12 <= argument.
//
    else
    {
      res = 0.0;

      if ( y <= frtbig )
      {
        res = c[6];
        ysq = y * y;
        for ( i = 0; i < 6; i++ )
        {
          res = res / ysq + c[i];
        }
      }
      res = res / y;
      corr = log ( y );
      res = res + sqrtpi - 0.5 * corr;
      res = res + y * ( corr - 1.0 );
    }
  }
//
//  Return for bad arguments.
//
  else
  {
    res = xinf;
  }

  return res;
}
//****************************************************************************80

double r8_gamma_pdf ( double beta, double alpha, double rval )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA_PDF evaluates the PDF of a gamma distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double BETA, the rate parameter.
//    0.0 < BETA.
//
//    Input, double ALPHA, the shape parameter.
//    0.0 < ALPHA.
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_GAMMA_PDF, the value of the PDF at RVAL.
//
{
  double temp;
  double value;

  if ( alpha <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_GAMMA_PDF - Fatal error!\n";
    cerr << "  Parameter ALPHA is not positive.\n";
    exit ( 1 );
  }

  if ( beta <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_GAMMA_PDF - Fatal error!\n";
    cerr << "  Parameter BETA is not positive.\n";
    exit ( 1 );
  }

  if ( rval <= 0.0 )
  {
    value = 0.0;
  }
  else
  {
    temp = alpha * log ( beta ) + ( alpha - 1.0 ) * log ( rval ) 
      - beta * rval - r8_gamma_log ( alpha );

    value = exp ( temp );
  }

  return value;
}
//****************************************************************************80

double r8_gamma_sample ( double a, double r )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA_SAMPLE generates a Gamma random deviate.
//
//  Discussion:
//
//    This procedure generates random deviates from the gamma distribution whose
//    density is (A^R)/Gamma(R) * X^(R-1) * Exp(-A*X)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN77 version by Barry Brown, James Lovato.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Joachim Ahrens, Ulrich Dieter,
//    Generating Gamma Variates by a Modified Rejection Technique,
//    Communications of the ACM,
//    Volume 25, Number 1, January 1982, pages 47-54.
//
//    Joachim Ahrens, Ulrich Dieter,
//    Computer Methods for Sampling from Gamma, Beta, Poisson and
//    Binomial Distributions,
//    Computing,
//    Volume 12, Number 3, September 1974, pages 223-246.
//
//  Parameters:
//
//    Input, double A, the location parameter.
//
//    Input, double R, the shape parameter.
//
//    Output, double R8_GAMMA_SAMPLE, a random deviate from the distribution.
//
{
  double value;

  value = r8_gamma_01_sample ( r ) / a;

  return value;
}
//****************************************************************************80

double r8_gamma_01_pdf ( double alpha, double rval )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA_01_PDF evaluates the PDF of a standard gamma distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 June 2013
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, double ALPHA, the shape parameter.
//    0.0 < ALPHA.
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_GAMMA_01_PDF, the value of the PDF at RVAL.
//
{
  double temp;
  double value;

  if ( alpha <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_GAMMA_01_PDF - Fatal error!\n";
    cerr << "  Parameter ALPHA is not positive.\n";
    exit ( 1 );
  }

  if ( rval <= 0.0 )
  {
    value = 0.0;
  }
  else
  {
    temp = ( alpha - 1.0 ) * log ( rval ) - rval - r8_gamma_log ( alpha );

    value = exp ( temp );
  }

  return value;
}
//****************************************************************************80

double r8_gamma_01_sample ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA_01_SAMPLE samples the standard Gamma distribution.
//
//  Discussion:
//
//    This procedure corresponds to algorithm GD in the reference.
//
//    pdf ( a; x ) = 1/gamma(a) * x^(a-1) * exp ( - x )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN77 version by Barry Brown, James Lovato.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Joachim Ahrens, Ulrich Dieter,
//    Generating Gamma Variates by a Modified Rejection Technique,
//    Communications of the ACM,
//    Volume 25, Number 1, January 1982, pages 47-54.
//
//  Parameters:
//
//    Input, double A, the shape parameter of the standard gamma
//    distribution.  0.0 < A.
//
//    Output, double R8_GAMMA_01_SAMPLE, a random deviate from the distribution.
//
{
  double a1 =  0.3333333;
  double a2 = -0.2500030;
  double a3 =  0.2000062;
  double a4 = -0.1662921;
  double a5 =  0.1423657;
  double a6 = -0.1367177;
  double a7 =  0.1233795;
  double b;
  double c;
  double d;
  double e;
  double e1 = 1.0;
  double e2 = 0.4999897;
  double e3 = 0.1668290;
  double e4 = 0.0407753;
  double e5 = 0.0102930;
  double p;
  double q;
  double q0;
  double q1 =  0.04166669;
  double q2 =  0.02083148;
  double q3 =  0.00801191;
  double q4 =  0.00144121;
  double q5 = -0.00007388;
  double q6 =  0.00024511;
  double q7 =  0.00024240;
  double r;
  double s;
  double s2;
  double si;
  double sqrt32 = 5.6568542494923801952;
  double t;
  double u;
  double v;
  double value;
  double w;
  double x;

  if ( 1.0 <= a )
  {
    s2 = a - 0.5;
    s = sqrt ( s2 );
    d = sqrt32 - 12.0 * s;
//
//  Immediate acceptance.
//
    t = r8_normal_01_sample ( );
    x = s + 0.5 * t;
    value = x * x;

    if ( 0.0 <= t )
    {
      return value;
    }
//
//  Squeeze acceptance.
//
    u = r8_uniform_01_sample ( );
    if ( d * u <= t * t * t )
    {
      return value;
    }

    r = 1.0 / a;
    q0 = (((((( 
        q7   * r 
      + q6 ) * r 
      + q5 ) * r 
      + q4 ) * r
      + q3 ) * r
      + q2 ) * r
      + q1 ) * r;
//
//  Approximation depending on size of parameter A.
//
    if ( 13.022 < a )
    {
      b = 1.77;
      si = 0.75;
      c = 0.1515 / s;
    }
    else if ( 3.686 < a )
    {
      b = 1.654 + 0.0076 * s2;
      si = 1.68 / s + 0.275;
      c = 0.062 / s + 0.024;
    }
    else
    {
      b = 0.463 + s + 0.178 * s2;
      si = 1.235;
      c = 0.195 / s - 0.079 + 0.16 * s;
    }
//
//  Quotient test.
//
    if ( 0.0 < x )
    {
      v = 0.5 * t / s;

      if ( 0.25 < fabs ( v ) )
      {
        q = q0 - s * t + 0.25 * t * t + 2.0 * s2 * log ( 1.0 + v );
      }
      else
      {
        q = q0 + 0.5 * t * t * (((((( 
            a7   * v 
          + a6 ) * v 
          + a5 ) * v 
          + a4 ) * v 
          + a3 ) * v 
          + a2 ) * v 
          + a1 ) * v;
      }

      if ( log ( 1.0 - u ) <= q )
      {
        return value;
      }
    }

    for ( ; ; )
    {
      e = r8_exponential_01_sample ( );
      u = 2.0 * r8_uniform_01_sample ( ) - 1.0;
 
      if ( 0.0 <= u )
      {
        t = b + fabs ( si * e );
      }
      else
      {
        t = b - fabs ( si * e );
      }
//
//  Possible rejection.
//
      if ( t < -0.7187449 )
      {
        continue;
      }
//
//  Calculate V and quotient Q.
//
      v = 0.5 * t / s;

      if ( 0.25 < fabs ( v ) )
      {
        q = q0 - s * t + 0.25 * t * t + 2.0 * s2 * log ( 1.0 + v );
      }
      else
      {
        q = q0 + 0.5 * t * t * (((((( 
            a7   * v 
          + a6 ) * v 
          + a5 ) * v 
          + a4 ) * v 
          + a3 ) * v 
          + a2 ) * v 
          + a1 ) * v;
      }
//
//  Hat acceptance.
//
      if ( q <= 0.0 )
      {
        continue;
      }

      if ( 0.5 < q )
      {
        w = exp ( q ) - 1.0;
      }
      else
      {
        w = (((( 
            e5   * q 
          + e4 ) * q 
          + e3 ) * q 
          + e2 ) * q 
          + e1 ) * q;
      }
//
//  May have to sample again.
//
      if ( c * fabs ( u ) <= w * exp ( e - 0.5 * t * t ) )
      {
        break;
      }
    }

    x = s + 0.5 * t;
    value = x * x;
  }
//
//  Method for A < 1.
//
  else if ( a < 1.0 )
  {
    b = 1.0 + 0.3678794 * a;

    for ( ; ; )
    {
      p = b * r8_uniform_01_sample ( );

      if ( p < 1.0 )
      {
        value = exp ( log ( p ) / a );
        if ( value <= r8_exponential_01_sample ( ) )
        {
          break;
        }
      }
      else
      {
        value = - log ( ( b - p ) / a );
        if ( ( 1.0 - a ) * log ( value ) <= r8_exponential_01_sample ( ) )
        {
          break;
        }
      }
    }
  }
  return value;
}
//****************************************************************************80

double r8_invchi_pdf ( double df, double rval )

//****************************************************************************80
//
//  Purpose:
//
//    R8_INVCHI_PDF evaluates the PDF of an inverse chi-squared distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double DF, the degrees of freedom.
//    0.0 < DF.
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_INVCHI_PDF, the value of the PDF at RVAL.
//
{
  double temp1;
  double temp2;
  double value; 

  if ( df <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_INVCHI_PDF - Fatal error!\n";
    cerr << "  Degrees of freedom must be positive.\n";
    exit ( 1 );
  }
      
  if ( rval <= 0.0 )
  {
    value = 0.0;
  }
  else
  {
    temp2 = df * 0.5;

    temp1 = - temp2 * log ( 2.0 ) - ( temp2 + 1.0 ) * log ( rval ) 
      - 0.5 / rval - r8_gamma_log ( temp2 );

    value = exp ( temp1 );
  }
  return value;
}
//****************************************************************************80

double r8_invchi_sample ( double df )

//****************************************************************************80
//
//  Purpose:
//
//    R8_INVCHI_SAMPLE samples the inverse chi-squared distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 May 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double DF, the degrees of freedom.
//    0.0 < DF.
//
//    Output, double R8_INVCHI_SAMPLE, the sample of the PDF.
//
{
  double a;
  double b;
  double value;

  a = 0.5;
  b = 0.5 * df;
  value = r8_gamma_sample ( a, b );

  if ( value != 0.0 )
  {
    value = 1.0 / value;
  }
  return value;
}
//****************************************************************************80

double r8_invgam_pdf ( double beta, double alpha, double rval )

//****************************************************************************80
//
//  Purpose:
//
//    R8_INVGAM_PDF evaluates the PDF of an inverse gamma distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double BETA, the rate parameter.
//    0.0 < BETA.
//
//    Input, double ALPHA, the shape parameter.
//    0.0 < ALPHA.
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_INVGAM_PDF, the value of the PDF at RVAL.
//
{
  double temp;
  double value;

  if ( alpha <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_INVGAM_PDF - Fatal error!\n";
    cerr << "  Parameter ALPHA is not positive.\n";
    exit ( 1 );
  }

  if ( beta <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_INVGAM_PDF - Fatal error!\n";
    cerr << "  Parameter BETA is not positive.\n";
    exit ( 1 );
  }

  if ( rval <= 0.0 )
  {
    value = 0.0;
  }
  else
  {
    temp = alpha * log ( beta ) - ( alpha + 1.0 ) * log ( rval ) 
      - beta / rval - r8_gamma_log ( alpha );

    value = exp ( temp );
  }

  return value;
}
//****************************************************************************80

double r8_invgam_sample ( double beta, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    R8_INVGAM_SAMPLE samples an inverse gamma distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 May 2013
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, double BETA, the rate parameter.
//    0.0 < BETA.
//
//    Input, double ALPHA, the shape parameter.
//    0.0 < ALPHA.
//
//    Output, double R8_INVGAM_SAMPLE, a sample of the PDF.
//
{
  double value;

  value = r8_gamma_sample ( beta, alpha );
  if ( value != 0.0 )
  {
    value = 1.0 / value;
  }
  return value;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = y;
  }
  else
  {
    value = x;
  }
  return value;
}
//****************************************************************************80

double r8_normal_pdf ( double av, double sd, double rval )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NORMAL_PDF evaluates the PDF of a normal distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double AV, the mean value.
//
//    Input, double SD, the standard deviation.
//    0.0 < SD.
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_NORMAL_PDF, the value of the PDF at RVAL.
//
{
  double pi = 3.141592653589793;
  double rtemp;
  double value;

  if ( sd <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_NORMAL_PDF - Fatal error!\n";
    cerr << "  Standard deviation must be positive.\n";
    exit ( 1 );
  }

  rtemp = ( rval - av ) * ( rval - av ) * 0.5 / ( sd * sd );

  value = exp ( - rtemp ) / sd / sqrt ( 2.0 * pi );

  return value;
}
//****************************************************************************80

double r8_normal_sample ( double av, double sd )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NORMAL_SAMPLE generates a normal random deviate.
//
//  Discussion:
//
//    This procedure generates a single random deviate from a normal 
//    distribution with mean AV, and standard deviation SD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN77 version by Barry Brown, James Lovato.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Joachim Ahrens, Ulrich Dieter,
//    Extensions of Forsythe's Method for Random
//    Sampling from the Normal Distribution,
//    Mathematics of Computation,
//    Volume 27, Number 124, October 1973, page 927-937.
//
//  Parameters:
//
//    Input, double AV, the mean.
//
//    Input, double SD, the standard deviation.
//
//    Output, double R8_NORMAL_SAMPLE, a random deviate from the distribution.
//
{
  double value;

  value = sd * r8_normal_01_sample ( ) + av;

  return value;
}
//****************************************************************************80

double r8_normal_01_pdf ( double rval )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NORMAL_01_PDF evaluates the PDF of a standard normal distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 June 2013
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_NORMAL_01_PDF, the value of the PDF at RVAL.
//
{
  double pi = 3.141592653589793;
  double value;

  value = exp ( - 0.5 * rval * rval ) / sqrt ( 2.0 * pi );

  return value;
}
//****************************************************************************80

double r8_normal_01_sample ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NORMAL_01_SAMPLE samples the standard normal probability distribution.
//
//  Discussion:
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//    The Box-Muller method is used, which is efficient, but
//    generates two values at a time.
//
//    Typically, we would use one value and save the other for the next call.
//    However, the fact that this function has saved memory makes it difficult
//    to correctly handle cases where we want to re-initialize the code,
//    or to run in parallel.  Therefore, we will instead use the first value
//    and DISCARD the second.
//
//    EFFICIENCY must defer to SIMPLICITY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_NORMAL_01_SAMPLE, a normally distributed random value.
//
{
  const double pi = 3.14159265358979323;
  double r1;
  double r2;
  double x;

  r1 = r8_uniform_01_sample ( );
  r2 = r8_uniform_01_sample ( );

  x = sqrt ( -2.0 * log ( r1 ) ) * cos ( 2.0 * pi * r2 );

  return x;
}
//****************************************************************************80

double r8_scinvchi_pdf ( double df, double s, double rval )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SCINVCHI_PDF: PDF for a scaled inverse chi-squared distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double DF, the degrees of freedom.
//    0.0 < DF.
//
//    Input, double S, the scale factor.
//    0.0 < S.
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_SCINVCHI_PDF, the value of the PDF at RVAL.
//    inverse-chi-square distribution.
//
{
  double temp1;
  double temp2;
  double value;

  if ( df <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_SCINVCHI_PDF - Fatal error!\n";
    cerr << "  Degrees of freedom must be positive.\n";
    exit ( 1 );
  }

  if ( s <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_SCINVCHI_PDF - Fatal error!\n";
    cerr << "  Scale parameter must be positive.\n";
    exit ( 1 );
  }

  if ( rval <= 0.0 )
  {
    value = 0.0;
  }
  else
  {
    temp2 = df * 0.5;

    temp1 = temp2 * log ( temp2 ) + temp2 * log ( s ) 
      - ( temp2 * s / rval ) 
      - ( temp2 + 1.0 ) * log ( rval ) - r8_gamma_log ( temp2 );

    value = exp ( temp1 );
  }

  return value;
}
//****************************************************************************80

double r8_scinvchi_sample ( double df, double s )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SCINVCHI_SAMPLE: sample a scaled inverse chi-squared distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 May 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double DF, the degrees of freedom.
//    0.0 < DF.
//
//    Input, double S, the scale factor.
//    0.0 < S.
//
//    Output, double R8_SCINVCHI_SAMPLE, a sample of the distribution
//
{
  double a;
  double b;
  double value;

  a = 0.5 * df * s;
  b = 0.5 * df;
  value = r8_gamma_sample ( a, b );
  if ( value != 0.0 )
  {
    value = 1.0 / value;
  }
  return value;
}
//****************************************************************************80

double r8_uniform_pdf ( double lower, double upper, double rval )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_PDF evaluates the PDF of a uniform distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double LOWER, UPPER, the lower and upper range limits.
//    LOWER < UPPER.
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_UNIFORM_PDF, the value of the PDF at RVAL.
//
{
  double value;

  if ( upper <= lower )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_PDF - Fatal error!\n";
    cerr << "  For uniform PDF, the lower limit must be\n";
    cerr << "  less than the upper limit\n";
    exit ( 1 );
  }

  if ( rval < lower )
  {
    value = 0.0;
  }
  else if ( rval <= upper )
  {
    value = 1.0 / ( upper - lower );
  }
  else
  {
    value = 0.0;
  }

  return value;
}
//****************************************************************************80

double r8_uniform_sample ( double low, double high )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_SAMPLE generates a uniform random deviate.
//
//  Discussion:
//
//    This procedure generates a real deviate uniformly distributed between
//    LOW and HIGH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN77 version by Barry Brown, James Lovato.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double LOW, HIGH, the lower and upper bounds.
//
//    Output, double R8_UNIFORM_SAMPLE, a random deviate from the distribution.
//
{
  double value;

  value = low + ( high - low ) * r8_uniform_01_sample ( );

  return value;
}
//****************************************************************************80

double r8_uniform_01_pdf ( double rval )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01_PDF evaluates the PDF of a standard uniform distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 June 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_UNIFORM_01_PDF, the value of the PDF at RVAL.
//
{
  double value;

  if ( rval < 0.0 )
  {
    value = 0.0;
  }
  else if ( rval <= 1.0 )
  {
    value = 1.0;
  }
  else
  {
    value = 0.0;
  }

  return value;
}
//****************************************************************************80

double r8_uniform_01_sample ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01_SAMPLE generates a uniform random deviate from [0,1].
//
//  Discussion:
//
//    This function should be the only way that the package accesses random
//    numbers.
//
//    Setting OPTION to 0 accesses the R8_UNI_01() function in RNGLIB,
//    for which there are versions in various languages, which should result
//    in the same values being returned.  This should be the only place
//    in this library that accesses a function from RNGLIB.
//
//    Setting OPTION to 1 in the C++ version calls the system
//    random number generator "rand()".
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Output, double R8_UNIFORM_01_SAMPLE, a random deviate.
//
{
  const int option = 0;
  double value;

  if ( option == 0 )
  {
    value = r8_uni_01 ( );
  }
  else
  {
    value = ( double ) rand ( ) / ( double ) RAND_MAX;
  }

  return value;
}
//****************************************************************************80

double *r8mat_mtv_new ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MTV_NEW multiplies a transposed matrix times a vector.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, double A[M,N], the M by N matrix.
//
//    Input, double X[M], the vector to be multiplied by A.
//
//    Output, double R8MAT_MTV_NEW[N], the product A'*X.
//
{
  int i;
  int j;
  double *y;

  y = new double[n];

  for ( j = 0; j < n; j++ )
  {
    y[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      y[j] = y[j] + a[i+j*m] * x[i];
    }
  }

  return y;
}
//****************************************************************************80

double *r8mat_mv_new ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MV_NEW multiplies a matrix times a vector.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, double A[M,N], the M by N matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8MAT_MV_NEW[M], the product A*X.
//
{
  int i;
  int j;
  double *y;

  y = new double[m];

  for ( i = 0; i < m; i++ )
  {
    y[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      y[i] = y[i] + a[i+j*m] * x[j];
    }
  }

  return y;
}
/******************************************************************************/

double r8mat_podet ( int n, double r[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PODET computes the determinant of a factored positive definite matrix.

  Discussion:

    This routine expects to receive R, the upper triangular factor of A,
    computed by R8MAT_POFAC, with the property that A = R' * R.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 June 2013

  Author:

    C version by John Burkardt.

  Reference:

    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979,
    ISBN13: 978-0-898711-72-1,
    LC: QA214.L56.

  Parameters:

    Input, int N, the order of the matrix.

    Input, double R[N*N], the Cholesky factor of A.

    Output, double R8MAT_PODET, the determinant of A.
*/
{
  double det;
  int i;

  det = 1.0;
  for ( i = 0; i < n; i++ )
  {
    det = det * r[i+i*n] * r[i+i*n];
  }
  return det;
}
//****************************************************************************80

double *r8mat_pofac ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_POFAC factors a real symmetric positive definite matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 June 2013
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch, 
//    Pete Stewart,
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the matrix to be factored.
//
//    Output, double R8MAT_POFAC[N*N], an upper triangular matrix such that
//    A = R'*R.
//
{
  double dot;
  int i;
  int j;
  int k;
  double *r;
  double s;
  double t;

  r = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i <= j; i++ )
    {
      r[i+j*n] = a[i+j*n];
    }
    for ( i = j + 1; i < n; i++ )
    {
      r[i+j*n] = 0.0;
    }
  }

  for ( j = 0; j < n; j++ )
  {
    s = 0.0;

    for ( k = 0; k < j; k++ )
    {
      dot = 0.0;
      for ( i = 0; i < k; i++ )
      {
        dot = dot + r[i+k*n] * r[i+j*n];
      }
      t = r[k+j*n] - dot;
      t = t / r[k+k*n];
      r[k+j*n] = t;
      s = s + t * t;
    }

    s = r[j+j*n] - s;

    if ( s < 0.0 )
    {
      cerr << "\n";
      cerr << "R8MAT_POFAC - Fatal error!\n";
      cerr << "  The matrix is not positive definite.\n";
      exit ( 1 );
    }

    if ( s == 0.0 )
    {
      cerr << "\n";
      cerr << "R8MAT_POFAC - Warning!\n";
      cerr << "  The matrix is not strictly positive definite.\n";
    }

    r[j+j*n] = sqrt ( s );
  }

  return r;
}
//****************************************************************************80

double *r8mat_poinv ( int n, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_POINV computes the inverse of a factored positive definite matrix.
//
//  Discussion:
//
//    This routine expects to receive R, the upper triangular factor of A,
//    computed by R8MAT_POFAC, with the property that A = R' * R.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 June 2013
//
//  Author:
//
//    C version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//  Parameters:
//
//    Input, int N, the order of the matrix A.
//
//    Input, double R[N*N], the Cholesky factor of A.
//
//    Input, double R8MAT_POINV[N*N], the inverse of A.
//
{
  double *b;
  int i;
  int j;
  int k;
  double t;

  b = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      b[i+j*n] = r[i+j*n];
    }
  }

  for ( k = 0; k < n; k++ )
  {
    b[k+k*n] = 1.0 / b[k+k*n];
    for ( i = 0; i < k; i++ )
    {
      b[i+k*n] = - b[i+k*n] * b[k+k*n];
    }
    for ( j = k + 1; j < n; j++ )
    {
      t = b[k+j*n];
      b[k+j*n] = 0.0;
      for ( i = 0; i <= k; i++ )
      {
        b[i+j*n] = b[i+j*n] + t * b[i+k*n];
      }
    }
  }
/*
  Form inverse(R) * (inverse(R))'.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( k = 0; k < j; k++ )
    {
      t = b[k+j*n];
      for ( i = 0; i <= k; i++ )
      {
        b[i+k*n] = b[i+k*n] + t * b[i+j*n];
      }
    }
    t = b[j+j*n];
    for ( i = 0; i <= j; i++ )
    {
      b[i+j*n] = b[i+j*n] * t;
    }
  }
  return b;
}
//****************************************************************************80

double *r8mat_upsol ( int n, double r[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UPSOL solves R * X = B, for an upper triangular matrix R.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 June 2013
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double R[N*N], the upper triangular matrix.
//
//    Input, double B[N], the right hand side.
//
//    Output, double R8MAT_UPSOL[N], the solution.
//
{
  int i;
  int j;
  double *x;

  x = new double[n];
  for ( j = 0; j < n; j++ )
  {
    x[j] = b[j];
  }
  return x;
}
//****************************************************************************80

double *r8mat_utsol ( int n, double r[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UTSOL solves R' * X = B for an upper triangular matrix R.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 June 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double R[N*N], the upper triangular matrix.
//
//    Input, double B[N], the right hand side.
//
//    Output, double X[N], the solution.
//
{
  int i;
  int j;
  double *x;

  x = new double[n];

  for ( j = 0; j < n; j++ )
  {
    x[j] = b[j];
  }

  x[0] = x[0] / r[0+0*n];

  for ( j = 1; j < n; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      x[j] = x[j] - r[i+j*n] * x[i];
    }
    x[j] = x[j] / r[j+j*n]; 
  }
  return x;
}
//****************************************************************************80

double r8vec_dot_product ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }
  return value;
}
//****************************************************************************80

double r8vec_multinormal_pdf ( int n, double mu[], double r[], double c_det, 
  double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MULTINORMAL_PDF evaluates a multivariate normal PDF.
//
//  Discussion:
//
//    PDF ( MU(1:N), C(1:N,1:N); X(1:N) ) = 
//      1 / ( 2 * pi ) ^ ( N / 2 ) * 1 / sqrt ( det ( C ) )
//      * exp ( - ( X - MU )' * inverse ( C ) * ( X - MU ) / 2 )
//
//    Here,
//
//      X is the argument vector of length N,
//      MU is the mean vector of length N,
//      C is an N by N positive definite symmetric covariance matrix.
//
//    The properties of C guarantee that it has an upper triangular
//    matrix R, the Cholesky factor, such that C = R' * R.  It is the
//    matrix R that is required by this routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double MU[N], the mean vector.
//
//    Input, double R[N*N], the upper triangular Cholesky
//    factor of the covariance matrix C.
//
//    Input, double C_DET, the determinant of the
//    covariance matrix C.
//
//    Input, double X[N], a sample of the distribution.
//
//    Output, double R8VEC_MULTINORMAL_PDF, the PDF evaluated
//    at X.
//
{
  double *b;
  int i;
  double pdf;
  double pi = 3.141592653589793;
  double xcx;
  double *y;
//
//  Compute:
//    inverse(R')*(x-mu) = y
//  by solving:
//    R'*y = x-mu
//
  b = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    b[i] = x[i] - mu[i];
  }
  y = r8mat_utsol ( n, r, b );
//
//  Compute:
//    (x-mu)' * inv(C)          * (x-mu)
//  = (x-mu)' * inv(R'*R)       * (x-mu)
//  = (x-mu)' * inv(R) * inv(R) * (x-mu)
//  = y' * y.
//
  xcx = r8vec_dot_product ( n, y, y );

  pdf = 1.0 / sqrt ( pow ( 2.0 * pi, n ) ) 
      * 1.0 / sqrt ( c_det ) 
      * exp ( - 0.5 * xcx );

  delete [] b;
  delete [] y;

  return pdf;
}
//****************************************************************************80

double *r8vec_multinormal_sample ( int n, double mu[], double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MULTINORMAL_SAMPLE samples a multivariate normal PDF.
//
//  Discussion:
//
//    PDF ( MU(1:N), C(1:N,1:N); X(1:N) ) = 
//      1 / ( 2 * pi ) ^ ( N / 2 ) * 1 / det ( C )
//      * exp ( - ( X - MU )' * inverse ( C ) * ( X - MU ) / 2 )
//
//    Here,
//
//      X is the argument vector of length N,
//      MU is the mean vector of length N,
//      C is an N by N positive definite symmetric covariance matrix.
//
//    The properties of C guarantee that it has an upper triangular
//    matrix R, the Cholesky factor, such that C = R' * R.  It is the
//    matrix R that is required by this routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 June 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double MU[N], the mean vector.
//
//    Input, double R[N*N], the upper triangular Cholesky
//    factor of the covariance matrix C.
//
//    Output, double R8VEC_MULTINORMAL_SAMPLE[N], a sample of the distribution.
//
{
  int i;
  int j;
  double *x;
  double *z;
//
//  Compute X = MU + R' * Z
//  where Z is a vector of standard normal variates.
//
  z = new double[n];
  for ( j = 0; j < n; j++ )
  {
    z[j] = r8_normal_01_sample ( );
  }

  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = mu[i];
    for ( j = 0; j <= i; j++ )
    {
      x[i] = x[i] + r[j+i*n] * z[j];
    }
  }

  delete [] z;

  return x;
}
