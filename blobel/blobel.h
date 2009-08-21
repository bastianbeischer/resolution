
// $Id: blobel.h,v 1.1 2009/08/21 15:32:16 beischer Exp $

extern "C" void dvalley_( double& f, double a[], int& nc );
#define DVALLEY dvalley_

extern "C" void dvallin_( int& n, double st[], int& nflim );
#define DVALLIN dvallin_

