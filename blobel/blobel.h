// $Id: blobel.h,v 1.2 2009/10/14 09:24:33 beischer Exp $

extern "C" void dvalley_( double& f, double a[], int& nc );
#define DVALLEY dvalley_

extern "C" void dvallin_( int& n, double st[], int& nflim );
#define DVALLIN dvallin_

