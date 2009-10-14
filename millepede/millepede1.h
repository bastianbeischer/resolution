// $Id: millepede1.h,v 1.1 2009/10/14 16:51:32 beischer Exp $

extern "C" void initgl_( int& nagb, int& nalc, int& nstdev, int& iprlim );
#define INITGL initgl_

extern "C" void zerloc_( float* dergb, float* derlc );
#define ZERLOC zerloc_

extern "C" void equloc_( float* dergb, float* derlc, float& rmeas, float& sigma );
#define EQULOC equloc_

extern "C" void fitloc_( );
#define FITLOC fitloc_

extern "C" void constf_( float* dercs, float& rhs );
#define CONSTF constf_

extern "C" void testmp_( int& iarg );
#define TESTMP testmp_

extern "C" void fitglo_( float* par );
#define FITGLO fitglo_

extern "C" void parsig_( int& i, float& value );
#define PARSIG parsig_

extern "C" void initun_( int& nIter, float& otherVal );
#define INITUN initun_
