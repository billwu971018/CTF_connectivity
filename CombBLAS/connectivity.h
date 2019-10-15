#include <mpi.h>

// These macros should be defined before stdint.h is included
#ifndef __STDC_CONSTANT_MACROS
#define __STDC_CONSTANT_MACROS
#endif
#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif
#include <stdint.h>

#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>
#include <cmath>
#include "CombBLAS/CombBLAS.h"


struct MAX_TIMES_SR
{
    typedef typename promote_trait<int,int>::T_promote T_promote;
    static T_promote id(){ return std::numeric_limits<T_promote>::max(); };
    static bool returnedSAID() { return false; }
    static MPI_Op mpi_op() { return MPI_MAX; };
    
    static T_promote add(const T_promote & arg1, const T_promote & arg2)
    {
        return arg1 * arg2;
    }
    
    static T_promote multiply(const int & arg1, const int & arg2)
    {
        return std::max(arg1, arg2);
    }
    
    static void axpy(const int a, const int & x, T_promote & y)
    {
        y = add(y, multiply(a, x));
    }
};

// Connectivity

FullyDistVec<int,int> hook_matrix(int n, SpParMat<int,int,SpDCCols<int,int>> & A);
FullyDistVec<int,int> supervertex_matrix(int n, SpParMat<int,int,SpDCCols<int,int>> & A, FullyDistVec<int,int> p, int sc2);

//Utilities
void shortcut(FullyDistVec<int,int> & p, FullyDistVec<int,int> & q, FullyDistVec<int,int> & rec_p, Vector<int> ** nonleaves=NULL, bool create_nonleaves=false);
