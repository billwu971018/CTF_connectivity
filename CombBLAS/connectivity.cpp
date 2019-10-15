#include "connectivity.h"



FullyDistVec<int,int> hook_matrix(int n, SpParMat<int,int,SpDCCols<int,int>> & A)
{
    FullyDistVec<int,int> p(A.getcommgrid());
    p.iota(A.getnrow(), 0);
    FullyDistVec<int,int> p_prev(n, 0);
    while(!(p == p_prev))
    {
    	p_prev = p;
    	FullyDistVec<int,int> q(n, 0);
        q = SpMV<MAX_TIMES_SR>(A, p);
        FullyDistVec<int,int> r(n, 0);
        p.EWiseApply(q, std::max<int>()); //max_vector
        r = p;
        FullyDistVec<int,int> pi(p);
        shortcut(p, p, p);
        while(!(pi == p)){
        	pi = p;
        	shortcut(p, p, p);
        }
    }
    return p;
}

FullyDistVec<int,int> supervertex_matrix(int n, SpParMat<int,int,SpDCCols<int,int>> & A, FullyDistVec<int,int> p, int sc2)
{

}
