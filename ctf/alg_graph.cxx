#include "alg_graph.h"

Int64Pair::Int64Pair(int64_t i1, int64_t i2) {
  this->i1 = i1;
  this->i2 = i2;
}

Int64Pair Int64Pair::swap() {
  return {this->i2, this->i1};
}

void mat_set(Matrix<int>* matrix, Int64Pair index, int value) {
  int64_t idx[1];
  idx[0] = index.i2 * matrix->nrow + index.i1;
  int fill[1];
  fill[0] = value;
  matrix->write(1, idx, fill);
}

Graph::Graph() {
  this->numVertices = 0;
  this->edges = new vector<Int64Pair>();
}

Matrix<int>* Graph::adjacencyMatrix(World* world, bool sparse) {
  auto attr = 0;
  if (sparse) {
    attr = SP;
  }
  auto A = new Matrix<int>(numVertices, numVertices,
      attr, *world, MAX_TIMES_SR);
  for (auto edge : *edges) {
    mat_set(A, edge);
    mat_set(A, edge.swap());
  }
  return A;
}

void init_pvector(Vector<int>* p)
{
  int64_t npairs;
  Pair<int> * loc_pairs;
  p->read_local(&npairs, &loc_pairs);
  for (int64_t i = 0; i < npairs; i++){
    loc_pairs[i].d = loc_pairs[i].k;
  }
  p->write(npairs, loc_pairs);
  delete [] loc_pairs;
}

// NOTE: can't use bool as return
template <typename dtype>
int64_t are_vectors_different(CTF::Vector<dtype> & A, CTF::Vector<dtype> & B)
{
  CTF::Scalar<int64_t> s;
  if (!A.is_sparse && !B.is_sparse){
    s[""] += CTF::Function<dtype,dtype,int64_t>([](dtype a, dtype b){ return a!=b; })(A["i"],B["i"]);
  } else {
    auto C = Vector<dtype>(A.len, SP*A.is_sparse, *A.wrld);
    C["i"] += A["i"];
    ((int64_t)-1)*C["i"] += B["i"];
    s[""] += CTF::Function<dtype,int64_t>([](dtype a){ return (int64_t)(a!=0); })(C["i"]);

  }
  return s.get_val();
}
template int64_t are_vectors_different<int>(CTF::Vector<int> & A, CTF::Vector<int> & B);

// p[i] = rec_p[q[i]]
// if create_nonleaves=true, computing non-leaf vertices in parent forest
void shortcut(Vector<int> & p, Vector<int> & q, Vector<int> & rec_p, Vector<int> ** nonleaves, bool create_nonleaves)
{
  TAU_FSTART(Unoptimized_shortcut);
  //Timer t_us("Unoptimized_shortcut");
  //t_us.start();
  int64_t npairs;
  Pair<int> * loc_pairs;
  if (q.is_sparse){
    //if we have updated only a subset of the vertices
    q.get_local_pairs(&npairs, &loc_pairs, true);
  } else {
    //if we have potentially updated all the vertices
    q.get_local_pairs(&npairs, &loc_pairs);
  }
  Pair<int> * remote_pairs = new Pair<int>[npairs];
  for (int64_t i=0; i<npairs; i++){
    remote_pairs[i].k = loc_pairs[i].d;
  }
  TAU_FSTART(Unoptimized_shortcut_rread);
  //Timer t_usr("Unoptimized_shortcut_rread");
  //t_usr.start();
  rec_p.read(npairs, remote_pairs); //obtains rec_p[q[i]]
  //t_usr.stop();
  TAU_FSTOP(Unoptimized_shortcut_rread);
  for (int64_t i=0; i<npairs; i++){
    loc_pairs[i].d = remote_pairs[i].d; //p[i] = rec_p[q[i]]
  }
  delete [] remote_pairs;
  TAU_FSTART(Unoptimized_shortcut_write);
  //Timer t_usw("Unoptimized_shortcut_write");
  //t_usw.start();
  p.write(npairs, loc_pairs); //enter data into p[i]
  TAU_FSTOP(Unoptimized_shortcut_write);
  //t_usw.stop();
  
  //prune out leaves
  if (create_nonleaves){
    TAU_FSTART(Unoptimized_shortcut_pruneleaves);
    //Timer t_usp("Unoptimzed_shortcut_pruneleaves");
    //t_usp.start();
    *nonleaves = new Vector<int>(p.len, *p.wrld, *p.sr);
    //set nonleaves[i] = max_j p[j], i.e. set nonleaves[i] = 1 if i has child, i.e. is nonleaf
    for (int64_t i=0; i<npairs; i++){
      loc_pairs[i].k = loc_pairs[i].d;
      loc_pairs[i].d = 1;
    }
    //FIXME: here and above potential optimization is to avoid duplicate queries to parent
    (*nonleaves)->write(npairs, loc_pairs);
    (*nonleaves)->operator[]("i") = (*nonleaves)->operator[]("i")*p["i"];
    (*nonleaves)->sparsify();
    TAU_FSTOP(Unoptimized_shortcut_pruneleaves);
    //t_usp.stop();
  }
   
  delete [] loc_pairs;
  TAU_FSTOP(Unoptimized_shortcut);
  //t_us.stop();
}

// p[i] = rec_p[q[i]]
// if create_nonleaves=true, computing non-leaf vertices in parent forest
void shortcut2(Vector<int> & p, Vector<int> & q, Vector<int> & rec_p, int64_t sc2, World * world, Vector<int> ** nonleaves, bool create_nonleaves)
{
  if (sc2 <= 0) { // run unoptimized shortcut
    shortcut(p, q, rec_p, nonleaves, create_nonleaves);
    return;
  }

  TAU_FSTART(Optimized_shortcut);
  //Timer t_os2("Optimized_shortcut2");
  //t_os2.start();
  
  int64_t rec_gf_npairs;
  Pair<int> * rec_p_loc_pairs;
  if (rec_p.is_sparse) {
    rec_p.get_local_pairs(&rec_gf_npairs, &rec_p_loc_pairs, true);
  } else {
    rec_p.get_local_pairs(&rec_gf_npairs, &rec_p_loc_pairs);
  }
  
  int64_t q_npairs;
  Pair<int> * q_loc_pairs;
  bool delete_p = true;
  if (&q == &rec_p) {
    q_npairs = rec_gf_npairs;
    q_loc_pairs = rec_p_loc_pairs;
    delete_p = false;
  } else {
    p["i"] += q["i"];
    if (q.is_sparse){
      //if we have updated only a subset of the vertices
      q.get_local_pairs(&q_npairs, &q_loc_pairs, true);
    } else {
      //if we have potentially updated all the vertices
      q.get_local_pairs(&q_npairs, &q_loc_pairs);
    }
  }
  
  int64_t global_roots_num;
  int64_t loc_roots_num;
  roots_num(rec_gf_npairs, rec_p_loc_pairs, &loc_roots_num, &global_roots_num, world);
  
  if (global_roots_num < sc2) {
    int global_roots[global_roots_num];
    roots(rec_gf_npairs, loc_roots_num, rec_p_loc_pairs, global_roots, world);
    
    int64_t * nontriv_loc_indices;
    int64_t loc_nontriv_num;
    create_nontriv_loc_indices(nontriv_loc_indices, &loc_nontriv_num, global_roots_num, global_roots, q_npairs, q_loc_pairs, world);

    Pair<int> * nontriv_loc_pairs = new Pair<int>[loc_nontriv_num];
    Pair<int> * remote_pairs = new Pair<int>[loc_nontriv_num];
    for (int64_t i = 0; i < loc_nontriv_num; i++) {
      int64_t nontriv_index = nontriv_loc_indices[i];
      nontriv_loc_pairs[i] = q_loc_pairs[nontriv_index];
      remote_pairs[i].k = q_loc_pairs[nontriv_index].d;
    }
  
    TAU_FSTART(Optimized_shortcut_read);
    //Timer t_osr("Optimized_shortcut2_read");
    //t_osr.start();
    rec_p.read(loc_nontriv_num, remote_pairs); //obtains rec_p[q[i]]
    //t_osr.stop();
    TAU_FSTOP(Optimized_shortcut_read);
    for(int64_t i = 0; i < loc_nontriv_num; i++) {
      nontriv_loc_pairs[i].d = remote_pairs[i].d;
    }
    
    for (int64_t i = 0; i < loc_nontriv_num; i++) { // update loc_pairs for create_nonleaves step
      int64_t nontriv_index = nontriv_loc_indices[i];
      q_loc_pairs[nontriv_index].d = remote_pairs[i].d; // p[i] = rec_p[q[i]]
    }
 
    TAU_FSTART(Optimized_shortcut_write); 
    //Timer t_osw("Optimized_shortcut2_write");
    //t_osw.start();
    p.write(loc_nontriv_num, nontriv_loc_pairs); //enter data into p[i]
    //t_osw.stop();
    TAU_FSTOP(Optimized_shortcut_write); 
    
    delete [] remote_pairs;
    delete [] nontriv_loc_pairs;
    delete [] nontriv_loc_indices;
  } else { // original shortcut
    Pair<int> * remote_pairs = new Pair<int>[q_npairs];
    for (int64_t i=0; i<q_npairs; i++) {
      remote_pairs[i].k = q_loc_pairs[i].d;
    }
    TAU_FSTART(Unoptimized_shortcut_Orread);
    //Timer t_uso("Unoptimized_shortcut_Orread");
    //t_uso.start();
    rec_p.read(q_npairs, remote_pairs); //obtains rec_p[q[i]]
    //t_uso.stop();
    TAU_FSTOP(Unoptimized_shortcut_Orread);
    for (int64_t i=0; i<q_npairs; i++){
      q_loc_pairs[i].d = remote_pairs[i].d; //p[i] = rec_p[q[i]]
    }
    TAU_FSTART(Unoptimized_shortcut_Owrite);
    //Timer t_usw("Unoptimized_shortcut_Owrite");
    //t_usw.start();
    p.write(q_npairs, q_loc_pairs); //enter data into p[i]
    //t_usw.stop();
    TAU_FSTOP(Unoptimized_shortcut_Owrite);

    delete [] remote_pairs;
  }
  
  //prune out leaves
  if (create_nonleaves) {
    TAU_FSTART(Unoptimized_shortcut_Opruneleaves);
    //Timer t_usoo("Unoptimized_shortcut_Opruneleaves");
    //t_usoo.start();
    *nonleaves = new Vector<int>(p.len, *p.wrld, *p.sr); //set nonleaves[i] = max_j p[j], i.e. set nonleaves[i] = 1 if i has child, i.e. is nonleaf
    for (int64_t i=0; i<q_npairs; i++){
      q_loc_pairs[i].k = q_loc_pairs[i].d;
      q_loc_pairs[i].d = 1;
    }
    //FIXME: here and above potential optimization is to avoid duplicate queries to parent
    (*nonleaves)->write(q_npairs, q_loc_pairs);
    (*nonleaves)->operator[]("i") = (*nonleaves)->operator[]("i")*p["i"];
    (*nonleaves)->sparsify();
    TAU_FSTOP(Unoptimized_shortcut_Opruneleaves);
    //t_usoo.stop();
  }
  TAU_FSTOP(Optimized_shortcut);
  //t_os2.stop();

  if (delete_p) delete [] q_loc_pairs;
  delete [] rec_p_loc_pairs;
}

void roots_num(int64_t npairs, Pair<int> * loc_pairs, int64_t * loc_roots_num, int64_t * global_roots_num, World * world) {
  *loc_roots_num = 0;
  for (int64_t i = 0; i < npairs; i++) {
    Pair<int> loc_pair = loc_pairs[i];
    if (loc_pair.d == loc_pair.k) {
      (*loc_roots_num)++;
    }
  }
    
  MPI_Allreduce(loc_roots_num, global_roots_num, 1, MPI_LONG_LONG, MPI_SUM, world->comm);
}

void roots(int64_t npairs, int64_t loc_roots_num, Pair<int> * loc_pairs, int * global_roots,  World * world) {
  int world_size;
  MPI_Comm_size(world->comm, &world_size);
 
  int loc_roots [loc_roots_num];
  int64_t j = 0;
  for (int64_t i=0; i<npairs; i++) {
    // same loop as roots_num but would introduce overhead
    Pair<int> loc_pair = loc_pairs[i];
    if (loc_pair.d == loc_pair.k) {
      loc_roots[j] = loc_pair.k;
      j++;
    }
  }

  int global_roots_nums [world_size];
  MPI_Allgather(&loc_roots_num, 1, MPI_INT, global_roots_nums, 1, MPI_INT, world->comm); // [3, 1, 2, 0, 4]

  // prefix sum
  int displs_roots [world_size];
  int64_t sum_roots = 0;
  for (int64_t i=0; i<world_size; i++) {
    displs_roots[i] = sum_roots;
    sum_roots += global_roots_nums[i];
  }

  MPI_Allgatherv(loc_roots, loc_roots_num, MPI_INT, global_roots, global_roots_nums, displs_roots, MPI_INT, world->comm); // [., ., ., ., ., ., ., ., ., ., .]?
}

void create_nontriv_loc_indices(int64_t *& nontriv_loc_indices, int64_t * loc_nontriv_num, int64_t global_roots_num, int * global_roots, int64_t q_npairs, Pair<int> * q_loc_pairs, World * world) {
  std::sort(global_roots, global_roots + global_roots_num);

  nontriv_loc_indices = new int64_t[q_npairs]; // wastes a bit of memory
  int64_t nontriv_index = 0;
  int64_t q_index = 0;
  bool end_roots = false;
  for (int64_t root_index = 0; root_index < global_roots_num && q_index < q_npairs; q_index++) { // construct nontrivial local indices O((n+m)log(n+m))
    if (q_loc_pairs[q_index].d < global_roots[root_index]) { // if a node's parent is not a root
      nontriv_loc_indices[nontriv_index] = q_index;
      nontriv_index++;
    }
    while (q_loc_pairs[q_index].d > global_roots[root_index]) {
      root_index++;
      if (root_index >= global_roots_num) { end_roots = true; break; }
       
      if (q_loc_pairs[q_index].d < global_roots[root_index]) { // if this node's parent is greater than previous root but less than current root
        nontriv_loc_indices[nontriv_index] = q_index;
        nontriv_index++;
      }
    }
    if (end_roots) { break; }
  }

  for (; q_index < q_npairs; q_index++) { // add nodes with parent greater than max root
    nontriv_loc_indices[nontriv_index] = q_index;
    nontriv_index++;
  }

  *loc_nontriv_num = nontriv_index;
}

// p[i] = rec_p[q[i]]
// read_local()
void shortcut3(Vector<int> & p, Vector<int> & q, Vector<int> & rec_p, Vector<int> & p_prev, MPI_Datatype &mpi_pkv, World * world)
{
  TAU_FSTART(Local_shortcut3);
  //Timer t_shortcut("Local_shortcut3");
  //t_shortcut.start();
  std::vector<struct parentkv> loc_cparents;
  // Find out the changed parents (locally)
  Pair<int> * pprs;
  Pair<int> * prev_pprs;
  int64_t npprs, prev_npprs;
  p.get_local_pairs(&npprs, &pprs);
  p_prev.get_local_pairs(&prev_npprs, &prev_pprs);
  assert(npprs == prev_npprs);
  for (int64_t i = 0; i < npprs; i++) {
    assert(pprs[i].k == prev_pprs[i].k);
    if (pprs[i].d != prev_pprs[i].d) {
      struct parentkv pkv;
      pkv.key = pprs[i].k;
      pkv.value = pprs[i].d;
      loc_cparents.push_back(pkv);
    }
  }
  int np = world->np;
  int *counts = new int[np];
  int nelements = (int)loc_cparents.size();
  MPI_Allgather(&nelements, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);
  int *disps = new int[np];
  for (int i = 0; i < np; i++)
    disps[i] = (i > 0) ? (disps[i-1] + counts[i-1]) : 0;

  struct parentkv *alldata;
  alldata = new struct parentkv[disps[np-1] + counts[np-1]];
  MPI_Allgatherv(loc_cparents.data(), nelements, mpi_pkv,
                  alldata, counts, disps, mpi_pkv, MPI_COMM_WORLD);
  
  // Build a map of the key value pair
  std::map<int64_t, int64_t> m_alldata;
  for(int i = 0; i < (disps[np-1] + counts[np-1]); i++) {
    m_alldata.insert({alldata[i].key, alldata[i].value});
  }
  delete [] counts;
  delete [] disps;
  delete [] alldata;
 
  // TODO: Make the shortcut3 function generic
  // recursive shortcut until the local parent vector is no longer updated
  while (1) {
    int64_t nf = 0;
    for (int64_t i = 0; i < npprs; i++) {
      std::map<int64_t, int64_t>::iterator it;
      it = m_alldata.find(pprs[i].d);
      if (it != m_alldata.end()) {
        // change the parent (shortcut)
        pprs[i].d = it->second;
      }
      else {
        nf++;
      }
    }
   if (nf == npprs) break; 
  }
  TAU_FSTART(Local_shortcut3_write);
  //Timer t_shortcut_write("Local_shortcut3_write");
  // TODO: can optimize here
  //t_shortcut_write.start();
  p.write(npprs, pprs);
  //t_shortcut_write.stop();
  TAU_FSTOP(Local_shortcut3_write);
  //t_shortcut.stop();
  TAU_FSTOP(Local_shortcut3);
}

// return B where B[i,j] = A[p[i],p[j]], or if P is P[i,j] = p[i], compute B = P^T A P
template<typename T>
Matrix<T>* PTAP(Matrix<T>* A, Vector<int>* p){
  TAU_FSTART(CONNECTIVITY_PTAP);
  //Timer t_cp("CONNECTIVITY_PTAP");
  //t_cp.start();
  int np = p->wrld->np;
  int64_t n = p->len;
  Pair<int> * pprs;
  int64_t npprs;
  //get local part of p
  p->get_local_pairs(&npprs, &pprs);
  assert((npprs <= (n+np-1)/np) && (npprs >= (n/np)));
  assert(A->ncol == n);
  assert(A->nrow == n);
  Pair<T> * A_prs;
  int64_t nprs;
  {
    //map matrix so rows are distributed as elements of p, ensures for each element of p, this process also owns the row of A (A1)
    Matrix<T> A1(n, n, "ij", Partition(1,&np)["i"], Idx_Partition(), SP*(A->is_sparse), *A->wrld, *A->sr);
    A1["ij"] = A->operator[]("ij");
    A1.get_local_pairs(&nprs, &A_prs, true);
    //use fact p and rows of A are distributed cyclically, to compute P^T * A
    for (int64_t i=0; i<nprs; i++){
      A_prs[i].k = (A_prs[i].k/n)*n + pprs[(A_prs[i].k%n)/np].d;
    }
  }
  {
    //map matrix so rows are distributed as elements of p, ensures for each element of p, this process also owns the column of A (A1)
    Matrix<T> A2(n, n, "ij", Partition(1,&np)["j"], Idx_Partition(), SP*(A->is_sparse), *A->wrld, *A->sr);
    //write in P^T A into A2
    A2.write(nprs, A_prs);
    delete [] A_prs;
    A2.get_local_pairs(&nprs, &A_prs, true);
    //use fact p and cols of A are distributed cyclically, to compute P^T A * P
    for (int64_t i=0; i<nprs; i++){
      A_prs[i].k = (A_prs[i].k%n) + pprs[(A_prs[i].k/n)/np].d*n;
    }
  }
  Matrix<T> * PTAP = new Matrix<T>(n, n, SP*(A->is_sparse), *A->wrld, *A->sr);
  PTAP->write(nprs, A_prs);
  delete [] A_prs;
  TAU_FSTOP(CONNECTIVITY_PTAP);
  //t_cp.stop();
  return PTAP;
}
template Matrix<int>* PTAP<int>(Matrix<int>* A, Vector<int>* p);

// bool does not work for some reason
/*
static Monoid<int> OR_STAR(
    1,
    [](int a, int b) { return a & b; },
    MPI_LAND);
*/
static Semiring<int> OR_STAR(
    1,
    [](int a, int b) { return a & b; },
    MPI_LAND,
    1,
    [](int a, int b) { return a * b; } // mult provided for correct accumulation with write
  );

Vector<int> * star_check(Vector<int> * p, Vector<int> * gf) {
  Vector<int> * star = new Vector<int>(p->len, *(p->wrld), OR_STAR);

  int64_t p_npairs;
  Pair<int> * p_loc_pairs;
  p->get_local_pairs(&p_npairs, &p_loc_pairs);

  int64_t gf_npairs = p_npairs;
  Pair<int> * gf_loc_pairs;
  if (gf == NULL) {
    gf_loc_pairs = new Pair<int>[gf_npairs];
    for (int64_t i = 0; i < gf_npairs; ++i) {
      gf_loc_pairs[i].k = p_loc_pairs[i].d;
    } 
    p->read(gf_npairs, gf_loc_pairs);
  } else {
    gf->get_local_pairs(&gf_npairs, &gf_loc_pairs);
  }
  assert(p_npairs == gf_npairs);

  // If F(i) =/= GF(i) then ST(i) <- FALSE and ST(GF(i)) <- FALSE
  // excludes vertices that have nontrivial grandparent or grandchild
  Pair<int> * nontriv_grandX = new Pair<int>[2 * gf_npairs];
  int64_t grandX_npairs = 0;
  for (int64_t i = 0; i < gf_npairs; ++i) {
    if (p_loc_pairs[i].d != gf_loc_pairs[i].d) {
      nontriv_grandX[grandX_npairs].k = p_loc_pairs[i].k;
      nontriv_grandX[grandX_npairs].d = 0;

      nontriv_grandX[grandX_npairs + 1].k = gf_loc_pairs[i].d;
      nontriv_grandX[grandX_npairs + 1].d = 0;

      grandX_npairs += 2;
    }
  }
  star->write(grandX_npairs, nontriv_grandX);

  // ST(i) <- ST(F(i))
  // excludes vertices that have nontrivial nephews
  Pair<int> * nontriv_nephews = new Pair<int>[gf_npairs];
  for (int64_t i = 0; i < gf_npairs; ++i) {
    nontriv_nephews[i].k = p_loc_pairs[i].d;
  }
  star->read(gf_npairs, nontriv_nephews);

  Pair<int> * updated_nephews = new Pair<int>[gf_npairs];
  for (int64_t i = 0; i < gf_npairs; ++i) {
    updated_nephews[i].k = p_loc_pairs[i].k;
    updated_nephews[i].d = nontriv_nephews[i].d;
  }
  star->write(gf_npairs, updated_nephews);

  delete [] updated_nephews;
  delete [] nontriv_nephews;
  delete [] nontriv_grandX;
  delete [] gf_loc_pairs;
  delete [] p_loc_pairs;

  return star;
}