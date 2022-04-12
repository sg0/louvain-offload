#pragma once
#ifndef DSPL_HPP
#define DSPL_HPP

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <list>
#include <numeric>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>

#include <omp.h>

#include "graph.hpp"
#include "utils.hpp"

struct Comm {
  GraphElem size;
  GraphWeight degree;

  Comm() : size(0), degree(0.0) {};
};

struct CommInfo {
    GraphElem community;
    GraphElem size;
    GraphWeight degree;
};

struct clmap_t {
  GraphElem f;
  GraphElem s;
};

#define CLMAP_MAX_NUM 32
#define COUNT_MAX_NUM 32

void sumVertexDegree(const Graph &g, std::vector<GraphWeight> &vDegree, std::vector<Comm> &localCinfo)
{
  const GraphElem nv = g.get_lnv();

#ifdef OMP_SCHEDULE_RUNTIME
#pragma omp parallel for default(none), shared(g, vDegree, localCinfo), schedule(runtime)
#else
#pragma omp parallel for default(none), shared(g, vDegree, localCinfo), firstprivate(nv) schedule(guided)
#endif
  for (GraphElem i = 0; i < nv; i++) {
    GraphElem e0, e1;
    GraphWeight tw = 0.0;

    g.edge_range(i, e0, e1);

    for (GraphElem k = e0; k < e1; k++) {
      const Edge &edge = g.get_edge(k);
      tw += edge.weight_;
    }

    vDegree[i] = tw;
   
    localCinfo[i].degree = tw;
    localCinfo[i].size = 1L;
  }
} // sumVertexDegree

GraphWeight calcConstantForSecondTerm(const std::vector<GraphWeight> &vDegree)
{
  GraphWeight localWeight = 0.0;

  const size_t vsz = vDegree.size();

#ifdef OMP_SCHEDULE_RUNTIME
#pragma omp parallel for default(none), shared(vDegree), reduction(+: localWeight) schedule(runtime)
#else
#pragma omp parallel for default(none), shared(vDegree), firstprivate(vsz), reduction(+: localWeight) schedule(static)
#endif  
  for (GraphElem i = 0; i < vsz; i++)
    localWeight += vDegree[i]; // Local reduction

  return (1.0 / static_cast<GraphWeight>(localWeight));
} // calcConstantForSecondTerm

void initComm(std::vector<GraphElem> &pastComm, std::vector<GraphElem> &currComm, const GraphElem base)
{
  const size_t csz = currComm.size();

#ifdef DEBUG_PRINTF  
  assert(csz == pastComm.size());
#endif

#ifdef OMP_SCHEDULE_RUNTIME
#pragma omp parallel for default(none), shared(pastComm, currComm), schedule(runtime)
#else
#pragma omp parallel for default(none), shared(pastComm, currComm), firstprivate(csz, base), schedule(static)
#endif
  for (GraphElem i = 0L; i < csz; i++) {
    pastComm[i] = i + base;
    currComm[i] = i + base;
  }
} // initComm

void initLouvain(const Graph &g, std::vector<GraphElem> &pastComm, 
        std::vector<GraphElem> &currComm, std::vector<GraphWeight> &vDegree, 
        std::vector<GraphWeight> &clusterWeight, std::vector<Comm> &localCinfo, 
        std::vector<Comm> &localCupdate, GraphWeight &constantForSecondTerm)
{
  const GraphElem nv = g.get_lnv();

  vDegree.resize(nv);
  pastComm.resize(nv);
  currComm.resize(nv);
  clusterWeight.resize(nv);
  localCinfo.resize(nv);
  localCupdate.resize(nv);
 
  distSumVertexDegree(dg, vDegree, localCinfo);
  constantForSecondTerm = distCalcConstantForSecondTerm(vDegree, gcomm);

  initComm(pastComm, currComm, base);
} // initLouvain

GraphElem getMaxIndex(clmap_t *clmap, int &clmap_size, GraphWeight *counter, int &counter_size,
			  const GraphWeight selfLoop, const Comm *localCinfo, const GraphWeight vDegree, 
        const GraphElem currSize, const GraphWeight currDegree, const GraphElem currComm,
			  const GraphWeight constant)
{
  clmap_t *storedAlready;
  GraphElem maxIndex = currComm;
  GraphWeight curGain = 0.0, maxGain = 0.0;
  GraphWeight eix = counter[0] - selfLoop;

  GraphWeight ax = currDegree - vDegree;
  GraphWeight eiy = 0.0, ay = 0.0;

  GraphElem maxSize = currSize; 
  GraphElem size = 0;

  storedAlready = clmap;

  do {
      if (currComm != storedAlready->f) {
          ay = localCinfo[storedAlready->f-base].degree;
          size = localCinfo[storedAlready->f - base].size;   

          if (storedAlready->s < counter_size)
            eiy = counter[storedAlready->s];

          curGain = 2.0 * (eiy - eix) - 2.0 * vDegree * (ay - ax) * constant;

          if ((curGain > maxGain) || ((curGain == maxGain) && (curGain != 0.0) && (storedAlready->f < maxIndex))) {
              maxGain = curGain;
              maxIndex = storedAlready->f;
              maxSize = size;
          }
      }
      storedAlready++;
  } while (storedAlready != clmap + clmap_size);

  if ((maxSize == 1) && (currSize == 1) && (maxIndex > currComm))
    maxIndex = currComm;

  return maxIndex;
} // distGetMaxIndex

GraphWeight buildLocalMapCounter(const GraphElem e0, const GraphElem e1, clmap_t *clmap, int &clmap_size, 
    GraphWeight *counter, int &counter_size, const Edge *edge_list, const GraphElem *currComm,
	  const GraphElem vertex)
{
  GraphElem numUniqueClusters = 1L;
  GraphWeight selfLoop = 0;
  clmap_t *storedAlready;

  for (GraphElem j = e0; j < e1; j++) {
        
    const Edge &edge = edge_list[j];
    const GraphElem &tail_ = edge.tail_;
    const GraphWeight &weight = edge.weight_;
    GraphElem tcomm;

    if (tail_ == vertex)
      selfLoop += weight;

    tcomm = currComm[tail_];

    storedAlready = clmap;
    for (int i = 0; i < clmap_size; i++, storedAlready++) {
      if (clmap[i].f == tcomm)
        break;
    }
    
    if (storedAlready != clmap + clmap_size && storedAlready->s < counter_size)
      counter[storedAlready->s] += weight;
    else {
        if (clmap_size < CLMAP_MAX_NUM) {
          clmap[clmap_size].f = tcomm;
          clmap[clmap_size].s = numUniqueClusters;
          clmap_size++;
        }
        if (counter_size < COUNT_MAX_NUM) {
          counter[counter_size] = weight;
          counter_size++;
        }
        numUniqueClusters++;
    }
  }

  return selfLoop;
} // buildLocalMapCounter

void execLouvainIteration(const GraphElem i, const GraphElem *edge_indices, const Edge *edge_list,
    const GraphElem *currComm, GraphElem *targetComm, const GraphWeight *vDegree, Comm *localCinfo, Comm *localCupdate,
    const GraphWeight constantForSecondTerm, GraphWeight *clusterWeight)
{
  GraphElem localTarget = -1;
  GraphElem e0, e1, selfLoop = 0;
  clmap_t clmap[CLMAP_MAX_NUM];
  int clmap_size = 0;
  GraphWeight counter[COUNT_MAX_NUM];
  int counter_size = 0;

  const GraphElem cc = currComm[i];
  GraphWeight ccDegree;
  GraphElem ccSize;  
  bool currCommIsLocal = false; 
  bool targetCommIsLocal = false;

  ccDegree=localCinfo[cc].degree;
  ccSize=localCinfo[cc].size;

  e0 = edge_indices[i];
  e1 = edge_indices[i+1];

  if (e0 != e1) {
    clmap[0].f = cc;
    clmap[0].s = 0;
    clmap_size++;
    counter[0] = 0.0;
    counter_size++;

    selfLoop =  buildLocalMapCounter(e0, e1, clmap, clmap_size, counter, 
        counter_size, edge_list, currComm, i);

    clusterWeight[i] += counter[0];

    localTarget = getMaxIndex(clmap, clmap_size, counter, counter_size, selfLoop, localCinfo,
                    vDegree[i], ccSize, ccDegree, cc, constantForSecondTerm);
  }
  else
    localTarget = cc;

  if ((localTarget != cc) && (localTarget != -1)) {
        
        #pragma omp atomic update
        localCupdate[localTarget].degree += vDegree[i];
        #pragma omp atomic update
        localCupdate[localTarget].size++;
        #pragma omp atomic update
        localCupdate[cc].degree -= vDegree[i];
        #pragma omp atomic update
        localCupdate[cc].size--;
     }	

#ifdef DEBUG_PRINTF  
  assert(localTarget != -1);
#endif
  targetComm[i] = localTarget;
} // execLouvainIteration

GraphWeight computeModularity(const Graph &g, Comm *localCinfo,
			     const GraphWeight *clusterWeight,
			     const GraphWeight constantForSecondTerm)
{
  const GraphElem nv = g.get_lnv();
  GraphWeight le_xx = 0.0, la2_x = 0.0;

#if defined(USE_OMP_OFFLOAD)
#pragma omp target teams distribute parallel for map(to: clusterWeight[0:nv], localCinfo[0:nv]) reduction(+: le_xx), reduction(+: la2_x)
#elif defined(OMP_SCHEDULE_RUNTIME)
#pragma omp parallel for default(none), shared(clusterWeight, localCinfo), \
  reduction(+: le_xx), reduction(+: la2_x) schedule(runtime)
#else
#pragma omp parallel for default(none), shared(clusterWeight, localCinfo), \
  reduction(+: le_xx), reduction(+: la2_x) schedule(static)
#endif
  for (GraphElem i = 0L; i < nv; i++) {
    le_xx += clusterWeight[i];
    la2_x += localCinfo[i].degree * localCinfo[i].degree; 
  } 

  GraphWeight currMod = (le_xx * constantForSecondTerm) - 
      (la2_x * constantForSecondTerm * constantForSecondTerm);
#ifdef DEBUG_PRINTF  
  std::cout << "[" << me << "]le_xx: " << le_xx << ", la2_x: " << la2_x << std::endl;
#endif

  return currMod;
} // computeModularity

void updateLocalCinfo(const GraphElem nv, Comm *localCinfo, const Comm *localCupdate)
{
#if defined(USE_OMP_OFFLOAD)
#pragma omp target teams distribute parallel for map(to                        \
                                                     : localCupdate [0:nv])    \
    map(tofrom                                                                 \
        : localCinfo [0:nv])
#elif defined(OMP_SCHEDULE_RUNTIME)
#pragma omp for schedule(runtime)
#else
#pragma omp for schedule(static)
#endif
    for (GraphElem i = 0L; i < nv; i++) {
        localCinfo[i].size += localCupdate[i].size;
        localCinfo[i].degree += localCupdate[i].degree;
    }
}

void cleanCWandCU(const GraphElem nv, GraphWeight *clusterWeight,
        Comm *localCupdate)
{
#if defined(USE_OMP_OFFLOAD)
#pragma omp target teams distribute parallel for map(from                      \
                                                     : clusterWeight [0:nv],   \
                                                       localCupdate [0:nv])
#elif defined(OMP_SCHEDULE_RUNTIME)
#pragma omp for schedule(runtime)
#else
#pragma omp for schedule(static)
#endif
    for (GraphElem i = 0L; i < nv; i++) {
        clusterWeight[i] = 0;
        localCupdate[i].degree = 0;
        localCupdate[i].size = 0;
    }
} // distCleanCWandCU

GraphWeight louvainMethod(const Graph &dg, const GraphWeight lower, const GraphWeight thresh, int &iters)
{
  std::vector<GraphElem> pastComm, currComm, targetComm;
  std::vector<GraphWeight> vDegree;
  std::vector<GraphWeight> clusterWeight;
  std::vector<Comm> localCinfo, localCupdate;
 
  std::unordered_map<GraphElem, GraphElem> remoteComm;
  std::map<GraphElem,Comm> remoteCinfo, remoteCupdate;
  
  const GraphElem nv = dg.get_lnv();

  GraphWeight constantForSecondTerm;
  GraphWeight prevMod = lower;
  GraphWeight currMod = -1.0;
  int numIters = 0;
  
  initLouvain(dg, pastComm, currComm, vDegree, clusterWeight, localCinfo, 
          localCupdate, constantForSecondTerm);
  targetComm.resize(nv);

#ifdef DEBUG_PRINTF  
  std::cout << "constantForSecondTerm: " << constantForSecondTerm << std::endl;
  std::cout << "Threshold: " << thresh << std::endl;
#endif

  const GraphElem *d_edge_indices = &dg.edge_indices_[0];
  const GraphElem *d_parts = &dg.parts_[0];
  const Edge *d_edge_list = &dg.edge_list_[0];
  GraphElem *d_currComm = &currComm[0];
  const GraphWeight *d_vDegree = &vDegree[0];
  GraphElem *d_targetComm = &targetComm[0];
  Comm *d_localCinfo = &localCinfo[0];
  Comm *d_localCupdate = &localCupdate[0];
  GraphWeight *d_clusterWeight = &clusterWeight[0];

  double size = sizeof(GraphElem) *
                    (dg.edge_indices_.size() + dg.parts_.size() + 2 * nv) +
                sizeof(Edge) * dg.edge_list_.size() +
                sizeof(GraphWeight) * 2 * nv + sizeof(Comm) * 2 * nv;
  double t_start = omp_get_wtime();

  // start Louvain iteration
  while(true) {
#ifdef DEBUG_PRINTF  
    double t2 = omp_get_wtime();
    std::cout << "Starting Louvain iteration: " << numIters << std::endl;
#endif
    numIters++;

#if defined(USE_OMP_OFFLOAD)
#else
#pragma omp parallel default(none), shared(clusterWeight, localCupdate, currComm, targetComm, \
        vDegree, localCinfo, remoteComm, pastComm, dg), \
    firstprivate(constantForSecondTerm)
#endif
    {
        cleanCWandCU(nv, d_clusterWeight, d_localCupdate);

#if defined(USE_OMP_OFFLOAD)
#pragma omp target teams distribute parallel for map(                          \
    to                                                                         \
    : d_edge_indices [0:dg.edge_indices_.size()],                              \
      d_parts [0:dg.parts_.size()], d_edge_list [0:dg.edge_list_.size()],      \
      d_currComm [0:nv], d_vDegree [0:nv], d_localCinfo [0:nv])                \
    map(from                                                                   \
        : d_targetComm [0:nv])                                                 \
        map(tofrom                                                             \
            : d_localCupdate [0:nv], d_clusterWeight [0:nv])
#elif defined(OMP_SCHEDULE_RUNTIME)
#pragma omp for schedule(runtime)
#else
#pragma omp for schedule(guided) 
#endif
        for (GraphElem i = 0; i < nv; i++) {
            execLouvainIteration(i, d_edge_indices, d_parts, d_edge_list, d_currComm, d_targetComm, d_vDegree, d_localCinfo, 
                    d_localCupdate, constantForSecondTerm, d_clusterWeight, me);
        }
    }

#if defined(USE_OMP_OFFLOAD)
#else
#pragma omp parallel default(none), shared(localCinfo, localCupdate)
#endif
    {
        updateLocalCinfo(nv, d_localCinfo, d_localCupdate);
    }

    // compute modularity
    currMod = computeModularity(dg, d_localCinfo, d_clusterWeight, constantForSecondTerm);

    // exit criteria
    if (currMod - prevMod < thresh)
        break;

    prevMod = currMod;
    if (prevMod < lower)
        prevMod = lower;

#ifdef OMP_SCHEDULE_RUNTIME
#pragma omp parallel for default(none) \
    shared(pastComm, d_currComm, d_targetComm) \
    schedule(runtime)
#else
#pragma omp parallel for default(none) \
    shared(pastComm, d_currComm, d_targetComm) firstprivate(nv) \
    schedule(static)
#endif
    for (GraphElem i = 0; i < nv; i++) {
        GraphElem tmp = pastComm[i];
        pastComm[i] = d_currComm[i];
        d_currComm[i] = d_targetComm[i];
        d_targetComm[i] = tmp;
    }
  } // end of Louvain iteration
  
  std::cout << "Louvain execution time: " << omp_get_wtime() - t_start << std::endl;


  iters = numIters;

  vDegree.clear();
  pastComm.clear();
  currComm.clear();
  targetComm.clear();
  clusterWeight.clear();
  localCinfo.clear();
  localCupdate.clear();
  
  return prevMod;
} // louvainMethod

#endif // __DSPL
