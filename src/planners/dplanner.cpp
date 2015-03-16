/*
 * Copyright (c) 2008, Maxim Likhachev
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of Pennsylvania nor the names of its
 *       contributors may be used to endorse or promote products derived from
 *       this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF DVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */



// Truncated


// SBPL Includes
// #include <sbpl/config.h>
#include <sbpl/planners/dplanner.h>
// #include <sbpl/utils/heap.h>



#include <iostream>
using namespace std;
#include <assert.h>
#include <set>
#include <algorithm>
// #include "../../sbpl/headers.h"


bool _state_print = false;
//bool _state_print = true;
bool use_suspension = true;
//bool use_suspension = false;
//#undef USE_HEUR
#define USE_HEUR 1

bool change_goal = false;
float _eps = 0.0;
int GS[10000];
int length = 0;
//-----------------------------------------------------------------------------------------------------

DPlanner::DPlanner(DiscreteSpaceInformation *environment, bool bForwardSearch) {
  environment_ = environment;

  bforwardsearch = bForwardSearch;

  bsearchuntilfirstsolution = false;
  searchexpands = 0;
  MaxMemoryCounter = 0;
#ifndef ROS
  const char *debug = "debug.txt";
#endif
  fDeb = SBPL_FOPEN(debug, "w");

  if (fDeb == NULL) {
    SBPL_ERROR("ERROR: could not open planner debug file\n");
    throw new SBPL_Exception();
  }

  SBPL_PRINTF("debug on\n");

  pSearchStateSpace_ = new DSearchStateSpace_t;

  //create the D planner
  if (CreateSearchStateSpace(pSearchStateSpace_) != 1) {
    SBPL_ERROR("ERROR: failed to create statespace\n");
    return;
  }

  //set the start and goal states
  if (InitializeSearchStateSpace(pSearchStateSpace_) != 1) {
    SBPL_ERROR("ERROR: failed to create statespace\n");
    return;
  }

  assert(pSearchStateSpace_->searchiteration == 0);
  num_of_expands_initial_solution = 0;
}

DPlanner::~DPlanner() {

  //delete the statespace
  DeleteSearchStateSpace(pSearchStateSpace_);
  delete pSearchStateSpace_;

  SBPL_FCLOSE(fDeb);
}


void DPlanner::Initialize_searchinfo(CMDPSTATE *state, DSearchStateSpace_t *pSearchStateSpace) {

  DState *searchstateinfo = (DState *)state->PlannerSpecificData;

  searchstateinfo->MDPstate = state;
  InitializeSearchStateInfo(searchstateinfo, pSearchStateSpace);
}


CMDPSTATE *DPlanner::CreateState(int stateID, DSearchStateSpace_t *pSearchStateSpace) {
  CMDPSTATE *state = NULL;

#if DEBUG

  if (environment_->StateID2IndexMapping[stateID][DMDP_STATEID2IND] != -1) {
    SBPL_ERROR("ERROR in CreateState: state already created\n");
    throw new SBPL_Exception();
  }

#endif

  //adds to the tail a state
  state = pSearchStateSpace->searchMDP.AddState(stateID);

  //remember the index of the state
  environment_->StateID2IndexMapping[stateID][DMDP_STATEID2IND] = pSearchStateSpace->searchMDP.StateArray.size() - 1;

#if DEBUG

  if (state != pSearchStateSpace->searchMDP.StateArray[environment_->StateID2IndexMapping[stateID][DMDP_STATEID2IND]]) {
    SBPL_ERROR("ERROR in CreateState: invalid state index\n");
    throw new SBPL_Exception();
  }

#endif


  //create search specific info
  state->PlannerSpecificData = (DState *)malloc(sizeof(DState));
  Initialize_searchinfo(state, pSearchStateSpace);
  MaxMemoryCounter += sizeof(DState);

  return state;

}


CMDPSTATE *DPlanner::GetState(int stateID, DSearchStateSpace_t *pSearchStateSpace) {

  if (stateID >= (int)environment_->StateID2IndexMapping.size()) {
    SBPL_ERROR("ERROR int GetState: stateID is invalid\n");
    throw new SBPL_Exception();
  }

  if (environment_->StateID2IndexMapping[stateID][DMDP_STATEID2IND] == -1)
    return CreateState(stateID, pSearchStateSpace);
  else
    return pSearchStateSpace->searchMDP.StateArray[environment_->StateID2IndexMapping[stateID][DMDP_STATEID2IND]];

}



//-----------------------------------------------------------------------------------------------------


CKey DPlanner::ComputeKey(DState *state) {
  CKey key;

  if (state->v >= state->g) {
    key.key[0] = state->g + state->h; //(int)(pSearchStateSpace_->eps*state->h);
    key.key[1] = 1; //state->g;
  } else {
    key.key[0] = state->v + state->h;
    key.key[1] = 0; //state->v;
  }

  return key;
}


int DPlanner::ComputeHeuristic(CMDPSTATE *MDPstate, DSearchStateSpace_t *pSearchStateSpace) {
  //compute heuristic for search
  if (bforwardsearch) {

#if MEM_CHECK == 1
    //int WasEn = DisableMemCheck();
#endif

    //forward search: heur = distance from state to searchgoal which is Goal DState
    int retv =  environment_->GetGoalHeuristic(MDPstate->id);

#if MEM_CHECK == 1
    //if (WasEn)
    //	EnableMemCheck();
#endif

    return retv;
  } else {
    //backward search: heur = distance from searchgoal to state
    return (int)(environment_->GetStartHeuristic(MDPstate->id));
  }

}


//initialization of a state
void DPlanner::InitializeSearchStateInfo(DState *state, DSearchStateSpace_t *pSearchStateSpace) {
  state->g = INFINITECOST;
  state->gp_iter = -1;
  state->gp_iter1 = -1;
  state->g_p = INFINITECOST;
  state->v = INFINITECOST;
  state->iterationclosed = -1;
  state->callnumberaccessed = pSearchStateSpace->callnumber;
  state->bestnextstate = NULL;
  state->costtobestnextstate = INFINITECOST;
  state->heapindex = 0;
  state->listelem[D_INCONS_LIST_ID] = NULL;
  state->numofexpands = 0;
  state->bestpredstate = NULL;
  state->path = NULL;
  //compute heuristics
#if USE_HEUR

  if (pSearchStateSpace->searchgoalstate != NULL)
    state->h = ComputeHeuristic(state->MDPstate, pSearchStateSpace);
  else
    state->h = 0;

#else
  state->h = 0;
#endif


}



//re-initialization of a state
void DPlanner::ReInitializeSearchStateInfo(DState *state, DSearchStateSpace_t *pSearchStateSpace) {
  state->g = INFINITECOST;
  state->gp_iter = -1;
  state->gp_iter1 = -1;
  state->g_p = INFINITECOST;
  state->v = INFINITECOST;
  state->iterationclosed = -1;
  state->callnumberaccessed = pSearchStateSpace->callnumber;
  state->bestnextstate = NULL;
  state->costtobestnextstate = INFINITECOST;
  state->heapindex = 0;
  state->listelem[D_INCONS_LIST_ID] = NULL;
  state->numofexpands = 0;
  state->bestpredstate = NULL;
  state->path = NULL;
  //compute heuristics
#if USE_HEUR

  if (pSearchStateSpace->searchgoalstate != NULL) {
    state->h = ComputeHeuristic(state->MDPstate, pSearchStateSpace);
  } else
    state->h = 0;

#else

  state->h = 0;

#endif


}



void DPlanner::DeleteSearchStateData(DState *state) {
  //no memory was allocated
  MaxMemoryCounter = 0;
  return;
}


void DPlanner::UpdateSetMembership(DState *state) {
  CKey key;

  if (state->v != state->g) {
    if ((state->gp_iter1 == 1209) && (state->v > state->g)) //Check whether the goal path has changed
      change_goal = true;

    key = ComputeKey(state);

    if (state->heapindex == 0) {
      pSearchStateSpace_->heap->insertheap(state, key);
    } else {
      pSearchStateSpace_->heap->updateheap(state, key);
    }
  } else {
    if (state->heapindex != 0) {
      pSearchStateSpace_->heap->deleteheap(state);
    }
  }
}


void DPlanner::SuspendPath(DSearchStateSpace_t *pSearchStateSpace, DState *s) {
  DState *searchstateinfo;
  CMDPSTATE *state = NULL;
  CMDPSTATE *goalstate = NULL;
  CMDPSTATE *startstate = NULL;
  startstate = s->MDPstate;
  goalstate = pSearchStateSpace->searchstartstate;
  state = startstate;
  int steps = 0;
  const int max_steps = 1000;

  if (s->path == NULL)
    s->path = new vector<CMDPSTATE *>;
  else
    s->path->clear();

  while (state->id != goalstate->id && steps < max_steps) {
    steps++;

    if (state->PlannerSpecificData == NULL) {
      assert(0);
      break;
    }

    searchstateinfo = (DState *)state->PlannerSpecificData;

    if ((searchstateinfo->listelem[D_INCONS_LIST_ID] != NULL) && (searchstateinfo != s)) {
      break;
    }

    /*if  ((searchstateinfo->iterationclosed == pSearchStateSpace->searchiteration) && (searchstateinfo!=s)) {
        break;
    }*/
    state = searchstateinfo->bestnextstate;

    if (state == NULL)
      break;

    s->path->push_back(state);
  }

  //assert(steps < max_steps);
}

bool DPlanner::RecomputegprimeStart(DSearchStateSpace_t *pSearchStateSpace, int bound) { //, DState* s)

  //set<CMDPSTATE*> v_state;
  DState *searchstateinfo;
  CMDPSTATE *state = NULL;

  CMDPSTATE *goalstate = NULL;
  CMDPSTATE *startstate = NULL;

  startstate = pSearchStateSpace->searchgoalstate;
  goalstate = pSearchStateSpace->searchstartstate;

  state = startstate;
  DState *s = (DState *) startstate->PlannerSpecificData;

  if (s->path == NULL)
    s->path = new vector<CMDPSTATE *>;
  else
    s->path->clear();

  int steps = 0;
  const int max_steps = 1000;
  unsigned int gcost = 0;

  while (state->id != goalstate->id && steps < max_steps) {
    searchstateinfo = (DState *)state->PlannerSpecificData;
    searchstateinfo->gp_iter1 = 1209;
    GS[steps] = state->id;
    /*if (searchstateinfo->iterationclosed == pSearchStateSpace->searchiteration) {
       gcost += _eps*searchstateinfo->g;
       break;
                } */
    steps++;

    if (use_suspension && searchstateinfo->listelem[D_INCONS_LIST_ID] != NULL) {
      gcost += searchstateinfo->g_p;
      assert(searchstateinfo->path != NULL);
      assert(!searchstateinfo->path->empty());

      break;
    }

    if (searchstateinfo->bestnextstate == NULL) {
      gcost = INFINITECOST;
      break;
    }

    gcost += searchstateinfo->costtobestnextstate;
    state = searchstateinfo->bestnextstate;
    s->path->push_back(state);
  }

  bool f = false;

  if ((steps >= max_steps) || (gcost >= INFINITECOST)) {
    gcost = INFINITECOST;
  } else
    f = true;

  s->g_p = gcost;
  return f;
}


bool seen = false;
void DPlanner::Recomputegprime(DSearchStateSpace_t *pSearchStateSpace, DState *s) {

  //set<CMDPSTATE*> v_state;
  DState *searchstateinfo;
  CMDPSTATE *state = NULL;

  CMDPSTATE *goalstate = NULL;
  CMDPSTATE *startstate = NULL;

  startstate = s->MDPstate;
  goalstate = pSearchStateSpace->searchstartstate;

  state = startstate;
  int steps = 0;
  const int max_steps = 1000;
  unsigned int gcost = 0;
  float p = _eps * (s->v + s->h);

  //float p = _eps*s->v;
  if (s->path == NULL)
    s->path = new vector<CMDPSTATE *>;
  else
    s->path->clear();

  while (state->id != goalstate->id && steps < max_steps) {
    steps++;
    searchstateinfo = (DState *)state->PlannerSpecificData;

    /*if (searchstateinfo->iterationclosed == pSearchStateSpace->searchiteration) {
       gcost += _eps*searchstateinfo->g;
       break;
                } */
    if (use_suspension && searchstateinfo->listelem[D_INCONS_LIST_ID] != NULL) {
      gcost += searchstateinfo->g_p;
      assert(searchstateinfo->path != NULL);
      assert(!searchstateinfo->path->empty());
      break;
    }

    if (searchstateinfo->bestnextstate == NULL) {
      gcost = INFINITECOST;
      break;
    }

    gcost += searchstateinfo->costtobestnextstate;

    if (gcost >= p) break;

    state = searchstateinfo->bestnextstate;
    s->path->push_back(state);
  }

  if ((gcost >= INFINITECOST) || (steps >= max_steps) || (gcost >= p)) {
    gcost = INFINITECOST;
    s->path->clear();
  }

  s->g_p = gcost;

}

void DPlanner::RecomputegprimeTest(DSearchStateSpace_t *pSearchStateSpace, DState *s) {
}

//used for backward search
void DPlanner::MarkPredsofOverconsState(DState *state, DSearchStateSpace_t *pSearchStateSpace) {
}


void DPlanner::Recomputegval(DState *node) {

  vector<nodeStub> *neighbours; //TODO: node should save it's neighbours
  CKey key;
  DState *searchpredstate;

  if (bforwardsearch)
// 		environment_->GetPreds(state->MDPstate->id, &searchpredsIDV, &costV);
    neighbours = environment_->GetPreds(node->MDPstate->id);
  else
    neighbours = environment_->GetSuccs(node->MDPstate->id);

//    environment_->GetSuccs(state->MDPstate->id, &searchpredsIDV, &costV);

  node->g = INFINITECOST;
  node->costtobestnextstate = INFINITECOST;
  node->bestnextstate = NULL;

  for (nodeStub n : *neighbours) {
    if (environment_->StateID2IndexMapping[n.id][DMDP_STATEID2IND] == -1)
      continue; //skip the states that do not exist - they can not be used to improve g-value anyway

    CMDPSTATE *predMDPState = GetState(n.id, pSearchStateSpace_);
    searchpredstate = (DState *)(predMDPState->PlannerSpecificData);
    /*if(predstate->callnumberaccessed != pSearchStateSpace->callnumber) {
    	ReInitializeSearchStateInfo(predstate, pSearchStateSpace);
    } */
    int max = searchpredstate->v;

    if (node->g > max + n.cost) {

      if (bforwardsearch) {
        node->g = max + n.cost;
        node->bestpredstate = predMDPState;
      } else {
        node->g = max + n.cost;
        node->bestnextstate = predMDPState;
        node->costtobestnextstate = n.cost;
      }
    }
  }
}


//used for backward search
void DPlanner::UpdatePredsofOverconsState(DState *node, DSearchStateSpace_t *pSearchStateSpace) {
  vector<nodeStub> *neighbours; //TODO: node should save it's neighbours
  CKey key;
  DState *predstate;

  neighbours = environment_->GetPreds(node->MDPstate->id);
  //printf("State -- %d, number of preds -- %d\n", state->MDPstate->StateID, (int)PredIDV.size());

  for (nodeStub n : *neighbours) {
    CMDPSTATE *PredMDPState = GetState(n.id, pSearchStateSpace);
    predstate = (DState *)(PredMDPState->PlannerSpecificData);
    /*if(predstate->callnumberaccessed != pSearchStateSpace->callnumber) {
    	ReInitializeSearchStateInfo(predstate, pSearchStateSpace);
    } */

    if ((use_suspension) && (predstate->listelem[D_INCONS_LIST_ID] != NULL))
      continue; // Skip suspended

    if ((use_suspension) && (predstate->iterationclosed == pSearchStateSpace->searchiteration))
      continue; // Skip suspended

    assert(node->v == node->g);
    int max = node->v;

    if (predstate->g > max + n.cost) { // && (predstate->listelem[D_INCONS_LIST_ID] == NULL))
      predstate->g = max + n.cost;
      predstate->bestnextstate = node->MDPstate;
      predstate->costtobestnextstate = n.cost;

      if (predstate->gp_iter != pSearchStateSpace->searchiteration) {
#if USE_HEUR

        if (pSearchStateSpace->searchgoalstate != NULL)
          predstate->h = ComputeHeuristic(predstate->MDPstate, pSearchStateSpace);
        else
          predstate->h = 0;

#else
        predstate->h = 0;
#endif
        predstate->gp_iter = pSearchStateSpace->searchiteration;
      }

      UpdateSetMembership(predstate);
    }
  }

}

//used for forward search
void DPlanner::UpdateSuccsofOverconsState(DState *state, DSearchStateSpace_t *pSearchStateSpace) {
}



//used for backward search
void DPlanner::UpdatePredsofUnderconsState(DState *state, DSearchStateSpace_t *pSearchStateSpace) {
  vector<nodeStub> *neighbours;
  CKey key;
  DState *predstate;

  neighbours = environment_->GetPreds(state->MDPstate->id);
  //printf("State -- %d, number of preds -- %d\n", state->MDPstate->StateID, (int)PredIDV.size());
  //iterate through predecessors of s
  vector <DState *> stv;

  for (nodeStub n : *neighbours) {
    CMDPSTATE *PredMDPState = GetState(n.id, pSearchStateSpace);
    predstate = (DState *)(PredMDPState->PlannerSpecificData);

    /*if(predstate->callnumberaccessed != pSearchStateSpace->callnumber)
      ReInitializeSearchStateInfo(predstate, pSearchStateSpace);*/
    if (predstate->bestnextstate == state->MDPstate) {
      if ((use_suspension) && (predstate->listelem[D_INCONS_LIST_ID] != NULL))
        continue; //DO NOTHING TO THE SUSPENDED

      if ((use_suspension) && (predstate->iterationclosed == pSearchStateSpace->searchiteration))
        continue; //DO NOTHING TO THE SUSPENDED

      Recomputegval(predstate);

      if (predstate->gp_iter != pSearchStateSpace->searchiteration) {
        if (pSearchStateSpace->searchgoalstate != NULL)
          predstate->h = ComputeHeuristic(predstate->MDPstate, pSearchStateSpace);
        else
          predstate->h = 0;

        predstate->gp_iter = pSearchStateSpace->searchiteration;
      }

      UpdateSetMembership(predstate);
    }

    if (_state_print) {
      SBPL_FPRINTF(stdout, "updated pred %d of undercons exp\n", predstate->MDPstate->id);
      PrintSearchState(predstate, stdout);
    }

  }
}
void DPlanner::UpdatePredsofUnderconsStateDelayed(DState *state, DSearchStateSpace_t *pSearchStateSpace) {
}




//used for forward search
void DPlanner::UpdateSuccsofUnderconsState(DState *state, DSearchStateSpace_t *pSearchStateSpace) {
}


int DPlanner::GetGVal(int StateID, DSearchStateSpace_t *pSearchStateSpace) {
  CMDPSTATE *cmdp_state = GetState(StateID, pSearchStateSpace);
  DState *state = (DState *)cmdp_state->PlannerSpecificData;
  return state->g;
}



//returns 1 if the solution is found, 0 if the solution does not exist and 2 if it ran out of time
int DPlanner::ComputePath(DSearchStateSpace_t *pSearchStateSpace, double MaxNumofSecs) {
  int expands;
  DState *state, *searchgoalstate;
  CKey key, minkey;
  CKey goalkey;
  expands = 0;

  if (pSearchStateSpace->searchgoalstate == NULL) {
    SBPL_ERROR("ERROR searching: no goal state is set\n");
    throw new SBPL_Exception();
  }

  //goal state
  searchgoalstate = (DState *)(pSearchStateSpace->searchgoalstate->PlannerSpecificData);
  goalkey = ComputeKey(searchgoalstate);
  //expand states until done
  minkey = pSearchStateSpace->heap->getminkeyheap();
  CKey oldkey = minkey;

  while (!pSearchStateSpace->heap->emptyheap() && minkey.key[0] < INFINITECOST && (goalkey > minkey || searchgoalstate->g > searchgoalstate->v) &&
         (clock() - TimeStarted) < MaxNumofSecs * (double)CLOCKS_PER_SEC) {

    state = (DState *)pSearchStateSpace->heap->deleteminheap();

    if (state->v == state->g) {
      printf("--------------------%zu\n", state->MDPstate->id);
      PrintSearchState(state, stdout);
      SBPL_ERROR("ERROR: consistent state is being expanded\n");
      //printf("--------------------\n");
      throw new SBPL_Exception();
    }

    //new expand
    expands++;
    state->numofexpands++;

    /*if (expands < 100) {
    	int wx = -1;
    	int wy = -1;
    	environment_->GetCoordFromState(state->MDPstate->id, wx, wy);
    EXP[wy][wx] = 3;
    }*/
    //printf("Id=%d, y = %d, x = %d\n", state->MDPstate->id, wy, wx);
    if (state->v > state->g) {
      //overconsistent expansion
      //if (_state_print) {

      state->v = state->g;
      state->iterationclosed = pSearchStateSpace->searchiteration;

      if (!bforwardsearch) {
        UpdatePredsofOverconsState(state, pSearchStateSpace);
      } else {
        UpdateSuccsofOverconsState(state, pSearchStateSpace);
      }
    } else {
      //underconsistent expansion
      //force the state to be overconsistent
      //      if ((_state_print) || (state->MDPstate->StateID==-4843731)) {

      state->v = INFINITECOST;
      UpdateSetMembership(state);

      if (!bforwardsearch) {
        UpdatePredsofUnderconsState(state, pSearchStateSpace);
      } else
        UpdateSuccsofUnderconsState(state, pSearchStateSpace);

    }

    //Storing the current g in g_last. Checking purpose
    //recompute minkey
    minkey = pSearchStateSpace->heap->getminkeyheap();
    //recompute goalkey if necessary
    goalkey = ComputeKey(searchgoalstate);

    if (_state_print)
      printf("GOALKEY --- %ld, %ld AND MINKEY -- %ld, %ld\n", goalkey.key[0], goalkey.key[1], minkey.key[0], minkey.key[1]);

    if (expands % 100000 == 0 && expands > 0) {
      SBPL_PRINTF("expands so far=%u\n", expands);
    }

  }

  //printf("searchgoalstate (g, v, g_p) -- %u, %u, %u\n", searchgoalstate->g, searchgoalstate->v, searchgoalstate->g_p);
  //assert (searchgoalstate->g == searchgoalstate->g_p);
  //printf("START OF COMPUTEPATH ------------\n\n");
  //PrintOpenList(pSearchStateSpace);
  //printf("\n\n");

  int retv = 1;

  if (searchgoalstate->g == INFINITECOST && pSearchStateSpace->heap->emptyheap()) {
    SBPL_PRINTF("solution does not exist: search exited because heap is empty\n");

#if DEBUG
    SBPL_FPRINTF(fDeb, "solution does not exist: search exited because heap is empty\n");
#endif

    retv = 0;
  } else if (!pSearchStateSpace->heap->emptyheap() && (goalkey > minkey || searchgoalstate->g > searchgoalstate->v)) {
    SBPL_PRINTF("search exited because it ran out of time\n");
#if DEBUG
    SBPL_FPRINTF(fDeb, "search exited because it ran out of time\n");
#endif
    retv = 2;
  } else if (searchgoalstate->g == INFINITECOST && !pSearchStateSpace->heap->emptyheap()) {
    SBPL_PRINTF("solution does not exist: search exited because all candidates for expansion have infinite heuristics\n");
#if DEBUG
    SBPL_FPRINTF(fDeb, "solution does not exist: search exited because all candidates for expansion have infinite heuristics\n");
#endif
    retv = 0;
  } else {
    SBPL_PRINTF("search exited with a solution for eps=%.3f\n", pSearchStateSpace->eps);
#if DEBUG
    SBPL_FPRINTF(fDeb, "search exited with a solution for eps=%.3f\n", pSearchStateSpace->eps);
#endif
    retv = 1;
  }

  //SBPL_FPRINTF(fDeb, "expanded=%d\n", expands);

  searchexpands += expands;
  assert(retv == 1);
  return retv;
}
int gp_count = 0;

int DPlanner::ComputePathEpsSave(DSearchStateSpace_t *pSearchStateSpace, double MaxNumofSecs) {
  printf("Compute path called with eps -- %0.2f\n", _eps);
  gp_count = 0;
  int expands;
  DState *state, *searchgoalstate;
  CKey key, minkey;
  CKey goalkey;
  expands = 0;

  if (pSearchStateSpace->searchgoalstate == NULL) {
    SBPL_ERROR("ERROR searching: no goal state is set\n");
    throw new SBPL_Exception();
  }

  searchgoalstate = (DState *)(pSearchStateSpace->searchgoalstate->PlannerSpecificData);

  if (_state_print) {
    printf("START OF COMPUTEPATH ------------\n\n");
    PrintOpenList(pSearchStateSpace);
    printf("\n\n");
  }

  goalkey = ComputeKey(searchgoalstate);
  minkey = pSearchStateSpace->heap->getminkeyheap();
  CKey oldkey = minkey;
  int visited_num = 0;

  //expand states until done
  bool hush = false;
  change_goal = true;

  while (!pSearchStateSpace->heap->emptyheap() && minkey.key[0] < INFINITECOST && (goalkey > minkey || ((searchgoalstate->g > searchgoalstate->v) && (searchgoalstate->listelem[D_INCONS_LIST_ID] == NULL)))) {
    state = (DState *)pSearchStateSpace->heap->deleteminheap();
    visited_num++;
    /*if((state->listelem[D_INCONS_LIST_ID] != NULL)) {
    	printf("--------------------%d\n", state->MDPstate->StateID);
    	PrintSearchState(state, stdout);
    	SBPL_ERROR("ERROR: suspended state is being expanded\n");
    	printf("--------------------\n");
    	throw new SBPL_Exception();
    }
    if (pSearchStateSpace->searchiteration == state->iterationclosed){
    	printf("--------------------%d\n", state->MDPstate->StateID);
    	PrintSearchState(state, stdout);
    	SBPL_ERROR("ERROR: closed state is being expanded\n");
    	printf("--------------------\n");
    	throw new SBPL_Exception();
                }
    if(state->v == state->g)
    {
    	printf("--------------------%d\n", state->MDPstate->StateID);
    	PrintSearchState(state, stdout);
    	SBPL_ERROR("ERROR: consistent state is being expanded\n");
    	printf("--------------------\n");
    	throw new SBPL_Exception();
    } */

    int kbound = 0;

    if (state->v > state->g) {
      kbound = (int)(_eps * (state->g + state->h));
    } else {
      kbound = (int)(_eps * (state->v + state->h));
    }

    bool better = false;
    if (change_goal) {
      for (int o = 0; o < 10000; o++)
        GS[o] = -1;

      better = RecomputegprimeStart(pSearchStateSpace, kbound);
      change_goal = false;
    }

    int goal_v = ((DState *) pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g_p ;

    if (goal_v <= kbound) {
      if (searchgoalstate->listelem[D_INCONS_LIST_ID] == NULL) {
        pSearchStateSpace->inconslist->insert(searchgoalstate, D_INCONS_LIST_ID);
      }

      UpdateSetMembership(state);
      //SuspendPath(pSearchStateSpace, state);
      hush = true;
      break;
    }

    if (state->v > state->g) {
      if (_state_print) {
        printf("OVERCONSISTENT EXPANSION :: Id=%zu, ", state->MDPstate->id);
        PrintSearchState(state, stdout);
      }

      state->v = state->g;
      expands++;
      state->numofexpands++;
      UpdatePredsofOverconsState(state, pSearchStateSpace);
      state->iterationclosed = pSearchStateSpace->searchiteration;

    } else {
      Recomputegprime(pSearchStateSpace, state);

      if ((state->g_p + state->h)  >= _eps * (state->v + state->h))  {
        if (_state_print) {
          printf("UNDERCONSISTENT EXPANSION :: Id=%zu, ", state->MDPstate->id);
          PrintSearchState(state, stdout);
        }

        state->v = INFINITECOST;
        expands++;
        state->numofexpands++;
        UpdateSetMembership(state);
        UpdatePredsofUnderconsState(state, pSearchStateSpace);
      } else {
        if (state->listelem[D_INCONS_LIST_ID] == NULL) {
          pSearchStateSpace->inconslist->insert(state, D_INCONS_LIST_ID);
        }

        //SuspendPath(pSearchStateSpace, state);
      }
    }

    goalkey = ComputeKey(searchgoalstate);
    minkey = pSearchStateSpace->heap->getminkeyheap();

    if (expands % 100000 == 0 && expands > 0) {
      SBPL_PRINTF("expands so far=%u\n", expands);
    }

  }

  int retv = 1;

  if (_state_print) {
    printf("END OF COMPUTEPATH ------------\n\n");
    PrintOpenList(pSearchStateSpace);
    printf("\n\n");
  }

  //searchgoalstate->g = searchgoalstate->g_p;
  //searchgoalstate->v = INFINITECOST;
  //UpdateSetMembership(searchgoalstate);
  //assert(searchgoalstate->heapindex == 0);
  /*if(_eps == 1.0) {
  	if (searchgoalstate->g != searchgoalstate->g_p)
  	printf("HUHA");
  } else {
  	searchgoalstate->g = searchgoalstate->g_p;
  	searchgoalstate->v = INFINITECOST;
  	//UpdateSetMembership(searchgoalstate);
  }*/


  if (hush) {
    searchexpands += expands;
    return 1;
  }

  if (searchgoalstate->g == INFINITECOST) {
    SBPL_PRINTF("solution does not exist: search exited heap empty\n");
    retv = 0;
  } else if (searchgoalstate->g == INFINITECOST && pSearchStateSpace->heap->emptyheap() && !hush) {
    SBPL_PRINTF("solution does not exist: search exited because heap is empty\n");

#if DEBUG
    SBPL_FPRINTF(fDeb, "solution does not exist: search exited because heap is empty\n");
#endif

    retv = 0;
  } else if (!pSearchStateSpace->heap->emptyheap() && (goalkey > minkey) && !hush) { // || searchgoalstate->g > searchgoalstate->v))
    SBPL_PRINTF("search exited because it ran out of time\n");
#if DEBUG
    SBPL_FPRINTF(fDeb, "search exited because it ran out of time\n");
#endif
    retv = 2;
  } else if (searchgoalstate->g == INFINITECOST && !pSearchStateSpace->heap->emptyheap() && !hush) {
    PrintOpenList(pSearchStateSpace);
    SBPL_PRINTF("solution does not exist: search exited because all candidates for expansion have infinite heuristics\n");
#if DEBUG
    SBPL_FPRINTF(fDeb, "solution does not exist: search exited because all candidates for expansion have infinite heuristics\n");
#endif
    retv = 0;
  } else {
    SBPL_PRINTF("search exited with a solution for eps=%.3f\n", pSearchStateSpace->eps);
#if DEBUG
    SBPL_FPRINTF(fDeb, "search exited with a solution for eps=%.3f\n", pSearchStateSpace->eps);
#endif
    retv = 1;
  }

  //SBPL_FPRINTF(fDeb, "expanded=%d\n", expands);
  searchexpands += expands;
  //assert ((retv == 1) || (retv == 101));
  return retv;
}






void DPlanner::BuildNewOPENList(DSearchStateSpace_t *pSearchStateSpace) {
  DState *state;
  CKey key;
  CHeap *pheap = pSearchStateSpace->heap;
  CList *pinconslist = pSearchStateSpace->inconslist;

  while (!pheap->emptyheap()) {
    //get the state
    state = (DState *)pheap->deleteminheap();
    pinconslist->insert(state, D_INCONS_LIST_ID); //predstate->listelem[D_INCONS_LIST_ID]
  }

  //move incons into open
  while (pinconslist->firstelement != NULL) {
    state = (DState *)pinconslist->firstelement->liststate;

    //compute f-value
    key = ComputeKey(state);

    //insert into OPEN
    if (state->heapindex == 0)
      pheap->insertheap(state, key);
    else
      pheap->updateheap(state, key); //should never happen, but sometimes it does - somewhere there is a bug TODO

    //remove from INCONS
    pinconslist->remove(state, D_INCONS_LIST_ID);
  }

  pSearchStateSpace->bRebuildOpenList = false;

}
void DPlanner::PrintOpenList(DSearchStateSpace_t *pSearchStateSpace) {
  DState *state;
  CKey key;
  CHeap *pheap = pSearchStateSpace->heap;
  vector<DState *> general;
  printf("Printing The Remaining Queue %d\n", pheap->currentsize);

  while (!pheap->emptyheap()) {
    //get the state
    state = (DState *)pheap->deleteminheap();
    PrintSearchState(state, stdout);
    general.push_back(state);
    //printf("what \n");
    //assert(pheap->emptyheap());
  }

  //move incons into open
  for (int pp = 0; pp < (int) general.size(); pp++) {
    state = general[pp];
    //compute f-value
    key = ComputeKey(state);

    //insert into OPEN
    if (state->heapindex == 0)
      pheap->insertheap(state, key);
    else
      pheap->updateheap(state, key); //should never happen, but sometimes it does - somewhere there is a bug TODO

    //remove from INCONS
  }

  printf("Queue Printing Ended\n");

}

void DPlanner::Reevaluatefvals(DSearchStateSpace_t *pSearchStateSpace) {
  CKey key;
  int i;
  CHeap *pheap = pSearchStateSpace->heap;

#if DEBUG
  SBPL_FPRINTF(fDeb, "re-computing heap priorities\n");
#endif

  //recompute priorities for states in OPEN and reorder it
  for (i = 1; i <= pheap->currentsize; ++i) {
    DState *state = (DState *)pheap->heap[i].heapstate;

    if (state->gp_iter != pSearchStateSpace->searchiteration) {
#if USE_HEUR

      if (pSearchStateSpace->searchgoalstate != NULL)
        state->h = ComputeHeuristic(state->MDPstate, pSearchStateSpace);
      else
        state->h = 0;

#else
      state->h = 0;
#endif
      state->gp_iter = pSearchStateSpace->searchiteration;
    }

    pheap->heap[i].key = ComputeKey(state);
  }

  pheap->makeheap();

  pSearchStateSpace->bReevaluatefvals = false;
}




//creates (allocates memory) search state space
//does not initialize search statespace
int DPlanner::CreateSearchStateSpace(DSearchStateSpace_t *pSearchStateSpace) {

  //create a heap
  pSearchStateSpace->heap = new CHeap;
  pSearchStateSpace->inconslist = new CList;
  MaxMemoryCounter += sizeof(CHeap);
  MaxMemoryCounter += sizeof(CList);

  pSearchStateSpace->searchgoalstate = NULL;
  pSearchStateSpace->searchstartstate = NULL;
  pSearchStateSpace->searchiteration = 0;

  searchexpands = 0;


  pSearchStateSpace->bReinitializeSearchStateSpace = false;

  return 1;
}

//deallocates memory used by SearchStateSpace
void DPlanner::DeleteSearchStateSpace(DSearchStateSpace_t *pSearchStateSpace) {
  if (pSearchStateSpace->heap != NULL) {
    pSearchStateSpace->heap->makeemptyheap();
    delete pSearchStateSpace->heap;
    pSearchStateSpace->heap = NULL;
  }

  if (pSearchStateSpace->inconslist != NULL) {
    pSearchStateSpace->inconslist->makeemptylist(D_INCONS_LIST_ID);
    delete pSearchStateSpace->inconslist;
    pSearchStateSpace->inconslist = NULL;
  }

  //delete the states themselves
  int iend = (int)pSearchStateSpace->searchMDP.StateArray.size();

  for (int i = 0; i < iend; i++) {
    CMDPSTATE *state = pSearchStateSpace->searchMDP.StateArray[i];
    DeleteSearchStateData((DState *)state->PlannerSpecificData);
    free(state->PlannerSpecificData); // allocated with malloc() on line 199 of revision 19485
    state->PlannerSpecificData = NULL;
  }

  pSearchStateSpace->searchMDP.Delete();
  environment_->StateID2IndexMapping.clear();
}



//reset properly search state space
//needs to be done before deleting states
int DPlanner::ResetSearchStateSpace(DSearchStateSpace_t *pSearchStateSpace) {
  pSearchStateSpace->heap->makeemptyheap();
  pSearchStateSpace->inconslist->makeemptylist(D_INCONS_LIST_ID);

  return 1;
}

//initialization before each search
void DPlanner::ReInitializeSearchStateSpace(DSearchStateSpace_t *pSearchStateSpace) {
  CKey key;

  //increase callnumber
  pSearchStateSpace->callnumber++;

  //reset iteration
  //pSearchStateSpace->searchiteration = 0;


#if DEBUG
  SBPL_FPRINTF(fDeb, "reinitializing search state-space (new call number=%d search iter=%d)\n",
               pSearchStateSpace->callnumber, pSearchStateSpace->searchiteration);
#endif



  pSearchStateSpace->heap->makeemptyheap();
  pSearchStateSpace->inconslist->makeemptylist(D_INCONS_LIST_ID);

  //reset
  pSearchStateSpace->eps = this->finitial_eps;
  pSearchStateSpace->eps_satisfied = INFINITECOST;

  //initialize start state
  DState *startstateinfo = (DState *)(pSearchStateSpace->searchstartstate->PlannerSpecificData);

  if (startstateinfo->callnumberaccessed != pSearchStateSpace->callnumber)
    ReInitializeSearchStateInfo(startstateinfo, pSearchStateSpace);

  startstateinfo->g = 0;

  //insert start state into the heap
  key = ComputeKey(startstateinfo);
  pSearchStateSpace->heap->insertheap(startstateinfo, key);

  pSearchStateSpace->bReinitializeSearchStateSpace = false;
  pSearchStateSpace->bReevaluatefvals = false;
  pSearchStateSpace->bRebuildOpenList = false;
}

//very first initialization
int DPlanner::InitializeSearchStateSpace(DSearchStateSpace_t *pSearchStateSpace) {

  if (pSearchStateSpace->heap->currentsize != 0 ||
      pSearchStateSpace->inconslist->currentsize != 0) {
    SBPL_ERROR("ERROR in InitializeSearchStateSpace: heap or list is not empty\n");
    throw new SBPL_Exception();
  }

  pSearchStateSpace->eps = this->finitial_eps;
  pSearchStateSpace->eps_satisfied = INFINITECOST;
  pSearchStateSpace->searchiteration = 0;
  pSearchStateSpace->callnumber = 0;
  pSearchStateSpace->bReevaluatefvals = false;
  pSearchStateSpace->bRebuildOpenList = false;


  //create and set the search start state
  pSearchStateSpace->searchgoalstate = NULL;
  //pSearchStateSpace->searchstartstate = GetState(SearchStartStateID, pSearchStateSpace);
  pSearchStateSpace->searchstartstate = NULL;


  pSearchStateSpace->bReinitializeSearchStateSpace = true;

  return 1;

}


int DPlanner::SetSearchGoalState(StateID SearchGoalStateID, DSearchStateSpace_t *pSearchStateSpace) {

  if (pSearchStateSpace->searchgoalstate == NULL ||
      pSearchStateSpace->searchgoalstate->id != SearchGoalStateID) {
    pSearchStateSpace->searchgoalstate = GetState(SearchGoalStateID, pSearchStateSpace);

    //current solution may be invalid
    pSearchStateSpace->eps_satisfied = INFINITECOST;
    pSearchStateSpace_->eps = this->finitial_eps;

    //recompute heuristic for the heap if heuristics is used
#if USE_HEUR
    int i;

    //TODO - should get rid of and instead use iteration to re-compute h-values online as needed
    for (i = 0; i < (int)pSearchStateSpace->searchMDP.StateArray.size(); i++) {
      CMDPSTATE *MDPstate = pSearchStateSpace->searchMDP.StateArray[i];
      DState *state = (DState *)MDPstate->PlannerSpecificData;
      state->h = ComputeHeuristic(MDPstate, pSearchStateSpace);
    }

#if DEBUG
    SBPL_PRINTF("re-evaluated heuristic values for %d states\n", i);
#endif

    pSearchStateSpace->bReevaluatefvals = true;
#endif
  }


  return 1;

}


int DPlanner::SetSearchStartState(StateID SearchStartStateID, DSearchStateSpace_t *pSearchStateSpace) {
  CMDPSTATE *MDPstate = GetState(SearchStartStateID, pSearchStateSpace);

  if (MDPstate !=  pSearchStateSpace->searchstartstate) {
    pSearchStateSpace->searchstartstate = MDPstate;
    pSearchStateSpace->bReinitializeSearchStateSpace = true;
    pSearchStateSpace->bRebuildOpenList = true;
  }

  return 1;

}



int DPlanner::ReconstructPath(DSearchStateSpace_t *pSearchStateSpace) {

  //nothing to do, if search is backward
  if (bforwardsearch) {

    CMDPSTATE *MDPstate = pSearchStateSpace->searchgoalstate;
    CMDPSTATE *PredMDPstate;
    DState *predstateinfo, *stateinfo;

    int steps = 0;
    const int max_steps = 100000;

    while (MDPstate != pSearchStateSpace->searchstartstate && steps < max_steps) {
      steps++;

      stateinfo = (DState *)MDPstate->PlannerSpecificData;

      if (stateinfo->g == INFINITECOST) {
        //SBPL_ERROR("ERROR in ReconstructPath: g of the state on the path is INFINITE\n");
        //throw new SBPL_Exception();
        return -1;
      }

      if (stateinfo->bestpredstate == NULL) {
        SBPL_ERROR("ERROR in ReconstructPath: bestpred is NULL\n");
        throw new SBPL_Exception();
      }

      //get the parent state
      PredMDPstate = stateinfo->bestpredstate;
      predstateinfo = (DState *)PredMDPstate->PlannerSpecificData;

      //set its best next info
      predstateinfo->bestnextstate = MDPstate;

      //check the decrease of g-values along the path
      if (predstateinfo->v >= stateinfo->g) {
        SBPL_ERROR("ERROR in ReconstructPath: g-values are non-decreasing\n");
        throw new SBPL_Exception();
      }

      //transition back
      MDPstate = PredMDPstate;
    }

    if (MDPstate != pSearchStateSpace->searchstartstate) {
      SBPL_ERROR("ERROR: Failed to reconstruct path (compute bestnextstate pointers): steps processed=%d\n", steps);
      return 0;
    }
  }

  return 1;
}


void DPlanner::PrintSearchState(DState *searchstateinfo, FILE *fOut) {

  CKey key = ComputeKey(searchstateinfo);
  SBPL_FPRINTF(fOut, "ID = %d, g=%d, v=%d, g_prime=%d, h=%d, key=[%d %d]\n",
               searchstateinfo->MDPstate->id, searchstateinfo->g, searchstateinfo->v, searchstateinfo->g_p, searchstateinfo->h,
               (int)key[0], (int)key[1]);

}

void DPlanner::PrintSearchPath(DSearchStateSpace_t *pSearchStateSpace, FILE *fOut) {
  DState *searchstateinfo;
  CMDPSTATE *state = pSearchStateSpace->searchgoalstate;
  CMDPSTATE *nextstate = NULL;

  if (fOut == NULL)
    fOut = stdout;

  int PathCost = ((DState *)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g;

  SBPL_FPRINTF(fOut, "Printing a path from state %d to the search start state %d\n",
               state->id, pSearchStateSpace->searchstartstate->id);
  SBPL_FPRINTF(fOut, "Path cost = %d:\n", PathCost);

  environment_->PrintState(state->id, true, fOut);

  int costFromStart = 0;
  int steps = 0;
  const int max_steps = 100000;

  while (state->id != pSearchStateSpace->searchstartstate->id && steps < max_steps) {
    steps++;

    SBPL_FPRINTF(fOut, "state %d ", state->id);

    if (state->PlannerSpecificData == NULL) {
      SBPL_FPRINTF(fOut, "path does not exist since search data does not exist\n");
      break;
    }

    searchstateinfo = (DState *)state->PlannerSpecificData;

    if (bforwardsearch)
      nextstate = searchstateinfo->bestpredstate;
    else
      nextstate = searchstateinfo->bestnextstate;

    if (nextstate == NULL) {
      SBPL_FPRINTF(fOut, "path does not exist since nextstate == NULL\n");
      break;
    }

    /*if(searchstateinfo->g == INFINITECOST)
    {
    SBPL_FPRINTF(fOut, "path does not exist since state->g == NULL\n");
    break;
    }
         */
    int costToGoal = PathCost - costFromStart;

    if (!bforwardsearch) {
      //otherwise this cost is not even set
      costFromStart += searchstateinfo->costtobestnextstate;
    }


#if DEBUG
    /* if(searchstateinfo->g > searchstateinfo->v){
       printf("State -- %d\n", state->StateID);
    SBPL_FPRINTF(fOut, "ERROR: underconsistent state %d is encountered\n", state->StateID);
    throw new SBPL_Exception();
     }

     if(!bforwardsearch) //otherwise this cost is not even set
    {
    if(nextstate->PlannerSpecificData != NULL && searchstateinfo->g < searchstateinfo->costtobestnextstate + ((DState*)(nextstate->PlannerSpecificData))->g)
     {
       SBPL_FPRINTF(fOut, "ERROR: g(source) < c(source,target) + g(target)\n");
       throw new SBPL_Exception();
     }
    }
    */
#endif

    //PrintSearchState(searchstateinfo, fOut);
    SBPL_FPRINTF(fOut, "-->state %d ctg = %d  ",
                 nextstate->id, costToGoal);

    state = nextstate;

    environment_->PrintState(state->id, true, fOut);

  }

  if (state->id != pSearchStateSpace->searchstartstate->id) {
    SBPL_ERROR("ERROR: Failed to printsearchpath, max_steps reached\n");
    return;
  }

}

int DPlanner::getHeurValue(DSearchStateSpace_t *pSearchStateSpace, int StateID) {
  CMDPSTATE *MDPstate = GetState(StateID, pSearchStateSpace);
  DState *searchstateinfo = (DState *)MDPstate->PlannerSpecificData;
  return searchstateinfo->h;
}

void DPlanner::FollowPath(DSearchStateSpace_t *pSearchStateSpace, DState *s, int &solcost, vector<StateID> &Ids) {

  //assert(s->listelem[D_INCONS_LIST_ID] != NULL);
  //printf("NEW FOLLOWPATH \n");
// PrintSearchState(s, stdout);
  if ((!s->path) || (s->path->empty())) {
    //printf("NO PATH FOUND\n");
    //PrintSearchState(s, stdout);
    return;
  }

  CMDPSTATE *cur = s->MDPstate;

  vector<nodeStub> *neighbours;
  for (int ii = 0; ii < (int) s->path->size(); ii++) {
    CMDPSTATE *ss = (*s->path)[ii];
    neighbours = environment_->GetSuccs(cur->id);
    uint actioncost = INFINITECOST;


    for (nodeStub n : *neighbours)
      if (n.id == ss->id && n.cost < actioncost)
        actioncost = n.cost;

    solcost += actioncost;
    Ids.push_back(ss->id);
    cur = ss;
    //PrintSearchState((DState*)ss->PlannerSpecificData, stdout);
  }

  FollowPath(pSearchStateSpace, (DState *)cur->PlannerSpecificData, solcost, Ids);
}
vector<StateID> DPlanner::GetSearchPathwithSus(DSearchStateSpace_t *pSearchStateSpace, int &solcost) {
  vector<StateID> wholePathIds;
  DState *searchstateinfo;
  CMDPSTATE *state = NULL;
  CMDPSTATE *goalstate = NULL;
  CMDPSTATE *startstate = NULL;

  if (bforwardsearch) {
    startstate = pSearchStateSpace->searchstartstate;
    goalstate = pSearchStateSpace->searchgoalstate;

    //reconstruct the path by setting bestnextstate pointers appropriately
    if (ReconstructPath(pSearchStateSpace) != 1) {
      solcost = INFINITECOST;
      return wholePathIds;
    }
  } else {
    startstate = pSearchStateSpace->searchgoalstate;
    goalstate = pSearchStateSpace->searchstartstate;
  }


  //PrintSearchPath(pSearchStateSpace, stdout);
#if DEBUG
  //PrintSearchPath(pSearchStateSpace, fDeb);
#endif


  state = startstate;

  wholePathIds.push_back(state->id);
  solcost = 0;

  FILE *fOut = stdout;

  if (fOut == NULL) {
    SBPL_ERROR("ERROR: could not open file\n");
    throw new SBPL_Exception();
  }

  int steps = 0;
  const int max_steps = 100000;

  while (state->id != goalstate->id && steps < max_steps) {
    steps++;

    if (state->PlannerSpecificData == NULL) {
      SBPL_FPRINTF(fOut, "path does not exist since search data does not exist\n");
      break;
    }

    searchstateinfo = (DState *)state->PlannerSpecificData;

    if (searchstateinfo->listelem[D_INCONS_LIST_ID] != NULL) {
      FollowPath(pSearchStateSpace, searchstateinfo, solcost, wholePathIds);
      return wholePathIds;
    }

    if (searchstateinfo->bestnextstate == NULL) {
      assert(0);
      SBPL_FPRINTF(fOut, "path does not exist since bestnextstate == NULL\n");
      break;
    }

    if (searchstateinfo->g == INFINITECOST) {
      SBPL_FPRINTF(fOut, "path does not exist since bestnextstate == NULL\n");
      assert(0);
      break;
    }

    vector<nodeStub> *neighbours;
    neighbours = environment_->GetSuccs(state->id);
    uint actioncost = INFINITECOST;

    for (nodeStub n : *neighbours)
      if (n.id == searchstateinfo->bestnextstate->id && n.cost < actioncost)
        actioncost = n.cost;


    if (actioncost == INFINITECOST) {
      printf("STATE (%zu), g --%d, g_p -- %d, v --%d, ratio1 -- %0.2f, ratio2 --%0.2f, next(%zu), cost -- %d\n", searchstateinfo->MDPstate->id, searchstateinfo->g, searchstateinfo->g_p, searchstateinfo->v, (searchstateinfo->g_p + 0.0) / searchstateinfo->g, (searchstateinfo->g_p + 0.0) / searchstateinfo->v, searchstateinfo->bestnextstate->id, actioncost);
      assert(0);
    }


    solcost += actioncost;

    if (_eps == 0.0) {
      if (searchstateinfo->v < searchstateinfo->g) {
        SBPL_ERROR("ERROR: underconsistent state on the path\n");

        //printf("%d\n",state->StateID);
        //PrintSearchState(searchstateinfo, stdout);
        //printf("UnderConsistency ratio --- %0.3f\n", (searchstateinfo->g + 0.0)/searchstateinfo->g_last);
        //PrintSearchState(searchstateinfo, fDeb);
        throw new SBPL_Exception();
      }
    } else {
      if (searchstateinfo->v < searchstateinfo->g) {
        if (_state_print) {
          SBPL_ERROR("ERROR: underconsistent state on the path\n");
          printf("%zu\n", state->id);
          PrintSearchState(searchstateinfo, stdout);
          RecomputegprimeTest(pSearchStateSpace, searchstateinfo);

          if (searchstateinfo->g_p > (_eps + 0.001)*searchstateinfo->g) {
            printf("THIS IS PROBLEM\n");
            //RecomputegprimeTest(pSearchStateSpace, searchstateinfo);
            assert(0);
          }
        }

        /*RecomputegprimeTest(pSearchStateSpace, searchstateinfo);
        if (searchstateinfo->g_p > (_eps+0.001)*searchstateinfo->g) {
        	printf("THIS IS PROBLEM\n");
        	//RecomputegprimeTest(pSearchStateSpace, searchstateinfo);
        	assert(0);
        } */
        //Recomputegprime(pSearchStateSpace, searchstateinfo);
        //assert(searchstateinfo->g_p <= (_eps + 0.01)*searchstateinfo->v); //printf("UnderConsistency ratio --- %0.3f\n", (searchstateinfo->g_p + 0.0)/searchstateinfo->g);
        //PrintSearchState(searchstateinfo, fDeb);
        //throw new SBPL_Exception();
      }


    }

    //
    //SBPL_FPRINTF(fDeb, "actioncost=%d between states %d and %d\n",
    //        actioncost, state->StateID, searchstateinfo->bestnextstate->StateID);
    //environment_->PrintState(state->StateID, false, fDeb);
    //environment_->PrintState(searchstateinfo->bestnextstate->StateID, false, fDeb);


    state = searchstateinfo->bestnextstate;

    wholePathIds.push_back(state->id);
  }

  //cin.get();
  if (state->id != goalstate->id) {
    SBPL_ERROR("ERROR: Failed to getsearchpath, steps processed=%d\n", steps);
    wholePathIds.clear();
    solcost = INFINITECOST;
    return wholePathIds;
  }

  return wholePathIds;
}


vector<StateID> DPlanner::GetSearchPath(DSearchStateSpace_t *pSearchStateSpace, int &solcost) {
  vector<StateID> wholePathIds;
  DState *searchstateinfo;
  CMDPSTATE *state = NULL;
  CMDPSTATE *goalstate = NULL;
  CMDPSTATE *startstate = NULL;

  if (bforwardsearch) {
    startstate = pSearchStateSpace->searchstartstate;
    goalstate = pSearchStateSpace->searchgoalstate;

    //reconstruct the path by setting bestnextstate pointers appropriately
    if (ReconstructPath(pSearchStateSpace) != 1) {
      solcost = INFINITECOST;
      return wholePathIds;
    }
  } else {
    startstate = pSearchStateSpace->searchgoalstate;
    goalstate = pSearchStateSpace->searchstartstate;
  }

  bool print_path = false;
  //if (pSearchStateSpace->searchiteration == 16)
//		print_path = true;
#if DEBUG
  //PrintSearchPath(pSearchStateSpace, fDeb);
#endif

  //PrintSearchPath(pSearchStateSpace, stdout);

  state = startstate;

  wholePathIds.push_back(state->id);
  solcost = 0;

  FILE *fOut = stdout;

  if (fOut == NULL) {
    SBPL_ERROR("ERROR: could not open file\n");
    throw new SBPL_Exception();
  }

  int steps = 0;
  const int max_steps = 100000;

  while (state->id != goalstate->id && steps < max_steps) {
    steps++;

    if (state->PlannerSpecificData == NULL) {
      SBPL_FPRINTF(fOut, "path does not exist since search data does not exist\n");
      break;
    }

    searchstateinfo = (DState *)state->PlannerSpecificData;

    if (searchstateinfo->bestnextstate == NULL) {
      assert(0);
      SBPL_FPRINTF(fOut, "path does not exist since bestnextstate == NULL\n");
      break;
    }

    if (searchstateinfo->g == INFINITECOST) {
      SBPL_FPRINTF(fOut, "path does not exist since bestnextstate == NULL\n");
      assert(0);
      break;
    }

    /*if(searchstateinfo->g_p == INFINITECOST)
     {
      SBPL_FPRINTF(fOut, "path does not exist since g_p == INFINITY\n");
      assert(0);
      break;
     }*/
    /*int wx = -1;
    int wy = -1;
    environment_->GetCoordFromState(state->StateID, wx, wy);
    //	  printf("Id=%d, y = %d, x = %d\n", state->StateID, wy, wx);
    EXP[wy][wx]= 5;
    */
    if (print_path)
      printf("STATE IN PATH (%zu), g -- %d, v -- %d, g_p -- %d, -- (%d,%d) -- ratio -- %0.3f\n", searchstateinfo->MDPstate->id, searchstateinfo->g, searchstateinfo->v,  searchstateinfo->g_p, searchstateinfo->iterationclosed, pSearchStateSpace->searchiteration, (searchstateinfo->g_p + 0.0) / searchstateinfo->g);

    vector<nodeStub> *neighbours;
    neighbours = environment_->GetSuccs(state->id);
    uint actioncost = INFINITECOST;

    for (nodeStub n : *neighbours)
      if (n.id == searchstateinfo->bestnextstate->id && n.cost < actioncost)
        actioncost = n.cost;

    if (actioncost == INFINITECOST) {
      printf("STATE (%zu), g --%d, g_p -- %d, v --%d, ratio1 -- %0.2f, ratio2 --%0.2f, next(%zu), cost -- %d\n", searchstateinfo->MDPstate->id, searchstateinfo->g, searchstateinfo->g_p, searchstateinfo->v, (searchstateinfo->g_p + 0.0) / searchstateinfo->g, (searchstateinfo->g_p + 0.0) / searchstateinfo->v, searchstateinfo->bestnextstate->id, actioncost);
      assert(0);
    }

    solcost += actioncost;

    //printf("SOL COST -- %d\n", solcost);
    if (_eps == 0.0) {
      if (searchstateinfo->v < searchstateinfo->g) {
        SBPL_ERROR("ERROR: underconsistent state on the path\n");

        //printf("%d\n",state->StateID);
        //PrintSearchState(searchstateinfo, stdout);
        //printf("UnderConsistency ratio --- %0.3f\n", (searchstateinfo->g + 0.0)/searchstateinfo->g_last);
        //PrintSearchState(searchstateinfo, fDeb);
        throw new SBPL_Exception();
      }
    } else {
      if (searchstateinfo->v < searchstateinfo->g) {
        if (_state_print) {
          SBPL_ERROR("ERROR: underconsistent state on the path\n");
          printf("%zu\n", state->id);
          PrintSearchState(searchstateinfo, stdout);
          RecomputegprimeTest(pSearchStateSpace, searchstateinfo);

          if (searchstateinfo->g_p > (_eps + 0.001)*searchstateinfo->g) {
            printf("THIS IS PROBLEM\n");
            //RecomputegprimeTest(pSearchStateSpace, searchstateinfo);
            assert(0);
          }
        }

        /*RecomputegprimeTest(pSearchStateSpace, searchstateinfo);
        if (searchstateinfo->g_p > (_eps+0.001)*searchstateinfo->g) {
        	printf("THIS IS PROBLEM\n");
        	//RecomputegprimeTest(pSearchStateSpace, searchstateinfo);
        	assert(0);
        } */
        //Recomputegprime(pSearchStateSpace, searchstateinfo);
        //assert(searchstateinfo->g_p <= (_eps + 0.01)*searchstateinfo->v); //printf("UnderConsistency ratio --- %0.3f\n", (searchstateinfo->g_p + 0.0)/searchstateinfo->g);
        //PrintSearchState(searchstateinfo, fDeb);
        //throw new SBPL_Exception();
      }


    }

    //
    //SBPL_FPRINTF(fDeb, "actioncost=%d between states %d and %d\n",
    //        actioncost, state->StateID, searchstateinfo->bestnextstate->StateID);
    //environment_->PrintState(state->StateID, false, fDeb);
    //environment_->PrintState(searchstateinfo->bestnextstate->StateID, false, fDeb);
    state = searchstateinfo->bestnextstate;
    wholePathIds.push_back(state->id);
  }

//cin.get();
  if (state->id != goalstate->id) {
    SBPL_ERROR("ERROR: Failed to getsearchpath, steps processed=%d\n", steps);
    wholePathIds.clear();
    solcost = INFINITECOST;
    return wholePathIds;
  }

  /*int wx = -1;
  int wy = -1;
  environment_->GetCoordFromState(goalstate->StateID, wx, wy);
  //printf("Id=%d, y = %d, x = %d\n", state->StateID, wy, wx);
  EXP[wy][wx]= 5;
  */
  return wholePathIds;
}


//Only applicable for backward search now
bool DPlanner::FindOrRaisePath(DSearchStateSpace_t *pSearchStateSpace) {
  return false;
}



bool DPlanner::Search(DSearchStateSpace_t *pSearchStateSpace, vector<StateID> &pathIds, int &PathCost, bool bFirstSolution, bool bOptimalSolution, double MaxNumofSecs) {
  CKey key;
  searchexpands = 0;
#if DEBUG
  SBPL_FPRINTF(fDeb, "new search call (call number=%d)\n", pSearchStateSpace->callnumber);
#endif

  if (pSearchStateSpace->bReinitializeSearchStateSpace == true) {
    //re-initialize state space
    ReInitializeSearchStateSpace(pSearchStateSpace);
  }

  MaxNumofSecs = INFINITECOST;

  //ensure heuristics are up-to-date
  environment_->EnsureHeuristicsUpdated((bforwardsearch == true));

  //the main loop of D*
  _eps = this->_eps;

  TimeStarted = clock();
  int prevexpands = 0;
  clock_t loop_time;
  int pop = -1;

  while ((clock() - TimeStarted) < MaxNumofSecs * (double)CLOCKS_PER_SEC) {
    loop_time = clock();
    //it will be a new search iteration
    //decrease eps for all subsequent iterations
    //pSearchStateSpace->searchiteration++;
    //printf("SEARCH ITERATION --- %d\n", pSearchStateSpace->searchiteration);
    //cin.get();
    //build a new open list by merging it with incons one
    /*if(pSearchStateSpace->bRebuildOpenList)
                  */
    //re-compute f-values if necessary and reorder the heap
    //BuildNewOPENList(pSearchStateSpace);
    //if(pSearchStateSpace->bReevaluatefvals)

    Reevaluatefvals(pSearchStateSpace);

    //improve or compute path
    bool goal_computed = false;

    if ((_eps == 0.0)) {
      //printf("EPS ___ %0.2f\n", _eps);
      if (ComputePath(pSearchStateSpace, MaxNumofSecs) == 1) {
        goal_computed = true;
      }
    } else {
      pop = -1;
      pop = ComputePathEpsSave(pSearchStateSpace, MaxNumofSecs);

      if ((pop  == 1) || (pop == 101)) {
        goal_computed = true;
      }
    }

    SBPL_PRINTF("expands=%d g(sstart)=%d\n", searchexpands - prevexpands,
                ((DState *)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g);
    prevexpands = searchexpands;

    //if just the first solution then we are done
    if (bFirstSolution || goal_computed)
      break;

    //no solution exists
    if (((DState *)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g == INFINITECOST)
      break;

    fflush(stdout);
  }

  double yy = (clock() - TimeStarted) / ((double)CLOCKS_PER_SEC);

#if DEBUG
  SBPL_FFLUSH(fDeb);
#endif

  PathCost = ((DState *)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g;
  MaxMemoryCounter += environment_->StateID2IndexMapping.size() * sizeof(int);

  SBPL_PRINTF("MaxMemoryCounter = %d\n", MaxMemoryCounter);

  int solcost = INFINITECOST;
  bool ret = false;

  if (PathCost == INFINITECOST) { // || pSearchStateSpace_->eps_satisfied == INFINITECOST)
    SBPL_PRINTF("could not find a solution\n");
    ret = false;
  } else {
    //SBPL_PRINTF("solution is found\n");
    if ((use_suspension) && (_eps != 0.0)) {
      pathIds = GetSearchPathwithSus(pSearchStateSpace, solcost);
      PathCost = solcost;
    }

    ret = true;
  }

  solcost = PathCost;
  SBPL_PRINTF("TDSTAR -- total expands this call = %d, planning time = %.3f secs, solution cost=%d\n",
              searchexpands, yy , solcost);

  return ret;

}


bool DPlanner::SearchDelayed(DSearchStateSpace_t *pSearchStateSpace, vector<int> &pathIds, int &PathCost, bool bFirstSolution, bool bOptimalSolution, double MaxNumofSecs) {
  return 1;

}


void DPlanner::Update_SearchSuccs_of_ChangedEdges(vector<StateID> const *statesIDV) {
  printf("Search iteration -- %d\n", pSearchStateSpace_->searchiteration);
  printf("Search callnumber -- %d\n", pSearchStateSpace_->callnumber);
  pSearchStateSpace_->searchiteration++;
  printf("Search iteration -- %d\n", pSearchStateSpace_->searchiteration);

  if (_state_print)
    SBPL_PRINTF("updating %d affected states\n", (unsigned int)statesIDV->size());


  //now we really do need to update it
  if ((_eps == 0.0)) {
    int numofstatesaffected = 0;

    for (auto stateID : *statesIDV) {

      //first check that the state exists (to avoid creation of additional states)
      if (environment_->StateID2IndexMapping[stateID][ADMDP_STATEID2IND] == -1)
        continue;

      //now get the state
      CMDPSTATE *state = GetState(stateID, pSearchStateSpace_);
      DState *searchstateinfo = (DState *)state->PlannerSpecificData;

      //now check that the state is not start state and was created after last search reset
      if (stateID != pSearchStateSpace_->searchstartstate->id) // && searchstateinfo->callnumberaccessed == pSearchStateSpace_->callnumber)
        //if(stateID != pSearchStateSpace_->searchstartstate->StateID  && searchstateinfo->callnumberaccessed == pSearchStateSpace_->callnumber)
      {
#if DEBUG
        SBPL_FPRINTF(fDeb, "updating affected state %d:\n", stateID);
        PrintSearchState(searchstateinfo, fDeb);
        SBPL_FPRINTF(fDeb, "\n");
#endif

        //now we really do need to update it
        Recomputegval(searchstateinfo);
        UpdateSetMembership(searchstateinfo);
        numofstatesaffected++;

#if DEBUG
        SBPL_FPRINTF(fDeb, "the state %d after update\n", stateID);
        PrintSearchState(searchstateinfo, fDeb);
        SBPL_FPRINTF(fDeb, "\n");
#endif

      }
    }

    return;
  }

  vector <DState *> stv;
  set<DState *> visited;
  int numofstatesaffected = 0;

  for (int pind = 0; pind < (int)statesIDV->size(); pind++) {
    StateID stateID = statesIDV->at(pind);

    if (environment_->StateID2IndexMapping[stateID][DMDP_STATEID2IND] == -1)
      continue;

    //now get the state
    CMDPSTATE *state = GetState(stateID, pSearchStateSpace_);
    DState *searchstateinfo = (DState *)state->PlannerSpecificData;
    assert(searchstateinfo);

    //now check that the state is not start state and was created after last search reset
    if ((stateID != pSearchStateSpace_->searchstartstate->id) && (visited.find(searchstateinfo) == visited.end()) && (searchstateinfo->callnumberaccessed == pSearchStateSpace_->callnumber)) {
#if DEBUG
      SBPL_FPRINTF(fDeb, "updating affected state %d:\n", stateID);
      PrintSearchState(searchstateinfo, fDeb);
      SBPL_FPRINTF(fDeb, "\n");
#endif

      //now we really do need to update it
      stv.push_back(searchstateinfo);
      visited.insert(searchstateinfo);
      /*Recomputegval(searchstateinfo);
      if (_eps == 0.0)
      	UpdateSetMembership(searchstateinfo);*/
    }

    numofstatesaffected++;

#if DEBUG
    SBPL_FPRINTF(fDeb, "the state %d after update\n", stateID);
    PrintSearchState(searchstateinfo, fDeb);
    SBPL_FPRINTF(fDeb, "\n");
#endif
  }

  CList *pinconslist = pSearchStateSpace_->inconslist;

  //move incons into open
  while (pinconslist->firstelement != NULL) {
    DState *state = (DState *)pinconslist->firstelement->liststate;
    //remove from INCONS
    pinconslist->remove(state, D_INCONS_LIST_ID);
    //assert(state->path != NULL);
    state->g_p = INFINITECOST;

    if (visited.find(state) == visited.end()) {
      stv.push_back(state);
      visited.insert(state);
    }

  }

  //std::sort(stv.begin(), stv.end(), my_g);
  for (int pind = 0; pind < (int)stv.size(); pind++) {
    DState *s = stv[pind];
    /*int g_old = s->g;
    if (_state_print) {
    SBPL_FPRINTF(stdout, "updating affected state %d:\n", s->MDPstate->StateID);
    PrintSearchState(s, stdout);
    //Recomputegprime(pSearchStateSpace_, s);
    SBPL_FPRINTF(stdout, "\n");
    }*/
    Recomputegval(s);
    UpdateSetMembership(s);
    /*if (_state_print) {
    	Recomputegprime(pSearchStateSpace_, s);
    	SBPL_FPRINTF(stdout, "the state %d after update\n", s->MDPstate->StateID);
    	PrintSearchState(s, stdout);
    	SBPL_FPRINTF(stdout, "\n");
    }*/

  }

  stv.clear();

  //TODO - check. I believe that there are cases when number of states generated is drastically smaller than the number of states really affected, which is a bug!
  if (_state_print)
    SBPL_PRINTF("%d states really affected (%d states generated total so far)\n", numofstatesaffected, (int)environment_->StateID2IndexMapping.size());
}




//-----------------------------Interface function-----------------------------------------------------


//returns 1 if found a solution, and 0 otherwise
int DPlanner::replan(double allocated_time_secs, vector<StateID> *solution_stateIDs_V, int *psolcost) {
  vector<StateID> pathIds;
  int PathCost = 0;
  bool bFound = false;
  *psolcost = 0;
  bool bOptimalSolution = false;

  SBPL_PRINTF("planner: replan called (bFirstSol=%d, bOptSol=%d)\n", bsearchuntilfirstsolution, bOptimalSolution);

  //plan for the first solution only
  if ((bFound = Search(pSearchStateSpace_, pathIds, PathCost, bsearchuntilfirstsolution, bOptimalSolution, allocated_time_secs)) == false) {
    SBPL_PRINTF("failed to find a solution\n");
  }

  //copy the solution
  *solution_stateIDs_V = pathIds;
  *psolcost = PathCost;

  return (int)bFound;

}

int DPlanner::set_goal(StateID goal_stateID) {

  SBPL_PRINTF("planner: setting goal to %d\n", goal_stateID);
  environment_->PrintState(goal_stateID, true, stdout);

  //it will be a new search iteration
  //pSearchStateSpace_->searchiteration++;
  pSearchStateSpace_->bRebuildOpenList = true; //is not really necessary for search goal changes

  if (bforwardsearch) {
    if (SetSearchGoalState(goal_stateID, pSearchStateSpace_) != 1) {
      SBPL_ERROR("ERROR: failed to set search goal state\n");
      return 0;
    }
  } else {
    if (SetSearchStartState(goal_stateID, pSearchStateSpace_) != 1) {
      SBPL_ERROR("ERROR: failed to set search start state\n");
      return 0;
    }
  }

  return 1;
}


int DPlanner::set_start(StateID start_stateID) {

  SBPL_PRINTF("planner: setting start to %d\n", start_stateID);
  environment_->PrintState(start_stateID, true, stdout);

  //it will be a new search iteration
  //pSearchStateSpace_->searchiteration++;
  pSearchStateSpace_->bRebuildOpenList = true;


  if (bforwardsearch) {
    if (SetSearchStartState(start_stateID, pSearchStateSpace_) != 1) {
      SBPL_ERROR("ERROR: failed to set search start state\n");
      return 0;
    }
  } else {
    if (SetSearchGoalState(start_stateID, pSearchStateSpace_) != 1) {
      SBPL_ERROR("ERROR: failed to set search goal state\n");
      return 0;
    }
  }

  return 1;

}
bool DPlanner::IsExpanded(StateID stateID) {
  //first check that the state exists (to avoid creation of additional states)
  if (environment_->StateID2IndexMapping[stateID][DMDP_STATEID2IND] == -1)
    return false;

  //continue;
  //now get the state
  CMDPSTATE *state = GetState(stateID, pSearchStateSpace_);
  DState *searchstateinfo = (DState *)state->PlannerSpecificData;

  //now check that the state is not start state and was created after last search reset
  if (searchstateinfo->iterationclosed == pSearchStateSpace_->searchiteration) {
    //State is created and expanded
    return true;
  }

  return false;
}

void DPlanner::update_succs_of_changededges(vector<StateID> *succstatesIDV) {
  SBPL_PRINTF("UpdateSuccs called on %d succs\n", (unsigned int)succstatesIDV->size());

  Update_SearchSuccs_of_ChangedEdges(succstatesIDV);
}

void DPlanner::update_preds_of_changededges(vector<StateID> *predstatesIDV) {
  SBPL_PRINTF("UpdatePreds called on %d preds\n", (unsigned int)predstatesIDV->size());

  Update_SearchSuccs_of_ChangedEdges(predstatesIDV);
}


int DPlanner::force_planning_from_scratch() {
  SBPL_PRINTF("planner: forceplanfromscratch set\n");

  pSearchStateSpace_->bReinitializeSearchStateSpace = true;

  return 1;
}


int DPlanner::set_search_mode(bool bSearchUntilFirstSolution) {
  SBPL_PRINTF("planner: search mode set to %d\n", bSearchUntilFirstSolution);

  bsearchuntilfirstsolution = bSearchUntilFirstSolution;

  return 1;
}


void DPlanner::costs_changed(StateChangeQuery const &stateChange) {
  if (pSearchStateSpace_->bReinitializeSearchStateSpace == true || pSearchStateSpace_->searchiteration == 0)
    return; //no processing if no search efforts anyway

  if (bforwardsearch)
    Update_SearchSuccs_of_ChangedEdges(stateChange.getSuccessors());
  else
    Update_SearchSuccs_of_ChangedEdges(stateChange.getPredecessors());
}


//---------------------------------------------------------------------------------------------------------

