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
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
using namespace std;

#include <sbpl/planners/adaptiveASTAR.h>

#include <iostream>
#include <stdlib.h>     /* malloc, free, rand */
#include <math.h>
#include <time.h>

// SBPL includes
#include <sbpl/utils/heap.h>
#include <sbpl/utils/list.h>
#include <sbpl/utils/key.h>


// REMOVE: quick fix
#define USE_XIAOXUN_STUFF 1
#define RUNS 10
#define SPEED 5
#define TARGET_MOVE_STEPS 2





// Create - Init - Destroy
// =======================


AAPlanner::AAPlanner(DiscreteSpaceInformation *environment, bool bSearchForward) {
    SBPL_DEBUG("call AAstar constructor_________________________\n");
    printf("AAPlanner(%p, bool)\n", environment);

    bforwardsearch = bSearchForward;
    environment_ = environment;
    bsearchuntilfirstsolution = true;       // Xiaoxun: should I change this to "true" ???
//  finitial_eps = AA_DEFAULT_INITIAL_EPS;
//  finitial_eps = 1.0;   // xiaoxun modify this for AA*

    searchexpands = 0;
    MaxMemoryCounter = 0;
    fDeb = fopen("AAstar_debug.txt", "w");

    SBPL_DEBUG("debug on\n");

    pSearchStateSpace_ = new AASearchStateSpace_t;

    // Create and initialize the search space
    if (CreateSearchStateSpace(pSearchStateSpace_) != 1) {
        SBPL_ERROR("ERROR: failed to create the state space\n");
        return;
    }
    else
        SBPL_DEBUG("Search Space created\n");

    if (InitializeSearchStateSpace(pSearchStateSpace_) != 1) {
        SBPL_ERROR("ERROR: failed to initialise the state space\n");
        return;
    }
    else
        SBPL_DEBUG("Search Space initialised\n");
}

AAPlanner::~AAPlanner() {
    if (pSearchStateSpace_ != NULL) {
        //delete the statespace
        DeleteSearchStateSpace(pSearchStateSpace_);
        delete pSearchStateSpace_;
    }

    fclose(fDeb);
}


void AAPlanner::Initialize_searchinfo(CMDPSTATE *state, AASearchStateSpace_t *pSearchStateSpace) {

    AAState *searchstateinfo = (AAState *)state->PlannerSpecificData;

    searchstateinfo->MDPstate = state;
    InitializeSearchStateInfo(searchstateinfo, pSearchStateSpace);
}



//Xiaoxun 1: CreateState()
CMDPSTATE *AAPlanner::CreateState(int stateID, AASearchStateSpace_t *pSearchStateSpace) {
    CMDPSTATE *state = NULL;

#if DEBUG

    if (environment_->StateID2IndexMapping[stateID][ARAMDP_STATEID2IND] != -1) {
        SBPL_ERROR("ERROR in CreateState: state already created\n");
        exit(1);
    }

#endif
    //adds to the tail a state
    state = pSearchStateSpace->searchMDP.AddState(stateID); // Xiaoxun 2:  searchMDP == vector<CMDPSTATE*> StateArray;

    //remember the index of the state
    environment_->StateID2IndexMapping[stateID][ARAMDP_STATEID2IND] = pSearchStateSpace->searchMDP.StateArray.size() - 1;

#if DEBUG

    if (state != pSearchStateSpace->searchMDP.StateArray[environment_->StateID2IndexMapping[stateID][ARAMDP_STATEID2IND]]) {
        SBPL_ERROR("ERROR in CreateState: invalid state index\n");
        exit(1);
    }

#endif

    //create search specific info
    state->PlannerSpecificData = (AAState *)malloc(sizeof(AAState));
    Initialize_searchinfo(state, pSearchStateSpace); // Xiaoxun 3:
    MaxMemoryCounter += sizeof(AAState);



    return state;
}

CMDPSTATE *AAPlanner::GetState(int stateID, AASearchStateSpace_t *pSearchStateSpace) {

    if (stateID >= (int)environment_->StateID2IndexMapping.size()) {
        SBPL_ERROR("ERROR int GetState: stateID %d is invalid\n", stateID);
        exit(1);
    }

    if (environment_->StateID2IndexMapping[stateID][ARAMDP_STATEID2IND] == -1) // Xiaoxun: ARAMDP_STATEID2IND == 0
        return CreateState(stateID, pSearchStateSpace);
    else
        return pSearchStateSpace->searchMDP.StateArray[environment_->StateID2IndexMapping[stateID][ARAMDP_STATEID2IND]];

}








// Search State handling
// =====================

int AAPlanner::ComputeHeuristic(CMDPSTATE *MDPstate, AASearchStateSpace_t *pSearchStateSpace) {

    //compute heuristic for search
    if (bforwardsearch) {
        //forward search: heur = distance from state to searchgoal which is Goal AAState
        int retv =  environment_->GetFromToHeuristic(MDPstate->StateID, pSearchStateSpace->searchgoalstate->StateID);
        return retv;
    } else {
        //backward search: heur = distance from state to searchstart
        int retv = environment_->GetFromToHeuristic(MDPstate->StateID, pSearchStateSpace->searchstartstate->StateID);
        return retv;
    }
}

//initialization of a state
void AAPlanner::InitializeSearchStateInfo(AAState *state, AASearchStateSpace_t *pSearchStateSpace) {
    state->g = INFINITECOST;
    state->iterationclosed = 0;
    state->bestnextstate = NULL;
    state->bestpredstate = NULL;
    state->costtobestnextstate = INFINITECOST;
    state->heapindex = 0;
/// state->numofexpands = 0;

//Xiaoxun added this:
    state->callnumberaccessed = pSearchStateSpace->callnumber;
    state->generated_iteration = 0;



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

//Xiaoxun 0: re-initialization of a state
void AAPlanner::ReInitializeSearchStateInfo(AAState *state, AASearchStateSpace_t *pSearchStateSpace) {
    state->g = INFINITECOST;
    state->iterationclosed = 0;   // xiaoxun 11111111111111111   ????
    state->heapindex = 0;         // xiaoxun 222222222222222222   ????
    state->costtobestnextstate = INFINITECOST;
    state->bestnextstate = NULL;
    state->bestpredstate = NULL;

    state->callnumberaccessed = pSearchStateSpace->callnumber;
    state->generated_iteration = pSearchStateSpace->searchiteration;

    //compute heuristics
#if USE_HEUR
    state->h = ComputeHeuristic(state->MDPstate, pSearchStateSpace);
#else
    state->h = 0;
#endif




    return;
}

void AAPlanner::DeleteSearchStateData(AAState *state) {
    //no memory was allocated
    MaxMemoryCounter = 0;
    return;
}







//TODO-debugmax - add obsthresh and other thresholds to other environments in 3dkin
int AAPlanner::GetGVal(int StateID, AASearchStateSpace_t *pSearchStateSpace) {
    CMDPSTATE *cmdp_state = GetState(StateID, pSearchStateSpace);
    AAState *state = (AAState *)cmdp_state->PlannerSpecificData;
    return state->g;
}


//Xiaoxun 0: function Improve Path()  ==  computeshortestpath() in grids.
//           returns 1 if the solution is found, 0 if the solution does not exist and 2 if it ran out of time
int AAPlanner::ImprovePath(AASearchStateSpace_t *pSearchStateSpace, double MaxNumofSecs) {
    int expands;
    AAState *state, *searchgoalstate, *searchstartstate;
    CKey key, minkey;
    CKey goalkey, startkey;

    expands = 0;

#ifdef ADAPTIVE_H
    pSearchStateSpace->keymod[pSearchStateSpace->searchiteration] = pSearchStateSpace->keymodifer;
#endif

    if (pSearchStateSpace->searchgoalstate == NULL) {
        SBPL_ERROR("ERROR searching: no goal state is set\n");
        exit(1);
    }

    if (pSearchStateSpace->searchstartstate == NULL) {
        SBPL_ERROR("ERROR searching: no start state is set\n");
        exit(1);
    }

#ifdef XIAOXUN_DEBUG
    SBPL_DEBUG("$$$$ call Improve Path , searchiteration == [%d], call = [%d]\n", pSearchStateSpace->searchiteration, pSearchStateSpace->callnumber);
    SBPL_FPRINTF(fDeb, "$$$$ call Improve Path , searchiteration == [%d], call = [%d]\n", pSearchStateSpace->searchiteration, pSearchStateSpace->callnumber);
#endif







#ifdef REVERSE  // search backward ____________(1)_____________

//Xiaoxun 1: set search goal state (=agent & start state ) and get the   (goalkey)
    searchstartstate = (AAState *)(pSearchStateSpace->searchstartstate->PlannerSpecificData);

    if (searchstartstate->generated_iteration != pSearchStateSpace->searchiteration)
        ReInitializeSearchStateInfo(searchstartstate, pSearchStateSpace);


    //set search goal key (= key of the agent)
    goalkey.key[0] = searchstartstate->g;
#ifdef TIE_LARGE_G
    goalkey.key[1] = searchstartstate->h;
#endif

    //expand states until done
    minkey = pSearchStateSpace->heap->getminkeyheap();    //Xiaoxun2: return the key of the Topheap()
    CKey oldkey = minkey;


    while (!pSearchStateSpace->heap->emptyheap() && minkey.key[0] < INFINITECOST && goalkey > minkey) {
        //Xiaoxun 2: pop out and get the state with the smallest f-value
        state = (AAState *)pSearchStateSpace->heap->deleteminheap();
        state->iterationclosed = pSearchStateSpace->searchiteration;   // expanded_iteration is set here
        expands++;


// Xiaoxun 3: expanding a state, here
///////////////////////////////////////////////////////////////

        UpdatePreds(state, pSearchStateSpace);    //search backwards

///////////////////////////////////////////////////////////////
        //recompute minkey
        minkey = pSearchStateSpace->heap->getminkeyheap();

        //recompute goalkey if necessary
        if (goalkey.key[0] != (int)searchstartstate->g) {
            //SBPL_DEBUG("re-computing goal key\n");      //recompute the goal key (heuristics should be zero)
            goalkey.key[0] = searchstartstate->g;
#ifdef TIE_LARGE_G
            goalkey.key[1] = searchstartstate->h;
#endif
        }
    } // end while


    int retv = 1;

    if (searchstartstate->g == INFINITECOST && pSearchStateSpace->heap->emptyheap()) {
//      SBPL_INFO("solution does not exist: search exited because heap is empty\n");
        retv = 0;
    } else {
        retv = 1;
#ifdef XIAOXUN_STATISTICS
        pSearchStateSpace->expansion_of_the_search[pSearchStateSpace->totalsearchiteration] = expands;
#endif

    }

    searchexpands += expands;





// ...
#else // forward search __________________(2)__________________________________________________________________

    //Xiaoxun 1: set goal state and get the   (goalkey)
    searchgoalstate = (AAState *)(pSearchStateSpace->searchgoalstate->PlannerSpecificData);

    if (searchgoalstate->generated_iteration != pSearchStateSpace->searchiteration)
        ReInitializeSearchStateInfo(searchgoalstate, pSearchStateSpace);

    //set goal key
    goalkey.key[0] = searchgoalstate->g;
#ifdef TIE_LARGE_G
    goalkey.key[1] = searchgoalstate->h;
#endif

    //expand states until done
    minkey = pSearchStateSpace->heap->getminkeyheap();    //Xiaoxun2: return the key of the Topheap()
    CKey oldkey = minkey;


    while (!pSearchStateSpace->heap->emptyheap() && minkey.key[0] < INFINITECOST && goalkey > minkey) {
        //Xiaoxun 2: pop out and get the state with the smallest f-value
        state = (AAState *)pSearchStateSpace->heap->deleteminheap();

#ifdef XIAOXUN_DEBUG

        if (pSearchStateSpace->callnumber == 1)
            SBPL_FPRINTF(fDeb, "expanding state(%4d): h=%d g=%u key=%u iterclosed=%d callnuma=%d expands=%d (g(goal)=%u)\n",
                    state->MDPstate->StateID, state->h, state->g, state->g + (int)(pSearchStateSpace->eps * state->h),
                    state->iterationclosed, state->callnumberaccessed, state->numofexpands, searchgoalstate->g);

//        SBPL_FPRINTF(fDeb, "expanding: ");
//        PrintSearchState(state, fDeb);
        fflush(fDeb);
#endif

        state->iterationclosed = pSearchStateSpace->searchiteration;   // expanded_iteration is set here
        expands++;


// Xiaoxun 3: expanding a state, here
///////////////////////////////////////////////////////////////
        UpdateSuccs(state, pSearchStateSpace);    //search forwards
///////////////////////////////////////////////////////////////

        //recompute minkey
        minkey = pSearchStateSpace->heap->getminkeyheap();

        //recompute goalkey if necessary
        if (goalkey.key[0] != (int)searchgoalstate->g) {
            //SBPL_DEBUG("re-computing goal key\n");      //recompute the goal key (heuristics should be zero)
            goalkey.key[0] = searchgoalstate->g;
#ifdef TIE_LARGE_G
            goalkey.key[1] = searchgoalstate->h;
#endif
        }

#ifdef XIAOXUN_DEBUG

        if (expands % 100000 == 0 && expands > 0)
            SBPL_DEBUG("expands so far=%u\n", expands);

#endif

    } // end while


#ifdef XIAOXUN_DEBUG
    SBPL_DEBUG("\n when [%d]-th search is done, goaloldkey == goal->g = %d, expands this search = %d\n", pSearchStateSpace->searchiteration, searchgoalstate->g, expands);
    SBPL_FPRINTF(fDeb, "\n when [%d]-th search is done, goaloldkey == goal->g = %d, expands this search = %d\n", pSearchStateSpace->searchiteration, searchgoalstate->g, expands);
#endif

    int retv = 1;

    if (searchgoalstate->g == INFINITECOST && pSearchStateSpace->heap->emptyheap()) {
//      SBPL_INFO("solution does not exist: search exited because heap is empty\n");
        retv = 0;
    } else if (!pSearchStateSpace->heap->emptyheap() && goalkey > minkey) {
//      SBPL_INFO("search exited because it ran out of time\n");
        retv = 2;
    } else if (searchgoalstate->g == INFINITECOST && !pSearchStateSpace->heap->emptyheap()) {
//      SBPL_INFO("solution does not exist: search exited because all candidates for expansion have infinite heuristics\n");
        retv = 0;
    } else {
//      SBPL_INFO("search exited with a solution for eps=%.3f\n", pSearchStateSpace->eps);
        retv = 1;
#ifdef XIAOXUN_STATISTICS
        pSearchStateSpace->expansion_of_the_search[pSearchStateSpace->totalsearchiteration] = expands;
#endif

    }

    //SBPL_FPRINTF(fDeb, "expanded=%d\n", expands);
    searchexpands += expands;

#endif


    return retv;
}








// Search Space handling
// =====================


//creates (allocates memory) search state space
//does not initialize search statespace
int AAPlanner::CreateSearchStateSpace(AASearchStateSpace_t *pSearchStateSpace) {

    // Create heap
    pSearchStateSpace->heap = new CHeap;
    pSearchStateSpace->inconslist = new CList;

    // Account memory
    MaxMemoryCounter += sizeof(CHeap);
    MaxMemoryCounter += sizeof(CList);

    // Reset start and goal
    pSearchStateSpace->searchgoalstate = NULL;
    pSearchStateSpace->searchstartstate = NULL;

    // Reset statistics
    searchexpands = 0;


    pSearchStateSpace->bReinitializeSearchStateSpace = false;

    return 1; // OK (not-so-wise 'success' code
}

//deallocates memory used by SearchStateSpace
void AAPlanner::DeleteSearchStateSpace(AASearchStateSpace_t *pSearchStateSpace) {
    if (pSearchStateSpace->heap != NULL) {
        pSearchStateSpace->heap->makeemptyheap();
        delete pSearchStateSpace->heap;
        pSearchStateSpace->heap = NULL;
    }

    if (pSearchStateSpace->inconslist != NULL) {
        pSearchStateSpace->inconslist->makeemptylist(AA_INCONS_LIST_ID);
        delete pSearchStateSpace->inconslist;
        pSearchStateSpace->inconslist = NULL;
    }




    //delete the states themselves
    int iend = (int)pSearchStateSpace->searchMDP.StateArray.size();

    for (int i = 0; i < iend; i++) {
        CMDPSTATE *state = pSearchStateSpace->searchMDP.StateArray[i];

        if (state != NULL && state->PlannerSpecificData != NULL) {
            DeleteSearchStateData((AAState *)state->PlannerSpecificData);
            free((AAState *)state->PlannerSpecificData);
            state->PlannerSpecificData = NULL;
        }
    }

    pSearchStateSpace->searchMDP.Delete();
}



//reset properly search state space
//needs to be done before deleting states
int AAPlanner::ResetSearchStateSpace(AASearchStateSpace_t *pSearchStateSpace) {
    pSearchStateSpace->heap->makeemptyheap();
    pSearchStateSpace->inconslist->makeemptylist(AA_INCONS_LIST_ID);

    return 1;
}

//initialization before each search
// Xiaoxun: note that: before EACH search
void AAPlanner::ReInitializeSearchStateSpace(AASearchStateSpace_t *pSearchStateSpace) {
    CKey key;

    //increase callnumber
    pSearchStateSpace->callnumber++;         // call number == the order of this test case
    pSearchStateSpace->searchiteration = 0;
#ifdef ADAPTIVE_H
    pSearchStateSpace->keymodifer = 0;
#endif

    pSearchStateSpace->bReinitializeSearchStateSpace = false;   // after R-initialize the search space, set it to (false)
    pSearchStateSpace->bReevaluatefvals = false;
}


//very first initialization
int AAPlanner::InitializeSearchStateSpace(AASearchStateSpace_t *pSearchStateSpace) {

//  if(pSearchStateSpace->heap->currentsize != 0 || pSearchStateSpace->inconslist->currentsize != 0)
    if (pSearchStateSpace->heap->currentsize != 0) {
        SBPL_ERROR("ERROR in InitializeSearchStateSpace: heap or list is not empty\n");
        exit(1);
    }


    pSearchStateSpace->searchiteration = 0;
    pSearchStateSpace->callnumber = 0;

    pSearchStateSpace->robot_steps = 0;
    pSearchStateSpace->target_steps = 0;
    pSearchStateSpace->total_robot_steps = 0;
    pSearchStateSpace->totalmovecost = 0;


// 2010.02.28
    pSearchStateSpace->hunter_movecost_testcase = 0;  // (1)
    pSearchStateSpace->target_movecost_testcase = 0;


    //create and set the search start state
    pSearchStateSpace->searchgoalstate = NULL;
    pSearchStateSpace->searchstartstate = NULL;
#ifdef ADAPTIVE_H
    pSearchStateSpace->keymodifer = 0;
#endif


    pSearchStateSpace->bReinitializeSearchStateSpace = true;   // (1) this is the only place to set it to be "true"


    return 1;
}


int AAPlanner::SetSearchGoalState(int SearchGoalStateID, AASearchStateSpace_t *pSearchStateSpace) {
//Xiaoxun1: if(1) ps->searchgoalstate has not been set,
//          if(2) goal has moved, so       ps->searchgoalstate->StateID != SearchGoalStateID
//
    if (pSearchStateSpace->searchgoalstate == NULL || pSearchStateSpace->searchgoalstate->StateID != SearchGoalStateID) {
        pSearchStateSpace->searchgoalstate = GetState(SearchGoalStateID, pSearchStateSpace);
        //recompute heuristic for the heap if heuristics is used
// Xiaoxun 1:  becasue the search goal state has changed, so the h-values need to be re-calculated
//             but why recalculated all states (generated + expanded)??????????

#if USE_HEUR

        for (int i = 0; i < (int)pSearchStateSpace->searchMDP.StateArray.size(); i++) {
            CMDPSTATE *MDPstate = pSearchStateSpace->searchMDP.StateArray[i];
            AAState *state = (AAState *)MDPstate->PlannerSpecificData;
            state->h = ComputeHeuristic(MDPstate, pSearchStateSpace);
        }

        pSearchStateSpace->bReevaluatefvals = true;
#endif
    }


    return 1;
}


int AAPlanner::SetSearchStartState(int SearchStartStateID, AASearchStateSpace_t *pSearchStateSpace) {
//Xiaoxun 1: this is to get the address of the (CMDPSTATE*) of the start state
    CMDPSTATE *MDPstate = GetState(SearchStartStateID, pSearchStateSpace); //Xiaoxun 1:  see P13

//Xiaoxun 2: if the start state changed, need to re-initialize it
    if (MDPstate !=  pSearchStateSpace->searchstartstate) {
        pSearchStateSpace->searchstartstate = MDPstate;
        pSearchStateSpace->bReinitializeSearchStateSpace = true;   // (2)
    }

    return 1;

}









// Solution handling
// =================

//Xiaoxun 0: to search the pointer of -->trace
int AAPlanner::ReconstructPath(AASearchStateSpace_t *pSearchStateSpace) {
    CMDPSTATE *MDPstate;
    CMDPSTATE *PredMDPstate;
    AAState *predstateinfo, *stateinfo;

    //nothing to do, if search is backward
    if (bforwardsearch) { // forward search for AA*
        //Xiaoxun added this if()
        if (pSearchStateSpace->temp_goal_state && pSearchStateSpace->searchgoalstate != pSearchStateSpace->temp_goal_state)
            MDPstate = pSearchStateSpace->temp_goal_state;
        else
            MDPstate = pSearchStateSpace->searchgoalstate;


#ifdef XIAOXUN_DEBUG
        SBPL_FPRINTF(fDeb, "reconstructing a path:\n");
#endif

        while (MDPstate != pSearchStateSpace->searchstartstate) {
            stateinfo = (AAState *)MDPstate->PlannerSpecificData;

#if DEBUG
            PrintSearchState(stateinfo, fDeb);
#endif

            if (stateinfo->g == INFINITECOST) {
                //SBPL_ERROR("ERROR in ReconstructPath: g of the state on the path is INFINITE\n");
                //exit(1);
                return -1;
            }

            if (stateinfo->bestpredstate == NULL) {
                SBPL_ERROR("ERROR in ReconstructPath: bestpred is NULL\n");
                exit(1);
            }

            //get the parent state
            PredMDPstate = stateinfo->bestpredstate;
            predstateinfo = (AAState *)PredMDPstate->PlannerSpecificData;

            //set its best next info
            predstateinfo->bestnextstate = MDPstate;

            //check the decrease of g-values along the path
            /*          if(predstateinfo->v >= stateinfo->g)
                        {
                            SBPL_ERROR("ERROR in ReconstructPath: g-values are non-decreasing\n");
                            PrintSearchState(predstateinfo, fDeb);
                            exit(1);
                        }
            */
            //transition back
            MDPstate = PredMDPstate;
        }
    } else { // backward search _____________________________________back____________________________________

        MDPstate = pSearchStateSpace->searchstartstate;

        while (MDPstate != pSearchStateSpace->searchgoalstate) {
            stateinfo = (AAState *)MDPstate->PlannerSpecificData;

            if (stateinfo->g == INFINITECOST) {
                SBPL_ERROR("ERROR in ReconstructPath: g of the state on the path is INFINITE\n");
                exit(1);
            }

            if (stateinfo->bestnextstate == NULL) {
                SBPL_ERROR("ERROR in ReconstructPath: bestpred is NULL\n");
                exit(1);
            }

            //get the parent state
            PredMDPstate = stateinfo->bestnextstate;
            predstateinfo = (AAState *)PredMDPstate->PlannerSpecificData;


            //check the decrease of g-values along the path
            if (predstateinfo->g >= stateinfo->g) {
                SBPL_ERROR("ERROR in ReconstructPath: g-values are non-decreasing\n");
                PrintSearchState(predstateinfo, fDeb);
                exit(1);
            }

            //transition back
            MDPstate = PredMDPstate;
        }


    }







#ifdef XIAOXUN_DEBUG
    SBPL_DEBUG("After reconstructing a path:\n");
    SBPL_FPRINTF(fDeb, "After reconstructing a path:\n");
#endif

    return 1;
}



void AAPlanner::PrintSearchPath(AASearchStateSpace_t *pSearchStateSpace, FILE *fOut) {
    AAState *searchstateinfo;
    CMDPSTATE *state;
    int goalID;
    int PathCost;

    if (bforwardsearch) {
        state  = pSearchStateSpace->searchstartstate;

        if (pSearchStateSpace->temp_goal_state && pSearchStateSpace->temp_goal_state != pSearchStateSpace->searchgoalstate)
            goalID = pSearchStateSpace->temp_goal_state->StateID;
        else
            goalID = pSearchStateSpace->searchgoalstate->StateID;

        PathCost = ((AAState *)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g;

    } else { // search backward _________________________________________
//      state = pSearchStateSpace->searchgoalstate;
        state = pSearchStateSpace->searchstartstate;
        goalID = pSearchStateSpace->searchstartstate->StateID;   // goalID = agent->StateID

        PathCost = ((AAState *)pSearchStateSpace->searchstartstate->PlannerSpecificData)->g;
    }

    if (fOut == NULL)
        fOut = stdout;

//  PathCost = ((AAState*)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g;


    SBPL_FPRINTF(fOut, "Printing a path from state %d to the goal state %d\n",
            state->StateID, pSearchStateSpace->searchgoalstate->StateID);
    SBPL_FPRINTF(fOut, "Path cost = %d:\n", PathCost);


    environment_->PrintState(state->StateID, false, fOut);

    int costFromStart = 0;

    while (state->StateID != goalID) {
        SBPL_FPRINTF(fOut, "state %d ", state->StateID);

        if (state->PlannerSpecificData == NULL) {
            SBPL_FPRINTF(fOut, "path does not exist since search data does not exist\n");
            break;
        }

        searchstateinfo = (AAState *)state->PlannerSpecificData;

        if (searchstateinfo->bestnextstate == NULL) {
            SBPL_FPRINTF(fOut, "path does not exist since bestnextstate == NULL  ????????????\n");
            break;
        }

        if (searchstateinfo->g == INFINITECOST) {
            SBPL_FPRINTF(fOut, "path does not exist since bestnextstate == NULL\n");
            break;
        }

        int costToGoal = PathCost - costFromStart;
        int transcost = searchstateinfo->g - ((AAState *)(searchstateinfo->bestnextstate->PlannerSpecificData))->g; //->v

        if (bforwardsearch)
            transcost = -transcost;

        costFromStart += transcost;

        SBPL_FPRINTF(fOut, "g=%d-->state %d, h = %d ctg = %d  ", searchstateinfo->g,
                searchstateinfo->bestnextstate->StateID, searchstateinfo->h, costToGoal);

        state = searchstateinfo->bestnextstate;

        environment_->PrintState(state->StateID, false, fOut);



    }
}

void AAPlanner::PrintSearchState(AAState *state, FILE *fOut) {
    int expansions = -1;
#if STATISTICS
    expansions = state->expansions;
#endif
    SBPL_FPRINTF(fOut, "state %d: h=%d g=%u iterc=%d callnuma=%d expands=%d heapind=%d \n",
            state->MDPstate->StateID, state->h, state->g,
            state->iterationclosed, state->callnumberaccessed, expansions, state->heapindex);
    environment_->PrintState(state->MDPstate->StateID, true, fOut);

}



int AAPlanner::getHeurValue(AASearchStateSpace_t *pSearchStateSpace, int StateID) {
    CMDPSTATE *MDPstate = GetState(StateID, pSearchStateSpace);
    AAState *searchstateinfo = (AAState *)MDPstate->PlannerSpecificData;
    return searchstateinfo->h;
}

//Xiaoxun 0: this is to get the <STATE_ID> of all cells along the path
vector<int> AAPlanner::GetSearchPath(AASearchStateSpace_t *pSearchStateSpace, int &solcost) {
    vector<int> SuccIDV;
    vector<int> CostV;
    vector<int> wholePathIds;
    AAState *searchstateinfo;
    CMDPSTATE *state = NULL;
    CMDPSTATE *goalstate = NULL;
    CMDPSTATE *startstate = NULL;

    if (bforwardsearch) { // search forward
        startstate = pSearchStateSpace->searchstartstate;
        goalstate = pSearchStateSpace->searchgoalstate;

        //reconstruct the path by setting -->bestnextstate pointers appropriately
        ReconstructPath(pSearchStateSpace);   //setting the (my)-->trace pointers
    } else {
//      startstate = pSearchStateSpace->searchgoalstate;
//      goalstate = pSearchStateSpace->searchstartstate;
        startstate = pSearchStateSpace->searchstartstate;
        goalstate = pSearchStateSpace->searchgoalstate;
    }


    state = startstate;

    wholePathIds.push_back(state->StateID);
    solcost = 0;

    FILE *fOut = stdout;

    while (state->StateID != goalstate->StateID) {
        if (state->PlannerSpecificData == NULL) {
            SBPL_FPRINTF(fOut, "path does not exist since search data does not exist\n");
            break;
        }

        searchstateinfo = (AAState *)state->PlannerSpecificData;

        if (searchstateinfo->bestnextstate == NULL) {
            SBPL_FPRINTF(fOut, "path does not exist since bestnextstate == NULL 2222222 HERE \n");
            break;
        }

        if (searchstateinfo->g == INFINITECOST) {
            SBPL_FPRINTF(fOut, "path does not exist since g == INFINITYCOST  22222 \n");
            break;
        }

//???????????????????????????????????????????????????????????????????????????????????????????????
#ifdef REVERSE
// temp empty


#else
        environment_->GetSuccs(state->StateID, &SuccIDV, &CostV);
        int actioncost = INFINITECOST;

        for (int i = 0; i < (int)SuccIDV.size(); i++) {
            if (SuccIDV.at(i) == searchstateinfo->bestnextstate->StateID) {
                actioncost = CostV.at(i);
                break; // Xiaoxun added this break;
            }
        }

        if (actioncost == INFINITECOST)
            SBPL_WARN("WARNING: actioncost = %d\n", actioncost);

        solcost += actioncost;
#endif

//SBPL_FPRINTF(fDeb, "actioncost=%d between states %d and %d\n",
//        actioncost, state->StateID, searchstateinfo->bestnextstate->StateID);
//environment_->PrintState(state->StateID, false, fDeb);
//environment_->PrintState(searchstateinfo->bestnextstate->StateID, false, fDeb);


#if DEBUG
        AAState *nextstateinfo = (AAState *)(searchstateinfo->bestnextstate->PlannerSpecificData);

        if (actioncost != abs((int)(searchstateinfo->g - nextstateinfo->g)) && pSearchStateSpace->eps_satisfied <= 1.001) {
            SBPL_FPRINTF(fDeb, "ERROR: actioncost=%d is not matching the difference in g-values of %d\n",
                    actioncost, abs((int)(searchstateinfo->g - nextstateinfo->g)));
            SBPL_ERROR("ERROR: actioncost=%d is not matching the difference in g-values of %d\n",
                   actioncost, abs((int)(searchstateinfo->g - nextstateinfo->g)));
            PrintSearchState(searchstateinfo, fDeb);
            PrintSearchState(nextstateinfo, fDeb);
        }

#endif

        state = searchstateinfo->bestnextstate;
        wholePathIds.push_back(state->StateID);

    }// end while


    return wholePathIds;
}







// Search
// ======


//Xiaoxun 0:     Search() == ComputeShortestPath()...... however, it contains function IMPROVE_PATH() inside.
bool AAPlanner::Search(AASearchStateSpace_t *pSearchStateSpace, vector<int> &pathIds, int &PathCost, bool bFirstSolution, bool bOptimalSolution, double MaxNumofSecs) {
    CKey key;
    TimeStarted = clock();
    searchexpands = 0;



    if (pSearchStateSpace->bReinitializeSearchStateSpace == true) {
//#ifdef XIAOXUN_DEBUG
        SBPL_DEBUG("call AAPlanner::search() function ____+++______ReinitializeSearchSpace \n");
        SBPL_FPRINTF(fDeb, "call AAPlanner::search() function ____+++______ReinitializeSearchSpace \n");
//#endif
        //re-initialize state space
        ReInitializeSearchStateSpace(pSearchStateSpace);
        ReInitializeNewSearch(pSearchStateSpace);      // insert the search start state into OPEN,  start->g= 0 etc...
    }





    MaxNumofSecs = INFINITECOST;
    int prevexpands = 0;
    clock_t loop_time;

    loop_time = clock();



//Xiaoxun added this for AA*  backward search:
#ifdef REVERSE
#ifdef ADAPTIVE_H
    pSearchStateSpace->previous_search_start_state = pSearchStateSpace->searchstartstate;
#endif
#endif



//call  improve_path()
    if (ImprovePath(pSearchStateSpace, MaxNumofSecs) == 1) {
#ifdef XIAOXUN_DEBUG
        SBPL_DEBUG("AA* improve_path is done, and path found.\n");
#endif
    } else { // Xiaoxun added
        SBPL_INFO("no path exists, exit. \n");
        exit(1);
    }


#ifdef XIAOXUN_DEBUG
    //print the solution cost and eps bound
    SBPL_DEBUG("***********************************************\n expands=%d g(searchgoal)=%d time=%.3f \n",
           searchexpands - prevexpands,
           ((AAState *)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g, double(clock() - loop_time) / CLOCKS_PER_SEC);
    SBPL_FPRINTF(fDeb, "***********************************************\n expands=%d g(searchgoal)=%d time=%.3f \n",
            searchexpands - prevexpands,
            ((AAState *)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g, double(clock() - loop_time) / CLOCKS_PER_SEC);
#endif

    prevexpands = searchexpands;
//  }// end while



#if DEBUG
    fflush(fDeb);
#endif



#ifdef REVERSE  // search backwards
    PathCost = ((AAState *)pSearchStateSpace->searchstartstate->PlannerSpecificData)->g;
#else           // search forwards
    PathCost = ((AAState *)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g;
#endif

    MaxMemoryCounter += environment_->StateID2IndexMapping.size() * sizeof(int);


#ifdef ADAPTIVE_H
    pSearchStateSpace->pathlength[pSearchStateSpace->searchiteration] = PathCost;  // for adaptive a*__________
#endif


    int solcost = INFINITECOST;
    bool ret = false;

    if (PathCost == INFINITECOST) {
        SBPL_INFO("AA* could not find a solution\n");
        ret = false;
    } else {
#ifdef XIAOXUN_DEBUG
        SBPL_DEBUG("[%3d]-th search is finished with a solution \n", pSearchStateSpace->searchiteration);
        SBPL_FPRINTF(fDeb, "[%3d]-th search is finished with a solution \n", pSearchStateSpace->searchiteration);
//      SBPL_INFO("solution is found\n");
#endif
        pathIds = GetSearchPath(pSearchStateSpace, solcost);   // get the STATE_ID of the cells along the path
        ret = true;


#ifdef XIAOXUN_DEBUG
        SBPL_DEBUG("[%3d]-th search, after tracing back a path \n", pSearchStateSpace->searchiteration);
        SBPL_FPRINTF(fDeb, "[%3d]-th search, after tracing back a path \n", pSearchStateSpace->searchiteration);
//      SBPL_INFO("solution is found\n");
#endif

    }


#ifdef XIAOXUN_DISPLAY
    SBPL_INFO("total state-expansion this search = %d, time = %.3f secs, solution cost=%d\n",
           searchexpands, (clock() - TimeStarted) / ((double)CLOCKS_PER_SEC), solcost);
#endif
    //SBPL_FPRINTF(fStat, "%d %d\n", searchexpands, solcost);

    return ret;
}








// Interface Functions
// ===================

//returns 1 if found a solution, and 0 otherwise
int AAPlanner::replan(double allocated_time_secs, vector<int> *solution_stateIDs_V) {
    int solcost;

    return replan(allocated_time_secs, solution_stateIDs_V, &solcost);
}





//Xiaoxun 0: return 1 if found a solution, and 0 otherwise
//           re-plan is to do multiple searches until the goal is caught
int AAPlanner::replan(double allocated_time_secs, vector<int> *solution_stateIDs_V, int *psolcost) {
    vector<int> pathIds;
    bool bFound = false;
    int PathCost;
    bool bFirstSolution = this->bsearchuntilfirstsolution;

    bool bOptimalSolution = true;   // Xiaoxun: for AA*, we are aiming at opt solution (2009.08.21)
    *psolcost = 0;
    bool robot_need_to_replan = true; // Xiaoxun added here
    double totaltime, time_per_search;


    totaltime = 0;
    srand(6);


    //Xiaoxun 1: ComputeShortestPath() is caled here......    //plan a path
//SBPL_DEBUG("_________AA planner: calling function replan() (bFirstSol=%d, bOptSol=%d)\n", bFirstSolution, bOptimalSolution);
    SBPL_DEBUG("_________AA planner: calling function replan() (bFirstSol=%d, bOptSol=%d)\n", bFirstSolution);









    for (int i = 0; i < RUNS; i++) {
#ifdef ADAPTIVE_H
        SBPL_INFO("Adaptive A* TEST CASE = [%d],  callnumber = [%d] \n", i, pSearchStateSpace_->callnumber);
#else
        SBPL_INFO("A* TEST CASE = [%d],  callnumber = [%d] \n", i, pSearchStateSpace_->callnumber);
#endif

        robot_need_to_replan = true;
        pSearchStateSpace_->temp_goal_state = pSearchStateSpace_->searchgoalstate;

        ////////////////////////////////////////////////////////////////////////////////////////
        if (pSearchStateSpace_->callnumber > 0) {
#ifdef XIAOXUN_DEBUG
            SBPL_DEBUG("\n\n +++++++++++++++++++++++Reset for the (%d)-th TEST CASE   ++++++++++++++++++\n\n", pSearchStateSpace_->callnumber);
            SBPL_FPRINTF(fDeb, "\n\n +++++++++++++++++++++++Reset for the (%d)-th TEST CASE   ++++++++++++++++++\n\n", pSearchStateSpace_->callnumber);
#endif
            Reset_New_Test_Case(pSearchStateSpace_);
        }

        ////////////////////////////////////////////////////////////////////////////////////////



        Case_TimeStarted = clock(); // the start time of the run (= n+searches until the target is caught)


        while (pSearchStateSpace_->searchgoalstate != pSearchStateSpace_->searchstartstate) {
            if (robot_need_to_replan == true) {

                SBPL_DEBUG("%d-th_s    ", pSearchStateSpace_->searchiteration);

                if ((bFound = Search(pSearchStateSpace_, pathIds, PathCost, bFirstSolution, bOptimalSolution, allocated_time_secs)) == false) {
                    SBPL_ERROR("AA* failed to find a solution after calling Search()... \n");
                    return (int)bFound;
                }

                robot_need_to_replan = false; // after compute path, set the flag to false.
                *solution_stateIDs_V = pathIds; //copy the solution
                *psolcost = PathCost;

                SBPL_DEBUG("  - Pathcost  == %d ", PathCost);



#ifdef XIAOXUN_DEBUG
                SBPL_DEBUG("\n\n ++++++++ 1 +++++++++++++++ The (%d)-th search is done  ++++++++++++++++++\n", pSearchStateSpace_->searchiteration);
                SBPL_FPRINTF(fDeb, "\n\n ++++++++ 1 +++++++++++++++ The (%d)-th search is done  ++++++++++++++++++\n\n", pSearchStateSpace_->searchiteration);
#endif
            }





            ////////////////////////////////////////////////
            //Xiaoxun 1:    Robot moves 1 step here       //
            ////////////////////////////////////////////////
            if (robot_need_to_replan == false) {
#ifdef XIAOXUN_DEBUG
                SBPL_DEBUG("1before robot moved, old start ID = [%d] \n", pSearchStateSpace_->searchstartstate->StateID);
                SBPL_FPRINTF(fDeb, "1before robot moved, old start ID = [%d] \n", pSearchStateSpace_->searchstartstate->StateID);
#endif
                start_moved(pSearchStateSpace_);
                pSearchStateSpace_->robot_steps++;  // ++ robot_steps
                pSearchStateSpace_->total_robot_steps++;  // ++ robot_steps

#ifdef XIAOXUN_DEBUG
                SBPL_DEBUG("2after robot moved, new start ID = [%d] \n", pSearchStateSpace_->searchstartstate->StateID);
                SBPL_FPRINTF(fDeb, "2after robot moved, new start ID = [%d] \n", pSearchStateSpace_->searchstartstate->StateID);
#endif
            }

            //  if(pSearchStateSpace_->searchstartstate == pSearchStateSpace_->searchgoalstate || pSearchStateSpace_->searchstartstate == pSearchStateSpace_->temp_goal_state)
            if (pSearchStateSpace_->searchstartstate == pSearchStateSpace_->temp_goal_state) {
                SBPL_INFO("**               robot caught the target                          **\n");
                break;
            }



//      if( (double) (pSearchStateSpace_->hunter_movecost_testcase * SPEED_RATIO) > (double) pSearchStateSpace_->target_movecost_testcase )
            //  if(pSearchStateSpace_->robot_steps % SPEED == 0)  // SPEED (robot : target) = 10 : 1
            if (pSearchStateSpace_->robot_steps % SPEED != 0) { // SPEED (robot : target) = 10 : 9
#ifdef XIAOXUN_DEBUG
                SBPL_DEBUG("call goal_mvoed_________________________________\n");
                SBPL_FPRINTF(fDeb, "call goal_mvoed_________________________________\n");
#endif
                // if retv == 0, goal stay in previous path, no need to re-plan
                // if retv == 1, move off the path,  need replan
                ////////////////////////////////////////////////
                //         Goal moves 1 step here             //
                ////////////////////////////////////////////////
                int retu_value = goal_moved(pSearchStateSpace_);         // ps->goal has NOT been reset here, ONLY the ps->temp_goal is set

                pSearchStateSpace_->target_steps++;  // ++ target_steps

                if (retu_value == 1)
                    robot_need_to_replan = true; //target moved off the path



#ifdef XIAOXUN_DEBUG
                SBPL_DEBUG("@@@@@@@@@  target moved to id--> [%d]\n", pSearchStateSpace_->temp_goal_state->StateID);
                SBPL_FPRINTF(fDeb, "@@@@@@@@@  target moved to id--> [%d]\n", pSearchStateSpace_->temp_goal_state->StateID);
#endif


                if (pSearchStateSpace_->searchstartstate == pSearchStateSpace_->temp_goal_state) {
                    SBPL_INFO("**               target moved to the robot                        **\n");
                    break;
                }


#ifdef XIAOXUN_DEBUG
                SBPL_DEBUG("\n  retu_value == [%d]\n", retu_value);
                SBPL_FPRINTF(fDeb, "\n  retu_value == [%d]\n", retu_value);

                if (retu_value == 1) {
                    SBPL_DEBUG("target moved off the path, so need replan \n");
                    SBPL_FPRINTF(fDeb, "target moved off the path, so need replan \n");
                } else {
                    SBPL_DEBUG("target stayed in the path, no re-plan needed \n");
                    SBPL_FPRINTF(fDeb, "target stayed in the path, no re-plan needed \n");
                }

                SBPL_DEBUG("finished goal_mvoed__________________________\n");
                SBPL_FPRINTF(fDeb, "finished goal_mvoed__________________________\n");
#endif


            }// end if (%SPEED != 0)


            if (robot_need_to_replan == true) {
#ifdef ADAPTIVE_H
#ifdef REVERSE
                pSearchStateSpace_->temp_start_state = pSearchStateSpace_->searchstartstate;  // buffer the current hunter state
                pSearchStateSpace_->searchstartstate = pSearchStateSpace_->previous_search_start_state;

                AAState *F_tempstart = (AAState *)pSearchStateSpace_->temp_start_state->PlannerSpecificData;
                AA_ReInitializeState(F_tempstart, pSearchStateSpace_);    // (1)
                pSearchStateSpace_->keymodifer = pSearchStateSpace_->keymodifer + F_tempstart->h;

                pSearchStateSpace_->searchstartstate = pSearchStateSpace_->temp_start_state;  // copy the current hunter state back


#else  // forward _______________________________

                AAState *F_tempgoal = (AAState *)pSearchStateSpace_->temp_goal_state->PlannerSpecificData;
                AA_ReInitializeState(F_tempgoal, pSearchStateSpace_);    // (1)
                pSearchStateSpace_->keymodifer = pSearchStateSpace_->keymodifer + F_tempgoal->h;
#endif
                //SBPL_DEBUG("F_tempgoal->h == %d \n", F_tempgoal->h);
#endif
                pSearchStateSpace_->searchgoalstate = pSearchStateSpace_->temp_goal_state;

                ReInitializeNewSearch(pSearchStateSpace_);
            }



        }// end while (start != goal)

        Case_TimeEnded = clock(); // the start time of the run (= n+searches until the target is caught)
        totaltime += (Case_TimeEnded - Case_TimeStarted) / CLOCKS_PER_SEC;
    }// end for (RUNS)





    time_per_search = (double)totaltime / (double)(pSearchStateSpace_->totalsearchiteration);



    unsigned long int total_expansion = 0;
    unsigned long int expansion_per_search = 0;
    long double variance_expa_persearch = 0;
    long double expa_SDOM;




    for (int j = 1; j <= pSearchStateSpace_->totalsearchiteration; j++) {
        total_expansion += pSearchStateSpace_->expansion_of_the_search[j];
        SBPL_FPRINTF(fDeb, "%d \n", pSearchStateSpace_->expansion_of_the_search[j]); // write the expansion of each search to a file
    }

    expansion_per_search = total_expansion / pSearchStateSpace_->totalsearchiteration;


    for (int j = 1; j <= pSearchStateSpace_->totalsearchiteration; j++)
        variance_expa_persearch += pow((pSearchStateSpace_->expansion_of_the_search[j] - expansion_per_search), 2);

    variance_expa_persearch = pow((variance_expa_persearch /  pSearchStateSpace_->totalsearchiteration), 0.5);

    expa_SDOM = variance_expa_persearch / (double)pow(pSearchStateSpace_->totalsearchiteration, 0.5);





/////////////////////////////////////////////////////////////////////////////////////
#ifdef ADAPTIVE_H
    SBPL_INFO(" Adaptive A*\n");
#else
    SBPL_INFO(" Pure A*\n");
#endif

#ifdef REVERSE
    SBPL_INFO("Backwards search\n");
#else
    SBPL_INFO("Forwards search\n");
#endif

    SBPL_INFO("Test case  RUNS        == [%10d]\n", RUNS);
    SBPL_INFO("Searches per test case == [%10d]\n", pSearchStateSpace_->totalsearchiteration / RUNS);
    SBPL_INFO("Robot total moves step == [%10d]\n", pSearchStateSpace_->total_robot_steps);
    SBPL_INFO("Robot moves step/ RUN  == [%10d]\n", (pSearchStateSpace_->total_robot_steps)  / RUNS);
    SBPL_INFO("Robot moves cost/ RUN  == [%10d]\n", (pSearchStateSpace_->totalmovecost)      / RUNS);
    SBPL_INFO("expansion / search     == [%10d]\n", expansion_per_search);
    SBPL_INFO("SDOM_expansion         == [%10d]\n", expa_SDOM);
    SBPL_INFO("Search total time      == [%10f]\n", totaltime);
    SBPL_INFO("Time per search        == [%10f]\n", time_per_search);

    SBPL_FPRINTF(fDeb, "Test case  RUNS        == [%10d]\n", RUNS);
    SBPL_FPRINTF(fDeb, "Searches per test case == [%10d]\n", pSearchStateSpace_->totalsearchiteration / RUNS);
    SBPL_FPRINTF(fDeb, "Robot moves step       == [%10d]\n", pSearchStateSpace_->total_robot_steps);
    SBPL_FPRINTF(fDeb, "Robot moves step/ RUN  == [%10d]\n", (pSearchStateSpace_->total_robot_steps)  / RUNS);
    SBPL_FPRINTF(fDeb, "Robot moves cost/ RUN  == [%10d]\n", (pSearchStateSpace_->totalmovecost)      / RUNS);
    SBPL_FPRINTF(fDeb, "expansion / search     == [%10d]\n", expansion_per_search);
    SBPL_FPRINTF(fDeb, "SDOM_expansion         == [%10d]\n", expa_SDOM);
    SBPL_FPRINTF(fDeb, "Search total time      == [%10f]\n", totaltime);
    SBPL_FPRINTF(fDeb, "Time per search        == [%10f]\n", time_per_search);

    return (int) bFound;
}


int AAPlanner::set_goal(int goal_stateID) {

    SBPL_DEBUG("planner: setting goal to %d\n", goal_stateID);
    environment_->PrintState(goal_stateID, true, stdout);

    if (bforwardsearch) { //Xiaoxun 1: if search forward
        if (SetSearchGoalState(goal_stateID, pSearchStateSpace_) != 1) {
            SBPL_ERROR("ERROR: failed to set search goal state\n");
            return 0;
        }
    } else { // reverse
        if (SetSearchGoalState(goal_stateID, pSearchStateSpace_) != 1) {
            SBPL_ERROR("ERROR: failed to set search goal state\n");
            return 0;
        }

//      SBPL_ERROR("___error__AA* can only search forwards____________\n");
//      exit(0);
    }

    return 1;
}

//Xiaoxun: set the start_state
int AAPlanner::set_start(int start_stateID) {
    SBPL_DEBUG("AA planner: setting start to %d\n", start_stateID);
    environment_->PrintState(start_stateID, true, stdout);


    if (SetSearchStartState(start_stateID, pSearchStateSpace_) != 1) {
        SBPL_ERROR("ERROR: failed to set search start state\n");
        return 0;
    }


    return 1;
}



void AAPlanner::costs_changed(StateChangeQuery const &stateChange) {


    pSearchStateSpace_->bReinitializeSearchStateSpace = true;


}

void AAPlanner::costs_changed() {

    pSearchStateSpace_->bReinitializeSearchStateSpace = true;

}



int AAPlanner::force_planning_from_scratch() {
    SBPL_DEBUG("planner: forceplanfromscratch set\n");

    pSearchStateSpace_->bReinitializeSearchStateSpace = true;

    return 1;
}


int AAPlanner::set_search_mode(bool bSearchUntilFirstSolution) { // set to true
// Xiaoxun 1: this is to set whether AA* will keep searching until the first solution is found
    SBPL_DEBUG("planner: search mode set to %d\n", bSearchUntilFirstSolution);

    bsearchuntilfirstsolution = bSearchUntilFirstSolution;

    return 1;
}


void AAPlanner::print_searchpath(FILE *fOut) {
    PrintSearchPath(pSearchStateSpace_, fOut);
}











// Xiaoxun's AA* functions
// =======================
#if USE_XIAOXUN_STUFF

//function 1:
// state is the oldsearchgoalstate
AAState *AAPlanner::ChooseOneSuccs_for_TargetMove(AAState *state, AASearchStateSpace_t *pSearchStateSpace) {
    //CKey key;

    // Get successors from environment
    vector<int> SuccIDV;
    vector<int> CostV;
    environment_->GetSuccs(state->MDPstate->StateID, &SuccIDV, &CostV);
    int successorCount = (int)SuccIDV.size();

    state->successors.clear();
    // state->successors.ensureCapacity(successorCount);
    // Populate node successor stubs
    for (int i=0; i<successorCount; i++) {
        neighborStub nS;
        nS.id = SuccIDV[i];
        nS.cost = CostV[i];
        state->successors.push_back(nS);
    }

//SBPL_DEBUG("2. in Choose_OneSuccs_for_TargetMove,  SuccIDV.size() == %d \n", SuccIDV.size());


// 2 generate a random number among all succerssors
    int rand_succ = rand() % successorCount ;

#ifdef XIAOXUN_DEBUG
    SBPL_DEBUG("3. in Choose_OneSuccs_for_TargetMove,  rand_succ == %d \n", rand_succ);
    SBPL_FPRINTF(fDeb, "new goal move cost == %d,    \n", state->SuccsCostV[rand_succ]);
#endif

    CMDPSTATE *SuccMDPState = GetState((state->successors[rand_succ]).id, pSearchStateSpace);
// 3 succstate will be the next goal state
    AAState *succstate = (AAState *)(SuccMDPState->PlannerSpecificData);



// post processing for generated empty state
    for(neighborStub nS : state->successors) {
        CMDPSTATE *SuccMDPState2 = GetState(nS.id, pSearchStateSpace);
        AAState *F_succstate = (AAState *)(SuccMDPState2->PlannerSpecificData);

        if (F_succstate->bestpredstate == NULL) {
            F_succstate->callnumberaccessed = 0;
#ifdef XIAOXUN_DEBUG
            SBPL_FPRINTF(fDeb, "post processing here once, for  Choose_OneSuccs_for_TargetMove, i == [%d] \n", i);
#endif
        }
    }




    return succstate;
}


//function 2:
int AAPlanner::ReSetSearchGoalState(int SearchGoalStateID, AASearchStateSpace_t *pSearchStateSpace) {
//Xiaoxun1: if(1) ps->searchgoalstate has not been set,
//          if(2) goal has moved, so       ps->searchgoalstate->StateID != SearchGoalStateID




    if (pSearchStateSpace->searchgoalstate == NULL || pSearchStateSpace->searchgoalstate->StateID != SearchGoalStateID) {
        pSearchStateSpace->searchgoalstate = GetState(SearchGoalStateID, pSearchStateSpace);

        //should be new search iteration
//      pSearchStateSpace->eps_satisfied = INFINITECOST;
//      pSearchStateSpace_->eps = this->finitial_eps;


        /*
        #if USE_HEUR
                for(int i = 0; i < (int)pSearchStateSpace->searchMDP.StateArray.size(); i++)
                {
                    CMDPSTATE* MDPstate = pSearchStateSpace->searchMDP.StateArray[i];
                    AAState* state = (AAState*)MDPstate->PlannerSpecificData;
                    state->h = ComputeHeuristic(MDPstate, pSearchStateSpace);
                }

                pSearchStateSpace->bReevaluatefvals = true;
        #endif
        */

    }

    return 1;
}



//function 3:
int AAPlanner::reset_goal(int goal_stateID) {
#ifdef XIAOXUN_DEBUG
    SBPL_DEBUG("planner: resetting goal to %d\n", goal_stateID);
//  environment_->PrintState(goal_stateID, true, stdout);
#endif

    if (bforwardsearch) { //Xiaoxun 1: if search forward
        if (ReSetSearchGoalState(goal_stateID, pSearchStateSpace_) != 1) {
            SBPL_ERROR("ERROR: failed to set search goal state\n");
            return 0;
        }
    } else { // reverse
        SBPL_ERROR("___error__AA* can only search forwards____________\n");
        exit(0);
    }

    return 1;
}


//function 4:
//Xiaoxun 0: to reset the ps->searchgoalstate
int AAPlanner::goal_moved(AASearchStateSpace_t *pSearchStateSpace) {
    AAState *oldsearchgoalstate, *newsearchgoalstate;
    int planned_move_steps = TARGET_MOVE_STEPS;
    int retv;
    int target_move_cost;

//Currently: only compute next move for the target
//TODO: plan 100 moves for the goal
//  for(i = 0; i < planned_move_steps; i++)
//  {
    oldsearchgoalstate = (AAState *)(pSearchStateSpace->searchgoalstate->PlannerSpecificData);

//Xiaoxun: step 1, choose a new state for goal to go
    newsearchgoalstate = ChooseOneSuccs_for_TargetMove(oldsearchgoalstate, pSearchStateSpace);


////////////////////////////////////////////////////////////////////////////////////
//2010.02.28
// FFS...

    target_move_cost = GetTargetMoveCost(oldsearchgoalstate, newsearchgoalstate, pSearchStateSpace);
    pSearchStateSpace->target_movecost_testcase += target_move_cost;

    SBPL_DEBUG("target move cost == %d \n", target_move_cost);
//      pSearchStateSpace->target_movecost_testcase += GetTargetMoveCost(oldsearchgoalstate, newsearchgoalstate, pSearchStateSpace);
////////////////////////////////////////////////////////////////////////////////////

    //temp_goal_state (not reset ps->goal yet)
    pSearchStateSpace->temp_goal_state = newsearchgoalstate->MDPstate;






#ifdef REVERSE  // backwards search

    if ((AAState *)newsearchgoalstate->bestnextstate != oldsearchgoalstate)
        retv = 1;   // moved off the path
    else
        retv = 0;   // stay in previous path

#else           // forward search

    if ((AAState *)oldsearchgoalstate->bestpredstate != newsearchgoalstate) // if the target moved off the previous path
        retv = 1;   // moved off the path
    else
        retv = 0;   // stay in previous path

#endif
//  }//end for (i)



// if retv == 0, stay in path, no need to re-plan
// if retv == 1, move off the path, need to further check
    return retv;
}


//function 5:
int AAPlanner::ReSetSearchStartState(int SearchStartStateID, AASearchStateSpace_t *pSearchStateSpace) {
    CMDPSTATE *MDPstate = GetState(SearchStartStateID, pSearchStateSpace); //Xiaoxun 1:  see P13

    if (MDPstate !=  pSearchStateSpace->searchstartstate) {
        pSearchStateSpace->searchstartstate = MDPstate;
//Xiaoxun 0: diable this because even start (= agent) moved, AA* may NOT need a new search.
///     pSearchStateSpace->bReinitializeSearchStateSpace = true;
    }

    return 1;

}


//function 6:
//Xiaoxun: set the start_state
int AAPlanner::reset_start(int start_stateID) {
#ifdef XIAOXUN_DEBUG
    SBPL_DEBUG("AA planner: re-setting start to STATEID = [%d]\n", start_stateID);
    environment_->PrintState(start_stateID, true, stdout);
#endif


    if (ReSetSearchStartState(start_stateID, pSearchStateSpace_) != 1) {
        SBPL_ERROR("ERROR: failed to set search start state\n");
        return 0;
    }

    return 1;
}


//function 7:
//Xiaoxun 0: to reset the ps->searchstartstate
int AAPlanner::start_moved(AASearchStateSpace_t *pSearchStateSpace) {
    int movecost;
    AAState *oldsearchstartstate, *newsearchstartstate;




    oldsearchstartstate = (AAState *)(pSearchStateSpace->searchstartstate->PlannerSpecificData);
    //2010.06.09 : search forward or backward are the same
    newsearchstartstate = (AAState *)(oldsearchstartstate->bestnextstate->PlannerSpecificData);

#ifdef XIAOXUN_STATISTICS
// old version before 2010
/// pSearchStateSpace->totalmovecost += newsearchstartstate->g - oldsearchstartstate->g;

// 2010.02.28


#ifdef REVERSE
//      movecost = GetMoveCost(newsearchstartstate, pSearchStateSpace);
    movecost = oldsearchstartstate->g - newsearchstartstate->g;
    pSearchStateSpace->totalmovecost +=            movecost;          // GetMoveCost(newsearchstartstate, pSearchStateSpace);
    pSearchStateSpace->hunter_movecost_testcase += movecost;          // GetMoveCost(newsearchstartstate, pSearchStateSpace);
#else

    movecost = newsearchstartstate->g - oldsearchstartstate->g;
//      movecost = GetMoveCost(oldsearchstartstate, pSearchStateSpace);
    pSearchStateSpace->totalmovecost +=            movecost;          // GetMoveCost(oldsearchstartstate, pSearchStateSpace);
    pSearchStateSpace->hunter_movecost_testcase += movecost;          // GetMoveCost(oldsearchstartstate, pSearchStateSpace);
#endif
#endif

    SBPL_DEBUG(" -hunter cost = %d ", movecost);

    reset_start(newsearchstartstate->MDPstate->StateID); // set new start state for the (ps) search space


#ifdef XIAOXUN_DEBUG
    SBPL_DEBUG("______robot moved, show the information of old and new start \n");
    environment_->PrintState(oldsearchstartstate->MDPstate->StateID, true, stdout);
    environment_->PrintState(newsearchstartstate->MDPstate->StateID, true, stdout);
#endif



    return 1;  // if robot moved 1 step
}


//function 8:
//Xiaoxun optimize this function here
void AAPlanner::UpdateSuccs(AAState *state, AASearchStateSpace_t *pSearchStateSpace) {
    CKey key;
    AAState *succstate;
    int cost;



    // Get successors from environment
    vector<int> SuccIDV;
    vector<int> CostV;
    environment_->GetSuccs(state->MDPstate->StateID, &SuccIDV, &CostV);

    int successorCount = SuccIDV.size();
    state->successors.clear();
    // TODO: state->successors.ensureCapacity(successorCount);
    // REVIEW: this probably can be done faster (at least by changing the env API)

    for(int i=0; i<successorCount; i++) {
        neighborStub nS;
        nS.id = SuccIDV[i];
        nS.cost = CostV[i];
        state->successors.push_back(nS);
    }


    for(neighborStub nS : state->successors) {
        CMDPSTATE *SuccMDPState = GetState(nS.id, pSearchStateSpace);
        succstate = (AAState *)(SuccMDPState->PlannerSpecificData);
        cost = nS.cost;

        if (succstate->iterationclosed != pSearchStateSpace->searchiteration) {
#ifdef ADAPTIVE_H // Adaptive A*
            AA_ReInitializeState(succstate, pSearchStateSpace);   // (2)
#else             // pure A*

            if (succstate->generated_iteration != pSearchStateSpace->searchiteration)  // calculate h-value(s) here
                ReInitializeSearchStateInfo(succstate, pSearchStateSpace);

#endif

            if (succstate->g > state->g + cost) {
                succstate->g = state->g + cost;
                succstate->bestpredstate = state->MDPstate;   //Xiaoxun 3: ->bestpredstate == (my)->searchtree

                key.key[0] = succstate->g + (int)(succstate->h);
#ifdef TIE_LARGE_G
                key.key[1] = succstate->h;
#endif

                if (succstate->heapindex != 0)
                    pSearchStateSpace->heap->updateheap(succstate, key);
                else
                    pSearchStateSpace->heap->insertheap(succstate, key);
            }
        }
    } //for (state->successors)



    return;
}


//function 9
//Xiaoxun optimize this function for AA*
void AAPlanner::UpdatePreds(AAState *state, AASearchStateSpace_t *pSearchStateSpace) {

///////////////    new version    ////////////////////////////////////////////
    CKey key;
    AAState *predstate;


    // Get predecessors from environment
    vector<int> PredIDV;
    vector<int> CostV;
    environment_->GetPreds(state->MDPstate->StateID, &PredIDV, &CostV);
    int predecessorCount = PredIDV.size();

    state->predecessors.clear();
    // TODO: state->predecessors.ensureCapacity(predecessorCount)

    for (int i=0; i<predecessorCount; i++) {
        neighborStub nS;
        nS.id = PredIDV[i];
        nS.cost = CostV[i];
        state->predecessors.push_back(nS);
    }

    int cost;
    for(neighborStub nS : state->predecessors) {
        CMDPSTATE *PredMDPState = GetState(nS.id, pSearchStateSpace);
        predstate = (AAState *)(PredMDPState->PlannerSpecificData);
        cost = nS.cost;

        if (predstate->iterationclosed != pSearchStateSpace->searchiteration) {
#ifdef ADAPTIVE_H // Adaptive A*
            AA_ReInitializeState(predstate, pSearchStateSpace);   // (3)
#else             // pure A*

            if (predstate->generated_iteration != pSearchStateSpace->searchiteration)  // calculate h-value(s) here
                ReInitializeSearchStateInfo(predstate, pSearchStateSpace);

#endif

            if (predstate->g > state->g + cost) {
                predstate->g = state->g + cost;
                predstate->bestnextstate = state->MDPstate;   //Xiaoxun 3: ->bestnextstate == (predecessor)->searchtree

                key.key[0] = predstate->g + (int)(predstate->h);
#ifdef TIE_LARGE_G
                key.key[1] = predstate->h;
#endif

                if (predstate->heapindex != 0)
                    pSearchStateSpace->heap->updateheap(predstate, key);
                else
                    pSearchStateSpace->heap->insertheap(predstate, key);
            }
        }
    } //for (sind)


    return;
}


//function 12:
//Xiaoxun:  reset for the new test case
void AAPlanner::Reset_New_Test_Case(AASearchStateSpace_t *pSearchStateSpace) {
    int newstart_id;
    int newgoal_id;



    srand(pSearchStateSpace->callnumber);

    //newstart_id = environment_->SetStart_NewTestCase(pSearchStateSpace->callnumber);
    //newgoal_id  = environment_->SetGoal_NewTestCase(pSearchStateSpace->callnumber);

    pSearchStateSpace->searchstartstate = GetState(newstart_id, pSearchStateSpace);
    pSearchStateSpace->searchgoalstate  = GetState(newgoal_id,  pSearchStateSpace);

    pSearchStateSpace->temp_goal_state = pSearchStateSpace->searchgoalstate ;
    pSearchStateSpace->robot_steps = 0;
#ifdef ADAPTIVE_H
    pSearchStateSpace->keymodifer = 0;
#endif

// 2010.02.28
    pSearchStateSpace->hunter_movecost_testcase = 0;   // (2)
    pSearchStateSpace->target_movecost_testcase = 0;




#ifdef XIAOXUN_DEBUG
    SBPL_DEBUG("When Rest_New_Test_Case_______________________\n");
    SBPL_FPRINTF(fDeb, "When Rest_New_Test_Case_______________________\n");
    environment_->PrintState(newstart_id, true, fDeb);
    environment_->PrintState(newgoal_id,  true, fDeb);
#endif



    for (int i = 0; i < (int)pSearchStateSpace->searchMDP.StateArray.size(); i++) {
        CMDPSTATE *MDPstate = pSearchStateSpace->searchMDP.StateArray[i];
        AAState *state = (AAState *)MDPstate->PlannerSpecificData;

        state->g = INFINITECOST;
        state->iterationclosed = 0;
        state->callnumberaccessed = 0;
        state->bestnextstate = NULL;
        state->bestpredstate = NULL;
        state->costtobestnextstate = INFINITECOST;
        state->heapindex = 0;
        state->generated_iteration = 0;
    }

    ReInitializeSearchStateSpace(pSearchStateSpace);  // do call_number++;
    ReInitializeNewSearch(pSearchStateSpace);         // set start->g = 0    and do     insertheap(start);
    return;
}


//Xiaoxun: this function is to prepare for the new search iteration
// (1) search iteration ++;
// (2) empty OPEN list
// (3) insert start into OPEN list
void AAPlanner::ReInitializeNewSearch(AASearchStateSpace_t *pSearchStateSpace) {
    CKey key;

#ifdef XIAOXUN_DEBUG
    SBPL_DEBUG("call ReInitialize_NewSearch, after [%d]-th search \n", pSearchStateSpace->searchiteration);
    SBPL_FPRINTF(fDeb, "call ReInitialize_NewSearch, after [%d]-th search \n", pSearchStateSpace->searchiteration);
#endif


    pSearchStateSpace->searchiteration++;              // search_iteration+++    (1)
    pSearchStateSpace->totalsearchiteration++;
    pSearchStateSpace->heap->makeemptyheap();          // empty OPEN             (2)





#ifdef REVERSE //________________________search backwards_____________________________________________________________
    //initialize search start state (= goal state)
    AAState *goalstateinfo = (AAState *)(pSearchStateSpace->searchgoalstate->PlannerSpecificData);
#ifdef ADAPTIVE_H
    AA_ReInitializeState(goalstateinfo, pSearchStateSpace);     // (4)
#else

    if (goalstateinfo->generated_iteration != pSearchStateSpace->searchiteration)
        ReInitializeSearchStateInfo(goalstateinfo, pSearchStateSpace);

#endif


    goalstateinfo->g = 0;                              // insert search start (= goal state) into OPEN (3)
    //insert start state into the heap
    key.key[0] = (int)(goalstateinfo->h);
#ifdef TIE_LARGE_G
    key.key[1] = goalstateinfo->h;
#endif
    pSearchStateSpace->heap->insertheap(goalstateinfo, key);

#else // ________________________________search forward________________________________________________________________
    //initialize start state
    AAState *startstateinfo = (AAState *)(pSearchStateSpace->searchstartstate->PlannerSpecificData);
#ifdef ADAPTIVE_H
    AA_ReInitializeState(startstateinfo, pSearchStateSpace);    // (5)
#else

    if (startstateinfo->generated_iteration != pSearchStateSpace->searchiteration)
        ReInitializeSearchStateInfo(startstateinfo, pSearchStateSpace);

#endif
    startstateinfo->g = 0;                             // insert start into OPEN (3)
    //insert start state into the heap
    key.key[0] = (int)(startstateinfo->h);
#ifdef TIE_LARGE_G
    key.key[1] = startstateinfo->h;
#endif
    pSearchStateSpace->heap->insertheap(startstateinfo, key);
#endif  // end search forward


    return;
}


void AAPlanner::AA_ReInitializeState(AAState *state, AASearchStateSpace_t *pSearchStateSpace) {
    int simple_h;


    if (state->generated_iteration == 0) {
        state->g = INFINITECOST;
        state->iterationclosed = 0;   // xiaoxun 11111111111111111   ????
        state->heapindex = 0;         // xiaoxun 222222222222222222   ????
        state->costtobestnextstate = INFINITECOST;
        state->bestnextstate = NULL;
        state->bestpredstate = NULL;
        state->callnumberaccessed = pSearchStateSpace->callnumber;
        state->generated_iteration = pSearchStateSpace->searchiteration;
        //compute heuristics
        state->h = ComputeHeuristic(state->MDPstate, pSearchStateSpace);
    } else if (state->generated_iteration != pSearchStateSpace->searchiteration) {
        if (state->g + state->h < pSearchStateSpace->pathlength[state->generated_iteration])
            state->h = pSearchStateSpace->pathlength[state->generated_iteration] - state->g;

        state->h = state->h - (pSearchStateSpace->keymodifer - pSearchStateSpace->keymod[state->generated_iteration]);
        simple_h = ComputeHeuristic(state->MDPstate, pSearchStateSpace);

        if (state->h < simple_h)
            state->h = simple_h;

        state->g = INFINITECOST;
        state->bestnextstate = NULL;
        state->bestpredstate = NULL;
        state->callnumberaccessed = pSearchStateSpace->callnumber;
        state->generated_iteration = pSearchStateSpace->searchiteration;
    }

    return;
}


//Xiaoxun 13: this is to get the action cost of the robot to move 1 step along the path
int AAPlanner::GetMoveCost(AAState *state, AASearchStateSpace_t *pSearchStateSpace) {
    int actioncost;
    AAState *predstate;


#ifdef REVERSE
    SBPL_DEBUG("state->pred_num  == %d \n", state->pred_num);


    for (int sind = 0; sind < state->pred_num; sind++) {
        CMDPSTATE *MDPstate = GetState(state->PredsIDV[sind], pSearchStateSpace); //Xiaoxun 1:  see P13
        predstate = (AAState *)MDPstate->PlannerSpecificData;

        if (predstate->bestnextstate  == state->MDPstate) {
            actioncost = state->PredsCostV[sind];
            break;
        }
    }// end (sind)

#else

    for(neighborStub nS : state->successors)
        if(nS.id == state->bestnextstate->StateID){
            actioncost = nS.cost;
            break;
        }

#endif

//  SBPL_DEBUG("actioncost == %d \n",actioncost);

    return actioncost;
}

//2010.02.28
//Xiaoxun 14: this is to get the action cost of the target to move 1 step randomly
int AAPlanner::GetTargetMoveCost(AAState *state, AAState *state2, AASearchStateSpace_t *pSearchStateSpace) {
    int actioncost;

    for(neighborStub nS : state->successors)
        if (nS.id == state2->MDPstate->StateID) {
            actioncost = nS.cost;
            break;
        }

    return actioncost;
}

#endif
