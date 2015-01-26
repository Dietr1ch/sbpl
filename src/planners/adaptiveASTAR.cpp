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
#include <assert.h>
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


// Adaptive A* State
// =================
AAState::AAState(AAPlanner *planner, AASpace *space) {
    resetSearchInfo();
    resetStatistics();

    callnumberaccessed = space->callnumber;
    generated_iteration = 0;

    computeH(planner, space);
}

AAState::~AAState() {
}

void AAState::reinitialize(AAPlanner* planner, AASpace* space){
    int simple_h;

    if (generated_iteration == 0) {

        g = INFINITECOST;
        iterationclosed = 0;
        heapindex = 0;
        costtobestnextstate = INFINITECOST;
        bestSuccState = NULL;
        bestpredstate = NULL;

        callnumberaccessed = space->callnumber;
        generated_iteration = space->searchiteration;

        h = planner->ComputeHeuristic(MDPstate, space);
    }
    else{
        if (generated_iteration != space->searchiteration) {
            if (g + h < space->pathlength[generated_iteration])
                h = space->pathlength[generated_iteration] - g;

            h = h - (space->keymodifer - space->keymod[generated_iteration]);
            simple_h = planner->ComputeHeuristic(MDPstate, space);
            if (h < simple_h)
                h = simple_h;

            g = INFINITECOST;
            bestSuccState = NULL;
            bestpredstate = NULL;

            // REVIEW
            callnumberaccessed = space->callnumber;
            generated_iteration = space->searchiteration;
        }
    }
}



inline
void
AAState::resetSearchInfo() {
    g = INFINITECOST;
    iterationclosed = 0;

    heapindex = 0;
    costtobestnextstate = INFINITECOST;
    bestSuccState = NULL;
    bestpredstate = NULL;
}

inline
void
AAState::resetStatistics() {
#if STATISTICS
    expansions = 0;
#endif
}



inline
void
AAState::computeH(AAPlanner* planner, AASpace* space){
#if USE_HEUR
    if (space->searchgoalstate != NULL)
        h = planner->ComputeHeuristic(MDPstate, space);
    else
        h = 0;
#else
    h = 0;
#endif
}






// Adaptive A* Space
// =================
AASpace::AASpace() {
}

AASpace::~AASpace(){
}







// Adaptive A* Planner
// ===================

// Create - Init - Destroy
// -----------------------
AAPlanner::AAPlanner(DiscreteSpaceInformation *environment, bool bSearchForward) {
    SBPL_DEBUG("call AAstar constructor_________________________\n");
    printf("AAPlanner(%p, bool)\n", (void*)environment);

    bforwardsearch = bSearchForward;
    environment_ = environment;
    bsearchuntilfirstsolution = true;       // Xiaoxun: should I change this to "true" ???
//  finitial_eps = AA_DEFAULT_INITIAL_EPS;
//  finitial_eps = 1.0;   // xiaoxun modify this for AA*

    searchexpands = 0;
    MaxMemoryCounter = 0;
    fDeb = fopen("AAstar_debug.txt", "w");

    SBPL_DEBUG("debug on\n");

    pSearchSpace = new AASpace();

    // Create and initialize the search space
    if (CreateSearchStateSpace(pSearchSpace) != 1) {
        SBPL_ERROR("ERROR: failed to create the state space\n");
        return;
    }
    else
        SBPL_DEBUG("Search Space created\n");

    if (InitializeSearchStateSpace(pSearchSpace) != 1) {
        SBPL_ERROR("ERROR: failed to initialise the state space\n");
        return;
    }
    else
        SBPL_DEBUG("Search Space initialised\n");
}

AAPlanner::~AAPlanner() {
    if (pSearchSpace != NULL) {
        //delete the statespace
        DeleteSearchStateSpace(pSearchSpace);
        delete pSearchSpace;
    }

    fclose(fDeb);
}



//Xiaoxun 1: CreateState()
CMDPSTATE*
AAPlanner::CreateState(int stateID, AASpace *space) {

#if DEBUG
    if (environment_->StateID2IndexMapping[stateID][ARAMDP_STATEID2IND] != -1) {
        SBPL_ERROR("ERROR in CreateState: state already created\n");
        exit(1);
    }
#endif


    //adds to the tail a state
    CMDPSTATE *state = space->searchMDP.AddState(stateID); // Xiaoxun 2:  searchMDP == vector<CMDPSTATE*> StateArray;

    //remember the index of the state
    environment_->StateID2IndexMapping[stateID][ARAMDP_STATEID2IND] = space->searchMDP.StateArray.size() - 1;


#if DEBUG
    if (state != space->searchMDP.StateArray[environment_->StateID2IndexMapping[stateID][ARAMDP_STATEID2IND]]) {
        SBPL_ERROR("ERROR in CreateState: invalid state index\n");
        exit(1);
    }
#endif

    // Create Adaptive A* State Info
    AAState *aaState = new AAState(this, space);
    aaState->MDPstate = state;
    state->PlannerSpecificData = aaState;


    return state;
}

CMDPSTATE*
AAPlanner::GetState(int stateID, AASpace *pSearchStateSpace) {

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
inline
int
AAPlanner::ComputeHeuristic(CMDPSTATE *MDPstate, AASpace* space) {
    //compute heuristic for search
    if (bforwardsearch) {
        //forward search: heur = distance from state to searchgoal which is Goal AAState
        int retv =  environment_->GetFromToHeuristic(MDPstate->StateID, space->searchgoalstate->StateID);
        return retv;
    } else {
        //backward search: heur = distance from state to searchstart
        int retv = environment_->GetFromToHeuristic(MDPstate->StateID, space->searchstartstate->StateID);
        return retv;
    }
}

//Xiaoxun 0: re-initialization of a state
void
AAPlanner::ReInitializeSearchStateInfo(AAState *state, AASpace *space) {
    state->g = INFINITECOST;
    state->iterationclosed = 0;   // xiaoxun 11111111111111111   ????
    state->heapindex = 0;         // xiaoxun 222222222222222222   ????
    state->costtobestnextstate = INFINITECOST;
    state->bestSuccState = NULL;
    state->bestpredstate = NULL;

    state->callnumberaccessed = space->callnumber;
    state->generated_iteration = space->searchiteration;

    // Compute heuristics
#if USE_HEUR
    state->h = ComputeHeuristic(state->MDPstate, space);
#else
    state->h = 0;
#endif

}

void
AAPlanner::DeleteSearchStateData(AAState *state) {
    // No memory was allocated
    MaxMemoryCounter = 0;
}




int
AAPlanner::GetGVal(int StateID, AASpace *space) {
    CMDPSTATE *cmdp_state = GetState(StateID, space);
    AAState *state = (AAState *)cmdp_state->PlannerSpecificData;
    return state->g;
}


//Xiaoxun 0: function Improve Path()  ==  computeshortestpath() in grids.
//           returns 1 if the solution is found, 0 if the solution does not exist and 2 if it ran out of time
int
AAPlanner::ImprovePath(AASpace *space, double MaxNumofSecs) {
    int expands;
    AAState *searchgoalstate;
    //AAState *searchstartstate;  // unused?
    CKey key, minkey;
    CKey goalkey, startkey;

    expands = 0;

#ifdef ADAPTIVE_H
    pSearchStateSpace->keymod[pSearchStateSpace->searchiteration] = pSearchStateSpace->keymodifer;
#endif

    if (space->searchgoalstate == NULL) {
        SBPL_ERROR("ERROR searching: no goal state is set\n");
        exit(1);
    }

    if (space->searchstartstate == NULL) {
        SBPL_ERROR("ERROR searching: no start state is set\n");
        exit(1);
    }

#ifdef XIAOXUN_DEBUG
    SBPL_DEBUG("$$$$ call Improve Path , searchiteration == [%d], call = [%d]\n", pSearchStateSpace->searchiteration, pSearchStateSpace->callnumber);
    SBPL_FPRINTF(fDeb, "$$$$ call Improve Path , searchiteration == [%d], call = [%d]\n", pSearchStateSpace->searchiteration, pSearchStateSpace->callnumber);
#endif







#ifdef REVERSE  // search backward ____________(1)_____________

//Xiaoxun 1: set search goal state (=agent & start state ) and get the   (goalkey)
    searchstartstate = (AAState *)(space->searchstartstate->PlannerSpecificData);

    if (searchstartstate->generated_iteration != space->searchiteration)
        ReInitializeSearchStateInfo(searchstartstate, pSearchStateSpace);


    //set search goal key (= key of the agent)
    goalkey.key[0] = searchstartstate->g;
#ifdef TIE_LARGE_G
    goalkey.key[1] = searchstartstate->h;
#endif

    //expand states until done
    minkey = pSearchStateSpace->heap->getminkeyheap();    //Xiaoxun2: return the key of the Topheap()
    CKey oldkey = minkey;


    while (!space->heap->emptyheap() && minkey.key[0] < INFINITECOST && goalkey > minkey) {
        //Xiaoxun 2: pop out and get the state with the smallest f-value
        AAState *state = (AAState*) space->heap->deleteminheap();
        state->iterationclosed = space->searchiteration;   // expanded_iteration is set here
        expands++;


// Xiaoxun 3: expanding a state, here
///////////////////////////////////////////////////////////////

        UpdatePreds(state, space);    //search backwards

///////////////////////////////////////////////////////////////
        //recompute minkey
        minkey = space->heap->getminkeyheap();

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

    if (searchstartstate->g == INFINITECOST && space->heap->emptyheap()) {
//      SBPL_INFO("solution does not exist: search exited because heap is empty\n");
        retv = 0;
    } else {
        retv = 1;
#ifdef XIAOXUN_STATISTICS
        space->expansion_of_the_search[space->totalsearchiteration] = expands;
#endif

    }

    searchexpands += expands;





// ...
#else // forward search __________________(2)__________________________________________________________________

    //Xiaoxun 1: set goal state and get the   (goalkey)
    searchgoalstate = (AAState*) (space->searchgoalstate->PlannerSpecificData);

    if (searchgoalstate->generated_iteration != space->searchiteration)
        ReInitializeSearchStateInfo(searchgoalstate, space);

    // Set goal key
    goalkey.key[0] = searchgoalstate->g;
    printf("Goal key: %ld\n", goalkey.key[0]);
#ifdef TIE_LARGE_G
    goalkey.key[1] = searchgoalstate->h;
#endif

    // Expand states until done
    minkey = space->heap->getminkeyheap();    //Xiaoxun2: return the key of the Topheap()
    CKey oldkey = minkey;


    while (!space->heap->emptyheap() && minkey.key[0] < INFINITECOST && goalkey > minkey) {
        //Xiaoxun 2: pop out and get the state with the smallest f-value
        // FIXME: AAState is invalid sometimes
        AAState *state = (AAState*) space->heap->deleteminheap();

#ifdef XIAOXUN_DEBUG
        if (pSearchStateSpace->callnumber == 1)
            SBPL_FPRINTF(fDeb, "expanding state(%4d): h=%d g=%u key=%u iterclosed=%d callnuma=%d expands=%d (g(goal)=%u)\n",
                    state->MDPstate->StateID, state->h, state->g, state->g + (int)(pSearchStateSpace->eps * state->h),
                    state->iterationclosed, state->callnumberaccessed, state->numofexpands, searchgoalstate->g);
//        SBPL_FPRINTF(fDeb, "expanding: ");
//        PrintSearchState(state, fDeb);
        fflush(fDeb);
#endif
        TRACE("Expanding [%5d] (f:%7.1f) (g:%7.1f) (h:%7.1f)\n",
                   state->MDPstate->StateID,
                   (state->g + state->h)/1000.0,
                   state->g/1000.0,
                   state->h/1000.0);

        state->iterationclosed = space->searchiteration;   // expanded_iteration is set here
        expands++;

// Xiaoxun 3: expanding a state, here
///////////////////////////////////////////////////////////////
        UpdateSuccs(state, space);    //search forwards
///////////////////////////////////////////////////////////////

        // Recompute minkey
        minkey = space->heap->getminkeyheap();

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

    if (searchgoalstate->g == INFINITECOST && space->heap->emptyheap()) {
//      SBPL_INFO("solution does not exist: search exited because heap is empty\n");
        retv = 0;
    } else if (!space->heap->emptyheap() && goalkey > minkey) {
//      SBPL_INFO("search exited because it ran out of time\n");
        retv = 2;
    } else if (searchgoalstate->g == INFINITECOST && !space->heap->emptyheap()) {
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
int
AAPlanner::CreateSearchStateSpace(AASpace *space) {
    // Create heap
    space->heap = new CHeap;
    space->inconslist = new CList;

    // Account memory
    MaxMemoryCounter += sizeof(CHeap);
    MaxMemoryCounter += sizeof(CList);

    // Reset start and goal
    space->searchgoalstate = NULL;
    space->searchstartstate = NULL;

    // Reset statistics
    searchexpands = 0;


    space->bReinitializeSearchStateSpace = false;

    return 1; // OK (not-so-wise 'success' code
}

//deallocates memory used by SearchStateSpace
void
AAPlanner::DeleteSearchStateSpace(AASpace *space) {
    if (space->heap != NULL) {
        space->heap->makeemptyheap();
        delete space->heap;
        space->heap = NULL;
    }

    if (space->inconslist != NULL) {
        space->inconslist->makeemptylist(AA_INCONS_LIST_ID);
        delete space->inconslist;
        space->inconslist = NULL;
    }


    //delete the states themselves
    int iend = (int)space->searchMDP.StateArray.size();
    for (int i=0; i<iend; i++) {
        CMDPSTATE *state = space->searchMDP.StateArray[i];

        if (state != NULL && state->PlannerSpecificData != NULL) {
            DeleteSearchStateData((AAState *)state->PlannerSpecificData);
            free((AAState *)state->PlannerSpecificData);
            state->PlannerSpecificData = NULL;
        }
    }

    space->searchMDP.Delete();
}



//reset properly search state space
//needs to be done before deleting states
int
AAPlanner::ResetSearchStateSpace(AASpace *pSearchStateSpace) {
    pSearchStateSpace->heap->makeemptyheap();
    pSearchStateSpace->inconslist->makeemptylist(AA_INCONS_LIST_ID);

    return 1;
}

//initialization before each search
// Xiaoxun: note that: before EACH search
void
AAPlanner::ReInitializeSearchStateSpace(AASpace *pSearchStateSpace) {
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
int
AAPlanner::InitializeSearchStateSpace(AASpace *space) {

//  if(pSearchStateSpace->heap->currentsize != 0 || pSearchStateSpace->inconslist->currentsize != 0)
    if (space->heap->currentsize != 0) {
        SBPL_ERROR("ERROR in InitializeSearchStateSpace: heap or list is not empty\n");
        exit(1);
    }


    space->searchiteration = 0;
    space->callnumber = 0;

    space->robot_steps = 0;
    space->target_steps = 0;
    space->total_robot_steps = 0;
    space->totalmovecost = 0;


// 2010.02.28
    space->hunter_movecost_testcase = 0;  // (1)
    space->target_movecost_testcase = 0;


    //create and set the search start state
    space->searchgoalstate = NULL;
    space->searchstartstate = NULL;
#ifdef ADAPTIVE_H
    pSearchStateSpace->keymodifer = 0;
#endif


    space->bReinitializeSearchStateSpace = true;   // (1) this is the only place to set it to be "true"

    return 1;
}


int
AAPlanner::SetSearchGoalState(int SearchGoalStateID, AASpace *pSearchStateSpace) {
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


int
AAPlanner::SetSearchStartState(int SearchStartStateID, AASpace *pSearchStateSpace) {
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
int
AAPlanner::ReconstructPath(AASpace *space) {
    CMDPSTATE *MDPstate;
    CMDPSTATE *PredMDPstate;
    AAState *predstateinfo, *stateinfo;

    //nothing to do, if search is backward
    if (bforwardsearch) { // forward search for AA*
        //Xiaoxun added this if()
        if (space->temp_goal_state && space->searchgoalstate != space->temp_goal_state)
            MDPstate = space->temp_goal_state;
        else
            MDPstate = space->searchgoalstate;


#ifdef XIAOXUN_DEBUG
        SBPL_FPRINTF(fDeb, "reconstructing a path:\n");
#endif

        while (MDPstate != space->searchstartstate) {
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
            predstateinfo->bestSuccState = MDPstate;

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

        MDPstate = space->searchstartstate;

        while (MDPstate != space->searchgoalstate) {
            stateinfo = (AAState *)MDPstate->PlannerSpecificData;

            if (stateinfo->g == INFINITECOST) {
                SBPL_ERROR("ERROR in ReconstructPath: g of the state on the path is INFINITE\n");
                exit(1);
            }

            if (stateinfo->bestSuccState == NULL) {
                SBPL_ERROR("ERROR in ReconstructPath: bestpred is NULL\n");
                exit(1);
            }

            //get the parent state
            PredMDPstate = stateinfo->bestSuccState;
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



void
AAPlanner::PrintSearchPath(AASpace *space, FILE *fOut) {
    AAState *searchstateinfo;
    CMDPSTATE *state;
    int goalID;
    int PathCost;

    if (bforwardsearch) {
        state  = space->searchstartstate;

        if (space->temp_goal_state && space->temp_goal_state != space->searchgoalstate)
            goalID = space->temp_goal_state->StateID;
        else
            goalID = space->searchgoalstate->StateID;

        PathCost = ((AAState *) space->searchgoalstate->PlannerSpecificData)->g;

    } else { // search backward _________________________________________
//      state = pSearchStateSpace->searchgoalstate;
        state = space->searchstartstate;
        goalID = space->searchstartstate->StateID;   // goalID = agent->StateID

        PathCost = ((AAState *)space->searchstartstate->PlannerSpecificData)->g;
    }

    if (fOut == NULL)
        fOut = stdout;

//  PathCost = ((AAState*)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g;


    SBPL_FPRINTF(fOut, "Printing a path from state %d to the goal state %d\n",
            state->StateID, space->searchgoalstate->StateID);
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

        if (searchstateinfo->bestSuccState == NULL) {
            SBPL_FPRINTF(fOut, "path does not exist since bestnextstate == NULL  ????????????\n");
            break;
        }

        if (searchstateinfo->g == INFINITECOST) {
            SBPL_FPRINTF(fOut, "path does not exist since bestnextstate == NULL\n");
            break;
        }

        int costToGoal = PathCost - costFromStart;
        int transcost = searchstateinfo->g - ((AAState *)(searchstateinfo->bestSuccState->PlannerSpecificData))->g; //->v

        if (bforwardsearch)
            transcost = -transcost;

        costFromStart += transcost;

        SBPL_FPRINTF(fOut, "g=%d-->state %d, h = %d ctg = %d  ", searchstateinfo->g,
                searchstateinfo->bestSuccState->StateID, searchstateinfo->h, costToGoal);

        state = searchstateinfo->bestSuccState;

        environment_->PrintState(state->StateID, false, fOut);



    }
}

void
AAPlanner::PrintSearchState(AAState *state, FILE *fOut) {
    int expansions = -1;
#if STATISTICS
    expansions = state->expansions;
#endif
    SBPL_FPRINTF(fOut, "state %d: h=%d g=%u iterc=%d callnuma=%d expands=%d heapind=%d \n",
            state->MDPstate->StateID, state->h, state->g,
            state->iterationclosed, state->callnumberaccessed, expansions, state->heapindex);
    environment_->PrintState(state->MDPstate->StateID, true, fOut);

}



int
AAPlanner::getHeurValue(AASpace *pSearchStateSpace, int StateID) {
    CMDPSTATE *MDPstate = GetState(StateID, pSearchStateSpace);
    AAState *searchstateinfo = (AAState *)MDPstate->PlannerSpecificData;
    return searchstateinfo->h;
}

//Xiaoxun 0: this is to get the <STATE_ID> of all cells along the path
vector<stateID>
AAPlanner::GetSearchPath(AASpace *space, int &solcost) {
    vector<stateID> SuccIDV;
    vector<int> CostV;
    vector<stateID> wholePathIds;
    AAState *searchstateinfo;
    CMDPSTATE *state = NULL;
    CMDPSTATE *goalstate = NULL;
    CMDPSTATE *startstate = NULL;

    if (bforwardsearch) { // search forward
        startstate = space->searchstartstate;
        goalstate = space->searchgoalstate;

        //reconstruct the path by setting -->bestnextstate pointers appropriately
        ReconstructPath(space);   //setting the (my)-->trace pointers
    } else {
//      startstate = space->searchgoalstate;
//      goalstate = space->searchstartstate;
        startstate = space->searchstartstate;
        goalstate = space->searchgoalstate;
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

        if (searchstateinfo->bestSuccState == NULL) {
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
            if (SuccIDV.at(i) == searchstateinfo->bestSuccState->StateID) {
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
        AAState *nextstateinfo = (AAState *)(searchstateinfo->bestSuccState->PlannerSpecificData);

        if (actioncost != abs((int)(searchstateinfo->g - nextstateinfo->g)) && space->eps_satisfied <= 1.001) {
            SBPL_FPRINTF(fDeb, "ERROR: actioncost=%d is not matching the difference in g-values of %d\n",
                    actioncost, abs((int)(searchstateinfo->g - nextstateinfo->g)));
            SBPL_ERROR("ERROR: actioncost=%d is not matching the difference in g-values of %d\n",
                   actioncost, abs((int)(searchstateinfo->g - nextstateinfo->g)));
            PrintSearchState(searchstateinfo, fDeb);
            PrintSearchState(nextstateinfo, fDeb);
        }

#endif

        state = searchstateinfo->bestSuccState;
        wholePathIds.push_back(state->StateID);

    }// end while


    return wholePathIds;
}







// Search
// ======


//Xiaoxun 0:     Search() == ComputeShortestPath()...... however, it contains function IMPROVE_PATH() inside.
bool
AAPlanner::Search(AASpace *space, vector<stateID> &pathIds, int &PathCost, bool bFirstSolution, bool bOptimalSolution, double MaxNumofSecs) {
    CKey key;
    TimeStarted = clock();
    searchexpands = 0;



    if (space->bReinitializeSearchStateSpace == true) {
//#ifdef XIAOXUN_DEBUG
        SBPL_DEBUG("call AAPlanner::search() function ____+++______ReinitializeSearchSpace \n");
        SBPL_FPRINTF(fDeb, "call AAPlanner::search() function ____+++______ReinitializeSearchSpace \n");
//#endif
        //re-initialize state space
        ReInitializeSearchStateSpace(space);
        ReInitializeNewSearch(space);      // insert the search start state into OPEN,  start->g= 0 etc...
    }





    MaxNumofSecs = INFINITECOST;
    int prevexpands = 0;  // UNUSED
    clock_t loop_time;

    loop_time = clock();



//Xiaoxun added this for AA*  backward search:
#ifdef REVERSE
#ifdef ADAPTIVE_H
    pSearchStateSpace->previous_search_start_state = pSearchStateSpace->searchstartstate;
#endif
#endif



//call  improve_path()
    if (ImprovePath(space, MaxNumofSecs) == 1) {
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
    PathCost = ((AAState *)space->searchstartstate->PlannerSpecificData)->g;
#else           // search forwards
    PathCost = ((AAState *)space->searchgoalstate->PlannerSpecificData)->g;
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
        pathIds = GetSearchPath(space, solcost);   // get the STATE_ID of the cells along the path
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
int AAPlanner::replan(double allocated_time_secs, vector<stateID> *solution_stateIDs_V) {
    int solcost;
    return replan(allocated_time_secs, solution_stateIDs_V, &solcost);
}


//Xiaoxun 0: return 1 if found a solution, and 0 otherwise
//           re-plan is to do multiple searches until the goal is caught
int
AAPlanner::replan(double allocated_time_secs, vector<stateID> *solution_stateIDs_V, int *psolcost) {
    vector<stateID> pathIds;
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
        SBPL_INFO("Adaptive A* TEST CASE = [%d],  callnumber = [%d] \n", i, pSearchSpace->callnumber);
#else
        SBPL_INFO("A* TEST CASE = [%d],  callnumber = [%d] \n", i, pSearchSpace->callnumber);
#endif

        robot_need_to_replan = true;
        pSearchSpace->temp_goal_state = pSearchSpace->searchgoalstate;

        ////////////////////////////////////////////////////////////////////////////////////////
        if (pSearchSpace->callnumber > 0) {
#ifdef XIAOXUN_DEBUG
            SBPL_DEBUG("\n\n +++++++++++++++++++++++Reset for the (%d)-th TEST CASE   ++++++++++++++++++\n\n", pSearchSpace->callnumber);
            SBPL_FPRINTF(fDeb, "\n\n +++++++++++++++++++++++Reset for the (%d)-th TEST CASE   ++++++++++++++++++\n\n", pSearchSpace->callnumber);
#endif
            // TODO: wtf happens here?
            //Reset_New_Test_Case(pSearchSpace);
            printf("--- Reset_New_Test_Case ---\n");
        }

        ////////////////////////////////////////////////////////////////////////////////////////



        Case_TimeStarted = clock(); // the start time of the run (= n+searches until the target is caught)


        while (pSearchSpace->searchgoalstate != pSearchSpace->searchstartstate) {
            if (robot_need_to_replan == true) {

                SBPL_DEBUG("%d-th_s    ", pSearchSpace->searchiteration);

                if ((bFound = Search(pSearchSpace, pathIds, PathCost, bFirstSolution, bOptimalSolution, allocated_time_secs)) == false) {
                    SBPL_ERROR("AA* failed to find a solution after calling Search()... \n");
                    return (int)bFound;
                }

                robot_need_to_replan = false; // after compute path, set the flag to false.
                *solution_stateIDs_V = pathIds; //copy the solution
                *psolcost = PathCost;

                SBPL_DEBUG("  - Pathcost  == %d ", PathCost);



#ifdef XIAOXUN_DEBUG
                SBPL_DEBUG("\n\n ++++++++ 1 +++++++++++++++ The (%d)-th search is done  ++++++++++++++++++\n", pSearchSpace->searchiteration);
                SBPL_FPRINTF(fDeb, "\n\n ++++++++ 1 +++++++++++++++ The (%d)-th search is done  ++++++++++++++++++\n\n", pSearchSpace->searchiteration);
#endif
            }





            ////////////////////////////////////////////////
            //Xiaoxun 1:    Robot moves 1 step here       //
            ////////////////////////////////////////////////
            if (robot_need_to_replan == false) {
#ifdef XIAOXUN_DEBUG
                SBPL_DEBUG("1before robot moved, old start ID = [%d] \n", pSearchSpace->searchstartstate->StateID);
                SBPL_FPRINTF(fDeb, "1before robot moved, old start ID = [%d] \n", pSearchSpace->searchstartstate->StateID);
#endif
                start_moved(pSearchSpace);
                pSearchSpace->robot_steps++;  // ++ robot_steps
                pSearchSpace->total_robot_steps++;  // ++ robot_steps

#ifdef XIAOXUN_DEBUG
                SBPL_DEBUG("2after robot moved, new start ID = [%d] \n", pSearchSpace->searchstartstate->StateID);
                SBPL_FPRINTF(fDeb, "2after robot moved, new start ID = [%d] \n", pSearchSpace->searchstartstate->StateID);
#endif
            }

            //  if(pSearchSpace->searchstartstate == pSearchSpace->searchgoalstate || pSearchSpace->searchstartstate == pSearchSpace->temp_goal_state)
            if (pSearchSpace->searchstartstate == pSearchSpace->temp_goal_state) {
                SBPL_INFO("**               robot caught the target                          **\n");
                break;
            }



//      if( (double) (pSearchSpace->hunter_movecost_testcase * SPEED_RATIO) > (double) pSearchSpace->target_movecost_testcase )
            //  if(pSearchSpace->robot_steps % SPEED == 0)  // SPEED (robot : target) = 10 : 1
            if (pSearchSpace->robot_steps % SPEED != 0) { // SPEED (robot : target) = 10 : 9
#ifdef XIAOXUN_DEBUG
                SBPL_DEBUG("call goal_mvoed_________________________________\n");
                SBPL_FPRINTF(fDeb, "call goal_mvoed_________________________________\n");
#endif
                // if retv == 0, goal stay in previous path, no need to re-plan
                // if retv == 1, move off the path,  need replan
                ////////////////////////////////////////////////
                //         Goal moves 1 step here             //
                ////////////////////////////////////////////////
                int retu_value = goal_moved(pSearchSpace);         // ps->goal has NOT been reset here, ONLY the ps->temp_goal is set

                pSearchSpace->target_steps++;  // ++ target_steps

                if (retu_value == 1)
                    robot_need_to_replan = true; //target moved off the path



#ifdef XIAOXUN_DEBUG
                SBPL_DEBUG("@@@@@@@@@  target moved to id--> [%d]\n", pSearchSpace->temp_goal_state->StateID);
                SBPL_FPRINTF(fDeb, "@@@@@@@@@  target moved to id--> [%d]\n", pSearchSpace->temp_goal_state->StateID);
#endif


                if (pSearchSpace->searchstartstate == pSearchSpace->temp_goal_state) {
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
                pSearchSpace->temp_start_state = pSearchSpace->searchstartstate;  // buffer the current hunter state
                pSearchSpace->searchstartstate = pSearchSpace->previous_search_start_state;

                AAState *F_tempstart = (AAState*) pSearchSpace->temp_start_state->PlannerSpecificData;
                AA_ReInitializeState(F_tempstart, pSearchSpace);    // (1)
                pSearchSpace->keymodifer = pSearchSpace->keymodifer + F_tempstart->h;

                pSearchSpace->searchstartstate = pSearchSpace->temp_start_state;  // copy the current hunter state back


#else  // forward _______________________________

                AAState *F_tempgoal = (AAState *)pSearchSpace->temp_goal_state->PlannerSpecificData;
                AA_ReInitializeState(F_tempgoal, pSearchSpace);    // (1)
                pSearchSpace->keymodifer = pSearchSpace->keymodifer + F_tempgoal->h;
#endif
                //SBPL_DEBUG("F_tempgoal->h == %d \n", F_tempgoal->h);
#endif
                pSearchSpace->searchgoalstate = pSearchSpace->temp_goal_state;

                ReInitializeNewSearch(pSearchSpace);
            }



        }// end while (start != goal)

        Case_TimeEnded = clock(); // the start time of the run (= n+searches until the target is caught)
        totaltime += (Case_TimeEnded - Case_TimeStarted) / CLOCKS_PER_SEC;
    }// end for (RUNS)





    time_per_search = (double)totaltime / (double)(pSearchSpace->totalsearchiteration);



    unsigned long int total_expansion = 0;
    unsigned long int expansion_per_search = 0;
    long double variance_expa_persearch = 0;
    long double expa_SDOM;




    for (long unsigned int j = 1; j <= pSearchSpace->totalsearchiteration; j++) {
        total_expansion += pSearchSpace->expansion_of_the_search[j];
        SBPL_FPRINTF(fDeb, "%d \n", pSearchSpace->expansion_of_the_search[j]); // write the expansion of each search to a file
    }

    expansion_per_search = total_expansion / pSearchSpace->totalsearchiteration;


    for (unsigned int j = 1; j <= pSearchSpace->totalsearchiteration; j++)
        variance_expa_persearch += pow((pSearchSpace->expansion_of_the_search[j] - expansion_per_search), 2);

    variance_expa_persearch = pow((variance_expa_persearch /  pSearchSpace->totalsearchiteration), 0.5);

    expa_SDOM = variance_expa_persearch / (double)pow(pSearchSpace->totalsearchiteration, 0.5);





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
    SBPL_INFO("Searches per test case == [%10d]\n", pSearchSpace->totalsearchiteration / RUNS);
    SBPL_INFO("Robot total moves step == [%10d]\n", pSearchSpace->total_robot_steps);
    SBPL_INFO("Robot moves step/ RUN  == [%10d]\n", (pSearchSpace->total_robot_steps)  / RUNS);
    SBPL_INFO("Robot moves cost/ RUN  == [%10d]\n", (pSearchSpace->totalmovecost)      / RUNS);
    SBPL_INFO("expansion / search     == [%10d]\n", expansion_per_search);
    SBPL_INFO("SDOM_expansion         == [%10d]\n", expa_SDOM);
    SBPL_INFO("Search total time      == [%10f]\n", totaltime);
    SBPL_INFO("Time per search        == [%10f]\n", time_per_search);

    SBPL_FPRINTF(fDeb, "Test case  RUNS        == [%10d]\n", RUNS);
    SBPL_FPRINTF(fDeb, "Searches per test case == [%10d]\n", pSearchSpace->totalsearchiteration / RUNS);
    SBPL_FPRINTF(fDeb, "Robot moves step       == [%10d]\n", pSearchSpace->total_robot_steps);
    SBPL_FPRINTF(fDeb, "Robot moves step/ RUN  == [%10d]\n", (pSearchSpace->total_robot_steps)  / RUNS);
    SBPL_FPRINTF(fDeb, "Robot moves cost/ RUN  == [%10d]\n", (pSearchSpace->totalmovecost)      / RUNS);
    SBPL_FPRINTF(fDeb, "expansion / search     == [%10d]\n", expansion_per_search);
    SBPL_FPRINTF(fDeb, "SDOM_expansion         == [%10d]\n", expa_SDOM);
    SBPL_FPRINTF(fDeb, "Search total time      == [%10f]\n", totaltime);
    SBPL_FPRINTF(fDeb, "Time per search        == [%10f]\n", time_per_search);

    return (int) bFound;
}


int
AAPlanner::set_goal(stateID goal_stateID) {

    SBPL_DEBUG("planner: setting goal to %d\n", goal_stateID);
    environment_->PrintState(goal_stateID, true, stdout);

    if (SetSearchGoalState(goal_stateID, pSearchSpace) != 1) {
        SBPL_ERROR("ERROR: failed to set search goal state\n");
        return 0;
    }

    return 1;
}

//Xiaoxun: set the start_state
int
AAPlanner::set_start(stateID start_stateID) {
    SBPL_DEBUG("AA planner: setting start to %d\n", start_stateID);
    environment_->PrintState(start_stateID, true, stdout);


    if (SetSearchStartState(start_stateID, pSearchSpace) != 1) {
        SBPL_ERROR("ERROR: failed to set search start state\n");
        return 0;
    }

    return 1;
}



void
AAPlanner::costs_changed(StateChangeQuery const &stateChange) {
    pSearchSpace->bReinitializeSearchStateSpace = true;
}

void
AAPlanner::costs_changed() {
    pSearchSpace->bReinitializeSearchStateSpace = true;
}



int
AAPlanner::force_planning_from_scratch() {
    SBPL_DEBUG("planner: forceplanfromscratch set\n");

    pSearchSpace->bReinitializeSearchStateSpace = true;

    return 1;
}


int
AAPlanner::set_search_mode(bool bSearchUntilFirstSolution) { // set to true
// Xiaoxun 1: this is to set whether AA* will keep searching until the first solution is found
    SBPL_DEBUG("planner: search mode set to %d\n", bSearchUntilFirstSolution);

    bsearchuntilfirstsolution = bSearchUntilFirstSolution;

    return 1;
}


void
AAPlanner::print_searchpath(FILE *fOut) {
    PrintSearchPath(pSearchSpace, fOut);
}











// Xiaoxun's AA* functions
// =======================
#if USE_XIAOXUN_STUFF

//function 1:
// state is the oldsearchgoalstate
AAState*
AAPlanner::ChooseOneSuccs_for_TargetMove(AAState *state, AASpace *pSearchStateSpace) {
    //CKey key;

    // Get successors from environment
    vector<stateID> SuccIDV;
    vector<int> CostV;
    environment_->GetSuccs(state->MDPstate->StateID, &SuccIDV, &CostV);
    int successorCount = (int)SuccIDV.size();

    state->successors.clear();
    state->successors.reserve(successorCount);
    // Populate node successor stubs
    for (int i=0; i<successorCount; i++) {
        nodeStub nS;
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
    for(nodeStub nS : state->successors) {
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
int
AAPlanner::ReSetSearchGoalState(int SearchGoalStateID, AASpace *pSearchStateSpace) {
//Xiaoxun1: if(1) ps->searchgoalstate has not been set,
//          if(2) goal has moved, so       ps->searchgoalstate->StateID != SearchGoalStateID




    if (pSearchStateSpace->searchgoalstate == NULL || pSearchStateSpace->searchgoalstate->StateID != SearchGoalStateID) {
        pSearchStateSpace->searchgoalstate = GetState(SearchGoalStateID, pSearchStateSpace);

        //should be new search iteration
//      pSearchStateSpace->eps_satisfied = INFINITECOST;
//      pSearchSpace->eps = this->finitial_eps;


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
int
AAPlanner::reset_goal(int goal_stateID) {
#ifdef XIAOXUN_DEBUG
    SBPL_DEBUG("planner: resetting goal to %d\n", goal_stateID);
//  environment_->PrintState(goal_stateID, true, stdout);
#endif

    if (bforwardsearch) { //Xiaoxun 1: if search forward
        if (ReSetSearchGoalState(goal_stateID, pSearchSpace) != 1) {
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
int
AAPlanner::goal_moved(AASpace *pSearchStateSpace) {
    AAState *oldsearchgoalstate, *newsearchgoalstate;
    //int planned_move_steps = TARGET_MOVE_STEPS;  // unused
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
int
AAPlanner::ReSetSearchStartState(int SearchStartStateID, AASpace *space) {
    CMDPSTATE *MDPstate = GetState(SearchStartStateID, space); //Xiaoxun 1:  see P13

    if (MDPstate != space->searchstartstate) {
        space->searchstartstate = MDPstate;
//Xiaoxun 0: diable this because even start (= agent) moved, AA* may NOT need a new search.
///     space->bReinitializeSearchStateSpace = true;
    }

    return 1;
}


//function 6:
//Xiaoxun: set the start_state
int
AAPlanner::reset_start(int start_stateID) {
#ifdef XIAOXUN_DEBUG
    SBPL_DEBUG("AA planner: re-setting start to STATEID = [%d]\n", start_stateID);
    environment_->PrintState(start_stateID, true, stdout);
#endif


    if (ReSetSearchStartState(start_stateID, this->pSearchSpace) != 1) {
        SBPL_ERROR("ERROR: failed to set search start state\n");
        return 0;
    }

    return 1;
}


//function 7:
//Xiaoxun 0: to reset the ps->searchstartstate
int
AAPlanner::start_moved(AASpace *space) {
    int movecost;
    AAState *oldsearchstartstate, *newsearchstartstate;




    oldsearchstartstate = (AAState *)(space->searchstartstate->PlannerSpecificData);
    //2010.06.09 : search forward or backward are the same
    newsearchstartstate = (AAState *)(oldsearchstartstate->bestSuccState->PlannerSpecificData);

#ifdef XIAOXUN_STATISTICS
// old version before 2010
/// space->totalmovecost += newsearchstartstate->g - oldsearchstartstate->g;

// 2010.02.28
#ifdef REVERSE
    // movecost = GetMoveCost(newsearchstartstate, space);
    movecost = oldsearchstartstate->g - newsearchstartstate->g;
    space->totalmovecost += movecost;  // GetMoveCost(newsearchstartstate, space);
    space->hunter_movecost_testcase += movecost;  // GetMoveCost(newsearchstartstate, space);
#else
    movecost = newsearchstartstate->g - oldsearchstartstate->g;
    //movecost = GetMoveCost(oldsearchstartstate, space);
    space->totalmovecost += movecost;  // GetMoveCost(oldsearchstartstate, space);
    space->hunter_movecost_testcase += movecost;  // GetMoveCost(oldsearchstartstate, space);
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
void
AAPlanner::UpdateSuccs(AAState *state, AASpace *space) {

    TRACE("Expanding [%d]", state->MDPstate->StateID);

    // Get successors from environment
    vector<stateID> SuccIDV;
    vector<int> CostV;
    environment_->GetSuccs(state->MDPstate->StateID, &SuccIDV, &CostV);
    // REVIEW: this probably can be done faster (at least by changing the env API)
    int successorCount = SuccIDV.size();
    TRACE("State node successors: %lu", state->successors.size());
    TRACE("State map  successors: %d", successorCount);
    state->successors.clear();
    state->successors.reserve(successorCount);
    for(int i=0; i<successorCount; i++) {
        nodeStub nS;
        nS.id = SuccIDV[i];
        nS.cost = CostV[i];
        state->successors.push_back(nS);
    }

    // Reach successors
    for(nodeStub nS : state->successors) {
        CMDPSTATE* SuccMDPState = GetState(nS.id, space);
        AAState* succstate = (AAState *)(SuccMDPState->PlannerSpecificData);
        printf("Reaching [%d]\n", succstate->MDPstate->StateID);

        if (succstate->iterationclosed != space->searchiteration) {
#ifdef ADAPTIVE_H // Adaptive A*
            AA_ReInitializeState(succstate, space);   // (2)
#else             // pure A*
            if (succstate->generated_iteration != space->searchiteration)  // calculate h-value(s) here
                ReInitializeSearchStateInfo(succstate, space);
#endif

            assert(nS.cost>=0);
            unsigned int newG = state->g + nS.cost;
            if (newG < succstate->g) {
                succstate->g = newG;
                succstate->bestpredstate = state->MDPstate;   //Xiaoxun 3: ->bestpredstate == (my)->searchtree

                CKey key;
                key.key[0] = succstate->g + succstate->h;
#ifdef TIE_LARGE_G
                key.key[1] = succstate->h;
#endif

                printf("Adding [%d] to the heap\n", succstate->MDPstate->StateID);
                if(succstate->MDPstate->StateID==31){
                    printf("--\n");
                    printf("succesors: %ld\n", succstate->successors.size());
                    printf("predecesors: %ld\n", succstate->predecessors.size());
                }
                if (succstate->heapindex != 0)
                    space->heap->updateheap(succstate, key);
                else
                    space->heap->insertheap(succstate, key);
            }
        }
    }
}


//function 9
//Xiaoxun optimize this function for AA*
void
AAPlanner::UpdatePreds(AAState *state, AASpace *space) {

///////////////    new version    ////////////////////////////////////////////
    CKey key;
    AAState *predstate;


    // Get predecessors from environment
    vector<stateID> PredIDV;
    vector<int> CostV;
    environment_->GetPreds(state->MDPstate->StateID, &PredIDV, &CostV);
    int predecessorCount = PredIDV.size();

    state->predecessors.clear();
    state->predecessors.reserve(predecessorCount);

    for (int i=0; i<predecessorCount; i++) {
        nodeStub nS;
        nS.id = PredIDV[i];
        nS.cost = CostV[i];
        state->predecessors.push_back(nS);
    }

    int cost;
    for(nodeStub nS : state->predecessors) {
        CMDPSTATE *PredMDPState = GetState(nS.id, space);
        predstate = (AAState *)(PredMDPState->PlannerSpecificData);
        cost = nS.cost;

        if (predstate->iterationclosed != space->searchiteration) {
#ifdef ADAPTIVE_H // Adaptive A*
            AA_ReInitializeState(predstate, space);   // (3)
#else             // pure A*

            if (predstate->generated_iteration != space->searchiteration)  // calculate h-value(s) here
                ReInitializeSearchStateInfo(predstate, space);
#endif

            if (predstate->g > state->g + cost) {
                predstate->g = state->g + cost;
                predstate->bestSuccState = state->MDPstate;   //Xiaoxun 3: ->bestSuccState == (predecessor)->searchtree

                key.key[0] = predstate->g + (int)(predstate->h);
#ifdef TIE_LARGE_G
                key.key[1] = predstate->h;
#endif

                if (predstate->heapindex != 0)
                    space->heap->updateheap(predstate, key);
                else
                    space->heap->insertheap(predstate, key);
            }
        }
    } //for (sind)


    return;
}


//Xiaoxun: this function is to prepare for the new search iteration
// (1) search iteration ++;
// (2) empty OPEN list
// (3) insert start into OPEN list
void
AAPlanner::ReInitializeNewSearch(AASpace *space) {
    CKey key;

#ifdef XIAOXUN_DEBUG
    SBPL_DEBUG("call ReInitialize_NewSearch, after [%d]-th search \n", space->searchiteration);
    SBPL_FPRINTF(fDeb, "call ReInitialize_NewSearch, after [%d]-th search \n", space->searchiteration);
#endif

    space->searchiteration++;              // search_iteration+++    (1)
    space->totalsearchiteration++;
    space->heap->makeemptyheap();          // empty OPEN             (2)


#ifdef REVERSE //________________________search backwards_____________________________________________________________
    //initialize search start state (= goal state)
    AAState *goalstateinfo = (AAState *)(space->searchgoalstate->PlannerSpecificData);
#ifdef ADAPTIVE_H
    AA_ReInitializeState(goalstateinfo, space);     // (4)
#else

    if (goalstateinfo->generated_iteration != space->searchiteration)
        ReInitializeSearchStateInfo(goalstateinfo, space);

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
    AAState *startstateinfo = (AAState *)(space->searchstartstate->PlannerSpecificData);
#ifdef ADAPTIVE_H
    AA_ReInitializeState(startstateinfo, pSearchStateSpace);    // (5)
#else

    if (startstateinfo->generated_iteration != space->searchiteration)
        ReInitializeSearchStateInfo(startstateinfo, space);

#endif
    startstateinfo->g = 0;                             // insert start into OPEN (3)
    //insert start state into the heap
    key.key[0] = (int)(startstateinfo->h);
#ifdef TIE_LARGE_G
    key.key[1] = startstateinfo->h;
#endif
    space->heap->insertheap(startstateinfo, key);
#endif  // end search forward


    return;
}


void
AAPlanner::AA_ReInitializeState(AAState *state, AASpace *space) {
}


//Xiaoxun 13: this is to get the action cost of the robot to move 1 step along the path
int
AAPlanner::GetMoveCost(AAState *state, AASpace *space) {
    int actioncost;


#ifdef REVERSE
    SBPL_DEBUG("state->pred_num  == %d \n", state->pred_num);

    for (int sind = 0; sind < state->pred_num; sind++) {
        CMDPSTATE *MDPstate = GetState(state->PredsIDV[sind], space); //Xiaoxun 1:  see P13
        AAState *predstate = (AAState *)MDPstate->PlannerSpecificData;

        if (predstate->bestSuccState  == state->MDPstate) {
            actioncost = state->PredsCostV[sind];
            break;
        }
    }// end (sind)
#else
    for(nodeStub nS : state->successors)
        if(nS.id == state->bestSuccState->StateID){
            actioncost = nS.cost;
            break;
        }
#endif

    //SBPL_DEBUG("actioncost == %d \n",actioncost);
    return actioncost;
}

//2010.02.28
//Xiaoxun 14: this is to get the action cost of the target to move 1 step randomly
int
AAPlanner::GetTargetMoveCost(AAState *state, AAState *state2, AASpace *space) {
    int actioncost;

    for(nodeStub nS : state->successors)
        if (nS.id == state2->MDPstate->StateID) {
            actioncost = nS.cost;
            break;
        }

    return actioncost;
}

#endif
