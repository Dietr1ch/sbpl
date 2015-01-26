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
#ifndef __AAPLANNER_H_
#define __AAPLANNER_H_

using namespace std;
#include <sbpl/planners/planner.h>

#include <time.h>

// SBPL includes
#include <sbpl/utils/mdp.h>
#include <sbpl/utils/mdpconfig.h>
#include <sbpl/discrete_space_information/environment.h>



class DiscreteSpaceInformation;


// Configuration
// =============

//control of EPS
#define AA_DEFAULT_INITIAL_EPS      5.0
#define AA_DECREASE_EPS    0.2
#define AA_FINAL_EPS        1.0
#define MAX_SEARCHES 100000

//---------------------

#define AA_INCONS_LIST_ID 0
//Xiaoxun added this
#define AA_REMOVE_LIST_ID 0

#define AA_MAGIC 'D'



class CMDP;
class CMDPSTATE;
class CMDPACTION;
class CHeap;
class CList;
//Xiaoxun added this
class RList;


class AAState;
class AASpace;
class AAPlanner;

/**
 * Adaptive A* Search State
 */
class AAState : public AbstractSearchState {

public:
    // MDP State
    // ---------
    CMDPSTATE *MDPstate;

    // Adaptive A* data
    // ----------------
    unsigned int g;
    unsigned int v;
    short unsigned int iterationclosed;
    short unsigned int callnumberaccessed;

    // Best next state
    // ---------------
    // REVIEW: use annonymous union {succs, preds}
    CMDPSTATE *bestSuccState;  // backward searches
    CMDPSTATE *bestpredstate;  //  forward searches
    unsigned int costtobestnextstate;
    int h;

    // Neighbourhood
    // -------------
    // REVIEW: use annonymous union {succs, preds}
    vector<nodeStub> successors;
    vector<nodeStub> predecessors;


    // Statistics
    // ----------
#if STATISTICS || 1
    short unsigned int expansions = 0;
    int succs_ever_generated;
    int preds_ever_generated;
#endif
    unsigned int generated_iteration;

    char _magic;

public:
    AAState(AAPlanner *planner, AASpace* space);
    ~AAState();

    void reinitialize(AAPlanner *planner, AASpace* space);
protected:
    inline void resetSearchInfo();
    inline void resetStatistics();
    inline void computeH(AAPlanner* planner, AASpace* space);
};



/**
 * Adaptive A* Search Space
 */
class AASpace {
public:
    double eps;
    double eps_satisfied;
    CHeap *heap;
    CList *inconslist;

    unsigned int searchiteration;             //Xiaoxun added this (2)
    unsigned long int totalsearchiteration;       //Xiaoxun added this (3)
    unsigned long int totalmovecost;


    unsigned long int hunter_movecost_testcase;  //Xiaoxun added this (13)
    unsigned long int target_movecost_testcase;  //Xiaoxun added this (14)

    short unsigned int callnumber;
    CMDPSTATE *searchgoalstate;
    CMDPSTATE *searchstartstate;
    CMDPSTATE *previous_search_start_state;   //Xiaoxun added this (3)
    CMDPSTATE *temp_start_state;              //Xiaoxun added this (4.1)
    CMDPSTATE *temp_goal_state;               //Xiaoxun added this (4.2)

    CMDP searchMDP;

    bool bReevaluatefvals;
    bool bReinitializeSearchStateSpace;
    bool bNewSearchIteration;

    unsigned int last_search_goal_key; // Xiaoxun added this (5)
    unsigned int new_search_goal_key;  // Xiaoxun added this (6)

    unsigned int robot_steps;         // Xiaoxun added this (7)
    unsigned int target_steps;        // Xiaoxun added this (8)
    unsigned int total_robot_steps;       // Xiaoxun added this (9)


    unsigned int expansion_of_the_search[MAX_SEARCHES];   // Xiaoxun added this (10)
    unsigned long int pathlength[MAX_SEARCHES];   // Xiaoxun added this (11)
    unsigned long int keymodifer;   // Xiaoxun added this (12)
    unsigned long int keymod[MAX_SEARCHES];

    AASpace();
    ~AASpace();
};



/**
 * Adaptive A* Planner
 */
class AAPlanner : public SBPLPlanner {

    public:
        int replan(double givenSeconds, vector<stateID> *solution_stateIDs_V);
        int replan(double givenSeconds, vector<stateID> *solution_stateIDs_V, int *solcost);

        int set_start(stateID startID);
        int set_goal(stateID goalID);

        void costs_changed(StateChangeQuery const &stateChange);
        void costs_changed();
        int force_planning_from_scratch();
        int set_search_mode(bool bSearchUntilFirstSolution);

        virtual double get_solution_eps() const {
            return pSearchSpace->eps_satisfied;
        };
        virtual int get_n_expands() const { return searchexpands; }
        virtual void set_initialsolution_eps(double initialsolution_eps) {finitial_eps = initialsolution_eps;};

        void print_searchpath(FILE *fOut);


        //constructors & destructors
        AAPlanner(DiscreteSpaceInformation *environment, bool bforwardsearch);
        ~AAPlanner();



    private:

        //member variables
        double finitial_eps;

        bool bforwardsearch; //if true, then search proceeds forward, otherwise backward

        bool bsearchuntilfirstsolution; //if true, then search until first solution only (see planner.h for search modes)

        AASpace *pSearchSpace;

        unsigned int searchexpands;
        int MaxMemoryCounter;
        clock_t TimeStarted;
//Xiaoxun added:
        clock_t Case_TimeStarted;
        clock_t Case_TimeEnded;

        FILE *fDeb;


        CMDPSTATE *CreateState(int stateID, AASpace *space);

        CMDPSTATE *GetState(int stateID, AASpace *space);

public:
        int ComputeHeuristic(CMDPSTATE *MDPstate, AASpace* space);
private:

        //initialization of a state
        void InitializeSearchStateInfo(AAState *state, AASpace *space);

        //re-initialization of a state
        void ReInitializeSearchStateInfo(AAState *state, AASpace *space);

        void DeleteSearchStateData(AAState *state);

        //used for backward search
        void UpdatePreds(AAState *state, AASpace *space);


        //used for forward search
        void UpdateSuccs(AAState *state, AASpace *space);

        int GetGVal(int StateID, AASpace *space);

        //returns 1 if the solution is found, 0 if the solution does not exist and 2 if it ran out of time
        int ImprovePath(AASpace *space, double MaxNumofSecs);

        void BuildNewOPENList(AASpace *space);

        void Reevaluatefvals(AASpace *space);

        //creates (allocates memory) search state space
        //does not initialize search statespace
        int CreateSearchStateSpace(AASpace *space);

        //deallocates memory used by SearchStateSpace
        void DeleteSearchStateSpace(AASpace *space);

        //debugging
        void PrintSearchState(AAState *state, FILE *fOut);


        //reset properly search state space
        //needs to be done before deleting states
        int ResetSearchStateSpace(AASpace *space);

        //initialization before each search
        void ReInitializeSearchStateSpace(AASpace *space);

        //very first initialization
        int InitializeSearchStateSpace(AASpace *space);

        int SetSearchGoalState(int SearchGoalStateID, AASpace *space);


        int SetSearchStartState(int SearchStartStateID, AASpace *space);

        //reconstruct path functions are only relevant for forward search
        int ReconstructPath(AASpace *space);


        void PrintSearchPath(AASpace *space, FILE *fOut);

        int getHeurValue(AASpace *space, int StateID);

        //get path
        vector<stateID> GetSearchPath(AASpace *space, int &solcost);

        bool Search(AASpace *space, vector<stateID> &pathIds, int &PathCost, bool bFirstSolution, bool bOptimalSolution, double MaxNumofSecs);





//Xiaoxun's AA* functions
//------------------------------------------------------------------------------
        AAState *ChooseOneSuccs_for_TargetMove(AAState *state, AASpace *space);

        int ReSetSearchGoalState(int SearchGoalStateID, AASpace *space);
        int reset_goal(int start_stateID);
// return 1 if target moved off the previous path, return 0 otherwise
        int goal_moved(AASpace *space);
        int ReSetSearchStartState(int SearchStartStateID, AASpace *space);
        int reset_start(int start_stateID);
// return 1 if robot moved one step
        int start_moved(AASpace *space);


        void Reset_New_Test_Case(AASpace *space);
        void ReInitializeNewSearch(AASpace *space);

        void AA_ReInitializeState(AAState *state, AASpace *space);

        int GetMoveCost(AAState *state, AASpace *space);
//2010.02.28
        int GetTargetMoveCost(AAState *state, AAState *state2, AASpace *space);
};


#endif
