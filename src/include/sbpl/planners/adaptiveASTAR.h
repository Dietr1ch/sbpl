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







class CMDP;
class CMDPSTATE;
class CMDPACTION;
class CHeap;
class CList;
//Xiaoxun added this
class RList;


struct neighborStub{
    int id;
    int cost;
};

//state structure
typedef class AASEARCHSTATEDATA : public AbstractSearchState {
    public:
        // MDP State
        CMDPSTATE *MDPstate;

        // Adaptive A* data
        // ----------------
        unsigned int v;
        unsigned int g;
        short unsigned int iterationclosed;
        short unsigned int callnumberaccessed;
        //best predecessor and the action from it, used only in forward searches
        CMDPSTATE *bestpredstate;
        //the next state if executing best action
        CMDPSTATE  *bestnextstate;
        unsigned int costtobestnextstate;
        int h;

        vector<neighborStub> successors;
        vector<neighborStub> predecessors;

#if STATISTICS || 1
        short unsigned int expansions = 0;
        int succs_ever_generated;
        int preds_ever_generated;
#endif

        unsigned int generated_iteration;

    public:
        AASEARCHSTATEDATA() {};
        ~AASEARCHSTATEDATA() {};
} AAState;




//statespace
typedef struct AASEARCHSTATESPACE {
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


} AASearchStateSpace_t;



//AA* planner
class AAPlanner : public SBPLPlanner {

    public:
        int replan(double allocated_time_secs, vector<int> *solution_stateIDs_V);
        int replan(double allocated_time_sec, vector<int> *solution_stateIDs_V, int *solcost);

        int set_goal(int goal_stateID);
        int set_start(int start_stateID);
        void costs_changed(StateChangeQuery const &stateChange);
        void costs_changed();
        int force_planning_from_scratch();
        int set_search_mode(bool bSearchUntilFirstSolution);

        virtual double get_solution_eps() const {return pSearchStateSpace_->eps_satisfied;};
        virtual int get_n_expands() const { return searchexpands; }
        virtual void set_initialsolution_eps(double initialsolution_eps) {finitial_eps = initialsolution_eps;};

        void print_searchpath(FILE *fOut);


        //constructors & destructors
        AAPlanner(DiscreteSpaceInformation *environment, bool bforwardsearch);
        ~AAPlanner();



    private:

        //member variables
        double finitial_eps;
        MDPConfig *MDPCfg_;

        bool bforwardsearch; //if true, then search proceeds forward, otherwise backward

        bool bsearchuntilfirstsolution; //if true, then search until first solution only (see planner.h for search modes)

        AASearchStateSpace_t *pSearchStateSpace_;

        unsigned int searchexpands;
        int MaxMemoryCounter;
        clock_t TimeStarted;
//Xiaoxun added:
        clock_t Case_TimeStarted;
        clock_t Case_TimeEnded;

        FILE *fDeb;


        //member functions
        void Initialize_searchinfo(CMDPSTATE *state, AASearchStateSpace_t *pSearchStateSpace);

        CMDPSTATE *CreateState(int stateID, AASearchStateSpace_t *pSearchStateSpace);

        CMDPSTATE *GetState(int stateID, AASearchStateSpace_t *pSearchStateSpace);

        int ComputeHeuristic(CMDPSTATE *MDPstate, AASearchStateSpace_t *pSearchStateSpace);

        //initialization of a state
        void InitializeSearchStateInfo(AAState *state, AASearchStateSpace_t *pSearchStateSpace);

        //re-initialization of a state
        void ReInitializeSearchStateInfo(AAState *state, AASearchStateSpace_t *pSearchStateSpace);

        void DeleteSearchStateData(AAState *state);

        //used for backward search
        void UpdatePreds(AAState *state, AASearchStateSpace_t *pSearchStateSpace);


        //used for forward search
        void UpdateSuccs(AAState *state, AASearchStateSpace_t *pSearchStateSpace);

        int GetGVal(int StateID, AASearchStateSpace_t *pSearchStateSpace);

        //returns 1 if the solution is found, 0 if the solution does not exist and 2 if it ran out of time
        int ImprovePath(AASearchStateSpace_t *pSearchStateSpace, double MaxNumofSecs);

        void BuildNewOPENList(AASearchStateSpace_t *pSearchStateSpace);

        void Reevaluatefvals(AASearchStateSpace_t *pSearchStateSpace);

        //creates (allocates memory) search state space
        //does not initialize search statespace
        int CreateSearchStateSpace(AASearchStateSpace_t *pSearchStateSpace);

        //deallocates memory used by SearchStateSpace
        void DeleteSearchStateSpace(AASearchStateSpace_t *pSearchStateSpace);

        //debugging
        void PrintSearchState(AAState *state, FILE *fOut);


        //reset properly search state space
        //needs to be done before deleting states
        int ResetSearchStateSpace(AASearchStateSpace_t *pSearchStateSpace);

        //initialization before each search
        void ReInitializeSearchStateSpace(AASearchStateSpace_t *pSearchStateSpace);

        //very first initialization
        int InitializeSearchStateSpace(AASearchStateSpace_t *pSearchStateSpace);

        int SetSearchGoalState(int SearchGoalStateID, AASearchStateSpace_t *pSearchStateSpace);


        int SetSearchStartState(int SearchStartStateID, AASearchStateSpace_t *pSearchStateSpace);

        //reconstruct path functions are only relevant for forward search
        int ReconstructPath(AASearchStateSpace_t *pSearchStateSpace);


        void PrintSearchPath(AASearchStateSpace_t *pSearchStateSpace, FILE *fOut);

        int getHeurValue(AASearchStateSpace_t *pSearchStateSpace, int StateID);

        //get path
        vector<int> GetSearchPath(AASearchStateSpace_t *pSearchStateSpace, int &solcost);

        bool Search(AASearchStateSpace_t *pSearchStateSpace, vector<int> &pathIds, int &PathCost, bool bFirstSolution, bool bOptimalSolution, double MaxNumofSecs);





//Xiaoxun's AA* functions
//------------------------------------------------------------------------------
        AAState *ChooseOneSuccs_for_TargetMove(AAState *state, AASearchStateSpace_t *pSearchStateSpace);

        int ReSetSearchGoalState(int SearchGoalStateID, AASearchStateSpace_t *pSearchStateSpace);
        int reset_goal(int start_stateID);
// return 1 if target moved off the previous path, return 0 otherwise
        int goal_moved(AASearchStateSpace_t *pSearchStateSpace);
        int ReSetSearchStartState(int SearchStartStateID, AASearchStateSpace_t *pSearchStateSpace);
        int reset_start(int start_stateID);
// return 1 if robot moved one step
        int start_moved(AASearchStateSpace_t *pSearchStateSpace);


        void Reset_New_Test_Case(AASearchStateSpace_t *pSearchStateSpace);
        void ReInitializeNewSearch(AASearchStateSpace_t *pSearchStateSpace);

        void AA_ReInitializeState(AAState *state, AASearchStateSpace_t *pSearchStateSpace);

        int GetMoveCost(AAState *state, AASearchStateSpace_t *pSearchStateSpace);
//2010.02.28
        int GetTargetMoveCost(AAState *state, AAState *state2, AASearchStateSpace_t *pSearchStateSpace);


};


#endif
