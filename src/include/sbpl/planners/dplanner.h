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

#pragma once
// SBPL includes
// -------------

// SBPL includes
// -------------

// Planner
#include <sbpl/planners/planner.h>
// MDP
#include <sbpl/utils/mdp.h>
// #include <sbpl/utils/mdpconfig.h>
#include <sbpl/discrete_space_information/environment.h>
// Heap
// #include <sbpl/utils/key.h>
#include <sbpl/utils/heap.h>
#include <sbpl/utils/list.h>



#define DMDP_STATEID2IND STATEID2IND_SLOT0

#define D_INCONS_LIST_ID 0
// class CMDP;
// class CMDPSTATE;
// class CMDPACTION;
// class CHeap;
// class CList;


// TODO: use typedef
#include<queue>
#define dq priority_queue<DState*, vector<DState*>, int (*)(DState*, DState*)>
//-------------------------------------------------------------


/** \brief D* Node */
typedef class DSEARCHSTATEDATA : public AbstractSearchState {
public:
  /** \brief the MDP state itself (pointer to the graph represented as MDPstates)
    */
  CMDPSTATE* MDPstate;
  /** \brief D* relevant data
    */
  unsigned int v;
  /** \brief D* relevant data
    */
  unsigned int g;
  /** \brief D* relevant data
    */
  unsigned int gp_iter;
  unsigned int gp_iter1;
  unsigned int g_p;
  short unsigned int iterationclosed;
  /** \brief D* relevant data
    */
  short unsigned int callnumberaccessed;
  short unsigned int numofexpands;
  /** \brief best predecessor and the action from it, used only in forward searches
    */
  CMDPSTATE *bestpredstate;
  DSEARCHSTATEDATA* father;
  /** \brief the next state if executing best action
    */
  CMDPSTATE  *bestnextstate;
  vector<CMDPSTATE *> *path;
  unsigned int costtobestnextstate;
  int h;


public:
  DSEARCHSTATEDATA() {};
  ~DSEARCHSTATEDATA() {
    if (path) {
      path->clear();
      delete path;
    }
    delete path;
  };
} DState;



/** \brief statespace of D* */
typedef struct DSEARCHSTATESPACE {
  double eps;
  double eps_satisfied;
  CHeap* heap;
  CList* inconslist;
  short unsigned int searchiteration;
  short unsigned int callnumber;
  CMDPSTATE* searchgoalstate;
  CMDPSTATE* searchstartstate;

  CMDP searchMDP;

  bool bReevaluatefvals;
  bool bReinitializeSearchStateSpace;
  bool bRebuildOpenList;

} DSearchStateSpace_t;



/** \brief Anytime D* search planner: all states are uniquely defined by stateIDs */
class DPlanner : public SBPLPlanner {

public:
  /** \brief replan a path within the allocated time, return the solution in the vector, also returns solution cost */
  int replan(double allocated_time_secs, vector<StateID>* solution_stateIDs_V, int* solcost);

  /** \brief set the goal state */
  int set_goal(StateID goal_stateID);

  /** \brief set the start state */
  int set_start(StateID start_stateID);

  /** \brief set a flag to get rid of the previous search efforts, release the memory and re-initialize the search, when the next replan is called */
  void set_eps (float eps) {_eps = eps;}
  int force_planning_from_scratch();

  /** \brief you can either search forwards or backwards */
  int set_search_mode(bool bSearchUntilFirstSolution);

  /** \brief inform the search about the new edge costs */
  void costs_changed(StateChangeQuery const & stateChange);

  /** \brief direct form of informing the search about the new edge costs
      \param succsIDV array of successors of changed edges
      \note this is used when the search is run forwards
  */
  void update_succs_of_changededges(vector<StateID>* succsIDV);

  /** \brief direct form of informing the search about the new edge costs
      \param predsIDV array of predecessors of changed edges
      \note this is used when the search is run backwards
  */
  void update_preds_of_changededges(vector<StateID>* predsIDV);

  /** \brief returns the number of states expanded so far */
  virtual int get_n_expands() const {
    return searchexpands;
  }

  int get_n_expands_init_solution() {
    return num_of_expands_initial_solution;
  };

  bool IsExpanded (StateID stateID);

  /** \brief constructor */
  DPlanner(DiscreteSpaceInformation* environment, bool bForwardSearch);

  /** \brief destructor */
  ~DPlanner();



private:

  //member variables
  double _eps, finitial_eps, finitial_eps_planning_time, final_eps_planning_time, final_eps;
  int num_of_expands_initial_solution;
  MDPConfig* MDPCfg_;

  bool bforwardsearch;
  bool bsearchuntilfirstsolution; //if true, then search until first solution (see planner.h for search modes)

  DSearchStateSpace_t* pSearchStateSpace_;

  unsigned int searchexpands;
  int MaxMemoryCounter;
  clock_t TimeStarted;
  FILE *fDeb;


  //member functions
  void Initialize_searchinfo(CMDPSTATE* state, DSearchStateSpace_t* pSearchStateSpace);

  CMDPSTATE* CreateState(int stateID, DSearchStateSpace_t* pSearchStateSpace);

  CMDPSTATE* GetState(int stateID, DSearchStateSpace_t* pSearchStateSpace);

  int ComputeHeuristic(CMDPSTATE* MDPstate, DSearchStateSpace_t* pSearchStateSpace);

  //initialization of a state
  void InitializeSearchStateInfo(DState* state, DSearchStateSpace_t* pSearchStateSpace);

  //re-initialization of a state
  void ReInitializeSearchStateInfo(DState* state, DSearchStateSpace_t* pSearchStateSpace);

  void DeleteSearchStateData(DState* state);


  //used for backward search
  void UpdatePredsofOverconsState(DState* state, DSearchStateSpace_t* pSearchStateSpace);
  void MarkPredsofOverconsState(DState* state, DSearchStateSpace_t* pSearchStateSpace);
  void UpdatePredsofUnderconsState(DState* state, DSearchStateSpace_t* pSearchStateSpace);
  void UpdatePredsofUnderconsStateDelayed(DState* state, DSearchStateSpace_t* pSearchStateSpace);
  void print_D(DSearchStateSpace_t* pSearchStateSpace, DState* s);

  //used for forward search
  void UpdateSuccsofOverconsState(DState* state, DSearchStateSpace_t* pSearchStateSpace);
  void UpdateSuccsofUnderconsState(DState* state, DSearchStateSpace_t* pSearchStateSpace);

  void UpdateSetMembership(DState* state);
  void UpdateSetMembershipOpt(DState* state);
  void UpdateSetMembershipUnder(DState* state);
  void UpdateSetMembershipLower(DState* state);
  void UpdateSetMembershipEps(DState* state);
  void Recomputegval(DState* state);
  void Recomputegprime(DSearchStateSpace_t* pSearchStateSpace, DState* state);
  bool RecomputegprimeStart (DSearchStateSpace_t* pSearchStateSpace, int bound);
  void RecomputegprimeTest(DSearchStateSpace_t* pSearchStateSpace, DState* state);
  void RecomputegvalDummy(DState* state);
  void RecomputegvalRecurring(DState* state);


  int GetGVal(int StateID, DSearchStateSpace_t* pSearchStateSpace);

  //returns 1 if the solution is found, 0 if the solution does not exist and 2 if it ran out of time
  int ComputePath(DSearchStateSpace_t* pSearchStateSpace, double MaxNumofSecs);
  int ComputePathAStar(DSearchStateSpace_t* pSearchStateSpace, double MaxNumofSecs);
  int ComputePathEps(DSearchStateSpace_t* pSearchStateSpace, double MaxNumofSecs);
  int ComputePathEpsOpt(DSearchStateSpace_t* pSearchStateSpace, double MaxNumofSecs);
  int ComputePathEpsSave(DSearchStateSpace_t* pSearchStateSpace, double MaxNumofSecs);
  int ComputePathDelayed(DSearchStateSpace_t* pSearchStateSpace, double MaxNumofSecs);
  void SuspendPath (DSearchStateSpace_t* pSearchStateSpace, DState* s);

  int ComputePathDijkstra (DSearchStateSpace_t* pSearchStateSpace, DState* s, double MaxNumofSecs);
  int GenerateChildren (DSearchStateSpace_t* pSearchStateSpace, DState* s,dq *open);
  void BuildNewOPENList(DSearchStateSpace_t* pSearchStateSpace);
  void PrintOpenList(DSearchStateSpace_t* pSearchStateSpace);

  void Reevaluatefvals(DSearchStateSpace_t* pSearchStateSpace);

  //creates (allocates memory) search state space
  //does not initialize search statespace
  int CreateSearchStateSpace(DSearchStateSpace_t* pSearchStateSpace);

  //deallocates memory used by SearchStateSpace
  void DeleteSearchStateSpace(DSearchStateSpace_t* pSearchStateSpace);


  //reset properly search state space
  //needs to be done before deleting states
  int ResetSearchStateSpace(DSearchStateSpace_t* pSearchStateSpace);

  //initialization before each search
  void ReInitializeSearchStateSpace(DSearchStateSpace_t* pSearchStateSpace);

  //very first initialization
  int InitializeSearchStateSpace(DSearchStateSpace_t* pSearchStateSpace);

  int SetSearchGoalState(StateID SearchGoalStateID, DSearchStateSpace_t* pSearchStateSpace);


  int SetSearchStartState(StateID SearchStartStateID, DSearchStateSpace_t* pSearchStateSpace);

  //reconstruct path functions are only relevant for forward search
  int ReconstructPath(DSearchStateSpace_t* pSearchStateSpace);


  void PrintSearchState(DState* searchstateinfo, FILE* fOut);
  void PrintSearchPath(DSearchStateSpace_t* pSearchStateSpace, FILE* fOut);

  int getHeurValue(DSearchStateSpace_t* pSearchStateSpace, int StateID);

  //get path
  vector<StateID> GetSearchPath(DSearchStateSpace_t* pSearchStateSpace, int& solcost);
  vector<StateID> GetSearchPathwithSus(DSearchStateSpace_t* pSearchStateSpace, int& solcost);
  void FollowPath (DSearchStateSpace_t* pSearchStateSpace, DState* s, int& solcost, vector <StateID>& Ids);
  bool FindOrRaisePath(DSearchStateSpace_t* pSearchStateSpace); //, int& solcost);


  bool Search(DSearchStateSpace_t* pSearchStateSpace, vector<StateID>& pathIds, int & PathCost, bool bFirstSolution, bool bOptimalSolution, double MaxNumofSecs);
  bool SearchDelayed(DSearchStateSpace_t* pSearchStateSpace, vector<int>& pathIds, int & PathCost, bool bFirstSolution, bool bOptimalSolution, double MaxNumofSecs);

  CKey ComputeKey(DState* state);

  void Update_SearchSuccs_of_ChangedEdges(vector<StateID> const* statesIDV);

};




/**
   \brief See comments in sbpl/src/planners/planner.h about the what and why
   of this class.
*/
/*class StateChangeQuery {
public:
  virtual ~StateChangeQuery() {}
  virtual std::vector<int> const * getPredecessors() const = 0;
  virtual std::vector<int> const * getSuccessors() const = 0;
};
*/
