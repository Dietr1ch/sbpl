// TODO: add updated license
#pragma once


using namespace std;

#include <time.h>

// SBPL includes
// -------------

// Planner
#include <sbpl/planners/planner.h>
// MDP
#include <sbpl/utils/mdp.h>
#include <sbpl/utils/mdpconfig.h>
#include <sbpl/discrete_space_information/environment.h>
// Heap
#include <sbpl/utils/key.h>
#include <sbpl/utils/heap.h>



// Forward Declarations
// --------------------
class ASTARSpace;
class ASTARPlanner;







/**
 * Adaptive A* Search State
 */
class ASTARNode : public AbstractSearchState {

// Data
// ====
public:
    // MDP State
    // ---------
    CMDPSTATE *MDPstate;
    stateID id() const;

    // A* data
    // -------
    uint g;
    uint h;
    /** g+h (Non-cached!!)*/
    uint f() const;

    short unsigned int iterationclosed;
    short unsigned int callnumberaccessed;

    // Best next state
    // ---------------
    union{
        CMDPSTATE *bestSuccessor;  // backward searches
        CMDPSTATE *bestPredecesor; //  forward searches
        CMDPSTATE *best;
    };
    unsigned int costToBestState;

    // Neighbourhood
    // -------------
    union{
        vector<nodeStub> successors;
        vector<nodeStub> predecessors;
    };


    // Statistics
    // ----------
#if STATISTICS || 1
    struct {
        uint expansions = 0;
        union{
            int successors;
            int predecessors;
        } generated;
    } stats;
#endif
    unsigned int generated_iteration;


// Functions
// =========
public:
    ASTARNode(stateID id, ASTARSpace* space);
    ~ASTARNode();

    void print(FILE *fOut);

protected:
    inline void resetSearchInfo();
    inline void resetStatistics();
    inline void computeH(ASTARPlanner* planner, ASTARSpace* space);
};







/**
 * Adaptive A* Search Space
 */
class ASTARSpace {
public:
    CMDP MDP;
    DiscreteSpaceInformation *problem;

    CHeap *heap;
    CMDPSTATE *startState;
    CMDPSTATE *currentState;
    CMDPSTATE *goalState;

    bool valid = true;

    uint searchiteration;             //Xiaoxun added this (2)
    unsigned long int totalsearchiteration;       //Xiaoxun added this (3)
    unsigned long int totalmovecost;

    short unsigned int callnumber;

    bool bReevaluatefvals;
    bool bReinitializeSearchStateSpace;
    bool bNewSearchIteration;

    /** Statistics */
    struct{
        /** Per-search-episode stuff */
        struct{
            /** Expansions made on this space */
            std::vector<uint> expansions;
            std::vector<uint> pathLength;
        } perSearch;
    } stats;



    // Object management
    // =================

    ASTARSpace(DiscreteSpaceInformation *problem);
    ~ASTARSpace();



    // Configuration
    // =============
    void setStart(stateID startID);
    void setGoal(stateID goalID);



    // Open list management
    // ====================

    /** Inserts or Updates (Upserts) a node in the open list */
    void upsertOpen(ASTARNode *node, CKey key);

    /** Inserts or Updates (Upserts) a node in the open list */
    void insertOpen(ASTARNode *node, CKey key) { upsertOpen(node, key); }
    /** Inserts or Updates (Upserts) a node in the open list */
    void updateOpen(ASTARNode *node, CKey key) { upsertOpen(node, key); }


    /** UNSAFE: inserts a node in open, NQA */
    void insertOpen_(ASTARNode *node, CKey key);
    /** UNSAFE: updates a node in open, NQA */
    void updateOpen_(ASTARNode *node, CKey key);

    /** gets the best node on the open 'list' */
    ASTARNode* popOpen();

    /** Checks wheter the open list is empty */
    bool openEmpty();



    // States acquisition
    // ==================

    /** Get an state */
    CMDPSTATE* getState(stateID id);
    /** Get a node */
    ASTARNode* getNode(stateID id);
    /** Get a node */
    ASTARNode* getNode(CMDPSTATE *mdpState);

    /** UNSAFE: Get an state */
    CMDPSTATE* getState_(stateID id);
    /** UNSAFE: Get an state */
    ASTARNode* getNode_(stateID id);

    /** Gets the starting state */
    ASTARNode* getStart();
    /** Gets the goal state */
    ASTARNode* getGoal();



    // Statistics
    // ==========
    void resetStatistics();

};







/**
 * Adaptive A* Planner
 */
class ASTARPlanner : public SBPLPlanner {

protected:
    ASTARSpace *space;

    // Configuration
    bool backwardSearch;
    bool firstSolutionOnly;

    /** Statistics */
    struct {
        uint expansions;
        uint iteration;
    } stats;

    /** Timing Data */
    struct{
        clock_t start;
        clock_t end;
        clock_t total;
    } time;

    FILE *dbg;


public:
    ASTARPlanner(DiscreteSpaceInformation *environment, bool backwardSearch=false);
    ~ASTARPlanner();

    inline int computeHeuristic(const CMDPSTATE &MDPstate);


protected:
    //used for backward search
    inline void updatePredecessors(const ASTARNode &node);
    inline void  reachPredecessor(const ASTARNode &node, const nodeStub &nS);
    //used for forward search
    inline void updateSuccessors(const ASTARNode &node);
    inline void  reachSuccessor(const ASTARNode &node, const nodeStub &nS);

    int GetGVal(int StateID);

    //1:Solution found. 0:No solution exists. 2: Out of time
    int ImprovePath(double MaxNumofSecs);

    void BuildNewOPENList();

    void Reevaluatefvals();

    int ReconstructPath();
    void printPath(FILE *fOut=stdout);

    int getHeurValue(int StateID);

    //get path
    vector<int> GetSearchPath(int &solcost);

    bool Search(vector<stateID> *pathIDs, int *cost, bool bFirstSolution, bool bOptimalSolution, double givenSeconds);



    // SBPL API
    // ========

public:
    /** Replan and drop the solution cost */
    int replan(double givenTime, vector<stateID>* pathIDs);
    /** Replan */
    int replan(double givenTime, vector<stateID>* pathIDs, int* cost);
    /** Replan (unimplemented) */
    int replan(vector<int>* pathIDs, ReplanParams params, int* cost);

    int set_start(int startID);
    int set_goal(int goalID);

    int force_planning_from_scratch();
    int set_search_mode(bool bSearchUntilFirstSolution);
    void costs_changed(StateChangeQuery const & stateChange);
};
