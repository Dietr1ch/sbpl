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
 * A* Node
 *
 * A* node contains the information of the search STATE (!=node)
 *  along information needed for the planner to search
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

    // Adaptive A* data
    // ----------------
    searchID lastUpdated;


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
        vector<nodeStub> *successors;
        vector<nodeStub> *predecessors;
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
    ASTARNode(stateID id, ASTARSpace* space, int initialH);
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
    // Problem data
    CMDP MDP;
    DiscreteSpaceInformation *problem;

    // Instance data
    CMDPSTATE *startState;
    CMDPSTATE *currentState;
    CMDPSTATE *goalState;

    // AdaptiveA* data
    searchID iteration;  // Search iteration identifier

    // Data
    CHeap *open;

    // Configuration
    bool backwardSearch;  // It should be on the planner, but eases heuristic calculation


    //unsigned long int totalsearchiteration;       //Xiaoxun added this (3)
    //unsigned long int totalmovecost;

    //bool bReevaluatefvals;
    //bool bReinitializeSearchStateSpace;
    //bool bNewSearchIteration;

    /** Statistics */
    struct{
        /** Per-search-episode stuff */
        struct{
            /** Expansions made on this space */
            vector<searchID> expansions;
            vector<searchID> pathLength;
        } perSearch;
    } stats;



    // Object management
    // =================

    ASTARSpace(DiscreteSpaceInformation *problem);
    ~ASTARSpace();

    /** Ensures Node values are updated up to this search iteration */
    inline void updateNode(ASTARNode &node);
    /** Computes the heuristic value from a State to the goal State */
    inline int computeHeuristic(const CMDPSTATE &origin);


    // Configuration
    // =============
    void setStart(stateID startID);
    void setGoal(stateID goalID);



    // Open list management
    // ====================

    /** Inserts or Updates (Upserts) a node in the open list */
    void upsertOpen(ASTARNode *node);


    /** UNSAFE: inserts a node in open, No questions asked */
    void insertOpen_(ASTARNode *node, CKey key);
    /** UNSAFE: updates a node in open, No questions asked */
    void updateOpen_(ASTARNode *node, CKey key);

    /** gets the best node on the open 'list' */
    ASTARNode* popOpen();

    /** Checks wheter the open list is empty */
    bool openEmpty();



    // Nodes acquisition
    // =================

    /** Get an updated node */
    ASTARNode* getNode(stateID id);
    /** Get a node without updating h (having h=0 if it's new) */
    ASTARNode* getNode0(stateID id);

    /** Get an updated node */
    ASTARNode* getNode(CMDPSTATE *mdpState);

    /** UNSAFE: Get an updated node */
    ASTARNode* getNode_(stateID id);

    /** Gets the updated starting node */
    ASTARNode* getStart();
    /** Gets the updated goal node */
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


protected:
    bool Search(vector<stateID> *pathIDs, int *cost, bool bFirstSolution, bool bOptimalSolution, double givenSeconds);

    typedef vector<nodeStub> (*neighborhoodFunction) (stateID id);
    typedef void                 (*reachingFunction) (const ASTARNode &node, const nodeStub &nS);

    //used for backward search
    inline void updatePredecessors(ASTARNode &node);
    inline void  reachPredecessor (const ASTARNode &node, const nodeStub &nS);
    //used for forward search
    inline void updateSuccessors(ASTARNode &node);
    inline void  reachSuccessor (const ASTARNode &node, const nodeStub &nS);

    //void reevaluateFVals();
    //int reconstructPath();
    //void printPath(FILE *fOut=stdout);
    //vector<int> getSearchPath(int &solcost);



    // SBPL API
    // ========

public:
    /** Replan and drop the solution cost */
    int replan(double givenTime, vector<stateID>* pathIDs);
    /** Replan */
    int replan(double givenTime, vector<stateID>* pathIDs, int* cost);
    /** Replan (unimplemented) */
    int replan(vector<stateID>* pathIDs, ReplanParams params, int* cost);

    int set_start(int startID);
    int set_goal(int goalID);

    int force_planning_from_scratch();
    int set_search_mode(bool bSearchUntilFirstSolution);
    void costs_changed(StateChangeQuery const & stateChange);
};
