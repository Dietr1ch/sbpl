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


// Configuration
// -------------
// Nodes consult their neighborhood once (+StateChangeQuery should give a boost)
#define ASTAR_SAVE_NEIGHBOURS 1





/**
 * A* Node
 *
 * A* node contains the information of the search STATE (!=node)
 *  along information needed for the planner to search
 */
class ASTARNode : public AbstractSearchState {

// Data
// ====

private:
public:
    // MDP State
    // ---------
    CMDPSTATE *MDPstate;
    StateID id() const;

    // A* data
    // -------
    uint _g;
    uint _h;
    uint g() const;
    uint h() const;
    /** g+h (Non-cached!!)*/
    uint f() const;

    // Adaptive A* data
    // ----------------
    SearchID lastUpdated;


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
#if ASTAR_SAVE_NEIGHBOURS
    union{
        vector<nodeStub> *successors;
        vector<nodeStub> *predecessors;
    };
#endif


    // Statistics
    // ----------
#if STATISTICS || 1
    struct {
        uint expansions = 0;
        union {
            int successors;
            int predecessors;
        } generated;
    } stats;
#endif
    unsigned int generated_iteration;


// Functions
// =========
public:
    ASTARNode(StateID id, ASTARSpace* space);
    ASTARNode(StateID id, ASTARSpace* space, int initialH);
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

    // Iterative Planner information
    // -----------------------------

    // Current search Iteration
    SearchID iteration;  // Search iteration identifier

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
            vector<uint> expansions;
            vector<uint> pathLength;
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
    void setStart(StateID startID);
    void setGoal(StateID goalID);



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
    ASTARNode* getNode(StateID id);
    /** Get a node without updating h (having h=0 if it's new) */
    ASTARNode* getNode0(StateID id);

    /** Get an updated node */
    ASTARNode* getNode(CMDPSTATE *mdpState);

    /** UNSAFE: Get an updated node */
    ASTARNode* getNode_(StateID id);

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
    bool Search(Path *pathIDs, int *cost, bool bFirstSolution, bool bOptimalSolution, double givenSeconds);

    typedef vector<nodeStub> (*neighborhoodFunction) (StateID id);
    typedef void                 (*reachingFunction) (const ASTARNode &node, const nodeStub &nS);

    /** Expand Node backwards
     * An expansion (backward) 'reaches' all the neighbors
     */
    inline void updatePredecessors(ASTARNode &node);
    /** Reach Node backwards
     * A (backwards) 'reach' reviews if a Node can reach the Node specified by
     *   the neighborStub in a better way than before, recording the improvement
     *   (if any)
     */
    inline void  reachPredecessor (const ASTARNode &node, const nodeStub &nS);
    /** Expand Node forward
     * An expansion (forward) 'reaches' all the neighbors
     */

    inline void updateSuccessors(ASTARNode &node);
    /** Reach Node forward
     * A (forward) 'reach' reviews if a Node can reach the Node specified by
     *   the neighborStub in a better way than before, recording the improvement
     *   (if any)
     */
    inline void  reachSuccessor (const ASTARNode &node, const nodeStub &nS);


    // TODO: implement these functions if needed
    //void reevaluateFVals();
    //int reconstructPath();
    //void printPath(FILE *fOut=stdout);
    //Path getSearchPath(int &solcost);



    // SBPL API
    // ========

public:
    /** Replan */
    int replan(double givenTime, Path* pathIDs, int* cost);
    /** Replan (unimplemented) */
    int replan(Path* pathIDs, ReplanParams params, int* cost);

    int set_start(StateID startID);
    int set_goal(StateID goalID);

    int force_planning_from_scratch();
    int set_search_mode(bool bSearchUntilFirstSolution);
    void costs_changed(StateChangeQuery const & stateChange);
};
