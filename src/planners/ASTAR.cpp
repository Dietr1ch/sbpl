#include <assert.h>
#include <time.h>
#include <stdlib.h>  // exit and (standard) failure codes

// SBPL Includes
#include <sbpl/config.h>
#include <sbpl/planners/ASTAR.h>
#include <sbpl/utils/heap.h>




// A* State
// ========
ASTARNode::ASTARNode(stateID id, ASTARSpace *space) : ASTARNode(id, space, 0) {
    h = space->computeHeuristic(*MDPstate);
}
ASTARNode::ASTARNode(stateID id, ASTARSpace *space, int initialH) {
//     SBPL_DEBUG("Creating           A*Node[%p] state[%d]", this, id);

    MDPstate = space->MDP.AddState(id);
    assert(MDPstate);
    assert(MDPstate->StateID);

    // Link state and node data
    MDPstate->PlannerSpecificData = this;

    // Register state
    int statesCount = space->MDP.StateArray.size();
    space->problem->StateID2IndexMapping[id][ASTARMDP_STATEID2IND] = statesCount-1;

    // Setup
    g = INFINITECOST;  // Not reached
    h = initialH;

    best = nullptr;
    bestSuccessor  = nullptr;
    bestPredecesor = nullptr;  // Just to be clear (only one exists)

    heapindex = 0;
    iterationclosed = 0;

    costToBestState = INFINITECOST;

}
ASTARNode::~ASTARNode(){
}

// Shortcuts
inline stateID ASTARNode::id() const {
#ifdef DEBUG
    if(!MDPstate)
        QUIT("Node @ %p is corrupted", (void*)this);
#endif
    assert(MDPstate->StateID);

    return MDPstate->StateID;
}
inline    uint ASTARNode::f()  const { return g+h; }








// A* Space
// ========
ASTARSpace::ASTARSpace(DiscreteSpaceInformation *problemDSI) {
    TRACE("Creating space with problem[%p]", (void*)problemDSI);
    assert(problemDSI);

    open = new CHeap;
    problem = problemDSI;

    startState = nullptr;
    goalState  = nullptr;

    // Asserts to remind whats true now (useless asserts)
    assert(stats.perSearch.expansions.empty());
    assert(stats.perSearch.pathLength.empty());
}

ASTARSpace::~ASTARSpace() {
}

inline
int
ASTARSpace::computeHeuristic(const CMDPSTATE &origin){
    auto target = backwardSearch ? startState
                                 : goalState;
    assert(target);
    return problem->GetFromToHeuristic(origin.StateID, target->StateID);
}

// Configuration
// -------------
inline
void
ASTARSpace::setStart(stateID startID) {
    SBPL_DEBUG("");
    SBPL_DEBUG("Setting start to state[%d]", startID);

    // Getting the A* node is an overhead (here), but will be generated later.
    auto newStartNode  = getNode0(startID);
    SBPL_DEBUG(" start = A*Node[%p]", newStartNode);
    auto newStartState = newStartNode->MDPstate;

    assert(newStartState);
    if(startState != newStartState){
        SBPL_DEBUG("Start state changed");
        startState = newStartState;
        // REVIEW: initialization
        valid = false;
    }
    SBPL_DEBUG("");
}
inline
void
ASTARSpace::setGoal(stateID goalID) {
    SBPL_DEBUG("");
    SBPL_DEBUG("Setting goal  to  state[%d]", goalID);

    // Getting the A* node is an overhead (here), but will be generated later.
    auto newGoalNode  = getNode0(goalID);
    SBPL_DEBUG("  goal = A*Node[%p]", newGoalNode);
    auto newGoalState = newGoalNode->MDPstate;

    assert(newGoalState);
    if(goalState != newGoalState){
        SBPL_DEBUG("Goal  state changed");
        goalState = newGoalState;

        // REVIEW: initialization
        for(CMDPSTATE *state : MDP.StateArray){
            assert(state);
            ASTARNode *node = (ASTARNode*) state->PlannerSpecificData;
            assert(node);
            node->h = 0;
            //TODO: recompute h
        }
    }
    SBPL_DEBUG("");
}




// Open list management
// --------------------
inline
void
ASTARSpace::upsertOpen(ASTARNode* node) {
    CKey k;
    k.key[0] = node->f();
    k.key[1] = node->h;

    assert(node);
    if(node->heapindex)
        updateOpen_(node, k);
    else
        insertOpen_(node, k);
}
inline
void
ASTARSpace::insertOpen_(ASTARNode *node, CKey key) {
    assert(node);
#ifdef DEBUG
    float f = key.key[0]/1000.0;
    float h = key.key[1]/1000.0;
    float g = f-h;
    TRACE("  open: Inserting  A*Node[%p] state[%4d] {%5.2f + %5.2f = %5.2f}",
               (void*)node, node->id(), g, h, f);
#endif
    open->insertheap(node, key);
}
inline
void
ASTARSpace::updateOpen_(ASTARNode* node, CKey key) {
    assert(node);
#ifdef DEBUG
    float f = key.key[0]/1000.0;
    float h = key.key[1]/1000.0;
    float g = f-h;
    TRACE("  open: Updating   A*Node[%p] state[%4d] (%5.2f + %5.2f = %5.2f)",
               (void*)node, node->id(), g, h, f);
#endif
    open->updateheap(node, key);
}
inline
ASTARNode*
ASTARSpace::popOpen() {
    auto node = (ASTARNode*) open->deleteminheap();
//     TRACE("  open: Retrieving A*Node[%p] state[%d]", node, node->id());
    return node;
}
inline
bool
ASTARSpace::openEmpty(){
    return open->emptyheap();
}




// Nodes management
// ----------------
inline
ASTARNode*
ASTARSpace::getNode(stateID id) {
    assert(id);

    if((size_t)id >= problem->StateID2IndexMapping.size())
        QUIT("State %d is invalid\n", id);
    if(problem->StateID2IndexMapping[id][ASTARMDP_STATEID2IND] == -1)
        return new ASTARNode(id, this);

    return getNode_(id);
}
inline
ASTARNode*
ASTARSpace::getNode0(stateID id) {
    assert(id);

    if((size_t)id >= problem->StateID2IndexMapping.size())
        QUIT("State %d is invalid\n", id);
    if(problem->StateID2IndexMapping[id][ASTARMDP_STATEID2IND] == -1)
        return new ASTARNode(id, this, 0);

    return getNode_(id);
}
inline
ASTARNode*
ASTARSpace::getNode_(stateID id) {
    assert(id);

    auto index = problem->StateID2IndexMapping[id][ASTARMDP_STATEID2IND];
    CMDPSTATE *cmdpState = MDP.StateArray[index];
    assert(cmdpState);
    auto node = (ASTARNode*) cmdpState->PlannerSpecificData;

    return node;
}

inline
ASTARNode*
ASTARSpace::getNode(CMDPSTATE* mdpState) {
    assert(mdpState);
    assert(mdpState->StateID);
    auto node = (ASTARNode*) mdpState->PlannerSpecificData;

#ifdef DEBUG
    assert(node);
    auto nodeFromID = getNode(mdpState->StateID);
    if(node!=nodeFromID)
        QUIT("\nCDMPState[%p] data unconsistent data:\n"
               "  state[%d]\n"
               "  cmdp planner data: A*Node[%p]\n"
               "  space id target:   A*Node[%p]",
              (void*)mdpState,
              mdpState->StateID,
              (void*)node,
              (void*)nodeFromID);
#endif

    return node;
}

// Shortcuts
inline ASTARNode* ASTARSpace::getStart()  { return  getNode(startState); }
inline ASTARNode* ASTARSpace::getGoal()   { return  getNode(goalState);  }





// Statistics
// ----------
inline
void
ASTARSpace::resetStatistics(){
}










// A* Planner
// ==========

// Object management
// -----------------
ASTARPlanner::ASTARPlanner(DiscreteSpaceInformation *environment, bool backwardSearch) {
    assert(environment);
    SBPL_DEBUG("Creating planner at %p\n", this);

    space = new ASTARSpace(environment);
    assert(space);

    SBPL_DEBUG("Planner created at %p\n", this);
}
ASTARPlanner::~ASTARPlanner(){
    SBPL_DEBUG("Planner at %p will be destroyed\n", this);
}




// Search
// ------
bool
ASTARPlanner::Search(vector<stateID> *pathIDs, int *cost, bool bFirstSolution, bool bOptimalSolution, double givenSeconds) {
    assert(space);  // redundant
    assert(space->startState);
    assert(space->goalState);
    assert(pathIDs->empty());


    // Initialization
    // ==============
    SBPL_DEBUG("\nPreparing A* loop");

    // Reset statistics (just in case)
    stats.iteration = 0;
    stats.expansions = 0;


    // Prepare goal node
    ASTARNode* goalNode = space->getGoal();
    SBPL_DEBUG("Looking for   A*Node[%p]", (void*)goalNode);
    goalNode->g = INFINITECOST;
    goalNode->h = 0;

    // Prepare starting node
    ASTARNode *startNode = space->getStart();
    SBPL_DEBUG("Starting from A*Node[%p]", (void*)startNode);
    startNode->g = 0;
    startNode->h = space->computeHeuristic(*startNode->MDPstate);

    space->insertOpen(startNode);


    // Expansion
    // =========
    SBPL_DEBUG("\nStarting A* loop");

    CKey bestKey;
    while(!space->openEmpty()){
        // Get the best node
        bestKey = space->open->getminkeyheap();
        ASTARNode *currentNode = space->popOpen();

        // Check the search should go on
        if(bestKey.key[0] >= INFINITECOST)         // Best node is not reachable
            break;

        // Check if the (current) best node can't improve the path
        if(goalNode->g < bestKey.key[0])
            break;

        // Expand the best node
        printf("\n\n");
        TRACE("Expanding state[%d]", currentNode->id());
        updateSuccessors(*currentNode);
    }
    SBPL_DEBUG("\nA*: Finished");


    // Open list information
    // =====================
    if(space->openEmpty())
        SBPL_DEBUG("A*: No more nodes to reach");

    if(bestKey.key[0] >= INFINITECOST)
        SBPL_DEBUG("A*: BestKey is infinite, nodes on open are not reachable");
    else
        SBPL_DEBUG("A*: BestKey (%5.2f, %5.2f)",
                   bestKey.key[0]/1000.0, bestKey.key[1]/1000.0);

    SBPL_DEBUG("Goal g: %5.2f", goalNode->g);


    return SBPL_OK;
}

inline
void
ASTARPlanner::updateSuccessors(const ASTARNode &node) {
    assert(space);
    assert(space->problem);

    stateID id = node.id();
    assert(id);

    // Get successors from environment
    vector<nodeStub> *successors = space->problem->GetSuccs(id);
    assert(successors);

    // Reach each node updating path&cost on improvements
    for(nodeStub nS : *successors)
        reachSuccessor(node, nS);
}

inline
void
ASTARPlanner::reachSuccessor(const ASTARNode &node, const nodeStub &nS) {
    assert(nS.id);
    assert(nS.cost>=0);
    auto *neigh = space->getNode(nS.id);

    // Reached from node and and action that costs nS.cost
    auto newG = node.g + nS.cost;
    assert(newG>=0);

    if (newG < neigh->g) {                            // A better path was found
//         SBPL_DEBUG("Path to state[%4d]     was enhanced (%5.2f > %5.2f)",
//                    nS.id, neigh->g/1000.0, newG/1000.0);

        // Use new path
        neigh->g = newG;  //Save new cost
        neigh->bestPredecesor = node.MDPstate;

        // Update this node on the open list
        space->upsertOpen(neigh);
    }
//     else
//         SBPL_DEBUG("Path to state[%4d] not was enhanced (%5.2f < %5.2f)",
//                    nS.id, neigh->g/1000.0, newG/1000.0);
}

inline
void
ASTARPlanner::updatePredecessors(const ASTARNode &node) {
    QUIT("%s", "Not yet implemented");
}

inline
void
ASTARPlanner::reachPredecessor(const ASTARNode &node, const nodeStub &nS) {
    QUIT("%s", "Not yet implemented");
}





// SBPL Planning API
// -----------------
inline
int
ASTARPlanner::replan(double givenSeconds, vector<int> *pathIDs) {
    assert(space);
    assert(space->problem);
    int cost;  // Drop the cost value
    return replan(givenSeconds, pathIDs, &cost);
}
int
ASTARPlanner::replan(double givenSeconds, vector<stateID> *pathIDs, int *cost) {

    // Response pointers should be ready
    assert(cost);
    assert(pathIDs);
    assert(pathIDs->empty());

    // Start timing
    time.start = clock();

    bool solutionFound;
    solutionFound = Search(pathIDs, cost, firstSolutionOnly, true, givenSeconds);

    if (!solutionFound) {
        TRACE("A* failed to find a solution (on %fs)...", givenSeconds/1000.0);
        return solutionFound;
    }

    // End timing
    time.end = clock();
    time.total = time.end - time.start;

    return SBPL_OK;
}
int
ASTARPlanner::replan(vector<int> *pathIDs, ReplanParams params, int *cost) {
    return SBPL_OK;
}


// SBPL Configuration API
// ----------------------
inline
int
ASTARPlanner::set_start(stateID startID) {
    space->setStart(startID);
    return SBPL_OK;
}
inline
int
ASTARPlanner::set_goal(stateID goalID) {
    space->setGoal(goalID);
    return SBPL_OK;
}

inline
int
ASTARPlanner::set_search_mode(bool firstSolutionOnly) {
    this->firstSolutionOnly = firstSolutionOnly;
    return SBPL_OK;
}
inline
int
ASTARPlanner::force_planning_from_scratch() {
    return SBPL_OK;
}


// SBPL Update API
// ---------------
inline
void
ASTARPlanner::costs_changed(StateChangeQuery const & stateChange) {
}
