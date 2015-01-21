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
    lastUpdated = space->iteration;
}
ASTARNode::ASTARNode(stateID id, ASTARSpace *space, int initialH) {
    MDPstate = space->MDP.AddState(id);
    assert(MDPstate);
    assert(MDPstate->StateID);

    // Link state and node data
    MDPstate->PlannerSpecificData = this;

    // Register state
    int statesCount = space->MDP.StateArray.size();
    space->problem->StateID2IndexMapping[id][ASTARMDP_STATEID2IND] = statesCount-1;

    // Setup
    g = INFINITECOST;  // Not reached yet
    h = initialH;

    best = nullptr;  // (==bestSuccessor==bestPredecesor, only one exists)

    heapindex = 0;
    lastUpdated = (searchID)0;
    successors = nullptr;

    costToBestState = INFINITECOST;
}
ASTARNode::~ASTARNode(){
    delete successors;
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

    backwardSearch = false;

    iteration = (searchID) 0;

    startState = nullptr;
    goalState  = nullptr;

    // Asserts to remind whats true now (useless asserts)
    assert(stats.perSearch.expansions.empty());
    assert(stats.perSearch.pathLength.empty());
}

ASTARSpace::~ASTARSpace() {
    SBPL_DEBUG("Space at %p will be destroyed\n", this);

    // Clear search instance (States and Nodes are cleared with the rest)
    // ---------------------
    startState = nullptr;
    goalState  = nullptr;


    // Delete heap (Clearing heap indexes)
    // -----------
    DELETE(open);


    // Delete states and nodes
    // -----------------------
    for(CMDPSTATE *state : MDP.StateArray){
        // Delete Node
        // REVIEW: heap index clear becomes useless (just here?)
        delete (ASTARNode*) state->PlannerSpecificData;
        state->PlannerSpecificData = nullptr;

        // Delete State
        delete state;
    }
    MDP.StateArray.clear();


    // Delete problem data
    // -------------------
    DELETE(problem);


    SBPL_DEBUG("Space at %p was destroyed\n", this);
}

inline
void
ASTARSpace::updateNode(ASTARNode &node) {
    if(node.lastUpdated!=iteration) {
        node.h = computeHeuristic(*node.MDPstate);
        node.lastUpdated = iteration;
    }
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
    // Getting the A* node is an overhead (here), but will be generated later.
    auto newStartNode  = getNode0(startID);
    SBPL_DEBUG(" start = A*Node[%p]", newStartNode);
    auto newStartState = newStartNode->MDPstate;

    assert(newStartState);
    if(startState != newStartState){
        startState = newStartState;

        // Heuristic remains consistent
    }
}
inline
void
ASTARSpace::setGoal(stateID goalID) {
    // Getting the A* node is an overhead (here), but will be generated later.
    auto newGoalNode  = getNode0(goalID);
    SBPL_DEBUG("  goal = A*Node[%p]", newGoalNode);
    auto newGoalState = newGoalNode->MDPstate;

    assert(newGoalState);
    if(goalState != newGoalState){
        goalState = newGoalState;

        // TODO: Use lazy update
        for(CMDPSTATE *state : MDP.StateArray){
            assert(state);
            ASTARNode *node = (ASTARNode*) state->PlannerSpecificData;
            assert(node);
            // REVIEW: compute h needed for the whole space?
            // REVIEW: use invalidation and re-compute on demand
            node->h = computeHeuristic(*state);
        }
    }
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
    assert(node->lastUpdated==iteration);
#ifdef DEBUG
    //float f = key.key[0]/1000.0;
    //float h = key.key[1]/1000.0;
    //float g = f-h;
    //SBPL_DEBUG("  open: Inserting  A*Node[%p] state[%4d] {%5.2f + %5.2f = %5.2f}",
               //(void*)node, node->id(), g, h, f);
#endif
    open->insertheap(node, key);
}
inline
void
ASTARSpace::updateOpen_(ASTARNode* node, CKey key) {
    assert(node);
    assert(node->lastUpdated==iteration);
#ifdef DEBUG
    //float f = key.key[0]/1000.0;
    //float h = key.key[1]/1000.0;
    //float g = f-h;
    //SBPL_DEBUG("  open: Updating   A*Node[%p] state[%4d] (%5.2f + %5.2f = %5.2f)",
               //(void*)node, node->id(), g, h, f);
#endif
    open->updateheap(node, key);
}
inline
ASTARNode*
ASTARSpace::popOpen() {
    auto node = (ASTARNode*) open->deleteminheap();
    assert(node);
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
    assert(node);

    updateNode(*node);

    return node;
}

inline
ASTARNode*
ASTARSpace::getNode(CMDPSTATE* mdpState) {
    assert(mdpState);
    assert(mdpState->StateID);
    auto node = (ASTARNode*) mdpState->PlannerSpecificData;
    assert(node);

#ifdef DEBUG
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


    updateNode(*node);

    return node;
}

// Shortcuts
inline ASTARNode* ASTARSpace::getStart()  {
    auto startNode = getNode(startState);
    assert(startNode);

    updateNode(*startNode);
    return startNode;
}
inline ASTARNode* ASTARSpace::getGoal()   {
    auto goalNode = getNode(goalState);
    assert(goalNode);

    goalNode->h = 0;
    return goalNode;
}





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
    SBPL_DEBUG("Creating planner at %p\n", this);
    assert(environment);

    space = new ASTARSpace(environment);
    assert(space);

    dbg = nullptr;

    SBPL_DEBUG("Planner created at %p\n", this);
}
ASTARPlanner::~ASTARPlanner(){
    SBPL_DEBUG("Planner at %p will be destroyed\n", this);
    DELETE(space);

    if(dbg){
        if(fclose(dbg))
            TRACE("File[%p] was not succesfully closed", (void*)dbg);
        dbg = nullptr;
    }

    SBPL_DEBUG("Planner at %p was destroyed\n", this);
}




// Search
// ------
bool
ASTARPlanner::Search(vector<stateID> *pathIDs, int *cost, bool bFirstSolution, bool bOptimalSolution, double givenSeconds) {
    assert(space);  // redundant
    assert(space->startState);
    assert(space->goalState);
    assert(pathIDs->empty());

    // TODO: allow backwardSearch

    // Initialization
    // ==============
    SBPL_DEBUG("\n  Preparing A* loop");

    // Reset statistics (just in case)
    stats.iteration = 0;
    stats.expansions = 0;


    // Prepare goal node
    ASTARNode* goalNode = space->getGoal();
    goalNode->g = INFINITECOST;


    // Prepare starting node
    ASTARNode *startNode = space->getStart();
    startNode->g = 0;

    space->upsertOpen(startNode);

    // Expansion
    // =========
    SBPL_DEBUG("\n  Starting A* loop");

    CKey bestKey;
    while(!space->openEmpty()){
        // Get the best node
        bestKey = space->open->getminkeyheap();
        ASTARNode *currentNode = space->popOpen();

        // Check the search should go on
        if(bestKey.key[0] >= INFINITECOST)         // Best node is not reachable
            break;

        // Check if the best node can't improve the path
        if(goalNode->g <= bestKey.key[0])  // Key[0] === f=g+h
            break;

        // Expand the best node
// #if DEBUG
//         float f = bestKey.key[0]/1000.0;
//         float h = bestKey.key[1]/1000.0;
//         float g = f-h;
//         TRACE("Expanding A*Node[%p] state[%4d] (%5.2f + %5.2f = %5.2f)",
//                (void*)currentNode, currentNode->id(), g, h, f);
// #endif
        updateSuccessors(*currentNode);
    }
    SBPL_DEBUG("\n  A*: Finished");


    // Open list information
    // =====================
    if(space->openEmpty())
        SBPL_DEBUG("    A*: No more nodes to reach");

    if(bestKey.key[0] >= INFINITECOST)
        SBPL_DEBUG("    Next bestKey is infinite, no more nodes are reachable");
    else
        SBPL_DEBUG("    Next bestKey is worse than goalKey, no more nodes are worth expanding");
    SBPL_DEBUG("    Goal g: %5.2f", goalNode->g/1000.0);

    return SBPL_OK;
}

inline
void
ASTARPlanner::updateSuccessors(ASTARNode &node) {
    assert(space);
    assert(space->problem);

    stateID id = node.id();
    assert(id);

    // Get successors from environment
    if(!node.successors)
        node.successors = space->problem->GetSuccs(id);
    assert(node.successors);

    // Reach each node updating path&cost on improvements
    for(nodeStub nS : *node.successors)
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
        // Use new path
        neigh->g = newG;  //Save new cost
        neigh->bestPredecesor = node.MDPstate;

        // Update this node on the open list
        space->upsertOpen(neigh);
    }
}

inline
void
ASTARPlanner::updatePredecessors(ASTARNode &node) {
    //TODO: Reuse updateSuccessors code. Find how to do it w/o performance hits
    assert(space);
    assert(space->problem);

    stateID id = node.id();
    assert(id);

    // Get predecessors from environment
    if(!node.predecessors)
        node.predecessors = space->problem->GetPreds(id);
    assert(node.predecessors);

    // Reach each node updating path&cost on improvements
    for(nodeStub nS : *node.predecessors)
        reachPredecessor(node, nS);

    QUIT("%s", "Not yet implemented");
}

inline
void
ASTARPlanner::reachPredecessor(const ASTARNode &node, const nodeStub &nS) {
    //TODO: reach predecessors
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
ASTARPlanner::replan(vector<stateID> *pathIDs, ReplanParams params, int *cost) {
    // TODO: replan with parameters (givenTime, repairTime, epsilon)
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
    // TODO: Drop search efforts
    return SBPL_OK;
}


// SBPL Update API
// ---------------
inline
void
ASTARPlanner::costs_changed(StateChangeQuery const & stateChange) {
    // TODO: update costs
}
