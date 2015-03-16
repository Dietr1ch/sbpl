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
#include <cmath>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <string>
#include <getopt.h>
#include <random>




// Macros (Maybe this should go somewhere else)
// ======


//Math (avoid ++, --, etc...)
#define MAX(a,b) ( (a)>(b) ? (a):(b) )
#define MIN(a,b) ( (a)<(b) ? (a):(b) )
#define CAP(x,m,M) ( (x) < (m) ? (m) : (MIN((x),(M))) )

//Loops
#define LOOP while(1)
#define UNTIL(exp) while(!(exp))
#define SKIP(exp) if(exp) continue;
#define ONLY(exp) if(!(exp)) continue;


using namespace std;

#include <sbpl/headers.h>







// Planners
// ========

// Planner names
#define PLANNER_STR_ADSTAR "adaptiveD*"
#define PLANNER_STR_ARASTAR "ARA*"
#define PLANNER_STR_PPCP "PPCP"
#define PLANNER_STR_RSTAR "R*"
#define PLANNER_STR_VI "vi"
#define PLANNER_STR_ANASTAR "ANA*"
#define PLANNER_STR_AASTAR "adaptiveA*"
#define PLANNER_STR_ASTAR "A*"
#define PLANNER_STR_D "D"

#define PLANNER_STR_INVALID "invalid"


// Planner enum
enum PlannerType {
    INVALID_PLANNER_TYPE = -1,
    PLANNER_TYPE_ADSTAR,
    PLANNER_TYPE_ARASTAR,
    PLANNER_TYPE_PPCP,
    PLANNER_TYPE_RSTAR,
    PLANNER_TYPE_VI,
    PLANNER_TYPE_ANASTAR,
    PLANNER_TYPE_AASTAR,
    PLANNER_TYPE_ASTAR,
    PLANNER_TYPE_D,

    NUM_PLANNER_TYPES
};

// Planner to string
string
PlannerTypeToStr(PlannerType plannerType) {

    switch (plannerType) {
        case PLANNER_TYPE_ADSTAR:
            return string(PLANNER_STR_ADSTAR);
        case PLANNER_TYPE_ARASTAR:
            return string(PLANNER_STR_ARASTAR);
        case PLANNER_TYPE_PPCP:
            return string(PLANNER_STR_PPCP);
        case PLANNER_TYPE_RSTAR:
            return string(PLANNER_STR_RSTAR);
        case PLANNER_TYPE_VI:
            return string(PLANNER_STR_VI);
        case PLANNER_TYPE_ANASTAR:
            return string(PLANNER_STR_ANASTAR);
        case PLANNER_TYPE_AASTAR:
            return string(PLANNER_STR_AASTAR);
        case PLANNER_TYPE_ASTAR:
            return string(PLANNER_STR_ASTAR);
        case PLANNER_TYPE_D:
            return string(PLANNER_STR_D);

        default:
            return string(PLANNER_STR_INVALID);
    }
}

// String to planner
PlannerType
StrToPlannerType(const char* str) {
    if (!strcmp(str, PLANNER_STR_ADSTAR))
        return PLANNER_TYPE_ADSTAR;
    if (!strcmp(str, PLANNER_STR_ARASTAR))
        return PLANNER_TYPE_ARASTAR;
    if (!strcmp(str, PLANNER_STR_PPCP))
        return PLANNER_TYPE_PPCP;
    if (!strcmp(str, PLANNER_STR_RSTAR))
        return PLANNER_TYPE_RSTAR;
    if (!strcmp(str, PLANNER_STR_VI))
        return PLANNER_TYPE_VI;
    if (!strcmp(str, PLANNER_STR_ANASTAR))
        return PLANNER_TYPE_ANASTAR;
    if (!strcmp(str, PLANNER_STR_AASTAR))
        return PLANNER_TYPE_AASTAR;
    if (!strcmp(str, PLANNER_STR_ASTAR))
        return PLANNER_TYPE_ASTAR;
    if (!strcmp(str, PLANNER_STR_D))
        return PLANNER_TYPE_D;

    return INVALID_PLANNER_TYPE;
}













// Environments
// ============
#define ENV_STR_2D "2D"
#define ENV_STR_2DUU "2DUU"
#define ENV_STR_ROBARM "RobotArm"
#define ENV_STR_XYTHETA "X-Y-Theta"
#define ENV_STR_XYTHETAMLEV "X-Y-Theta_mlev"

#define ENV_STR_INVALID "Invalid"

// Environment types
enum EnvironmentType {
    INVALID_ENV_TYPE = -1,
    ENV_TYPE_2D,
    ENV_TYPE_2DUU,
    ENV_TYPE_ROBARM,
    ENV_TYPE_XYTHETA,
    ENV_TYPE_XYTHETAMLEV,

    NUM_ENV_TYPES
};

// Environment to string
inline
string
EnvironmentTypeToStr(EnvironmentType environmentType) {
    switch (environmentType) {
    case ENV_TYPE_2D:
        return string(ENV_STR_2D);
    case ENV_TYPE_2DUU:
        return string(ENV_STR_2DUU);
    case ENV_TYPE_ROBARM:
        return string(ENV_STR_ROBARM);
    case ENV_TYPE_XYTHETA:
        return string(ENV_STR_XYTHETA);
    case ENV_TYPE_XYTHETAMLEV:
        return string(ENV_STR_XYTHETAMLEV);

    default:
        return string(ENV_STR_INVALID);
    }
}

// String to Environment
inline
EnvironmentType
StrToEnvironmentType(const char* str) {
    if (!strcmp(str, ENV_STR_2D))
        return ENV_TYPE_2D;
    if (!strcmp(str, ENV_STR_2DUU))
        return ENV_TYPE_2DUU;
    if (!strcmp(str, ENV_STR_ROBARM))
        return ENV_TYPE_ROBARM;
    if (!strcmp(str, ENV_STR_XYTHETA))
        return ENV_TYPE_XYTHETA;
    if (!strcmp(str, ENV_STR_XYTHETAMLEV))
        return ENV_TYPE_XYTHETAMLEV;


    return INVALID_ENV_TYPE;
}




enum
MainResultType {
    INVALID_MAIN_RESULT = -1,

    MAIN_RESULT_SUCCESS = 0,
    MAIN_RESULT_FAILURE = 1,
    MAIN_RESULT_INSUFFICIENT_ARGS = 2,
    MAIN_RESULT_INCORRECT_OPTIONS = 3,
    MAIN_RESULT_UNSUPPORTED_ENV = 4,

    NUM_MAIN_RESULTS
};
























// Command line arguments parsing
// ==============================


// Settings
// --------
// Planner
char *planner_type     = nullptr;
char *planner_settings = nullptr;
PlannerType selectedPlanner = PLANNER_TYPE_D;
SBPLPlanner *planner = nullptr;

// Random environments (defaults to no random environments)
char *env_type     = nullptr;
char *env_settings = nullptr;
EnvironmentType selectedEnvironment = ENV_TYPE_2D;
DiscreteSpaceInformation *environment = nullptr;

int env_random_start = 0;
int env_random_count = 0;
int env_random_end;  // start+count
int env_id;

// Runs (Defaults to test run 0 only)
int run_start = 0;
int run_count = 1;
int run_end;  // start+count
int run_id;

// Runs (Defaults to test run 0 only)
int maxTries = 100;
int it_start = 0;
int it_count = 0;
int it_end;  // start+count
int it_id;

// Output text
char *output_raw_fileName = nullptr;
bool  output_print_raw    = false;
// Output sql
char *output_sql_fileName = nullptr;
bool  output_print_sql    = false;

#define DEFAULT_TIME 1000
double givenSeconds = DEFAULT_TIME;
double givenTime = DEFAULT_TIME;
bool backwardSearch = true;


/**
 * Prints the help  **flies away**
 */
//void printHelp(option *opts[]){
void
printHelp(){
    // TODO: print help (from option[]?)
    cout << "" << endl;
}

/**
 * Parses all settings.
 */
int
parseInput(int argc, char **argv){


    // Available options
    // =================
    static struct option long_options[] = {
    //  {name, (no_argument/required_argument/optional_argument), int *flag, val}

        // Planner settings
        {"planner",           optional_argument, 0, 'p'},
        {"planner-settings",  optional_argument, 0, 'P'},

        // Random environment type
        {"env-type",          optional_argument, 0, 'e'},
        {"env-settings",      optional_argument, 0, 'E'},
        // Environments to generate
        {"env-random-count",  optional_argument, 0, 'n'},
        {"env-random-start",  optional_argument, 0, 'N'},


        // Search Algorithm
        {"alg", optional_argument, 0, 'a'},


        // Runs to test on each environment
        {"run-count", optional_argument, 0, 'r'},
        {"run-start", optional_argument, 0, 'R'},

        // Runs to test on each environment
        {"iteration-count", optional_argument, 0, 'i'},
        {"iteration-start", optional_argument, 0, 'I'},

        // Output files
        {"csv-file", optional_argument, 0, 'o'},
        {"csv-file", optional_argument, 0, 'O'},
        {"sql-file", optional_argument, 0, 's'},
        {"sql-file", optional_argument, 0, 'S'},


        // Help
        {"help", optional_argument, 0, 'h'},

        {0, 0, 0, 0}  // End of arguments
    };


    // Parsing loop
    // ============
    opterr = 0;
    int getopt_ret, option_index;
    char *tail;
    LOOP {
        getopt_ret = getopt_long(
            argc,
            argv,
            "p::P::e::E:::n:N:a:r:R:i:I:o::O::s::S::",
            long_options,
            &option_index);

        if (getopt_ret == -1)
            break;

        switch (getopt_ret) {
            // Help
            case 'h':
                printHelp();
                exit(0);

            // Random environments
            case 'p':
                if (!optarg || !strlen(optarg))
                    break;
                planner_type = optarg;
                selectedPlanner = StrToPlannerType(planner_type);
                break;
            case 'P':
                if (!optarg || !strlen(optarg))
                    break;
                planner_settings = optarg;
                break;
            case 'e':
                if (!optarg || !strlen(optarg))
                    break;
                env_type = optarg;
                selectedEnvironment = StrToEnvironmentType(env_type);
                break;
            case 'E':
                if (!optarg || !strlen(optarg))
                    break;
                env_settings = optarg;
                break;


            // Random Environments
            case 'n':
                env_random_count = strtol(optarg, &tail, 0);
                env_random_count = MAX(env_random_count, 0);
                break;
            case 'N':
                env_random_start = strtol(optarg, &tail, 0);
                break;

            // Runs to test on each environment
            case 'r':
                run_count = strtol(optarg, &tail, 0);
                run_count = MAX(run_count, 1);
                break;
            case 'R':
                run_start = strtol(optarg, &tail, 0);
                break;
            // Replanning itearations to test on each run
            case 'i':
                it_count = strtol(optarg, &tail, 0);
                it_count = MAX(it_count, 1);
                break;
            case 'I':
                it_start = strtol(optarg, &tail, 0);
                break;



             // Output files
            case 'o':
            case 'O':
                if (!optarg || !strlen(optarg))
                    break;
                output_raw_fileName = optarg;
                output_print_raw = 1;
                break;
            case 's':
            case 'S':
                if (!optarg || !strlen(optarg))
                    break;
                output_sql_fileName = optarg;
                output_print_sql = 1;
                break;

            case 0:
                break;
            case '?':
            default:
                printf("Unknown option\n");
                return 10;
        }
    }

    // Post Parsing
    // ============
    env_random_end = env_random_start + env_random_count;
    run_end = run_start + run_count;
    it_end = it_start + it_count;

    return 0;
}















// Planning environment handlers
// =============================

/*******************************************************************************
 * plan2d
 * @brief An example of planning in two-dimensional space without motion primitives.
 *
 * @param plannerType The type of planner to be used in this example
 * @param envCfgFilename The environment config file. See
 *                       sbpl/env_examples/nav2d/ for examples
 * @return 1 if the planner successfully found a solution; 0 otherwise
 *******************************************************************************/
int
plan2d(PlannerType plannerType, char* envCfgFilename, bool forwardSearch) {
    int bRet = 0;
    double allocated_time_secs = 100.0; // in seconds
    double initialEpsilon = 3.0;
    MDPConfig MDPCfg;
    bool bsearchuntilfirstsolution = false;
    bool bforwardsearch = forwardSearch;

    // Initialize Environment (should be called before initializing anything else)
    EnvironmentNAV2D environment_nav2D;
    if (!environment_nav2D.InitializeEnv(envCfgFilename)) {
        printf("ERROR: InitializeEnv failed\n");
        throw new SBPL_Exception();
    }

    // Initialize MDP Info
    if (!environment_nav2D.InitializeMDPCfg(&MDPCfg)) {
        printf("ERROR: InitializeMDPCfg failed\n");
        throw new SBPL_Exception();
    }

    // plan a path
    Path solution_stateIDs_V;

    SBPLPlanner* planner = NULL;
    switch (plannerType) {
//TODO add new planners
    case PLANNER_TYPE_ARASTAR:
        printf("Initializing ARAPlanner...\n");
        planner = new ARAPlanner(&environment_nav2D, bforwardsearch);
        break;
    case PLANNER_TYPE_ADSTAR:
        printf("Initializing ADPlanner...\n");
        planner = new ADPlanner(&environment_nav2D, bforwardsearch);
        break;
    case PLANNER_TYPE_RSTAR:
        printf("Initializing RSTARPlanner...\n");
        planner = new RSTARPlanner(&environment_nav2D, bforwardsearch);
        break;
    case PLANNER_TYPE_ANASTAR:
        printf("Initializing anaPlanner...\n");
        planner = new anaPlanner(&environment_nav2D, bforwardsearch);
        break;

    default:
        printf("Invalid planner type\n");
        break;
    }

    // set search mode
    planner->set_search_mode(bsearchuntilfirstsolution);

    if (planner->set_start(MDPCfg.startstateid) == 0) {
        printf("ERROR: failed to set start state\n");
        throw new SBPL_Exception();
    }

    if (planner->set_goal(MDPCfg.goalstateid) == 0) {
        printf("ERROR: failed to set goal state\n");
        throw new SBPL_Exception();
    }

    planner->set_initialsolution_eps(initialEpsilon);

    printf("start planning...\n");
    bRet = planner->replan(allocated_time_secs, &solution_stateIDs_V);
    printf("done planning\n");
    cout << "size of solution=" << solution_stateIDs_V.size() << endl;

    environment_nav2D.PrintTimeStat(stdout);

    const char* sol = "sol.txt";
    FILE* fSol = fopen(sol, "w");
    if (fSol == NULL) {
        printf("ERROR: could not open solution file\n");
        throw new SBPL_Exception();
    }
    for (unsigned int i = 0; i < solution_stateIDs_V.size(); i++) {
        environment_nav2D.PrintState(solution_stateIDs_V[i], false, fSol);
    }
    fclose(fSol);

    environment_nav2D.PrintTimeStat(stdout);

    //print a path
    if (bRet) {
        //print the solution
        printf("Solution is found\n");
    }
    else {
        printf("Solution does not exist\n");
    }

    fflush(NULL);

    delete planner;

    return bRet;
}

/*******************************************************************************
 * plan2duu
 * @brief An examle of planning in two-dimensional space and under uncertainty.
 *
 * @param plannerType The type of planner to be used in this example
 * @param envCfgFilename The environment config file. See
 *                       sbpl/env_examples/nav2duu/ for examples
 * @return 1 if the planner successfully found a solution; 0 otherwise
 *******************************************************************************/
int
plan2duu(PlannerType plannerType, char* envCfgFilename) {
    int bRet = 0;
    double allocated_time_secs = 10.0; //in seconds
    MDPConfig MDPCfg;

    //Initialize Environment (should be called before initializing anything else)
    EnvironmentNAV2DUU environment_nav2Duu;
    if (!environment_nav2Duu.InitializeEnv(envCfgFilename)) {
        printf("ERROR: InitializeEnv failed\n");
        throw new SBPL_Exception();
    }

    //Initialize MDP Info
    if (!environment_nav2Duu.InitializeMDPCfg(&MDPCfg)) {
        printf("ERROR: InitializeMDPCfg failed\n");
        throw new SBPL_Exception();
    }

    //create the planner
//TODO add new planners
    PPCPPlanner planner(&environment_nav2Duu, environment_nav2Duu.SizeofCreatedEnv(), environment_nav2Duu.SizeofH());

    //set start and goal
    if (planner.set_start(MDPCfg.startstateid) == 0) {
        printf("ERROR: failed to set start state\n");
        throw new SBPL_Exception();
    }
    if (planner.set_goal(MDPCfg.goalstateid) == 0) {
        printf("ERROR: failed to set goal state\n");
        throw new SBPL_Exception();
    }

    printf("start planning...\n");
    float ExpectedCost, ProbofReachGoal;
    vector<sbpl_PolicyStatewithBinaryh_t> SolutionPolicy;
    bRet = planner.replan(allocated_time_secs, &SolutionPolicy, &ExpectedCost, &ProbofReachGoal);
    printf("done planning\n");

    if (bRet) {
        //print the solution
        printf("Solution is found: exp. cost=%f probofreachgoal=%f\n", ExpectedCost, ProbofReachGoal);
    }
    else {
        printf("Solution does not exist\n");
    }

    fflush(NULL);

    return bRet;
}

/*******************************************************************************
 * planxythetalat
 * @brief An example of planning in three-dimensional space (x,y,theta)
 *
 * @param plannerType The type of planner to be used in this example
 * @param envCfgFilename The environment config file. See
 *                       sbpl/env_examples/nav3d/ for examples
 * @param motPrimFilename The motion primitives file. See
 *                        sbpl/matlab/mprim/ for examples
 * @return 1 if the planner successfully found a solution; 0 otherwise
 *******************************************************************************/
int
planxythetalat(PlannerType plannerType, char* envCfgFilename, char* motPrimFilename, bool forwardSearch) {
    int bRet = 0;
    double allocated_time_secs = 10.0; // in seconds
    double initialEpsilon = 3.0;
    MDPConfig MDPCfg;
    bool bsearchuntilfirstsolution = false;
    bool bforwardsearch = forwardSearch;

    // set the perimeter of the robot (it is given with 0,0,0 robot ref. point for which planning is done)
    vector<sbpl_2Dpt_t> perimeterptsV;
    sbpl_2Dpt_t pt_m;
    double halfwidth = 0.01; //0.3;
    double halflength = 0.01; //0.45;
    pt_m.x = -halflength;
    pt_m.y = -halfwidth;
    perimeterptsV.push_back(pt_m);
    pt_m.x = halflength;
    pt_m.y = -halfwidth;
    perimeterptsV.push_back(pt_m);
    pt_m.x = halflength;
    pt_m.y = halfwidth;
    perimeterptsV.push_back(pt_m);
    pt_m.x = -halflength;
    pt_m.y = halfwidth;
    perimeterptsV.push_back(pt_m);

    // clear the footprint
    perimeterptsV.clear();

    // Initialize Environment (should be called before initializing anything else)
    EnvironmentNAVXYTHETALAT environment_navxythetalat;

    if (!environment_navxythetalat.InitializeEnv(envCfgFilename, perimeterptsV, motPrimFilename)) {
        printf("ERROR: InitializeEnv failed\n");
        throw new SBPL_Exception();
    }

    // Initialize MDP Info
    if (!environment_navxythetalat.InitializeMDPCfg(&MDPCfg)) {
        printf("ERROR: InitializeMDPCfg failed\n");
        throw new SBPL_Exception();
    }

    // plan a path
    Path solution_stateIDs_V;

    SBPLPlanner* planner = NULL;
    switch (plannerType) {
//TODO add new planners
    case PLANNER_TYPE_ARASTAR:
        printf("Initializing ARAPlanner...\n");
        planner = new ARAPlanner(&environment_navxythetalat, bforwardsearch);
        break;
    case PLANNER_TYPE_ADSTAR:
        printf("Initializing ADPlanner...\n");
        planner = new ADPlanner(&environment_navxythetalat, bforwardsearch);
        break;
    case PLANNER_TYPE_RSTAR:
        printf("Invalid configuration: xytheta environment does not support rstar planner...\n");
        return 0;
    case PLANNER_TYPE_ANASTAR:
        printf("Initializing anaPlanner...\n");
        planner = new anaPlanner(&environment_navxythetalat, bforwardsearch);
        break;

    default:
        printf("Invalid planner type\n");
        break;
    }

    // set planner properties
    if (planner->set_start(MDPCfg.startstateid) == 0) {
        printf("ERROR: failed to set start state\n");
        throw new SBPL_Exception();
    }
    if (planner->set_goal(MDPCfg.goalstateid) == 0) {
        printf("ERROR: failed to set goal state\n");
        throw new SBPL_Exception();
    }
    planner->set_initialsolution_eps(initialEpsilon);
    planner->set_search_mode(bsearchuntilfirstsolution);

    // plan
    printf("start planning...\n");
    bRet = planner->replan(allocated_time_secs, &solution_stateIDs_V);
    printf("done planning\n");
    printf("size of solution=%d\n", (unsigned int)solution_stateIDs_V.size());

    environment_navxythetalat.PrintTimeStat(stdout);

    // write solution to sol.txt
    const char* sol = "sol.txt";
    FILE* fSol = fopen(sol, "w");
    if (fSol == NULL) {
        printf("ERROR: could not open solution file\n");
        throw new SBPL_Exception();
    }

    // write the discrete solution to file
    //	for (StateID = 0; i < solution_stateIDs_V.size(); i++) {
    //		int x;
    //		int y;
    //		int theta;
    //		environment_navxythetalat.GetCoordFromState(solution_stateIDs_V[i], x, y, theta);
    //
    //		fprintf(fSol, "%d %d %d\t\t%.3f %.3f %.3f\n", x, y, theta,
    //              DISCXY2CONT(x, 0.1), DISCXY2CONT(y, 0.1), DiscTheta2Cont(theta, 16));
    //	}

    // write the continuous solution to file
    vector<sbpl_xy_theta_pt_t> xythetaPath;
    environment_navxythetalat.ConvertStateIDPathintoXYThetaPath(&solution_stateIDs_V, &xythetaPath);
    printf("solution size=%d\n", (unsigned int)xythetaPath.size());
    for (unsigned int i = 0; i < xythetaPath.size(); i++) {
        fprintf(fSol, "%.3f %.3f %.3f\n", xythetaPath.at(i).x, xythetaPath.at(i).y, xythetaPath.at(i).theta);
    }
    fclose(fSol);

    environment_navxythetalat.PrintTimeStat(stdout);

    // print a path
    if (bRet) {
        // print the solution
        printf("Solution is found\n");
    }
    else {
        printf("Solution does not exist\n");
    }

    fflush(NULL);

    delete planner;

    return bRet;
}

/*******************************************************************************
 * planxythetamlevlat
 * @brief An example of planning with a multiple-level (x,y,theta) lattice
 *
 * @desc An example of planning with a multiple-level x,y,theta lattice (for
 *       example, a base and upper body). There are no additional degrees of
 *       freedom but each level may have a different footprint and should be
 *       checked against a corresponding costmap at its height. Useful for doing
 *       navigation of a tall ground robot operating in a cluttered 3d map.
 * @param plannerType The type of planner to be used in this example
 * @param envCfgFilename The environment config file. See
 *                       sbpl/env_examples/nav3d/ for examples
 * @param motPrimFilename The motion primitives file. See
 *                        sbpl/matlab/mprim/ for examples
 * @return 1 if the planner successfully found a solution; 0 otherwise
 *******************************************************************************/
int
planxythetamlevlat(PlannerType plannerType, char* envCfgFilename, char* motPrimFilename, bool forwardSearch) {
    int bRet = 0;

    double allocated_time_secs = 10.0; //in seconds
    double initialEpsilon = 3.0;
    MDPConfig MDPCfg;
    bool bsearchuntilfirstsolution = false;
    bool bforwardsearch = forwardSearch;

    //set the perimeter of the robot (it is given with 0,0,0 robot ref. point for which planning is done)
    //this is for the default level - base level
    vector<sbpl_2Dpt_t> perimeterptsV;
    sbpl_2Dpt_t pt_m;
    double halfwidth = 0.02; //0.3;
    double halflength = 0.02; //0.45;
    pt_m.x = -halflength;
    pt_m.y = -halfwidth;
    perimeterptsV.push_back(pt_m);
    pt_m.x = halflength;
    pt_m.y = -halfwidth;
    perimeterptsV.push_back(pt_m);
    pt_m.x = halflength;
    pt_m.y = halfwidth;
    perimeterptsV.push_back(pt_m);
    pt_m.x = -halflength;
    pt_m.y = halfwidth;
    perimeterptsV.push_back(pt_m);

    //Initialize Environment (should be called before initializing anything else)
    EnvironmentNAVXYTHETAMLEVLAT environment_navxythetalat;

    if (!environment_navxythetalat.InitializeEnv(envCfgFilename, perimeterptsV, motPrimFilename)) {
        printf("ERROR: InitializeEnv failed\n");
        throw new SBPL_Exception();
    }

    //this is for the second level - upper body level
    vector<sbpl_2Dpt_t> perimeterptsVV[2];
    perimeterptsVV[0].clear();
    halfwidth = 0.02;
    halflength = 0.02;
    pt_m.x = -halflength;
    pt_m.y = -halfwidth;
    perimeterptsVV[0].push_back(pt_m);
    pt_m.x = halflength;
    pt_m.y = -halfwidth;
    perimeterptsVV[0].push_back(pt_m);
    pt_m.x = halflength;
    pt_m.y = halfwidth;
    perimeterptsVV[0].push_back(pt_m);
    pt_m.x = -halflength;
    pt_m.y = halfwidth;
    perimeterptsVV[0].push_back(pt_m);

    //	perimeterptsV.clear();
    //	perimeterptsVV[0].clear();

    //enable the second level
    int numofaddlevels = 1;
    printf("Number of additional levels = %d\n", numofaddlevels);
    unsigned char cost_inscribed_thresh_addlevels[2]; //size should be at least numofaddlevels
    unsigned char cost_possibly_circumscribed_thresh_addlevels[2]; //size should be at least numofaddlevels
    //no costs are indicative of whether a cell is within inner circle
    cost_inscribed_thresh_addlevels[0] = 255;
    //no costs are indicative of whether a cell is within outer circle
    cost_possibly_circumscribed_thresh_addlevels[0] = 0;
    //no costs are indicative of whether a cell is within inner circle
    cost_inscribed_thresh_addlevels[1] = 255;
    //no costs are indicative of whether a cell is within outer circle
    cost_possibly_circumscribed_thresh_addlevels[1] = 0;
    if (!environment_navxythetalat.InitializeAdditionalLevels(numofaddlevels, perimeterptsVV,
                                                              cost_inscribed_thresh_addlevels,
                                                              cost_possibly_circumscribed_thresh_addlevels))
    {
        printf("ERROR: InitializeAdditionalLevels failed with numofaddlevels=%d\n", numofaddlevels);
        throw new SBPL_Exception();
    }

    //set the map for the second level (index parameter for the additional levels and is zero based)
    //for this example, we pass in the same map as the map for the base. In general, it can be a totally different map
    //as it corresponds to a different level
    //NOTE: this map has to have costs set correctly with respect to inner and outer radii of the robot
    //if the second level of the robot has these radii different than at the base level, then costs
    //should reflect this
    //(see explanation for cost_possibly_circumscribed_thresh and
    //cost_inscribed_thresh parameters in environment_navxythetalat.h file)
    int addlevind = 0;
    if (!environment_navxythetalat.Set2DMapforAddLev(
            (const unsigned char**)(environment_navxythetalat.GetEnvNavConfig()->Grid2D), addlevind))
    {
        printf("ERROR: Setting Map for the Additional Level failed with level %d\n", addlevind);
        throw new SBPL_Exception();
    }

    //Initialize MDP Info
    if (!environment_navxythetalat.InitializeMDPCfg(&MDPCfg)) {
        printf("ERROR: InitializeMDPCfg failed\n");
        throw new SBPL_Exception();
    }

    //plan a path
    Path solution_stateIDs_V;

    SBPLPlanner* planner = NULL;
    switch (plannerType) {
//TODO add new planners
    case PLANNER_TYPE_ARASTAR:
        printf("Initializing ARAPlanner...\n");
        planner = new ARAPlanner(&environment_navxythetalat, bforwardsearch);
        break;
    case PLANNER_TYPE_ADSTAR:
        printf("Initializing ADPlanner...\n");
        planner = new ADPlanner(&environment_navxythetalat, bforwardsearch);
        break;
    case PLANNER_TYPE_RSTAR:
        printf("Invalid configuration: xytheta environment does not support rstar planner...\n");
        return 0;
    case PLANNER_TYPE_ANASTAR:
        printf("Initializing anaPlanner...\n");
        planner = new anaPlanner(&environment_navxythetalat, bforwardsearch);
        break;

    default:
        printf("Invalid planner type\n");
        break;
    }

    if (planner->set_start(MDPCfg.startstateid) == 0) {
        printf("ERROR: failed to set start state\n");
        throw new SBPL_Exception();
    }

    if (planner->set_goal(MDPCfg.goalstateid) == 0) {
        printf("ERROR: failed to set goal state\n");
        throw new SBPL_Exception();
    }
    planner->set_initialsolution_eps(initialEpsilon);

    //set search mode
    planner->set_search_mode(bsearchuntilfirstsolution);

    printf("start planning...\n");
    bRet = planner->replan(allocated_time_secs, &solution_stateIDs_V);
    printf("done planning\n");
    cout << "size of solution=" << solution_stateIDs_V.size() << endl;

    environment_navxythetalat.PrintTimeStat(stdout);

    const char* sol = "sol.txt";
    FILE* fSol = fopen(sol, "w");
    if (fSol == NULL) {
        printf("ERROR: could not open solution file\n");
        throw new SBPL_Exception();
    }
    vector<sbpl_xy_theta_pt_t> xythetaPath;
    environment_navxythetalat.ConvertStateIDPathintoXYThetaPath(&solution_stateIDs_V, &xythetaPath);
    printf("solution size=%d\n", (unsigned int)xythetaPath.size());
    for (unsigned int i = 0; i < xythetaPath.size(); i++) {
        fprintf(fSol, "%.3f %.3f %.3f\n", xythetaPath.at(i).x, xythetaPath.at(i).y, xythetaPath.at(i).theta);
    }
    fclose(fSol);

    environment_navxythetalat.PrintTimeStat(stdout);

    //print a path
    if (bRet) {
        //print the solution
        printf("Solution is found\n");
    }
    else {
        printf("Solution does not exist\n");
    }

    fflush(NULL);

    delete planner;

    return bRet;
}

/*******************************************************************************
 * planandnavigate2d
 * @brief An example simulation of how a robot would use two-dimensional
 *        environment planning.
 *
 * @param planner The planner to be used in this example
 * @param envCfgFilename The environment config file. See
 *                       sbpl/env_examples/nav2d/ for examples
 *******************************************************************************/
int
planandnavigate2d(PlannerType plannerType, char* envCfgFilename) {
    SBPL_DEBUG("Starting 2D planning");

    double allocated_time_secs_foreachplan = 0.2; //in seconds
    MDPConfig MDPCfg;
    EnvironmentNAV2D environment_nav2D;
    EnvironmentNAV2D trueenvironment_nav2D;
    int size_x = -1, size_y = -1;
    int startx = 0, starty = 0;
    int goalx = -1, goaly = -1;
    const char* sol = "sol.txt";
    FILE* fSol = fopen(sol, "w");
    if (fSol == NULL) {
        printf("ERROR: could not open solution file\n");
        throw new SBPL_Exception();
    }
    //int dx[8] = {-1, -1, -1,  0,  0,  1,  1,  1};
    //int dy[8] = {-1,  0,  1, -1,  1, -1,  0,  1};
    bool bPrint = false;
    int x, y;
    Path preds_of_changededgesIDV;
    vector<nav2dcell_t> changedcellsV;
    nav2dcell_t nav2dcell;
    unsigned char obsthresh = 0;
    srand(0);
    int plantime_over1secs = 0, plantime_over0p5secs = 0, plantime_over0p1secs = 0, plantime_over0p05secs = 0,
        plantime_below0p05secs = 0;

    //set parameters - should be done before initialization
    if (!trueenvironment_nav2D.SetEnvParameter("is16connected", 1)) {
        printf("ERROR: failed to set parameters\n");
        throw new SBPL_Exception();
    }
    if (!environment_nav2D.SetEnvParameter("is16connected", 1)) {
        printf("ERROR: failed to set parameters\n");
        throw new SBPL_Exception();
    }

    //initialize true map and robot map
    if (!trueenvironment_nav2D.InitializeEnv(envCfgFilename)) {
        printf("ERROR: InitializeEnv failed\n");
        throw new SBPL_Exception();
    }
    trueenvironment_nav2D.GetEnvParms(&size_x, &size_y, &startx, &starty, &goalx, &goaly, &obsthresh);
    unsigned char* map = (unsigned char*)calloc(size_x * size_y, sizeof(unsigned char));

    //print the map
    if (bPrint) printf("true map:\n");
    for (y = 0; bPrint && y < size_y; y++) {
        for (x = 0; x < size_x; x++) {
            printf("%d ", (int)trueenvironment_nav2D.IsObstacle(x, y));
        }
        printf("\n");
    }
    if (bPrint) printf("System Pause (return=%d)\n", system("pause"));

    //Initialize Environment (should be called before initializing anything else)
    if (!environment_nav2D.InitializeEnv(size_x, size_y, map, startx, starty, goalx, goaly, obsthresh)) {
        printf("ERROR: InitializeEnv failed\n");
        throw new SBPL_Exception();
    }

    //Initialize MDP Info
    if (!environment_nav2D.InitializeMDPCfg(&MDPCfg)) {
        printf("ERROR: InitializeMDPCfg failed\n");
        throw new SBPL_Exception();
    }

    //create a planner
    Path solution_stateIDs_V;
    bool bforwardsearch = false;

    SBPLPlanner* planner = NULL;
    switch (plannerType) {
//TODO add new planners
    case PLANNER_TYPE_ARASTAR:
        printf("Initializing ARAPlanner...\n");
        planner = new ARAPlanner(&environment_nav2D, bforwardsearch);
        break;
    case PLANNER_TYPE_ADSTAR:
        printf("Initializing ADPlanner...\n");
        planner = new ADPlanner(&environment_nav2D, bforwardsearch);
        break;
    case PLANNER_TYPE_RSTAR:
        printf("Initializing RSTARPlanner...\n");
        planner = new RSTARPlanner(&environment_nav2D, bforwardsearch);
        break;
    case PLANNER_TYPE_ANASTAR:
        printf("Initializing anaPlanner...\n");
        planner = new anaPlanner(&environment_nav2D, bforwardsearch);
        break;

    default:
        printf("Invalid planner type\n");
        break;
    }

    planner->set_initialsolution_eps(2.0);

    //set the start and goal configurations
    if (planner->set_start(MDPCfg.startstateid) == 0) {
        printf("ERROR: failed to set start state\n");
        throw new SBPL_Exception();
    }
    if (planner->set_goal(MDPCfg.goalstateid) == 0) {
        printf("ERROR: failed to set goal state\n");
        throw new SBPL_Exception();
    }

    //set search mode
    planner->set_search_mode(false);

    //now comes the main loop
    int goalthresh = 0;
    while (abs(startx - goalx) > goalthresh || abs(starty - goaly) > goalthresh) {

        //simulate sensor data update
        bool bChanges = false;
        preds_of_changededgesIDV.clear();
        changedcellsV.clear();

        int dX = 0;
        int dY = 0;
        for (dX = -2; dX <= 2; dX++) {
            for (dY = -2; dY <= 2; dY++) {
                int x = startx + dX;
                int y = starty + dY;
                if (x < 0 || x >= size_x || y < 0 || y >= size_y) {
                    continue;
                }
                int index = x + y * size_x;
                unsigned char truecost = trueenvironment_nav2D.GetMapCost(x, y);
                if (map[index] != truecost) {
                    map[index] = truecost;
                    environment_nav2D.UpdateCost(x, y, map[index]);
                    printf("setting cost[%d][%d] to %d\n", x, y, map[index]);
                    bChanges = true;
                    // store the changed cells
                    nav2dcell.x = x;
                    nav2dcell.y = y;
                    changedcellsV.push_back(nav2dcell);
                }
            }
        }

        double TimeStarted = clock();

        if (bChanges) {
            if (dynamic_cast<ARAPlanner*> (planner) != NULL) {
                ((ARAPlanner*)planner)->costs_changed(); //use by ARA* planner (non-incremental)
            }
            else if (dynamic_cast<ADPlanner*> (planner) != NULL) {
                //get the affected states
                environment_nav2D.GetPredsofChangedEdges(&changedcellsV, &preds_of_changededgesIDV);
                //let know the incremental planner about them
                ((ADPlanner*)planner)->update_preds_of_changededges(&preds_of_changededgesIDV);
            }
        }
        //planner.force_planning_from_scratch();

        fprintf(fSol, "%d %d ", startx, starty);

        //plan a path
        bool bPlanExists = false;
        while (bPlanExists == false) {
            printf("new planning...\n");
            bPlanExists = (planner->replan(allocated_time_secs_foreachplan, &solution_stateIDs_V) == 1);
            printf("done with the solution of size=%d\n", (unsigned int)solution_stateIDs_V.size());
            environment_nav2D.PrintTimeStat(stdout);
            if (bPlanExists == false) {
                throw new SBPL_Exception();
            }

            //for(unsigned int i = 0; i < solution_stateIDs_V.size(); i++) {
            //environment_nav2D.PrintState(solution_stateIDs_V[i], true, fSol);
            //}
            //fprintf(fSol, "*********\n");
        }

        double plantime_secs = (clock() - TimeStarted) / ((double)CLOCKS_PER_SEC);
        fprintf(fSol, "%.5f %.5f\n", plantime_secs, planner->get_solution_eps());
        fflush(fSol);
        if (plantime_secs > 1.0)
            plantime_over1secs++;
        else if (plantime_secs > 0.5)
            plantime_over0p5secs++;
        else if (plantime_secs > 0.1)
            plantime_over0p1secs++;
        else if (plantime_secs > 0.05)
            plantime_over0p05secs++;
        else
            plantime_below0p05secs++;

        //print the map
        int startindex = startx + starty * size_x;
        int goalindex = goalx + goaly * size_x;
        for (y = 0; bPrint && y < size_y; y++) {
            for (x = 0; x < size_x; x++) {
                int index = x + y * size_x;

                //check to see if it is on the path
                bool bOnthePath = false;
                for (int j = 1; j < (int)solution_stateIDs_V.size(); j++) {
                    int newx, newy;
                    environment_nav2D.GetCoordFromState(solution_stateIDs_V[j], newx, newy);
                    if (x == newx && y == newy) bOnthePath = true;
                }

                if (index != startindex && index != goalindex && !bOnthePath)
                    printf("%3d ", map[index]);
                else if (index == startindex)
                    printf("  R ");
                else if (index == goalindex)
                    printf("  G ");
                else if (bOnthePath)
                    printf("  * ");
                else
                    printf("  ? ");
            }
            printf("\n");
        }
        if (bPrint) printf("System Pause (return=%d)\n", system("pause"));

        //move along the path
        if (bPlanExists && (int)solution_stateIDs_V.size() > 1) {
            //get coord of the successor
            int newx, newy;
            environment_nav2D.GetCoordFromState(solution_stateIDs_V[1], newx, newy);

            if (trueenvironment_nav2D.GetMapCost(newx, newy) >= obsthresh) {
                printf("ERROR: robot is commanded to move into an obstacle\n");
                throw new SBPL_Exception();
            }

            //move
            printf("moving from %d %d to %d %d\n", startx, starty, newx, newy);
            startx = newx;
            starty = newy;

            //update the environment
            environment_nav2D.SetStart(startx, starty);

            //update the planner
            if (planner->set_start(solution_stateIDs_V[1]) == 0) {
                printf("ERROR: failed to update robot pose in the planner\n");
                throw new SBPL_Exception();
            }
        }

    }

    //print stats
    printf("stats: plantimes over 1 secs=%d; over 0.5 secs=%d; over 0.1 secs=%d; "
           "over 0.05 secs=%d; below 0.05 secs=%d\n",
           plantime_over1secs, plantime_over0p5secs, plantime_over0p1secs, plantime_over0p05secs,
           plantime_below0p05secs);
    fprintf(fSol, "stats: plantimes over 1 secs=%d; over 0.5; secs=%d; over 0.1 secs=%d; "
            "over 0.05 secs=%d; below 0.05 secs=%d\n",
            plantime_over1secs, plantime_over0p5secs, plantime_over0p1secs, plantime_over0p05secs,
            plantime_below0p05secs);

    if (bPrint) printf("System Pause (return=%d)\n", system("pause"));

    fflush(NULL);
    fclose(fSol);

    delete planner;

    return 1;
}

/*******************************************************************************
 * planandnavigatexythetalat
 * @brief An example simulation of how a robot would use (x,y,theta) lattice
 *        planning.
 *
 * @param envCfgFilename The environment config file. See
 *                       sbpl/env_examples/nav3d/ for examples
 *******************************************************************************/
int
planandnavigatexythetalat(PlannerType plannerType, char* envCfgFilename, char* motPrimFilename, bool forwardSearch) {
    double allocated_time_secs_foreachplan = 10.0; // in seconds
    double initialEpsilon = 3.0;
    bool bsearchuntilfirstsolution = false;
    bool bforwardsearch = forwardSearch;

    double goaltol_x = 0.001, goaltol_y = 0.001, goaltol_theta = 0.001;

    bool bPrint = false;
    bool bPrintMap = true;

    EnvironmentNAVXYTHETALAT environment_navxythetalat;
    EnvironmentNAVXYTHETALAT trueenvironment_navxythetalat;

    // set the perimeter of the robot
    // it is given with 0, 0, 0 robot ref. point for which planning is done.
    vector<sbpl_2Dpt_t> perimeterptsV;
    sbpl_2Dpt_t pt_m;
    double halfwidth = 0.01;
    double halflength = 0.01;
    pt_m.x = -halflength;
    pt_m.y = -halfwidth;
    perimeterptsV.push_back(pt_m);
    pt_m.x = halflength;
    pt_m.y = -halfwidth;
    perimeterptsV.push_back(pt_m);
    pt_m.x = halflength;
    pt_m.y = halfwidth;
    perimeterptsV.push_back(pt_m);
    pt_m.x = -halflength;
    pt_m.y = halfwidth;
    perimeterptsV.push_back(pt_m);

    //	perimeterptsV.clear();

    // initialize true map from the environment file without perimeter or motion primitives
    if (!trueenvironment_navxythetalat.InitializeEnv(envCfgFilename)) {
        printf("ERROR: InitializeEnv failed\n");
        throw new SBPL_Exception();
    }

    // environment parameters
    int size_x = -1, size_y = -1, num_thetas = -1;
    double startx = -1, starty = -1, starttheta = -1;
    double goalx = -1, goaly = -1, goaltheta = -1;
    double cellsize_m = 0.0, nominalvel_mpersecs = 0.0, timetoturn45degsinplace_secs = 0.0;
    unsigned char obsthresh = 0;
    vector<SBPL_xytheta_mprimitive> motionprimitiveV;

    // additional environment parameters
    unsigned char costinscribed_thresh = 0;
    unsigned char costcircum_thresh = 0;

    // get environment parameters from the true environment
    trueenvironment_navxythetalat.GetEnvParms(&size_x, &size_y, &num_thetas, &startx, &starty, &starttheta, &goalx,
                                              &goaly, &goaltheta, &cellsize_m, &nominalvel_mpersecs,
                                              &timetoturn45degsinplace_secs, &obsthresh, &motionprimitiveV);
    costinscribed_thresh = trueenvironment_navxythetalat.GetEnvParameter("cost_inscribed_thresh");
    costcircum_thresh = trueenvironment_navxythetalat.GetEnvParameter("cost_possibly_circumscribed_thresh");

    // print the map
    if (bPrintMap) {
        printf("true map:\n");
        for (int y = 0; y < size_y; y++) {
            for (int x = 0; x < size_x; x++) {
                printf("%3d ", trueenvironment_navxythetalat.GetMapCost(x, y));
            }
            printf("\n");
        }
        printf("System Pause (return=%d)\n", system("pause"));
    }

    // create an empty map
    unsigned char* map = new unsigned char[size_x * size_y];
    for (int i = 0; i < size_x * size_y; i++) {
        map[i] = 0;
    }

    // check the start and goal obtained from the true environment
    printf("start: %f %f %f, goal: %f %f %f\n", startx, starty, starttheta, goalx, goaly, goaltheta);

    // set robot environment parameters (should be done before initialize function is called)
    if (!environment_navxythetalat.SetEnvParameter("cost_inscribed_thresh", costinscribed_thresh)) {
        printf("ERROR: failed to set parameters\n");
        throw new SBPL_Exception();
    }
    if (!environment_navxythetalat.SetEnvParameter("cost_possibly_circumscribed_thresh", costcircum_thresh)) {
        printf("ERROR: failed to set parameters\n");
        throw new SBPL_Exception();
    }

    // initialize environment (should be called before initializing anything else)
    EnvNAVXYTHETALAT_InitParms params;
    params.startx = startx;
    params.starty = starty;
    params.starttheta = starttheta;
    params.goalx = goalx;
    params.goaly = goaly;
    params.goaltheta = goaltheta;
    params.goaltol_x = goaltol_x;
    params.goaltol_y = goaltol_y;
    params.goaltol_theta = goaltol_theta;
    params.mapdata = map;
    params.numThetas = num_thetas;

    bool envInitialized = environment_navxythetalat.InitializeEnv(size_x, size_y, perimeterptsV, cellsize_m,
                                                                  nominalvel_mpersecs, timetoturn45degsinplace_secs,
                                                                  obsthresh, motPrimFilename, params);

    if (!envInitialized) {
        printf("ERROR: InitializeEnv failed\n");
        throw new SBPL_Exception();
    }

    // set start and goal states
    environment_navxythetalat.SetStart(startx, starty, starttheta);
    environment_navxythetalat.SetGoal(goalx, goaly, goaltheta);

    MDPConfig MDPCfg;

    // initialize MDP info
    if (!environment_navxythetalat.InitializeMDPCfg(&MDPCfg)) {
        printf("ERROR: InitializeMDPCfg failed\n");
        throw new SBPL_Exception();
    }

    // create a planner
    Path solution_stateIDs_V;

    SBPLPlanner* planner = NULL;
    switch (plannerType) {
//TODO add new planners
    case PLANNER_TYPE_ARASTAR:
        printf("Initializing ARAPlanner...\n");
        planner = new ARAPlanner(&environment_navxythetalat, bforwardsearch);
        break;
    case PLANNER_TYPE_ADSTAR:
        printf("Initializing ADPlanner...\n");
        planner = new ADPlanner(&environment_navxythetalat, bforwardsearch);
        break;
    case PLANNER_TYPE_RSTAR:
        printf("Invalid configuration: xytheta environment does not support rstar planner...\n");
        return 0;
    case PLANNER_TYPE_ANASTAR:
        printf("Initializing anaPlanner...\n");
        planner = new anaPlanner(&environment_navxythetalat, bforwardsearch);
        break;

    default:
        printf("Invalid planner type\n");
        break;
    }

    // set the start and goal states for the planner and other search variables
    if (planner->set_start(MDPCfg.startstateid) == 0) {
        printf("ERROR: failed to set start state\n");
        throw new SBPL_Exception();
    }
    if (planner->set_goal(MDPCfg.goalstateid) == 0) {
        printf("ERROR: failed to set goal state\n");
        throw new SBPL_Exception();
    }
    planner->set_initialsolution_eps(initialEpsilon);
    planner->set_search_mode(bsearchuntilfirstsolution);

    // compute sensing as a square surrounding the robot with length twice that of the
    // longest motion primitive
    vector<sbpl_2Dcell_t> sensecells;
    double maxMotPrimLengthSquared = 0.0;
    double maxMotPrimLength = 0.0;
    const EnvNAVXYTHETALATConfig_t* cfg = environment_navxythetalat.GetEnvNavConfig();
    for (int i = 0; i < (int)cfg->mprimV.size(); i++) {
        const SBPL_xytheta_mprimitive& mprim = cfg->mprimV.at(i);
        int dx = mprim.endcell.x;
        int dy = mprim.endcell.y;
        if (dx * dx + dy * dy > maxMotPrimLengthSquared) {
            cout << "Found a longer motion primitive with dx = " << dx
                 << " and dy = " << dy
                 << " from starttheta = " << (int)mprim.starttheta_c
                 << endl;
            maxMotPrimLengthSquared = dx * dx + dy * dy;
        }
    }
    maxMotPrimLength = sqrt((double)maxMotPrimLengthSquared);
    cout << "Maximum motion primitive length: " << maxMotPrimLength << endl;

    int sensingRange = (int)ceil(maxMotPrimLength);
    for (int x = -sensingRange; x <= sensingRange; x++) {
        for (int y = -sensingRange; y <= sensingRange; y++) {
            sensecells.push_back(sbpl_2Dcell_t(x, y));
        }
    }

    // create a file to hold the solution vector
    const char* sol = "sol.txt";
    FILE* fSol = fopen(sol, "w");
    if (fSol == NULL) {
        printf("ERROR: could not open solution file\n");
        throw new SBPL_Exception();
    }

    // print the goal pose
    int goalx_c = CONTXY2DISC(goalx, cellsize_m);
    int goaly_c = CONTXY2DISC(goaly, cellsize_m);
    int goaltheta_c = ContTheta2Disc(goaltheta, num_thetas);
    printf("goal_c: %d %d %d\n", goalx_c, goaly_c, goaltheta_c);

    Path preds_of_changededgesIDV;
    vector<nav2dcell_t> changedcellsV;
    nav2dcell_t nav2dcell;
    vector<sbpl_xy_theta_pt_t> xythetaPath;

    // now comes the main loop
    while (fabs(startx - goalx) > goaltol_x || fabs(starty - goaly) > goaltol_y || fabs(starttheta - goaltheta)
        > goaltol_theta) {
        //simulate sensor data update
        bool bChanges = false;
        preds_of_changededgesIDV.clear();
        changedcellsV.clear();

        // simulate sensing the cells
        for (int i = 0; i < (int)sensecells.size(); i++) {
            int x = CONTXY2DISC(startx, cellsize_m) + sensecells.at(i).x;
            int y = CONTXY2DISC(starty, cellsize_m) + sensecells.at(i).y;

            // ignore if outside the map
            if (x < 0 || x >= size_x || y < 0 || y >= size_y) {
                continue;
            }

            int index = x + y * size_x;
            unsigned char truecost = trueenvironment_navxythetalat.GetMapCost(x, y);
            // update the cell if we haven't seen it before
            if (map[index] != truecost) {
                map[index] = truecost;
                environment_navxythetalat.UpdateCost(x, y, map[index]);
                printf("setting cost[%d][%d] to %d\n", x, y, map[index]);
                bChanges = true;
                // store the changed cells
                nav2dcell.x = x;
                nav2dcell.y = y;
                changedcellsV.push_back(nav2dcell);
            }
        }

        double TimeStarted = clock();

        // if necessary notify the planner of changes to costmap
        if (bChanges) {
            if (dynamic_cast<ARAPlanner*> (planner) != NULL) {
                ((ARAPlanner*)planner)->costs_changed(); //use by ARA* planner (non-incremental)
            }
            else if (dynamic_cast<ADPlanner*> (planner) != NULL) {
                // get the affected states
                environment_navxythetalat.GetPredsofChangedEdges(&changedcellsV, &preds_of_changededgesIDV);
                // let know the incremental planner about them
                //use by AD* planner (incremental)
                ((ADPlanner*)planner)->update_preds_of_changededges(&preds_of_changededgesIDV);
                printf("%d states were affected\n", (int)preds_of_changededgesIDV.size());
            }
        }

        int startx_c = CONTXY2DISC(startx, cellsize_m);
        int starty_c = CONTXY2DISC(starty, cellsize_m);
        int starttheta_c = ContTheta2Disc(starttheta, num_thetas);

        // plan a path
        bool bPlanExists = false;

        printf("new planning...\n");
        bPlanExists = (planner->replan(allocated_time_secs_foreachplan, &solution_stateIDs_V) == 1);
        printf("done with the solution of size=%d and sol. eps=%f\n", (unsigned int)solution_stateIDs_V.size(),
               planner->get_solution_eps());
        environment_navxythetalat.PrintTimeStat(stdout);

        // write the solution to sol.txt
        fprintf(fSol, "plan time=%.5f eps=%.2f\n", (clock() - TimeStarted) / ((double)CLOCKS_PER_SEC),
                planner->get_solution_eps());
        fflush(fSol);

        xythetaPath.clear();
        environment_navxythetalat.ConvertStateIDPathintoXYThetaPath(&solution_stateIDs_V, &xythetaPath);
        printf("actual path (with intermediate poses) size=%d\n", (unsigned int)xythetaPath.size());
        for (unsigned int i = 0; i < xythetaPath.size(); i++) {
            fprintf(fSol, "%.3f %.3f %.3f\n", xythetaPath.at(i).x, xythetaPath.at(i).y, xythetaPath.at(i).theta);
        }
        fprintf(fSol, "*********\n");

        for (int j = 1; j < (int)solution_stateIDs_V.size(); j++) {
            int newx, newy, newtheta = 0;
            environment_navxythetalat.GetCoordFromState(solution_stateIDs_V[j], newx, newy, newtheta);
            fprintf(fSol, "%d %d %d\n", newx, newy, newtheta);
        }
        fflush(fSol);

        // print the map (robot's view of the world and current plan)
        int startindex = startx_c + starty_c * size_x;
        int goalindex = goalx_c + goaly_c * size_x;
        for (int y = 0; bPrintMap && y < size_y; y++) {
            for (int x = 0; x < size_x; x++) {
                int index = x + y * size_x;
                int cost = map[index];
                cost = environment_navxythetalat.GetMapCost(x, y);

                // check to see if it is on the path
                bool bOnthePath = false;
                for (int j = 1; j < (int)solution_stateIDs_V.size(); j++) {
                    int newx, newy, newtheta = 0;
                    environment_navxythetalat.GetCoordFromState(solution_stateIDs_V[j], newx, newy, newtheta);
                    if (x == newx && y == newy) bOnthePath = true;
                }

                if (index != startindex && index != goalindex && !bOnthePath) {
                    printf("%3d ", cost);
                }
                else if (index == startindex) {
                    printf("  X ");
                }
                else if (index == goalindex) {
                    printf("  G ");
                }
                else if (bOnthePath) {
                    printf("  * ");
                }
                else {
                    printf("? ");
                }
            }
            printf("\n");
        }

        // move along the path
        if (bPlanExists && (int)xythetaPath.size() > 1) {
            //get coord of the successor
            int newx, newy, newtheta;

            // move until we move into the end of motion primitive
            environment_navxythetalat.GetCoordFromState(solution_stateIDs_V[1], newx, newy, newtheta);

            printf("moving from %d %d %d to %d %d %d\n", startx_c, starty_c, starttheta_c, newx, newy, newtheta);

            // this check is weak since true configuration does not know the actual perimeter of the robot
            if (!trueenvironment_navxythetalat.IsValidConfiguration(newx, newy, newtheta)) {
                printf("ERROR: robot is commanded to move into an invalid configuration "
                       "according to true environment\n");
                throw new SBPL_Exception();
            }

            // move
            startx = DISCXY2CONT(newx, cellsize_m);
            starty = DISCXY2CONT(newy, cellsize_m);
            starttheta = DiscTheta2Cont(newtheta, num_thetas);

            // update the environment
            int newstartstateID = environment_navxythetalat.SetStart(startx, starty, starttheta);

            // update the planner
            if (planner->set_start(newstartstateID) == 0) {
                printf("ERROR: failed to update robot pose in the planner\n");
                throw new SBPL_Exception();
            }
        }
        else {
            printf("No move is made\n");
        }

        if (bPrint) {
            printf("System Pause (return=%d)\n", system("pause"));
        }
    }

    printf("goal reached!\n");

    fflush(NULL);
    fclose(fSol);

    delete map;
    delete planner;

    return 1;
}

/*******************************************************************************
 * planrobarm - An example of planning a robot arm with six degrees-of-freedom.
 *
 * @param envCfgFilename The environment config file. See
 *                       sbpl/env_examples/robarm for examples
 * @return 1 if the planner successfully found a solution; 0 otherwise
 *******************************************************************************/
int
planrobarm(PlannerType plannerType, char* envCfgFilename, bool forwardSearch) {
    int bRet = 0;
    double allocated_time_secs = 5.0; //in seconds
    MDPConfig MDPCfg;
    bool bforwardsearch = forwardSearch;

    //Initialize Environment (should be called before initializing anything else)
    EnvironmentROBARM environment_robarm;
    if (!environment_robarm.InitializeEnv(envCfgFilename)) {
        printf("ERROR: InitializeEnv failed\n");
        throw new SBPL_Exception();
    }

    //Initialize MDP Info
    if (!environment_robarm.InitializeMDPCfg(&MDPCfg)) {
        printf("ERROR: InitializeMDPCfg failed\n");
        throw new SBPL_Exception();
    }

    //srand(1);

    //plan a path
    Path solution_stateIDs_V;

    SBPLPlanner* planner = NULL;
    switch (plannerType) {
//TODO add new planners
    case PLANNER_TYPE_ARASTAR:
        printf("Initializing ARAPlanner...\n");
        planner = new ARAPlanner(&environment_robarm, bforwardsearch);
        break;
    case PLANNER_TYPE_ADSTAR:
        printf("Initializing ADPlanner...\n");
        planner = new ADPlanner(&environment_robarm, bforwardsearch);
        break;
    case PLANNER_TYPE_RSTAR:
        printf("Initializing RSTARPlanner...\n");
        planner = new RSTARPlanner(&environment_robarm, bforwardsearch);
        break;
    case PLANNER_TYPE_ANASTAR:
        printf("Initializing anaPlanner...\n");
        planner = new anaPlanner(&environment_robarm, bforwardsearch);
        break;

    default:
        printf("Invalid planner type\n");
        break;
    }

    if (planner->set_start(MDPCfg.startstateid) == 0) {
        printf("ERROR: failed to set start state\n");
        throw new SBPL_Exception();
    }

    if (planner->set_goal(MDPCfg.goalstateid) == 0) {
        printf("ERROR: failed to set goal state\n");
        throw new SBPL_Exception();
    }

    printf("start planning...\n");
    bRet = planner->replan(allocated_time_secs, &solution_stateIDs_V);
    printf("done planning\n");
    cout << "size of solution=" << solution_stateIDs_V.size() << endl;

    const char* sol = "sol.txt";
    FILE* fSol = fopen(sol, "w");
    if (fSol == NULL) {
        printf("ERROR: could not open solution file\n");
        throw new SBPL_Exception();
    }
    for (unsigned int i = 0; i < solution_stateIDs_V.size(); i++) {
        environment_robarm.PrintState(solution_stateIDs_V[i], true, fSol);
    }
    fclose(fSol);

    //print a path
    if (bRet) {
        //print the solution
        printf("Solution is found\n");
    }
    else
        printf("Solution does not exist\n");

    fflush(NULL);

    delete planner;

    return bRet;
}






SBPLPlanner*
getPlanner(PlannerType plannerType) {
    assert(environment);

    bool bforwardsearch = !backwardSearch;

    switch (plannerType) {
        case INVALID_PLANNER_TYPE:
            printf("ERROR: tried to build an INVALID type planner\n");
            throw new SBPL_Exception();
            return nullptr;

        case PLANNER_TYPE_ARASTAR:
            return new ARAPlanner(environment, bforwardsearch);
        case PLANNER_TYPE_ADSTAR:
            return new ADPlanner(environment, bforwardsearch);
        case PLANNER_TYPE_RSTAR:
            return new RSTARPlanner(environment, bforwardsearch);
        case PLANNER_TYPE_ANASTAR:
            return new anaPlanner(environment, bforwardsearch);
        case PLANNER_TYPE_AASTAR:
            return new AAPlanner(environment, bforwardsearch);
        case PLANNER_TYPE_ASTAR:
            return new ASTARPlanner(environment, backwardSearch);
        case PLANNER_TYPE_D:
            return new DPlanner(environment, bforwardsearch);

        //TODO add new planners
        case PLANNER_TYPE_PPCP:
        case PLANNER_TYPE_VI:
        default:
            printf("Invalid planner type\n");
            throw new SBPL_Exception();
            return nullptr;
    }
}

DiscreteSpaceInformation*
getSearchSpace(EnvironmentType environmentType) {
    switch(environmentType){
        case INVALID_ENV_TYPE:
            printf("ERROR: tried to build an INVALID type environment\n");
            throw new SBPL_Exception();
            return nullptr;

        case ENV_TYPE_2D:
            return new EnvironmentNAV2D();
        case ENV_TYPE_2DUU:
        case ENV_TYPE_ROBARM:
        case ENV_TYPE_XYTHETA:
        case ENV_TYPE_XYTHETAMLEV:
        default:
            printf("Environment type not yet exposed\n");
            throw new SBPL_Exception();
            return nullptr;
    }
}

void
setStartGoal(MDPConfig &instancePair) {
    // TODO: Handle instance change
    assert(planner);
    if (planner->set_start(instancePair.startstateid) == 0) {
        printf("ERROR: failed to set start state\n");
        throw new SBPL_Exception();
    }
    if (planner->set_goal(instancePair.goalstateid) == 0) {
        printf("ERROR: failed to set goal state\n");
        throw new SBPL_Exception();
    }
}


void
loadEnvironment(DiscreteSpaceInformation *space, char* path){

    if (!space->InitializeEnv(path)) {
        printf("ERROR: InitializeEnv failed\n");
        throw new SBPL_Exception();
    }
}



/**
 * Generates a randomized seed based on run and iteration numbers
 * \note The generation of random numbers is too important to be left to chance
 */
Seed
runItseed(){
    Seed n;
    n = (Seed) env_id;

    n += n<<7;
    n ^= n<<11;

    n += (Seed) run_id;

    n ^= n<<5;
    n += n<<3;

    n += (Seed) it_id;

    n ^= n<<15;
    n += n<<2;
    return n;
}



void
testEnvironment() {
    assert(planner);
    assert(environment);

    MDPConfig instancePair;
    it_id = it_start;
    Seed seed = runItseed();
    if(environment->generateRandomProblem(&instancePair, seed, maxTries)){
        // Config
        // TODO: allow testing multiple planner configurations (weight, lookahead, etc..)

        // Compute first plan
        // ------------------
        printf("\n\nInstance solving (Env:%d, Run:%d, It:%d)->(%u)",
               env_id, run_id, it_id, seed);
        printf("\n================\n");
        // Set instance
        setStartGoal(instancePair);
        // Plan (1st time)
        Path path;
        planner->replan(givenTime, &path);

        for(it_id=it_start+1; it_id<=it_end; it_id++){
            // Replanning
            // ----------
            seed = runItseed();
            printf("\nReplanning setup (Env:%d, Run:%d, It:%d)->(%u)",
                   env_id, run_id, it_id, seed);
            printf("\n----------------\n");
            if(environment->generateRandomStart(&instancePair, seed, maxTries)){
                // Modify the start
                setStartGoal(instancePair);
                environment->modifyEnvironment(seed, 0.1);
                // Replan
                planner->replan(givenTime, &path);
            }
            else{
                // Stop at this iteration
                printf("  Stopping at iteration %d on run %d on random environment %d."
                       " No new instance could be generated\n\n",
                       it_id,
                       run_id,
                       env_id
              );
                return;
            }
        }
    }
    else
        printf("  Skipping run %d on random environment %d."
               " No problem instance could be generated\n",
               run_id,
               env_id
              );
}



// Main
// ====
int
main(int argc, char *argv[]) {

    if(parseInput(argc, argv))
        exit(1);

    cout << "Options:" << endl;

    cout << "  Environment:         " << EnvironmentTypeToStr(selectedEnvironment) << endl;
    cout << "  Planner:             " << PlannerTypeToStr(selectedPlanner) << endl;
    cout << "  Random environments: " << "[" << env_random_start << ", " << env_random_end << ")" << endl;
    cout << "  Runs:                " << "[" << run_start << ", " << run_end << ")" << endl;
    cout << "  Iteration:           " << "[" << it_start << ", " << it_end << "]" << endl;


    // Benchmark
    // =========

    if(selectedEnvironment == INVALID_ENV_TYPE) {
        printf("Invalid environment selected\n");
        exit(2);
     }

    // Random environments
    // -------------------
    for(env_id=env_random_start; env_id<env_random_end; env_id++){
        printf("\n\n\n\n\n");
        // Initialize environment
        // ----------------------
        // TODO: Use _env_settings_
        printf("Building Environment[%d]...\n", env_id);
        assert(!environment);
        environment = getSearchSpace(selectedEnvironment);
        assert(environment);
        environment->generateRandomEnvironment((Seed)env_id);

        // Initialize planner (for this environment :/).
        // ------------------
        // TODO: planners should be able to swap environments and cfgs
        // TODO: Use _planner_settings_
        printf("Building Planner for Env[%d]...\n", env_id);
        assert(!planner);
        planner = getPlanner(selectedPlanner);
        assert(planner);

        for(run_id=run_start; run_id<run_end; run_id++)
            testEnvironment();

        DELETE(planner);
        DELETE(environment);
    }

    // Given environments
    // ------------------
    // TODO: same loop on given environments
}
