// CPLEX
#include <ilcplex/ilocplex.h>

#include <iostream>
#include <stdlib.h>

using namespace std;

// CPLEX STL macros
ILOSTLBEGIN

// Subtour elimination callback (lazy constraint)
// This callback inspects the current fractional/integer solution of the y variables
// and detects cycles (subtours) in each bus route. When a subtour (cycle that
// does not include the final mandatory station) is found, it adds a cut to
// eliminate that subtour: sum(y_{k->j} for (k,j) in cycle) <= |cycle| - 1
// Arguments:
//  y        : 3D array y[i][j][k] indicating arc j->k used by bus i
//  Stations : total number of stations
//  nBuses   : number of buses
//  N        : number of mandatory stops (first N entries are mandatory)
//  C        : number of passengers (not used in callback but kept for signature)
ILOLAZYCONSTRAINTCALLBACK5(SubtourEliminationCallback, IloArray<IloArray<IloBoolVarArray>>, y, int, Stations, int, nBuses, int, N, int, C) {
    int i, j, k, p, n;
    vector<int> sta; // list of optional/mandatory stations yet to explore

    // Iterate over each bus and detect subtours in its current route
    for (i = 0; i < nBuses; i++) {
        // determine the stations visited in this bus according to current solution
        n = N; // start counting with mandatory stops
        sta.clear();
        // include mandatory stops 1..N-1 (indexing assumes 0 is origin)
        for (j = 1; j < N; j++) {
            sta.push_back(j);
        }
        // check optional stops (those indexed >= N) and add ones used by this bus
        for (j = N; j < Stations; j++) {
            for (p = 0; p < Stations; p++) {
                if (((int)getValue(y[i][j][p] + 0.001) == 1)) {
                    sta.push_back(j);
                    n++;
                    break;
                }
            }
        }

        // Arrays used for DFS-like traversal to identify cycles
        bool *seen = new bool[Stations];
        int *visited = new int[n + 1];

        memset(seen, false, sizeof(bool) * Stations);
        memset(visited, 0, sizeof(int) * n + 1);
        int length = 0;
        k = 0;          // start from station 0 (origin)
        seen[0] = true; // mark origin visited

        // Try to trace paths following arcs that are 'selected' in current solution
        while (true) {
            for (j = 0; j < Stations; j++) {
                if (k != j) {
                    // if arc k->j is active in current solution
                    if (((int)getValue(y[i][k][j] + 0.001) == 1)) {
                        // record visiting this station
                        visited[++length] = j;
                        if (seen[j]) {
                            // Found a repeated node -> a cycle (subtour) candidate
                            break;
                        }
                        seen[j] = true;
                        k = j;  // continue from this new node
                        j = -1; // reset inner loop to start scanning from 0

                        // remove visited station from the 'sta' list (stations to explore)
                        for (p = 0; p < (int)sta.size(); p++) {
                            if (sta[p] == k) {
                                sta.erase(sta.begin() + p);
                                break;
                            }
                        }
                    }
                }
            }

            // heuristics to decide whether to continue exploring or reset search
            if (length < n - 1 && visited[0] != visited[length]) {
                // if there are still stations not covered, restart search from first remaining
                if (n == N && length >= N - 1) {
                    break;
                }
                k = sta[0];
                length = 0;
                memset(seen, false, sizeof(bool) * Stations);
                memset(visited, 0, sizeof(int) * n + 1);
                visited[length] = k;
                seen[k] = true;
            }
            else if (length < n && length > 0 || sta.size() == 0) {
                // finished exploring a connected component or no more stations left
                break;
            }
            else if (length == 0) {
                // advance starting node if no progress was made
                if (k < n - 1) {
                    k++;
                }
                else {
                    length = n; // mark finished
                    break;
                }
                length = 0;
                memset(seen, false, sizeof(bool) * Stations);
                memset(visited, 0, sizeof(int) * n + 1);
                visited[length] = k;
                seen[k] = true;
            }
        }

        // If a subtour (cycle) smaller than the full set was detected, add a cut
        int l0 = length;
        if (length > 0 && length < n - 1) {
            length = l0;
            IloExpr clique(getEnv());
            int firstjobinloop = visited[length];
            int newlen = 0;
            while (true) {
                j = visited[length];
                k = visited[--length];
                // accumulate y variables that form the cycle
                clique += y[i][k][j];
                newlen++;
                if (k == firstjobinloop) {
                    break;
                }
            }
            // add the subtour elimination cut: sum(y in cycle) <= |cycle| - 1
            add(clique <= newlen - 1).end();
            clique.end();
        }
        delete[] seen;
        delete[] visited;
    }

    return;
}

int main() {
    // Instance counter used only for naming the instance file (useful when
    // generating multiple instances in a single run). Here it's fixed to 0.
    int instance = 0;
    ofstream inst("data/output/Instance_S" + to_string(instance) + ".txt");
    inst << "---------- Weight factors of the objective function -------- \n" << endl;

    // ------------------------------------------------------------------
    // Objective weights: these scale the different components of the cost
    // c1: bus travel + dwell time, c2: passenger walking, c3: arrival deviation
    // ------------------------------------------------------------------
    float c1 = 0.25f; // weight for bus travel time
    inst << "c1: " << c1 << " \t (For travel-time of the buses)" << endl;
    float c2 = 0.35f; // weight for passenger walking time
    inst << "c2: " << c2 << " \t (For walking time of the passengers)" << endl;
    float c3 = 1 - c1 - c2; // remaining weight for arrival time deviation
    inst << "c3: " << c3 << " \t (For the absolute difference in desired arrival time and actual arrival time of the passengers)" << endl;
    float lvsea = 1; // linearization/weight for early/late (used in waiting term)

    inst << "\n---------------------- Parameters -------------------------- \n" << endl;

    // -------------------- Problem parameters -----------------------
    // These are example (default) values. You can change them or load
    // from a configuration file if desired.
    const int nBuses = 2; // number of available buses
    inst << "Number of buses: " << nBuses << endl;
    const int N = 5; // number of mandatory stations (these must be visited in order)
    inst << "Number of mandatory bus stops: " << N << endl;
    const int M = 3; // number of optional stops in each inter-mandatory cluster
    inst << "Number optional bus stops per cluster: " << M << " \n --> One cluster between each mandatory stop: " + to_string((N - 1) * M) + " optional stops in total" << endl;
    const int Stations = (N - 1) * M + N; // total number of stations (mandatory + optional)
    inst << "Total number of bus stops: " << Stations << endl;
    const int C = 20; // number of passenger requests
    inst << "Number of passenger requests: " << C << endl;

    // -------------------- Physical / time parameters ----------------
    const int bCapacity = 15; // bus capacity (passengers)
    inst << "Bus capacity: " << bCapacity << " passengers" << endl;
    const float pspeed = 1.0f;         // passenger walking speed (m/s)
    const float bspeed = 40.0f / 3.6f; // bus speed (40 km/h converted to m/s)
    const int delta = 30;              // per-arc acceleration/deceleration overhead (s)
    inst << "Acceleration and deceleration time: " << delta << " seconds" << endl;
    const int tau = 5; // dwell time per boarding/alighting passenger (s)
    inst << "Dwell time per passenger at a bus stop: " << tau << " seconds" << endl;
    const int d = 20 * 60; // maximum allowed walking time (seconds)
    inst << "Maximum walking time for any passenger: " << d << " seconds" << endl;
    const int d_time1 = 15 * 60; // max early arrival tolerance (s)
    const int d_time2 = 5 * 60;  // max late arrival tolerance (s)
    inst << "Maximum amount of time a passenger can arrive too early: " << d_time1 << " seconds" << endl;
    inst << "Maximum amount of time a passenger can arrive too late: " << d_time2 << " seconds" << endl;

    // -------------------- Read input data (format expectations) -------
    // - data/input/passengers.txt : C rows of "x y" coordinates (meters or km consistent with speed)
    // - data/input/mandatory.txt : N rows of mandatory stop coordinates
    // - data/input/optional.txt  : (N-1)*M rows with optional stop coordinates
    // - data/input/arrivals.txt  : C desired arrival times (seconds)
    double passengers[C][2];
    ifstream filep("data/input/passengers.txt"); // passenger coordinates
    int i = 0;
    while (i < C) {
        // each line: two floating point values (x y)
        filep >> passengers[i][0] >> passengers[i][1];
        i++;
    }

    double mandatory[N][2];
    ifstream filem("data/input/mandatory.txt");
    i = 0;
    while (i < N) {
        filem >> mandatory[i][0] >> mandatory[i][1];
        i++;
    }

    double optional[(N - 1) * M][2];
    ifstream fileo("data/input/optional.txt");
    i = 0;
    while (i < (N - 1) * M) {
        fileo >> optional[i][0] >> optional[i][1];
        i++;
    }

    // Arrival times (desired) for passengers, in seconds
    double arrivals[C];
    ifstream filea("data/input/arrivals.txt");
    inst << endl
         << "Desired arrival times of the passengers in seconds: " << endl;
    i = 0;
    while (i < C) {
        filea >> arrivals[i];
        inst << "p_" << i + 1 << ": \t" << arrivals[i] << endl;
        i++;
    }

    // -------------------- Compute travel times -----------------------
    // traveltimep[p][j] : walking time (seconds) for passenger p to station j
    // traveltimes[j][k] : bus travel time (seconds) between station j and k
    double traveltimep[C][Stations];
    double traveltimes[Stations][Stations];

    // Compute walking times using Manhattan distance scaled to seconds.
    for (int i = 0; i < C; i++) {
        for (int j = 0; j < N; j++) {
            traveltimep[i][j] = d / 20.0 / 60.0 * (abs(passengers[i][0] - mandatory[j][0]) + abs(passengers[i][1] - mandatory[j][1])) * 1000 / pspeed;
        }
        for (int j = N; j < Stations; j++) {
            traveltimep[i][j] = d / 20.0 / 60.0 * (abs(passengers[i][0] - optional[j - N][0]) + abs(passengers[i][1] - optional[j - N][1])) * 1000 / pspeed;
        }
    }

    // Pretty-print walking matrix (filtered by threshold/min heuristics)
    inst << "\nWalking time between passengers and bus stops in seconds: "
         << "\npassengers correspond with the rows, bus stops correspond with the columns \nthe mandatory stops are listed first, then the optional stops are listed" << endl;
    for (int i = 0; i < C; i++) {
        double minp = 10000000000;
        for (int k = 0; k < N; k++) {
            if (minp > traveltimep[i][k]) {
                minp = traveltimep[i][k];
            }
        }
        for (int j = 0; j < Stations; j++) {
            // show '/' for stations that are too far or dominated by a closer mandatory stop
            if (traveltimep[i][j] > d || (j > N && traveltimep[i][j] > minp)) {
                inst << "/\t";
            }
            else {
                inst << int(traveltimep[i][j]) << "\t";
            }
        }
        inst << endl;
    }

    // Bus travel times (Manhattan) between all station pairs
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            traveltimes[i][j] = (abs(mandatory[i][0] - mandatory[j][0]) + abs(mandatory[i][1] - mandatory[j][1])) * 1000 / bspeed;
        }
        for (int j = N; j < Stations; j++) {
            traveltimes[i][j] = (abs(mandatory[i][0] - optional[j - N][0]) + abs(mandatory[i][1] - optional[j - N][1])) * 1000 / bspeed;
        }
    }
    for (int i = N; i < Stations; i++) {
        for (int j = 0; j < N; j++) {
            traveltimes[i][j] = (abs(optional[i - N][0] - mandatory[j][0]) + abs(optional[i - N][1] - mandatory[j][1])) * 1000 / bspeed;
        }
        for (int j = N; j < Stations; j++) {
            traveltimes[i][j] = (abs(optional[i - N][0] - optional[j - N][0]) + abs(optional[i - N][1] - optional[j - N][1])) * 1000 / bspeed;
        }
    }
    inst << "\nTravel time between bus stops in seconds: " << endl;
    for (int i = 0; i < Stations; i++) {
        for (int j = 0; j < Stations; j++) {
            inst << int(traveltimes[i][j]) << "\t";
        }
        inst << endl;
    }
    inst.close();

    // -------------------- Build index sets (simple integer arrays) ------
    // These arrays are convenience iterables later used with range-based for loops.
    int I[nBuses];
    int J[Stations];
    int F[N];
    int O[(N - 1) * M];
    int P[C];
    for (int i = 0; i < nBuses; i++) {
        I[i] = i;
    }
    for (int i = 0; i < Stations; i++) {
        J[i] = i;
    }
    for (int i = 0; i < N; i++) {
        F[i] = i;
    }
    for (int i = N; i < Stations; i++) {
        O[i - N] = i;
    }
    for (int i = 0; i < C; i++) {
        P[i] = i;
    }

    double elapsed_time;
    clock_t start_time;
    start_time = clock();

    // -------------------- Big-M and time bounds -------------------------
    // Earliest and latest requested arrival times are used to compute a
    // safe upper bound (Big-M) for time linking constraints below.
    double EA = 10000000000, LA = -1;
    for (int p = 0; p < C; p++) {
        if (EA > arrivals[p]) {
            EA = arrivals[p];
        }
        if (LA < arrivals[p]) {
            LA = arrivals[p];
        }
    }
    // M_1 is a conservative upper bound on time differences used in logical constraints
    int M_1 = (int)(LA + d_time2 - EA + d_time1) + 1;

    // *************************************************|| MODEL || ***************************************************************
    // Create CPLEX environment and model container
    IloEnv env;
    IloModel model(env);

    // -------------------- Decision variables ---------------------------
    // x[p][i][j] : binary, passenger p assigned to bus i at station j
    // y[i][j][k] : binary, bus i travels from station j to k
    // D[i]       : continuous, bus start time
    // A[p]       : continuous, actual arrival time for passenger p
    // q_early/q_late : non-negative slacks for early/late arrival
    IloArray<IloArray<IloBoolVarArray>> x(env);
    for (int p = 0; p < C; p++) {
        IloArray<IloBoolVarArray> x1(env);
        for (int i = 0; i < nBuses; i++) {
            x1.add(IloBoolVarArray(env, Stations));
        }
        x.add(x1);
    }

    IloArray<IloArray<IloBoolVarArray>> y(env);
    for (int i = 0; i < nBuses; i++) {
        IloArray<IloBoolVarArray> y1(env);
        for (int j = 0; j < Stations; j++) {
            y1.add(IloBoolVarArray(env, Stations));
        }
        y.add(y1);
    }

    IloNumVarArray D(env); // bus start times
    for (int i = 0; i < nBuses; i++) {
        D.add(IloNumVar(env, 0, INFINITY));
    }

    IloNumVarArray A(env); // passenger actual arrival times
    for (int i = 0; i < C; i++) {
        A.add(IloNumVar(env, 0, INFINITY));
    }
    IloNumVarArray q_early(env);
    for (int i = 0; i < C; i++) {
        q_early.add(IloNumVar(env, 0, INFINITY));
    }
    IloNumVarArray q_late(env);
    for (int i = 0; i < C; i++) {
        q_late.add(IloNumVar(env, 0, INFINITY));
    }

    // -------------------- Objective construction -----------------------
    // TravelTime: bus travel and dwell times
    IloExpr TravelTime(env);
    for (int i = 0; i < nBuses; i++) {
        for (int j = 0; j < Stations; j++) {
            for (int k = 0; k < Stations; k++) {
                TravelTime += y[i][j][k] * (traveltimes[j][k] + delta);
            }
            for (int p = 0; p < C; p++) {
                TravelTime += x[p][i][j] * tau; // dwell contributions
            }
        }
    }

    // WalkingTime: passenger walking time to assigned station
    IloExpr WalkingTime(env);
    for (int p = 0; p < C; p++) {
        for (int i = 0; i < nBuses; i++) {
            for (int j = 0; j < Stations; j++) {
                WalkingTime += x[p][i][j] * traveltimep[p][j];
            }
        }
    }

    // WaitingTime: penalty for arriving early/late (linearized using q_early/q_late)
    IloExpr WaitingTime(env);
    for (int p = 0; p < C; p++) {
        WaitingTime += lvsea * q_late[p] + q_early[p];
    }

    // Full objective: weighted sum
    model.add(IloMinimize(env, c1 * TravelTime + c2 * WalkingTime + c3 * WaitingTime));
    TravelTime.end();
    WaitingTime.end();
    WalkingTime.end();

    // -------------------- CONSTRAINTS ------------------------
    
    // ROUTING ******

    // one arc at most going out the optional stops
    for (const auto &i : I) {
        for (const auto &j : O) {
            IloExpr sumY1(env);
            for (int k = 0; k < Stations; k++) {
                sumY1 += y[i][j][k];
            }
            // before clearing sumConst expr, add it to model
            model.add(sumY1 <= 1);
            sumY1.end();
        }
    }

    // exactly one arc going out of the mandatory stops (or one going in for the last stop)
    for (const auto &i : I) {
        for (const auto &j : F) {
            if (j == N - 1) {
                IloExpr sumY1(env);
                for (int k = 0; k < Stations; k++) {
                    sumY1 += y[i][k][j];
                }
                model.add(sumY1 == 1);
                sumY1.end();
            }
            else {
                IloExpr sumY2(env);
                for (int k = 0; k < Stations; k++) {
                    sumY2 += y[i][j][k];
                }
                model.add(sumY2 == 1);
                sumY2.end();
            }
        }
    }

    // No loops
    for (const auto &i : I) {
        for (const auto &j : J) {
            model.add(y[i][j][j] == 0);
        }
    }

    ///*
    // edges
    for (const auto &i : I) {
        IloExpr sumY1(env);
        IloExpr sumY2(env);
        for (const auto &j : J) {
            sumY1 += y[i][N - 1][j];
            sumY2 += y[i][j][0];
        }
        model.add(sumY1 == 0);
        model.add(sumY2 == 0);
        sumY1.end();
        sumY2.end();
    }
    //*/

    // What goes in, goes out (except for boundary nodes)
    for (const auto &i : I) {
        for (const auto &j : J) {
            if (j != 0 && j != N - 1) {
                IloExpr sumY1(env);
                IloExpr sumY2(env);
                for (int k = 0; k < Stations; k++) {
                    sumY1 += y[i][j][k];
                    sumY2 += y[i][k][j];
                }
                model.add(sumY1 == sumY2);
                sumY1.end();
                sumY2.end();
            }
        }
    }

    // ASSIGNMENT ******

    // All passengers need a bus
    for (const auto &p : P) {
        IloExpr sumX1(env);
        for (const auto &i : I) {
            for (const auto &j : J) {
                sumX1 += x[p][i][j];
            }
        }
        model.add(sumX1 == 1);
        sumX1.end();
    }

    // DEFINITONS *******

    // Binding constraints
    for (const auto &i : I) {
        for (const auto &j : O) {
            IloExpr sumX(env);
            IloExpr sumY(env);
            for (const auto &p : P) {
                sumX += x[p][i][j];
            }

            for (int l = 0; l < Stations; l++) {
                sumY += y[i][j][l];
            }

            model.add(sumX <= bCapacity * sumY);
            // model.add(sumY <= bCapacity * sumX);
            sumX.end();
            sumY.end();
        }
    }

    // Defining Ap variable
    for (const auto &i : I) {
        for (const auto &p : P) {
            IloExpr sumX1(env);
            IloExpr sumX2(env);
            IloExpr sumY(env);

            for (const auto &j : J) {
                sumX2 += x[p][i][j];
                for (const auto &k : J) {
                    sumY += (traveltimes[j][k] + delta) * y[i][j][k];
                }

                for (const auto &c : P) {
                    sumX1 += x[c][i][j] * tau;
                }
            }
            model.add(sumY + sumX1 + D[i] - A[p] <= (1 - sumX2) * M_1);
            model.add(-sumY - sumX1 - D[i] + A[p] <= (1 - sumX2) * M_1);
            sumX1.end();
            sumX2.end();
            sumY.end();
        }
    }

    // CAPACITIES *****

    // Walking threshold
    for (const auto &p : P) {
        IloExpr sumX1(env);
        for (const auto &i : I) {
            for (const auto &j : J) {
                sumX1 += (traveltimep[p][j]) * x[p][i][j];
            }
        }
        double mindist = 10000000000;
        for (const auto &j : J) {
            if (traveltimep[p][j] < mindist) {
                mindist = traveltimep[p][j];
            }
        }

        // model.add(sumX1 >= mindist);
        model.add(sumX1 <= d);
        sumX1.end();
    }

    /*
    for (const auto& p : P) {
            model.add( A[p]<= arrivals[p]);
    }
    //*/

    //  Arrival time
    for (const auto &p : P) {
        model.add(q_early[p] <= d_time1);
        model.add(q_late[p] <= d_time2);
        model.add(arrivals[p] - A[p] + q_late[p] - q_early[p] == 0);
    }

    // Bus capacity
    for (const auto &i : I) {
        IloExpr sumX(env);
        for (const auto &j : J) {
            for (const auto &p : P) {
                sumX += x[p][i][j];
            }
        }
        model.add(sumX <= bCapacity);
        sumX.end();
    }

    // SOLVE ----------------------------------------------------
    // Configure and run the CPLEX solver. A few solver params are shown
    // commented-out below for easy tuning (gap, emphasis, time limit, etc.).
    IloCplex cplex(model);
    // Suppress solver output if desired:
    // cplex.setOut(env.getNullStream());
    // Stop earlier by setting acceptable optimality gap (e.g., 1%):
    // cplex.setParam(IloCplex::EpGap, 0.01);
    // Change solver focus/emphasis (3 is aggressive towards feasible solutions):
    // cplex.setParam(IloCplex::MIPEmphasis, 3);
    // Set number of threads (adjust to your machine):
    cplex.setParam(cplex.Threads, 12);

    // Install lazy constraint callback for subtour elimination. This callback
    // dynamically adds constraints during the branch-and-cut search to remove
    // infeasible subtours found in intermediate solutions.
    cplex.use(SubtourEliminationCallback(env, y, Stations, nBuses, N, C));

    // Optional time limit (seconds):
    // cplex.setParam(cplex.TiLim,3610);

    // Solve the model (blocking call). The callback will be invoked
    // internally by CPLEX when nodes are explored.
    cplex.solve();

    // Record elapsed CPU time for profiling/logging
    elapsed_time = (double)(clock() - start_time) / CLK_TCK;
    std::cout << "Computational time: " << elapsed_time << " seconds\n" << endl;

    // Retrieve objective value after solve
    double objval = cplex.getObjValue();

    // ------------------ Extract solution values ------------------
    // Local containers to store the (binary/integer) solution values read from CPLEX.
    // xsol[p][i][j] : passenger p assigned to bus i at station j (0/1)
    int xsol[C][nBuses][Stations];
    ofstream txt_xsol("data/output/xsol.txt");

    // ysol[i][j][k] : bus i travels from station j to station k (0/1)
    int ysol[nBuses][Stations][Stations];
    ofstream txt_ysol("data/output/ysol.txt");

    // Dsol: bus start times, Asol: passenger actual arrival times
    double Dsol[nBuses];
    ofstream txt_Dsol("data/output/Dsol.txt");
    double Asol[C];
    ofstream txt_Asol("data/output/Asol.txt");

    // Populate the local arrays by querying CPLEX. A small numeric tolerance
    // is used when converting floating point solver values to integers.
    for (int i = 0; i < nBuses; i++) {
        // continuous start time
        Dsol[i] = cplex.getValue(D[i]);
        txt_Dsol << Dsol[i] << endl;

        txt_xsol << "Bus " << i << endl;
        txt_ysol << "Bus " << i << endl;
        for (int j = 0; j < Stations; j++) {
            for (int k = 0; k < Stations; k++) {
                // routing decisions
                ysol[i][j][k] = (int)cplex.getValue(y[i][j][k] + 0.001);
                txt_ysol << ysol[i][j][k] << " ";
            }
            for (int p = 0; p < C; p++) {
                // passenger assignments
                xsol[p][i][j] = (int)cplex.getValue(x[p][i][j] + 0.001);
                txt_xsol << xsol[p][i][j] << " ";
            }
            txt_xsol << endl;
            txt_ysol << endl;
        }
        txt_xsol << "end" << endl;
        txt_ysol << "end" << endl;
    }

    // Passenger arrival times
    for (int p = 0; p < C; p++) {
        Asol[p] = cplex.getValue(A[p]);
        txt_Asol << Asol[p] << endl;
    }

    // Close files used to dump solver outputs
    txt_xsol.close();
    txt_ysol.close();
    txt_Asol.close();
    txt_Dsol.close();

    env.end();

    cout << "\n********************************\nSolution is : " << objval <<  "s\n********************************\n" << endl;

    cout << " Bus starting times in seconds -------- " << endl;

    for (int i = 0; i < nBuses; i++) {
        cout << "Bus " << i << "--> " << Dsol[i] <<  "s" << endl;
    }

    cout << endl;

    cout << " Early times in seconds  --------- " << endl;
    double AD = 0;
    for (int p = 0; p < C; p++) {
        cout << "Passeger " << p << " --> " << arrivals[p] - Asol[p] << "s" <<  endl;
        AD += abs(arrivals[p] - Asol[p]);
    }

    // ------------------ Write visits (timetable) ------------------
    // For each bus, follow the selected arcs (ysol) starting from node 0
    // until reaching the final mandatory stop (index N-1). At each station
    // record the time when the bus arrives.
    ofstream txt_visits("data/output/visits.txt");
    int j = 0;
    double travel = 0, travelp = 0;
    for (int i = 0; i < nBuses; i++) {
        j = 0;            // start at station 0
        travel = Dsol[i]; // initial time is the bus start time
        txt_visits << "Bus " << i + 1 << endl;
        // follow the route until the last mandatory station
        while (j != N - 1) {
            for (int k = 0; k < Stations; k++) {
                if (ysol[i][j][k] == 1) {
                    // record visit: station index and arrival time
                    txt_visits << j << " " << travel << endl;
                    // add travel time for arc and dwell times from passengers
                    travel += (traveltimes[j][k] + delta);
                    for (int p = 0; p < C; p++) {
                        travel += xsol[p][i][j] * tau;
                    }
                    j = k; // advance to next station
                }
            }
        }
        // final station visit
        txt_visits << j << " " << travel << endl;
        txt_visits << "end" << endl;
    }
    txt_visits.close();

    // ------------------ Walking output ------------------
    // For each passenger, find the assigned station and write the walking time
    double W = 0;
    ofstream txt_walking("data/output/walking.txt");
    for (int p = 0; p < C; p++) {
        for (int j = 0; j < Stations; j++) {
            for (int i = 0; i < nBuses; i++) {
                if (xsol[p][i][j] == 1) {
                    txt_walking << j << " " << traveltimep[p][j] << endl;
                    W += traveltimep[p][j];
                }
            }
        }
    }
    txt_walking.close();

    // ------------------ Aggregate passenger trip times ------------------
    // travelp is the total bus route time experienced by passengers. For each
    // passenger, find their boarding station and sum the bus travel and dwell
    // times from that station to the final mandatory stop.
    for (int p = 0; p < C; p++) {
        for (int i = 0; i < nBuses; i++) {
            for (int l = 0; l < Stations; l++) {
                if (xsol[p][i][l] == 1) {
                    j = l;
                    travel = 0;
                    while (j != N - 1) {
                        for (int k = 0; k < Stations; k++) {
                            if (ysol[i][j][k] == 1) {
                                travel += (traveltimes[j][k] + delta);
                                for (int p = 0; p < C; p++) {
                                    travel += xsol[p][i][j] * tau;
                                }
                                j = k;
                            }
                        }
                    }
                }
            }
        }
        travelp += travel;
    }

    // Print aggregated metrics: absolute deviation (AD), total passenger travel time (T), total walking (W)
    cout << "\n----------------------\n- Absolute deviation (AD): " << AD <<  "s" << endl;
    cout << "- Total passenger travel time (T): " << travelp <<  "s" << endl;
    cout << "- Total walking (W): " << W <<  "s" << endl;
}