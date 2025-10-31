// CPLEX
#include <ilcplex/ilocplex.h>

#include <iostream>
#include <stdlib.h>

void Gen(int N, int M, int C);
using namespace std;

// for cplex
ILOSTLBEGIN
ILOLAZYCONSTRAINTCALLBACK5(SubtourEliminationCallback, IloArray<IloArray<IloBoolVarArray>>, y, int, Stations, int, nBuses, int, N, int, C) {
    int i, j, k, p, n;
    vector<int> sta;

    for (i = 0; i < nBuses; i++) {
        // std::cout << "---------------------------- bus " << i + 1 << endl;
        // determine the stations visited in each bus
        n = N;
        sta.clear();
        for (j = 1; j < N; j++) {
            sta.push_back(j);
        }
        for (j = N; j < Stations; j++) {
            for (p = 0; p < Stations; p++) {
                if (((int)getValue(y[i][j][p] + 0.001) == 1)) {
                    sta.push_back(j);
                    n++;
                    break;
                }
            }
        }

        // cout << "----------------------------------------  n: " << n << endl;
        /*
        cout << " y-variable \n";
        for (j = 0; j < Stations; j++) {
                for (k = 0; k < Stations; k++) {
                        cout << (int)getValue(y[i][j][k] + 0.001) << " ";
                }
                cout << endl;
        }
        //*/

        bool *seen = new bool[Stations];
        int *visited = new int[n + 1];

        memset(seen, false, sizeof(bool) * Stations);
        memset(visited, 0, sizeof(int) * n + 1);
        int length = 0;
        k = 0;
        seen[0] = true;

        while (true) {
            for (j = 0; j < Stations; j++) {
                if (k != j) {
                    // cout << k << ", " << j;
                    // cout << " y=" <<(int)getValue(y[i][k][j] + 0.001) << endl;
                    if (((int)getValue(y[i][k][j] + 0.001) == 1)) {
                        // cout << "  ----------------------------- " << k << "-->" << j << endl;
                        visited[++length] = j;
                        if (seen[j]) {
                            // cout << "break \n";
                            break;
                        }
                        seen[j] = true;
                        k = j;
                        j = -1;

                        for (p = 0; p < sta.size(); p++) {

                            if (sta[p] == k) {
                                /*
                                for (int p1 = 0; p1 < sta.size(); p1++) {
                                        cout << sta[p1] << " ";
                                }
                                cout << endl;
                                cout << " erased " << sta[p] << endl;
                                */
                                sta.erase(sta.begin() + p);
                                /*
                                for (int p1 = 0; p1 < sta.size(); p1++) {
                                        cout << sta[p1] << " ";
                                }
                                cout << endl;
                                cout << " erased " << sta[p] << endl;
                                */
                                break;
                            }
                        }
                    }
                }
            }
            // cout << " done ----- \n";
            if (length < n - 1 && /*k == N - 1*/ visited[0] != visited[length]) {
                // cout << " choose 1 ----- \n";
                // cout << "length: " << length << endl;
                if (n == N && length >= N - 1) {
                    break;
                }
                k = sta[0];
                // cout << k << endl;
                length = 0;
                memset(seen, false, sizeof(bool) * Stations);
                memset(visited, 0, sizeof(int) * n + 1);
                visited[length] = k;
                seen[k] = true;
            }
            else if (length < n && length > 0 || sta.size() == 0) {
                break;
            }
            else if (length == 0) {
                // cout << " choose 2 ----- \n";
                if (k < n - 1) {
                    k++;
                }
                else {
                    length = n;
                    break;
                }
                length = 0;
                memset(seen, false, sizeof(bool) * Stations);
                memset(visited, 0, sizeof(int) * n + 1);
                visited[length] = k;
                seen[k] = true;
            }
        }
        // cout << "length: " <<length << endl;
        int l0 = length;
        if (length > 0 && length < n - 1) {
            // for (int ii = 0; ii < nBuses; ii++) {
            length = l0;
            IloExpr clique(getEnv());
            int firstjobinloop = visited[length];
            int newlen = 0;
            while (true) {
                j = visited[length];
                k = visited[--length];
                // cout << "y_{ " << i << ", "<< k << ", " << j<< "}";
                clique += y[i][k][j];
                newlen++;
                if (k == firstjobinloop) {
                    break;
                }
                // cout << " + ";
            }
            // cout << " <=  " << newlen - 1 << endl;
            add(clique <= newlen - 1).end();
            clique.end();
            //}
        }
        delete[] seen;
        delete[] visited;
        // exit(0);
    }

    return;
}

int main() {
    int instance = 0;
    ofstream inst("data/output/Instance_S" + to_string(instance) + ".txt");
    inst << "---------- Weight factors of the objective function -------- \n"<< endl;
    // WEIGHT FACTORS--------------------------------------------------------------------------
    float c1 = 0.25f;
    inst << "c1: " << c1 << " \t (For travel-time of the buses)" << endl;
    float c2 = 0.35f;
    inst << "c2: " << c2 << " \t (For walking time of the passengers)" << endl;
    float c3 = 1 - c1 - c2;
    inst << "c3: " << c3 << " \t (For the absolute difference in desired arrival time and actual arrival time of the passengers)" << endl;
    float lvsea = 1;
    inst << "\n---------------------- Parameters -------------------------- \n" << endl;
    // Define parameters-----------------------------------------------------------------------
    const int nBuses = 2; // amount of buses available
    inst << "Number of buses: " << nBuses << endl;
    const int N = 5; // number of mandatory stations
    inst << "Number of mandatory bus stops: " << N << endl;
    const int M = 3; // number of stations in cluster
    inst << "Number optional bus stops per cluster: " << M << " \n --> One cluster between each mandatory stop: " + to_string((N - 1) * M) + " optional stops in total" << endl;
    const int Stations = (N - 1) * M + N; // amount of Stations
    inst << "Total number of bus stops: " << Stations << endl;
    const int C = 20; // number of clients in opt horizon
    inst << "Number of passenger requests: " << C << endl;

    const int bCapacity = 15; // Bus capcity
    inst << "Bus capacity: " << bCapacity << " passengers" << endl;
    const float pspeed = 1.0f;         // passengers speed in meter per scond
    const float bspeed = 40.0f / 3.6f; // bus speed in m/s
    const int delta = 30;              // acceleration and deceleration time  in seconds
    inst << "Acceleration and deceleration time: " << delta << " seconds" << endl;
    const int tau = 5; // dwell time coeficient in seconds
    inst << "Dwell time per passenger at a bus stop: " << tau << " seconds" << endl;
    const int d = 20 * 60; // threshold of individual walking time in sec
    inst << "Maximum walking time for any passenger: " << d << " seconds" << endl;
    const int d_time1 = 15 * 60;
    const int d_time2 = 5 * 60;
    inst << "Maximum amount of time a passenger can arrive too early: " << d_time1 << " seconds" << endl;
    inst << "Maximum amount of time a passenger can arrive too late: " << d_time2 << " seconds" << endl;
    // const int M0 = 10000; // Big M

    // Gen(N, M, C); // to generate locations

    // Read in locations
    double passengers[C][2];
    ifstream filep("data/input/passengers.txt"); // Passengers
    int i = 0;
    while (i < C) {
        filep >> passengers[i][0] >> passengers[i][1]; // extracts 2 floating point values seperated by whitespace
        i++;
    }

    double mandatory[N][2]; // mandatory Stations
    ifstream filem("data/input/mandatory.txt");
    i = 0;
    while (i < N) {
        filem >> mandatory[i][0] >> mandatory[i][1]; // extracts 2 floating point values seperated by whitespace
        i++;
    }

    double optional[(N - 1) * M][2]; // optinal stations
    ifstream fileo("data/input/optional.txt");
    i = 0;
    while (i < (N - 1) * M) {
        fileo >> optional[i][0] >> optional[i][1]; // extracts 2 floating point values seperated by whitespace
        i++;
    }

    // Arrival times of the passengers
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

    // calculate travel time using manhattan distance
    double traveltimep[C][Stations];        // travel times of people between passangers and stations
    double traveltimes[Stations][Stations]; // travel times of buses between stations
    for (int i = 0; i < C; i++) {
        for (int j = 0; j < N; j++) {
            traveltimep[i][j] = d / 20.0 / 60.0 * (abs(passengers[i][0] - mandatory[j][0]) + abs(passengers[i][1] - mandatory[j][1])) * 1000 / pspeed;
        }

        for (int j = N; j < Stations; j++) {
            traveltimep[i][j] = d / 20.0 / 60.0 * (abs(passengers[i][0] - optional[j - N][0]) + abs(passengers[i][1] - optional[j - N][1])) * 1000 / pspeed;
        }
    }

    inst << "\nWalking time between passengers and bus stops in seconds: " << "\npassengers correspond with the rows, bus stops correspond with the columns \nthe mandatory stops are listed first, then the optional stops are listed" << endl;
    for (int i = 0; i < C; i++) {
        double minp = 10000000000;
        for (int k = 0; k < N; k++) {
            if (minp > traveltimep[i][k]) {
                minp = traveltimep[i][k];
            }
        }
        for (int j = 0; j < Stations; j++) {
            if (traveltimep[i][j] > d || (j > N && traveltimep[i][j] > minp)) {
                inst << "/\t";
            }
            else {
                inst << int(traveltimep[i][j]) << "\t";
            }
        }
        inst << endl;
    }

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
    // exit(0);

    // define sets----------------------------------------------------------------------------------------------
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

    // earliest and latest arrival times
    double EA = 10000000000, LA = -1;
    for (int p = 0; p < C; p++) {
        if (EA > arrivals[p]) {
            EA = arrivals[p];
        }
        if (LA < arrivals[p]) {
            LA = arrivals[p];
        }
    }

    int M_1 = (int)(LA + d_time2 - EA + d_time1) + 1;

    //*************************************************|| MODEL || ***********************************************************************
    IloEnv env;
    IloModel model(env);

    // VARIBLES-----------------------------------------------------------
    IloArray<IloArray<IloBoolVarArray>> x(env); // x variable
    for (int p = 0; p < C; p++) {
        IloArray<IloBoolVarArray> x1(env);
        for (int i = 0; i < nBuses; i++) {
            x1.add(IloBoolVarArray(env, Stations));
        }
        x.add(x1);
    }

    IloArray<IloArray<IloBoolVarArray>> y(env); // y variable
    for (int i = 0; i < nBuses; i++) {
        IloArray<IloBoolVarArray> y1(env);
        for (int j = 0; j < Stations; j++) {
            y1.add(IloBoolVarArray(env, Stations));
        }
        y.add(y1);
    }

    IloNumVarArray D(env);             // Di variable
    for (int i = 0; i < nBuses; i++) { // D
        D.add(IloNumVar(env, 0, INFINITY));
    }

    IloNumVarArray A(env); // Ap variable
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

    // OBJECTIVE -------------------------------------------------------
    IloExpr TravelTime(env);
    for (int i = 0; i < nBuses; i++) {
        for (int j = 0; j < Stations; j++) {
            for (int k = 0; k < Stations; k++) {
                TravelTime += y[i][j][k] * (traveltimes[j][k] + delta);
            }
            for (int p = 0; p < C; p++) {
                TravelTime += x[p][i][j] * tau;
            }
        }
    }

    IloExpr WalkingTime(env);
    for (int p = 0; p < C; p++) {
        for (int i = 0; i < nBuses; i++) {
            for (int j = 0; j < Stations; j++) {
                WalkingTime += x[p][i][j] * (traveltimep[p][j]);
            }
        }
    }

    IloExpr WaitingTime(env);
    for (int p = 0; p < C; p++) {
        WaitingTime += lvsea * q_late[p] + q_early[p];
    }

    model.add(IloMinimize(env, c1 * TravelTime + c2 * WalkingTime + c3 * WaitingTime));
    TravelTime.end();
    WaitingTime.end();
    WalkingTime.end();

    // CONSTRAINTS-------------------------------------------------------

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

    /*
    //Subtour elimination
    int  count = pow(2, Stations); // number of subsets
    for (int i = 0; i < nBuses; i++) {
            for (int s = 0; s < count; s++) {// This loop will generate a subset
                    IloExpr Ysum(env);
                    int S = 0;
                    for (int j = 0; j < Stations; j++) {
                            // This if condition will check if jth bit in binary representation of  s  is set or not
                            // if the value of (i & (1 << j)) is greater than 0 , include y var in the current subset
                            // otherwise exclude it
                            if ((s & (1 << j)) > 0) {
                                    S++;
                                    for (int k = 0; k < Stations; k++) {
                                            if ((s & (1 << k)) > 0)
                                                    Ysum += y[i][j][k];
                                    }
                            }
                    }
                    if (2<= S && S<=Stations-1) {
                            model.add(Ysum <= S - 1);
                    }
                    Ysum.end();
            }
    }
    /*/

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

    /*

    int minbus = nBuses * bCapacity;
    int minpass = C;
    for (const auto& i : I) {
            if (minpass - bCapacity > 0) {
                    minbus -= bCapacity;
                    minpass -= bCapacity;
            }
            else {
                    minbus = i;
                    break;
            }
    }
    for (i = 0; i < minbus; i++) {
            IloExpr sumX(env);
            for (const auto& j : J) {
                    for (const auto& p : P) {
                            sumX += x[p][i][j];
                    }
            }
            model.add(sumX >= minpass);
            sumX.end();
    }
    */

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

    // SOLVE----------------------------------------------------
    IloCplex cplex(model);
    // cplex.setOut(env.getNullStream());
    // cplex.setParam(IloCplex::EpGap, 0.01);
    // cplex.setParam(IloCplex::MIPEmphasis, 3);
    cplex.setParam(cplex.Threads, 12);
    cplex.use(SubtourEliminationCallback(env, y, Stations, nBuses, N, C));
    // cplex.setParam(cplex.TiLim,3610);
    cplex.solve();
    elapsed_time = (double)(clock() - start_time) / CLK_TCK;
    std::cout << "Computational time: " << elapsed_time << " seconds\n" << endl;

	// Get objective function value 
    double objval = cplex.getObjValue();

    int xsol[C][nBuses][Stations];
    ofstream txt_xsol("data/output/xsol.txt");
    int ysol[nBuses][Stations][Stations];
    ofstream txt_ysol("data/output/ysol.txt");
    double Dsol[nBuses];
    ofstream txt_Dsol("data/output/Dsol.txt");
    double Asol[C];
    ofstream txt_Asol("data/output/Asol.txt");

    for (int i = 0; i < nBuses; i++) {
        Dsol[i] = cplex.getValue(D[i]);
        txt_Dsol << Dsol[i] << endl;

        txt_xsol << "Bus " << i << endl;
        txt_ysol << "Bus " << i << endl;
        for (int j = 0; j < Stations; j++) {
            for (int k = 0; k < Stations; k++) {
                ysol[i][j][k] = (int)cplex.getValue(y[i][j][k] + 0.001);
                txt_ysol << ysol[i][j][k] << " ";
            }
            for (int p = 0; p < C; p++) {
                xsol[p][i][j] = (int)cplex.getValue(x[p][i][j] + 0.001);
                txt_xsol << xsol[p][i][j] << " ";
            }
            txt_xsol << endl;
            txt_ysol << endl;
        }
        txt_xsol << "end" << endl;
        txt_ysol << "end" << endl;
    }

    for (int p = 0; p < C; p++) {
        Asol[p] = cplex.getValue(A[p]);
        txt_Asol << Asol[p] << endl;
    }

    txt_xsol.close();
    txt_ysol.close();
    txt_Asol.close();
    txt_Dsol.close();

    env.end();

    cout << endl;

    cout << " Solution is : " << objval << endl;

    cout << " Bus starting times in seconds -------- " << endl;

    for (int i = 0; i < nBuses; i++) {
        cout << "Bus " << i << "--> " << Dsol[i] << endl;
    }

    cout << endl;

    cout << " Early times in seconds  --------- " << endl;
    double AD = 0;
    for (int p = 0; p < C; p++) {
        cout << "Passeger " << p << " --> " << arrivals[p] - Asol[p] << endl;
        AD += abs(arrivals[p] - Asol[p]);
    }

    ofstream txt_visits("data/output/visits.txt");
    int j = 0;
    double travel = 0, travelp = 0;
    for (int i = 0; i < nBuses; i++) {
        j = 0;
        travel = Dsol[i];
        txt_visits << "Bus " << i + 1 << endl;
        while (j != N - 1) {
            for (int k = 0; k < Stations; k++) {
                if (ysol[i][j][k] == 1) {
                    txt_visits << j << " " << travel << endl;
                    travel += (traveltimes[j][k] + delta);
                    for (int p = 0; p < C; p++) {
                        travel += xsol[p][i][j] * tau;
                    }
                    j = k;
                }
            }
        }
        txt_visits << j << " " << travel << endl;
        cout << j << endl;
        txt_visits << "end" << endl;
    }
    txt_visits.close();

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
    cout << " AD: " << AD << endl;
    cout << " T: " << travelp << endl;
    cout << " W: " << W << endl;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started:
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
