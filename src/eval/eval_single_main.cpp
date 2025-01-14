#include <fstream>
#include <iomanip>

#include "other_methods/interface.h"

// Constants
// not needed by most tree construction methods (except Bonn's Algorithm)
// for evaluation by Elmore delay
double dbuPerMicron = 2000;
double unitResistance = 0.0012675;  // Ohm/dbu
double unitCapacitance = 8e-20;     // Farad/dbu
double driverResistance = 25.35;    // Ohm
void SetUnitRCSing();

// From argument parsing
string netFile = "../toys/toy1.net";
int pinNum = 10, seed = 0;
double eps = 0.1;
bool ParseArgs(int argc, char **argv);

int main(int argc, char **argv) {
    printlog("================================================================================");
    printlog("                                   Start ...                                    ");
    printlog("================================================================================");

    // Prepare
    if (!ParseArgs(argc, argv)) return 1;
    SetUnitRCSing();
    printlog("Epsilon is", eps);

    // Read or create a net
    salt::Net net;
    if (!netFile.empty()) {
        net.Read(netFile);
    } else {
        printlog("Generate a random net (seed =", seed, ", #pins =", pinNum, ")");
        printlog();
        srand(seed);
        net.RanInit(0, pinNum);
        net.Write("");
    }

    // vector<Method> methods = {Method::FLUTE,
    //                           Method::KRY,
    //                           Method::PD,
    //                           Method::SALT_R0,
    //                           Method::SALT_R1,
    //                           Method::SALT_R2,
    //                           Method::SALT_R3,
    //                           Method::DelayTree_R1,
    //                           Method::DelayTree_R2,
    //                           Method::Critical};
    vector<Method> methods = {Method::DelayTree_R1, Method::DelayTree_R2, Method::SALT_R3, Method::SALT_R2};
    printlog("method norWL maxStretch avgStretch maxNorDelay avgNorDelay time");
    for (auto method : methods) {
        salt::Tree tree;
        GetATree(net, tree, method, eps);
        salt::CompleteEval eval(driverResistance, tree);
        printGreenLog("\t============ Method = ", method._to_string(), " ============");
        printlog("Wire length is", eval.wireLength);
        printlog("Max path length is", eval.maxPathLength);
        printlog("Avg path length is", eval.avgPathLength);
        printlog("Max stretch (shallowness) is", eval.maxStretch);
        printlog("Avg stretch is", eval.avgStretch);
        printlog("maxNorDelay is", eval.maxNorDelay);
        printlog("avgNorDelay is", eval.avgNorDelay);

        // Output tree to .tree file
        tree.Write(method._to_string());
    }

    printlog("================================================================================");
    printlog("                                      Done ...                                  ");
    printlog("================================================================================");

    return 0;
}

bool ParseArgs(int argc, char **argv) {
    for (int i = 1; i < argc; ++i) {
        if (string(argv[i]) == "-net")
            netFile = string(argv[++i]);
        else if (string(argv[i]) == "-pin")
            pinNum = atoi(argv[++i]);
        else if (string(argv[i]) == "-seed")
            seed = atoi(argv[++i]);
        else if (string(argv[i]) == "-eps")
            eps = atof(argv[++i]);
        else {
            cerr << "Unknown parameter: " << argv[i] << endl;
            cerr << "Usage 1: " << argv[0] << " -net <.net> [-eps <epsilon>]" << endl;
            cerr << "Usage 2: " << argv[0] << " -pin <pin_num> [-seed <rand_seed> -eps <epsilon>]" << endl;
            return false;
        }
    }

    return true;
}

void SetUnitRCSing() {
    printlog("dbu_per_micron :", dbuPerMicron);
    printlog("unit_resistance :", unitResistance, "Ohm/dbu");
    printlog("unit_capacitance :", unitCapacitance, "Farad/dbu");
    printlog("driver_resistance :", driverResistance, "Ohm");
    printlog();

    salt::ElmoreDelayEval::unitRes = unitResistance;
    salt::ElmoreDelayEval::unitCap = unitCapacitance;
    salt::BonnBuilder::unitCap = unitCapacitance;
}
