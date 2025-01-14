#include <fstream>
#include <sstream>

#include "other_methods/interface.h"

// From argument parsing
vector<salt::Net> nets;
double driverResistance;  // TODO: net-specific

// For classified stat
// the kth bucket is [pinClasses[k], pinClasses[k+1])
// 0 is for safety, in case of nets with # pins < 4
vector<int> pinClasses = {0, 4, 8, 16, 32};

bool ReadNets(const string& netFileName);
void EvalAllMethods(const string& evalFileNameSuffix);
vector<salt::CompleteStat> EvalAllNets(Method method, double eps);

int main(int argc, char** argv) {
    printlog("================================================================================");
    printlog("                                   Start ...                                    ");
    printlog("================================================================================");

    // if (!ReadNets("../iccad15_nets/iccad15_nets/superblue1.nets")) return 1;
    // if (!ReadNets("../iccad15_nets/iccad15_nets/superblue3.nets")) return 1;
    // if (!ReadNets("../iccad15_nets/iccad15_nets/superblue4.nets")) return 1;
    // if (!ReadNets("../iccad15_nets/iccad15_nets/superblue5.nets")) return 1;
    // if (!ReadNets("../iccad15_nets/iccad15_nets/superblue7.nets")) return 1;
    // if (!ReadNets("../iccad15_nets/iccad15_nets/superblue10.nets")) return 1;
    // if (!ReadNets("../iccad15_nets/iccad15_nets/superblue16.nets")) return 1;
    if (!ReadNets("../iccad15_nets/iccad15_nets/superblue18.nets")) return 1;

    // Evaluate all methods
    EvalAllMethods("v2");

    printlog("================================================================================");
    printlog("                                      Done ...                                  ");
    printlog("================================================================================");

    return 0;
}

bool ReadNets(const string& netFileName) {
    cout << netFileName << endl;
    ifstream netFile(netFileName);

    // 1. Skip the header
    string line;
    do {
        getline(netFile, line);
    } while (!line.empty() && line[0] == '#');

    // 2. Read unit RC
    // 2.1 check keyword PARAMETERS
    string buf;
    netFile >> buf;
    if (buf != "PARAMETERS") {
        cerr << "Cannot find keyword PARAMETERS" << endl;
        return false;
    }
    getline(netFile, buf);
    getline(netFile, buf);  // skip an empty line
    // 2.2 read
    vector<string> targetNames = {"dbu_per_micron", "unit_resistance", "unit_capacitance", "driver_resistance"};
    vector<string> values(targetNames.size());
    string name, colon;
    for (unsigned i = 0; i < targetNames.size(); ++i) {
        getline(netFile, buf);
        printlog(buf);
        istringstream iss(buf);
        iss >> name >> colon >> values[i];
        if (name != targetNames[i]) {
            cerr << "Parameter " << i + 1 << " should be " << targetNames[i] << endl;
            return false;
        }
    }
    double unitResistance = stod(values[1]);
    double unitCapacitance = stod(values[2]);
    driverResistance = stod(values[3]);
    // 2.3 set
    salt::ElmoreDelayEval::unitRes = unitResistance;
    salt::ElmoreDelayEval::unitCap = unitCapacitance;
    salt::BonnBuilder::unitCap = unitCapacitance;

    // 3. Read nets
    do {
        nets.emplace_back();
    } while (nets.back().Read(netFile));
    nets.pop_back();

    printlog("# nets is", nets.size());
    return true;
}

void EvalAllMethods(const string& evalFileNameSuffix) {
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
    vector<Method> methods = {Method::DelayTree_R1, Method::DelayTree_R2};

    for (auto method : methods) {
        vector<vector<salt::CompleteStat>> allStats;  // <eps, pin_class> -> stat
        allStats.push_back(EvalAllNets(method, 0));
        if (method != +Method::FLUTE && method != +Method::RSA && method != +Method::MST) {
            // double epsMax = 5.4;  // 13个点
            double epsMax = 0.08;  // 3个点
            for (double eps = 0.05; eps <= epsMax; eps *= 1.5) {
                allStats.push_back(EvalAllNets(method, eps));
            }
        }

        // write eval stat
        for (size_t j = 0; j < (pinClasses.size() + 1); ++j) {
            string pinClass = (j == pinClasses.size()) ? "_all" : "_pin" + to_string(pinClasses[j]);
            ofstream ofs("../result/" + string(method._to_string()) + "_" + evalFileNameSuffix + pinClass);
            ofs << "# norWL maxStretch avgStretch maxNorDelay avgNorDelay cnt eps avgTime" << endl;
            for (const auto& stat : allStats) {
                const auto& s = stat[j];
                ofs << s.norWL << " " << s.maxStretch << " " << s.avgStretch << " " << s.norPathLength << " "
                    << s.maxNorDelay << " " << s.avgNorDelay << " " << s.cnt << " " << s.eps << " " << s.time << endl;
            }
            ofs << endl;
            ofs << "# norWL maxStretch avgStretch maxNorDelay avgNorDelay cnt eps avgTime Score" << endl;
            for (const auto& stat : allStats) {
                const auto& s = stat[j];
                ofs << s.norWL << " " << s.maxStretch << " " << s.avgStretch << " " << s.norPathLength << " "
                    << s.maxNorDelay << " " << s.avgNorDelay << " " << s.cnt << " " << s.eps << " " << s.time << " "
                    << s.score << endl;
            }

            ofs.close();
        }
    }
}

vector<salt::CompleteStat> EvalAllNets(Method method, double eps) {
    vector<salt::CompleteStat> stats(pinClasses.size() + 1);  // the last is the total summary
    // printlog("Method = ", method, ", eps = ", eps);
    printGreenLog("Method = ", method, ", eps = ", eps);

    int id = 0;
    int total_net_num = nets.size();
    for (const auto& net : nets) {
        id++;

        if (id % (total_net_num / 10) == 0) {
            int percent = static_cast<int>((id + 1.0) / total_net_num * 10000);
            log() << " " << percent / 100 << "." << std::setw(2) << std::setfill('0') << percent % 100 << "% "
                  << " id = " << id << ", total_net_num = " << total_net_num << std::endl;
        }

        salt::Tree tree;
        // route
        utils::timer routeTimer;
        if (method == +Method::DelayTree_R1 || method == +Method::DelayTree_R2 || method == +Method::SALT_R3 ||
            method == +Method::SALT_R2) {
            GetATree(net, tree, method, eps, false);
            double routeTime = routeTimer.elapsed();
            // check
            tree.QuickCheck();
            // eval and stat
            salt::CompleteEval eval(driverResistance, tree);
            int i = pinClasses.size() - 1;
            while (pinClasses[i] > net.pins.size()) {
                --i;
            }
            stats[i].Inc(eval, routeTime);
            stats.back().Inc(eval, routeTime);  // total
            continue;
        }
        GetATree(net, tree, method, eps, false);
        double routeTime = routeTimer.elapsed();
        // check
        tree.QuickCheck();
        // eval and stat
        salt::CompleteEval eval(driverResistance, tree);
        int i = pinClasses.size() - 1;
        while (pinClasses[i] > net.pins.size()) {
            --i;
        }
        stats[i].Inc(eval, routeTime);
        stats.back().Inc(eval, routeTime);  // total
    }

    for (auto& stat : stats) {
        stat.eps = eps;
        stat.Avg();
    }
    return stats;
}