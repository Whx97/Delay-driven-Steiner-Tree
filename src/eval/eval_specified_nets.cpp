#include "other_methods/interface.h"
#include "salt/base/eval.h"  // for salt::WireLengthEval
#include "salt/salt.h"

#include <iostream>

double dbuPerMicron = 2000;
double unitResistance = 0.0012675;  // Ohm/dbu
double unitCapacitance = 8e-20;     // Farad/dbu
double driverResistance = 25.35;    // Ohm
void SetUnitRCSing();

int pinNum = 10, seed = 0;
double eps = 1;

pair<salt::CompleteEval, salt::CompleteEval> run_net(string netFile, double eps) {
    printlog("================================================================================");
    printlog("                                   Start ...                                    ");
    printlog("================================================================================");

    SetUnitRCSing();

    salt::Net net;
    net.Read(netFile);
    printlog("Run SALT algorithm on net", net.name, "with", net.pins.size(), "pins using epsilon =", eps);

    salt::Tree tree;
    utils::timer time;

    salt::DelayTreeBuilder DTB(net.pins.size());
    DTB.Run(net, tree, eps);

    // run SALT
    salt::Net net_salt;
    net_salt.Read(netFile);
    salt::Tree tree_salt;
    salt::SaltBuilder saltB;
    saltB.Run(net_salt, tree_salt, eps - 1);

    // Report
    printlog("Tree topology is as follows:");
    // tree.Write("SELF");
    // tree_salt.Write("SALT");

    salt::CompleteEval eval(driverResistance, tree);
    salt::CompleteEval eval_salt(driverResistance, tree_salt);

    printlog("================================================================================");
    printlog("                                      Done ...                                  ");
    printlog("================================================================================");

    return make_pair(eval, eval_salt);
}

void SetUnitRCSing() {
    printlog("dbu_per_micron :", dbuPerMicron);
    printlog("unit_resistance :", unitResistance, "Ohm/dbu");
    printlog("unit_capacitance :", unitCapacitance, "Farad/dbu");
    printlog("driver_resistance :", driverResistance, "Ohm");
    printlog();

    salt::ElmoreDelayEval::unitRes = unitResistance;
    salt::ElmoreDelayEval::unitCap = unitCapacitance;
}

int string_get_int(string s) {
    std::string tempNumberString = "";

    for (char c : s) {
        if (isdigit(c)) {
            tempNumberString += c;
        }
    }

    if (!tempNumberString.empty()) {
        int extractedNumber = std::stoi(tempNumberString);
        return extractedNumber;
    }

    return 0;
}

int main(int argc, char **argv) {
    vector<string> targetNames = {
        "../toys/toy1.net",
        "../toys/test/test_17.net",
        "../toys/test/test_19.net",
        "../toys/test/test_26.net",
        "../toys/test/test_28.net",
        "../toys/test/test_30.net",
        "../toys/test/test_31.net",
        "../toys/test/test_33.net",
    };

    vector<int> net_num{};
    vector<pair<salt::CompleteEval, salt::CompleteEval>> trees;
    for (auto net_path : targetNames) {
        auto tree = run_net(net_path, 1);
        trees.push_back(tree);
        net_num.push_back(string_get_int(net_path));
    }
    for (int i = 0; i < static_cast<int>(trees.size()); i++) {
        auto tree = trees[i];
        auto eval = tree.first;
        auto eval_salt = tree.second;

        COUT_BLUE_START;
        cout << "\t\t" << targetNames[i] << endl;
        COUT_RED_START;
        std::cout << std::left << std::setw(25) << " " << std::setw(15) << "SELF" << std::setw(15) << "SALT"
                  << std::setw(15) << "delta" << std::endl;

        std::cout << std::left << std::setw(25) << "Wire length is " << std::setw(15) << eval.wireLength
                  << std::setw(15) << eval_salt.wireLength << std::setw(15) << showpos
                  << eval.wireLength - eval_salt.wireLength << endl;
        std::cout << std::left << std::setw(25) << "Max path length is" << std::setw(15) << noshowpos
                  << eval.maxPathLength << std::setw(15) << eval_salt.maxPathLength << std::setw(15) << showpos
                  << eval.maxPathLength - eval_salt.maxPathLength << endl;
        std::cout << std::left << std::setw(25) << "Avg path length is" << std::setw(15) << noshowpos
                  << eval.avgPathLength << std::setw(15) << eval_salt.avgPathLength << std::setw(15) << showpos
                  << eval.avgPathLength - eval_salt.avgPathLength << endl;
        std::cout << std::left << std::setw(25) << "shallowness is" << std::setw(15) << noshowpos << eval.maxStretch
                  << std::setw(15) << eval_salt.maxStretch << std::setw(15) << showpos
                  << eval.maxStretch - eval_salt.maxStretch << endl;
        std::cout << std::left << std::setw(25) << "Avg stretch is" << std::setw(15) << noshowpos << eval.avgStretch
                  << std::setw(15) << eval_salt.avgStretch << std::setw(15) << showpos
                  << eval.avgStretch - eval_salt.avgStretch << endl;
        std::cout << std::left << std::setw(25) << "maxNorDelay is" << std::setw(15) << noshowpos << eval.maxNorDelay
                  << std::setw(15) << eval_salt.maxNorDelay << endl;
        std::cout << std::left << std::setw(25) << "avgNorDelay is" << std::setw(15) << noshowpos << eval.avgNorDelay
                  << std::setw(15) << eval_salt.avgNorDelay << endl;
        std::cout << "\033[0m";
    }

    COUT_GREEN_START;
    bool show_salt = 0;
    if (show_salt) {
        cout << std::left << "SALT\n";
    }
    cout << std::left << std::setw(15) << " " << std::setw(15) << " " << std::setw(15) << "   SELF" << std::setw(15)
         << " " << std::setw(15) << " " << std::right << "|   " << endl;
    cout << std::left << std::setw(15) << "pin number" << std::setw(15) << "Wire length" << std::setw(15)
         << "shallowness" << std::setw(15) << "maxNorDelay" << std::setw(15) << "avgNorDelay" << std::right << "|   "
         << endl;
    for (int i = 0; i < static_cast<int>(trees.size()); i++) {
        auto tree = trees[i];
        auto eval = tree.first;
        auto eval_salt = tree.second;
        if (show_salt) {
            eval = eval_salt;
        }

        auto show = [](int net_num, salt::CompleteEval eval) {
            cout << std::left << std::setw(15) << net_num << std::setw(15) << eval.wireLength << std::setw(15)
                 << eval.maxStretch << std::setw(15) << eval.maxNorDelay << std::setw(15) << eval.avgNorDelay
                 << std::right << "|   ";
        };

        show(net_num[i], eval);
        show(net_num[i], eval_salt);
        cout << endl;
    }
    COUT_COLOR_END;
}