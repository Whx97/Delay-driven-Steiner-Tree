#include <fstream>
#include <sstream>

#include "other_methods/interface.h"

using namespace std;

vector<salt::Net> nets;
double driverResistance;  // TODO: net-specific
bool ReadNets(const string& netFileName);
bool WriteNets(const string& netFileName);

// 产生不同slack的nets
void benchmark() {
    // vector<string> cases = {"superblue1",
    //                         "superblue3",
    //                         "superblue4",
    //                         "superblue5",
    //                         "superblue7",
    //                         "superblue10",
    //                         "superblue16",
    //                         "superblue18"};
    vector<string> cases = {"superblue18"};
    string out_file_path = "../iccad15_nets/bench/";
    for (auto cas : cases) {
        string in_file_path = "../iccad15_nets/iccad15_nets/" + cas + ".nets";
        string out_file_name = out_file_path + cas + "_slack.nets";
        if (!ReadNets(in_file_path)) return;
        WriteNets(out_file_name);
        nets.clear();
    }
}

bool ReadNets(const string& netFileName) {
    ifstream netFile(netFileName);

    // 0. Check the file
    if (!netFile.is_open()) {
        cerr << "ERROR: Cannot open file \"" << netFileName << "\"" << endl;
        return false;
    }

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

    // 3. Read nets
    do {
        nets.emplace_back();
    } while (nets.back().Read(netFile));
    nets.pop_back();

    printlog("# nets is", nets.size());
    return true;
}

void writeNet(ostream& os, const salt::Net& net) {
    // header
    string header = to_string(net.id) + " " + net.name + " " + to_string(net.pins.size());
    // string header = name + " " + to_string(pins.size());
    header += " -cap -slack";
    os << "Net " << header << endl;

    // pins
    auto pins = net.pins;
    for (const auto& pin : pins) {
        os << pin->id << " " << pin->loc.x << " " << pin->loc.y;
        os << " " << pin->cap;
        // 随机分配一个slack，值域为[-4,4]
        int slack = rand() % 9 - 4;
        if (slack >= 0) {
            slack = 0;
        }
        os << " " << slack;

        os << endl;
    }
    os << endl;
}

bool WriteNets(const string& netFileName) {
    ofstream netFile(netFileName);

#include <chrono>
#include <iomanip>
#include <sstream>

    // 1. Write the header
    netFile << "# Network Configuration File" << endl;
    netFile << "# Routing benchmark generated from ICCAD 2015 benchmark superblue4" << endl;
    // 获取当前时间点
    auto now = std::chrono::system_clock::now();
    // 转换为时间结构
    std::time_t now_time_t = std::chrono::system_clock::to_time_t(now);
    std::tm now_tm = *std::localtime(&now_time_t);

    std::stringstream ss;
    ss << std::put_time(&now_tm, "%Y-%m-%d %H:%M:%S");

    netFile << "# Date    : " << ss.str() << endl;
    netFile << "# Note    : Length unit is dbu\n" << endl;

    // 2. Write unit RC
    // 2.1 Write the keyword PARAMETERS
    netFile << "PARAMETERS" << endl;
    netFile << endl;  // Write an empty line

    netFile << "dbu_per_micron : 2000" << endl;
    netFile << "unit_resistance : 0.0012675 Ohm/dbu" << endl;
    netFile << "unit_capacitance : 8e-20 Farad/dbu" << endl;
    netFile << "driver_resistance : 25.35 Ohm" << endl;

    netFile << endl;
    netFile << "NETS\n" << endl;

    // 3. Write nets
    for (auto& net : nets) {
        writeNet(netFile, net);
    }

    printlog("# nets written:", nets.size());
    return true;
}

int main() {
    printlog("================================================================================");
    printlog("                               Benchmark Start ...                              ");
    printlog("================================================================================");
    // 产生具有不同slack的nets
    benchmark();

    printlog("================================================================================");
    printlog("                                      Done ...                                  ");
    printlog("================================================================================");
    return 0;
}