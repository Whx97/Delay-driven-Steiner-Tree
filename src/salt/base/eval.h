#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>

#include "tree.h"
namespace salt {

// Used to record whether WL changes after optimizing PL or not
using IsConstWL = bool;
// Change in PL divided by change in WL
using DeltaPL_WL = double;
// The change in max delay divided by the change in WL
using DeltaMaxDelay_WL = double;
// avg change in delay divided by change in WL
using DeltaAvgDelay_WL = double;

class WireLengthEvalBase {
public:
    DTYPE wireLength;

    WireLengthEvalBase() = default;
    void Update(const Tree& tree);
    WireLengthEvalBase(const Tree& tree) { Update(tree); }
};

class WireLengthEval : public WireLengthEvalBase {
public:
    DTYPE maxPathLength;
    double avgPathLength;
    double norPathLength;  // avgPathLength / avgShortestPathLength
    double maxStretch;     // max{pathLength / shortestPathLength}
    double avgStretch;     // avg{pathLength / shortestPathLength}
    int max_path_length_id = -1;
    vector<int> nodeLevel;

    WireLengthEval() = default;
    void Update(const Tree& tree);
    WireLengthEval(const Tree& tree) { Update(tree); }
};

inline ostream& operator<<(ostream& os, const WireLengthEval& eval) {
    os << " wl=" << eval.wireLength << " mp=" << eval.maxPathLength << " ap=" << eval.avgPathLength
       << " ms=" << eval.maxStretch << " as=" << eval.avgStretch;
    return os;
}

//********************************************************************************

class ElmoreDelayEval {
public:
    static double unitRes;
    static double unitCap;

    DTYPE fluteWL = -1;
    double maxDelay;
    double avgDelay;
    double sumDelay;
    double sumNorDelay;
    double maxNorDelay;
    double avgNorDelay;
    int max_delay_node_id = -1;
    int max_delay_node_level = -1;
    vector<double> nodeDelay;
    vector<double> nodeNorDelay;
    vector<double> low_bound_delay;

    double score = 0;

    ElmoreDelayEval() {}
    void Calc(double rd, Tree& tree);                           // Not all pins
    void Update(double rd, Tree& tree, bool normalize = true);  // tree node id will be updated
    ElmoreDelayEval(double rd, Tree& tree, bool normalize = true) { Update(rd, tree, normalize); }

    double _maxLb = 0.0;  // For each net, it is recalculated each time, so cache it, and the value is only constant for one net
    ElmoreDelayEval(double rd, Tree& tree, double maxLb) {
        assert(rd > 0);
        assert(unitRes > 0 && unitCap > 0);
        int numNodes = tree.ObtainNodes().size();
        maxDelay = avgDelay = maxNorDelay = avgNorDelay = sumDelay = sumNorDelay = 0;
        auto delay = GetDelay(rd, tree, numNodes);  // delay for all tree nodes
        nodeDelay = delay;
        tree.PreOrder([&](const shared_ptr<TreeNode>& node) {
            if (!node->pin || node == tree.source) return;
            sumDelay += delay[node->id];
        });
        sumNorDelay = sumDelay / maxLb;

        nodeNorDelay.resize(numNodes, 0);
        for (int i = 0; i < numNodes; i++) {
            nodeNorDelay[i] = nodeDelay[i] / maxLb;
        }
    }

private:
    vector<double> GetDelay(double rd, const Tree& tree, int numNode);
    vector<double> GetDelayLB(double rd, const Tree& tree);
};

inline ostream& operator<<(ostream& os, const ElmoreDelayEval& eval) {
    os << " md=" << eval.maxDelay << " ad=" << eval.avgDelay << " mnd=" << eval.maxNorDelay
       << " and=" << eval.avgNorDelay;
    return os;
}

//********************************************************************************

class CompleteEval : public WireLengthEval, public ElmoreDelayEval {
public:
    double norWL;

    CompleteEval() = default;
    void Update(double rd, Tree& tree) {
        WireLengthEval::Update(tree);
        ElmoreDelayEval::Update(rd, tree);
        norWL = double(wireLength) / fluteWL;
    }

    CompleteEval(double rd, Tree& tree) { Update(rd, tree); }
};

class CompleteStat {
public:
    double norWL = 0, maxStretch = 0, avgStretch = 0, norPathLength = 0, maxNorDelay = 0, avgNorDelay = 0;

    int unchanged_topo_num = 0;
    int num_const_wl_net = 0;
    double NorDeltaPL_WL = 0.0;
    double NorDeltaMaxDelay_WL = 0.0;
    double NorDeltaAvgDelay_WL = 0.0;
    long long int flute_WL = 0;
    long long int optPL_WL = 0;
    double bound_WL = 0.0;  // bound_WL = optPL_WL / flute_WL
    double delta_WL = 0.0;  // delta_WL = optPL_WL - flute_WL
    double delta_NorWL = 0.0;
    double delta_NorWL_reduce = 0.0;  // WL shortened not added to the calculation

    long long int deltaPL = 0.0;  // delta max PL
    double deltaNorPL = 0.0;      // Norm. delta max PL
    double deltaMaxDelay = 0.0;   // Norm. delta Max Delay
    double deltaAvgDelay = 0.0;   // Norm. delta Avg Delay

    vector<double> deltaMaxDelay_vec;
    vector<double> deltaAvgDelay_vec;

    long long int WL_reduce_num = 0;
    double bound_WL_reduce = 0.0;  // WL shortened not added to the calculation

    double delta_max_delay_same_wl = 0.0;
    double delta_avg_delay_same_wl = 0.0;
    double delta_pl_same_wl = 0.0;

    vector<double> WL;
    int cnt = 0;
    double eps;
    double time = 0;
    double score = 0;  // delay score
    void Inc(const CompleteEval& eval, double runtime = 0.0) {
        WL.push_back(eval.norWL);
        ++cnt;
        norWL += eval.norWL;
        maxStretch += eval.maxStretch;
        avgStretch += eval.avgStretch;
        norPathLength += eval.norPathLength;
        maxNorDelay += eval.maxNorDelay;
        avgNorDelay += eval.avgNorDelay;
        time += runtime;
        score += eval.score;
    }
    void Avg() {
        norWL /= cnt;
        maxStretch /= cnt;
        avgStretch /= cnt;
        norPathLength /= cnt;
        maxNorDelay /= cnt;
        avgNorDelay /= cnt;
        time /= cnt;
        score /= cnt;

        if (cnt > 0) {
            flute_WL /= cnt;
            optPL_WL /= cnt;
            bound_WL /= cnt;
            delta_WL /= cnt;
            deltaPL /= cnt;

            delta_NorWL /= cnt;
            deltaNorPL /= cnt;
            deltaMaxDelay /= cnt;
            deltaAvgDelay /= cnt;
        }

        if ((cnt - num_const_wl_net) > 0) {
            NorDeltaPL_WL /= (cnt - num_const_wl_net);
            NorDeltaMaxDelay_WL /= (cnt - num_const_wl_net);
            NorDeltaAvgDelay_WL /= (cnt - num_const_wl_net);
        }

        if ((cnt - WL_reduce_num) > 0) {
            bound_WL_reduce /= (cnt - WL_reduce_num);
            delta_NorWL_reduce /= (cnt - WL_reduce_num);
        }

        if (num_const_wl_net > 0) {
            delta_max_delay_same_wl /= num_const_wl_net;
            delta_avg_delay_same_wl /= num_const_wl_net;
            delta_pl_same_wl /= num_const_wl_net;
        }
    }
};

class TimeStat {
public:
    double time = 0;
    int cnt = 0;
    double change_topo_time = 0;
    int change_topo_cnt = 0;
    double recover_topo_time = 0;
    int recover_topo_cnt = 0;
    double calc_delay_time = 0;
    int calc_delay_cnt = 0;

    void Inc_change_topo(double runtime) {
        change_topo_time += runtime;
        ++change_topo_cnt;
    }
    void Inc_recover_topo(double runtime) {
        recover_topo_time += runtime;
        ++recover_topo_cnt;
    }
    void Inc_calc_delay(double runtime) {
        calc_delay_time += runtime;
        ++calc_delay_cnt;
    }
    void Avg() { time /= cnt; }
};

}  // namespace salt