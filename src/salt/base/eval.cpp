#include "eval.h"

#include <algorithm>

#include "flute.h"

namespace salt {

void WireLengthEvalBase::Update(const Tree& tree) {
    wireLength = 0;
    tree.PostOrder([&](const shared_ptr<TreeNode>& node) {
        if (node->parent) {
            wireLength += node->WireToParent();
        }
    });
}

void WireLengthEval::Update(const Tree& tree) {
    // wirelength
    WireLengthEvalBase::Update(tree);
    // path
    vector<DTYPE> pathLength(tree.net->pins.size());
    nodeLevel.resize(tree.net->pins.size());
    function<void(const shared_ptr<TreeNode>&, DTYPE)> traverse = [&](const shared_ptr<TreeNode>& node, DTYPE curDist) {
        if (node->pin) {
            pathLength[node->pin->id] = curDist;
            nodeLevel[node->pin->id] = node->level;
        }
        for (auto c : node->children) traverse(c, curDist + c->WireToParent());
    };
    traverse(tree.source, 0);
    maxPathLength = 0;
    double totalPathLength = 0;
    double totalShortestPathLength = 0;
    maxStretch = 0;
    avgStretch = 0;
    for (auto p : tree.net->pins) {
        if (p->IsSink()) {
            DTYPE pl = pathLength[p->id], sp = Dist(tree.source->loc, p->loc);
            double stretch = double(pl) / sp;
            if (pl > maxPathLength) {
                maxPathLength = pl;
                max_path_length_id = p->id;
            }
            totalPathLength += pl;
            totalShortestPathLength += sp;
            if (stretch > maxStretch) {
                maxStretch = stretch;
            }
            avgStretch += stretch;
        }
    }
    auto numSink = tree.net->pins.size() - 1;
    avgPathLength = totalPathLength / numSink;
    norPathLength = totalPathLength / totalShortestPathLength;
    avgStretch /= numSink;
}

//********************************************************************************

double ElmoreDelayEval::unitRes = -1;
double ElmoreDelayEval::unitCap = -1;

void ElmoreDelayEval::Calc(double rd, Tree& tree) {
    assert(rd > 0);
    assert(unitRes > 0 && unitCap > 0);

    int numNodes = 5 * tree.net->pins.size();
    maxDelay = avgDelay = maxNorDelay = avgNorDelay = sumDelay = sumNorDelay = 0;

    auto delay = GetDelay(rd, tree, numNodes);  // delay for all tree nodes
    tree.PreOrder([&](const shared_ptr<TreeNode>& node) {
        if (!node->pin || node == tree.source) return;
        maxDelay = max(maxDelay, delay[node->id]);
        avgDelay += delay[node->id];
        sumDelay += delay[node->id];
    });
}

void ElmoreDelayEval::Update(double rd, Tree& tree, bool normalize) {
    assert(rd > 0);
    assert(unitRes > 0 && unitCap > 0);
    int numPins = tree.net->pins.size();
    int numNodes = tree.UpdateId();
    maxDelay = avgDelay = maxNorDelay = avgNorDelay = sumDelay = sumNorDelay = 0;
    nodeDelay.resize(numNodes, 0);
    nodeNorDelay.resize(numNodes, 0);

    auto delay = GetDelay(rd, tree, numNodes);  // delay for all tree nodes
    nodeDelay = delay;                          // nodeDelay[node->id] = delay[node->id];

    std::vector<double> node_slack;
    node_slack.resize(numNodes, 0);
    tree.PreOrder([&](const shared_ptr<TreeNode>& node) {
        if (!node->pin || node == tree.source) return;
        maxDelay = max(maxDelay, delay[node->id]);
        avgDelay += delay[node->id];
        sumDelay += delay[node->id];
        max_delay_node_id = maxDelay == delay[node->id] ? node->id : max_delay_node_id;
        max_delay_node_level = maxDelay == delay[node->id] ? node->level : max_delay_node_level;
        node_slack[node->id] = node->pin->slack;
    });
    avgDelay /= (numPins - 1);

    if (!normalize) return;

    auto lb = GetDelayLB(rd, tree);  // delay lb for all pins, 0 is source
    low_bound_delay = lb;
    auto maxLb = *max_element(lb.begin(), lb.end());
    _maxLb = maxLb;

    for (int i = 0; i < numNodes; i++) {
        if (i >= numPins || i == tree.source.get()->id) continue;
        // nodeNorDelay[i] = nodeDelay[i] / maxLb;
        nodeNorDelay[i] = nodeDelay[i] / lb[i];

        double delay_score = (nodeNorDelay[i] - node_slack[i]);
        delay_score = max(0.0, delay_score);
        score += delay_score;
    }
    maxNorDelay = maxDelay / maxLb;
    avgNorDelay = avgDelay / maxLb;
    sumNorDelay = sumDelay / maxLb;
}

vector<double> ElmoreDelayEval::GetDelay(double rd, const Tree& tree, int numNode) {
    // get node cap by post-order traversal
    vector<double> cap(numNode, 0);
    tree.PostOrder([&](const shared_ptr<TreeNode>& node) {
        if (node->pin && node != tree.source) cap[node->id] = node->pin->cap;
        for (auto c : node->children) {
            cap[node->id] += cap[c->id];
            cap[node->id] += c->WireToParent() * unitCap;
        }
    });

    // get delay by post-order traversal
    vector<double> delay(numNode, 0);
    tree.PreOrder([&](const shared_ptr<TreeNode>& node) {
        if (node == tree.source)
            delay[node->id] = rd * cap[node->id];
        else {
            double dist = node->WireToParent();
            delay[node->id] = dist * unitRes * (0.5 * dist * unitCap + cap[node->id]) + delay[node->parent->id];
        }
    });
    return delay;
}

vector<double> ElmoreDelayEval::GetDelayLB(double rd, const Tree& tree) {
    vector<double> lb(tree.net->pins.size(), 0);

    // call flute and get smt
    Tree flute;
    FluteBuilder fluteB;
    fluteB.Run(*tree.net, flute);
    WireLengthEvalBase wl(flute);
    fluteWL = wl.wireLength;

    double totalcap = 0;
    for (auto p : tree.net->pins) totalcap += p->cap;

    double lb_sd = rd * (fluteWL * unitCap + totalcap);
    for (auto pin : tree.net->pins) {
        if (pin->IsSource()) continue;
        double dist = Dist(tree.source->loc, pin->loc);
        lb[pin->id] = dist * unitRes * (0.5 * dist * unitCap + pin->cap) + lb_sd;
    }

    return lb;
}

}  // namespace salt