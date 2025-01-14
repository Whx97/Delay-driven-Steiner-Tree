#include "critical.h"

#include <boost/functional/hash.hpp>
#include <unordered_map>

#include "../salt/base/flute/flute.h"

namespace salt {

void CriticalBuilder::Run(const Net& net, Tree& tree) {
    // load LUT
    static bool once = false;
    if (!once) {
        flute::readLUT();
        once = true;
    }

    // Obtain flute tree
    flute::Tree fluteTree;
    fluteTree.branch = nullptr;
    int d = net.pins.size();
    assert(d <= MAXD);
    int x[MAXD], y[MAXD];

    vector<shared_ptr<Pin>> vio_pin;
    vector<shared_ptr<Pin>> unvio_pin;

    for (size_t i = 0; i < d; ++i) {
        if (net.pins[i]->id == net.source()->id) {
            continue;
        }
        if (net.pins[i]->slack < 0.0) {
            vio_pin.push_back(net.pins[i]);
            // continue;
        } else {
            unvio_pin.push_back(net.pins[i]);
        }
    }

    // unvio_pin need add source node
    unvio_pin.push_back(net.source());
    int unvio_pin_num = unvio_pin.size();
    if (unvio_pin_num > 1) {
        for (size_t i = 0; i < unvio_pin_num; ++i) {
            x[i] = unvio_pin[i]->loc.x;
            y[i] = unvio_pin[i]->loc.y;
        }
        if (fluteTree.branch) free(fluteTree.branch);  // is it complete for mem leak?
        fluteTree = flute::flute(unvio_pin_num, x, y, ACCURACY);
    }

    // Build adjacency list
    unordered_map<pair<DTYPE, DTYPE>, shared_ptr<salt::TreeNode>, boost::hash<pair<DTYPE, DTYPE>>> key2node;
    unordered_map<pair<DTYPE, DTYPE>, shared_ptr<salt::TreeNode>, boost::hash<pair<DTYPE, DTYPE>>> key2steinernode;
    for (auto p : net.pins) {
        key2node[{p->loc.x, p->loc.y}] = make_shared<salt::TreeNode>(p);
    }
    auto& t = fluteTree;

    auto FindOrCreate = [&](DTYPE x, DTYPE y) {
        auto it = key2node.find({x, y});
        if (it == key2node.end()) {
            shared_ptr<salt::TreeNode> node = make_shared<salt::TreeNode>(x, y);
            key2node[{x, y}] = node;
            return node;
        } else
            return it->second;
    };
    auto SteinerFineOrCreate = [&](DTYPE x, DTYPE y) {
        auto it = key2steinernode.find({x, y});
        if (it == key2node.end()) {
            shared_ptr<salt::TreeNode> node = make_shared<salt::TreeNode>(x, y);
            key2steinernode[{x, y}] = node;
            return node;
        } else
            return it->second;
    };

    auto Find = [&](DTYPE x, DTYPE y) {
        auto it = key2node.find({x, y});
        assert(it != key2node.end());
        return it->second;
    };

    // printtree(t);

    for (int i = 0; i < 2 * t.deg - 2; i++) {
        if (unvio_pin_num < 2) break;
        int j = t.branch[i].n;
        if (i >= t.deg) {
            if (t.branch[i].x == t.branch[j].x && t.branch[i].y == t.branch[j].y) continue;
        }

        shared_ptr<salt::TreeNode> n1;
        shared_ptr<salt::TreeNode> n2;

        if (i < t.deg) {
            n1 = Find(t.branch[i].x, t.branch[i].y);
        } else {
            n1 = SteinerFineOrCreate(t.branch[i].x, t.branch[i].y);
        }

        if (j < t.deg) {
            n2 = Find(t.branch[j].x, t.branch[j].y);
        } else {
            n2 = SteinerFineOrCreate(t.branch[j].x, t.branch[j].y);
        }

        n1->children.push_back(n2);
        n2->children.push_back(n1);
    }

    // Reverse parent-child orders
    tree.source = key2node[{net.source()->loc.x, net.source()->loc.y}];

    for (auto sink : vio_pin) {
        auto node = Find(sink.get()->loc.x, sink.get()->loc.y);
        node->children.clear();
        node->parent = nullptr;
        node->children.push_back(tree.source);
        tree.source->children.push_back(node);
    }
    tree.SetParentFromUndirectedAdjList();
    tree.net = &net;

    free(fluteTree.branch);
}

}  // namespace salt