#pragma once

#include <algorithm>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>
#include <vector>

#include "../salt/base/flute.h"
#include "../salt/base/rsa.h"
#include "../salt/refine/refine.h"
#include "salt/base/tree.h"
#include "salt/utils/utils.h"

#include <boost/functional/hash.hpp>
#include <limits>
#include <unordered_map>

#include "../salt/base/eval.h"
#include "../salt/base/flute/flute.h"
#include "util.h"

// #define CRITICAL_SINK

namespace salt {

class DelayTreeBuilder {
public:
    DelayTreeBuilder(int pin_num)
        : _pin_num(pin_num), _length_to_source(2 * pin_num), _length_Mah_source(2 * pin_num) {}
    double run_flute(const salt::Net& net, salt::Tree& saltTree);
    void Run(const Net& net, Tree& tree, double eps, int refineLevel = 2);

private:
    void Init(Tree& minTree, shared_ptr<Pin> srcP);

    bool Relax(const shared_ptr<TreeNode>& u, const shared_ptr<TreeNode>& v);  // from u to v
    void DFS(const shared_ptr<TreeNode>& mstNode, const shared_ptr<TreeNode>& slNode, double eps);

    void update_length(shared_ptr<TreeNode> node);

    vector<pair<shared_ptr<TreeNode>, int>> find_all_stretch_node(shared_ptr<salt::TreeNode> source, Tree& tree);

    // Determine which points need to be fixed based on slack
    vector<pair<shared_ptr<TreeNode>, double>> find_node_be_fix_PL(shared_ptr<salt::TreeNode> source, Tree& tree);

    bool PathLengthOpt(Tree& tree, shared_ptr<salt::TreeNode> node, double eps, DTYPE best_WL);

    vector<shared_ptr<TreeNode>> find_path(shared_ptr<TreeNode> s_node, shared_ptr<TreeNode> e_node);

    void flip_order(shared_ptr<TreeNode> p_node, shared_ptr<TreeNode> c_node, shared_ptr<TreeNode> c_node_father);

    shared_ptr<salt::TreeNode> find_ancestor_closest_source(shared_ptr<salt::TreeNode> node, Tree& tree);

    shared_ptr<salt::TreeNode> find_father_closest_source_across_path(shared_ptr<salt::TreeNode> node, Tree& tree);

    vector<DTYPE> _shortestDists;
    vector<DTYPE> _curDists;
    vector<shared_ptr<TreeNode>> _slNodes;  // nodes of the shallow-light tree
    shared_ptr<TreeNode> _slSrc;            // source node of the shallow-light tree

    int _pin_num;
    DTYPE _best_WL;
    DTYPE _current_WL;
    shared_ptr<TreeNode> _source;
    vector<int> _length_to_source;
    vector<int> _length_Mah_source;

    double _eps;

    int _make_steiner_id;
};
}  // namespace salt