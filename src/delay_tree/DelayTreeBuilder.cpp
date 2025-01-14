#include "DelayTreeBuilder.h"

#include <cmath>

namespace salt {
bool use_salt_class = 1;

vector<shared_ptr<TreeNode>> get_children(shared_ptr<TreeNode> node) { return node->children; }

void DelayTreeBuilder::Init(Tree& minTree, shared_ptr<Pin> srcP) {
    auto mtNodes = minTree.ObtainNodes();
    _slNodes.resize(mtNodes.size());
    _shortestDists.resize(mtNodes.size());
    _curDists.resize(mtNodes.size());
    for (auto mtN : mtNodes) {
        _slNodes[mtN->id] = mtN;
        _shortestDists[mtN->id] = Dist(mtN->loc, srcP->loc);
        _curDists[mtN->id] = (std::numeric_limits<DTYPE>::max)();
    }
    _curDists[srcP->id] = 0;
    _slSrc = _slNodes[srcP->id];
}

vector<pair<shared_ptr<TreeNode>, int>> Break_points_1;
bool DelayTreeBuilder::Relax(const shared_ptr<TreeNode>& u, const shared_ptr<TreeNode>& v) {
    DTYPE newDist = _curDists[u->id] + Dist(u->loc, v->loc);
    if (_curDists[v->id] > newDist) {
        _curDists[v->id] = newDist;
        return true;
    } else if (_curDists[v->id] == newDist && Dist(u->loc, v->loc) < v->WireToParentChecked()) {
        return true;
    } else
        return false;
}

void DelayTreeBuilder::DFS(const shared_ptr<TreeNode>& smtNode, const shared_ptr<TreeNode>& slNode, double eps) {
    if (smtNode->pin && _curDists[slNode->id] > (1 + eps) * _shortestDists[slNode->id]) {
        _curDists[slNode->id] = _shortestDists[slNode->id];
        Break_points_1.push_back(make_pair(slNode, _length_to_source[slNode->id] - _length_Mah_source[slNode->id]));
    }
    for (auto c : smtNode->children) {
        Relax(slNode, _slNodes[c->id]);
        DFS(c, _slNodes[c->id], eps);
        Relax(_slNodes[c->id], slNode);
    }
}

double flute_critical(const salt::Net& net, salt::Tree& saltTree) {
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
        if (net.pins[i]->slack < 0.0) {
            vio_pin.push_back(net.pins[i]);
            // continue;
        } else {
            unvio_pin.push_back(net.pins[i]);
        }
    }

    int unvio_pin_num = unvio_pin.size();
    for (size_t i = 0; i < unvio_pin_num; ++i) {
        x[i] = unvio_pin[i]->loc.x;
        y[i] = unvio_pin[i]->loc.y;
    }
    if (fluteTree.branch) free(fluteTree.branch);  // is it complete for mem leak?
    fluteTree = flute::flute(unvio_pin_num, x, y, ACCURACY);

    // Build adjacency list
    unordered_map<pair<DTYPE, DTYPE>, shared_ptr<salt::TreeNode>, boost::hash<pair<DTYPE, DTYPE>>> key2node;
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

    for (int i = 0; i < 2 * t.deg - 2; i++) {
        int j = t.branch[i].n;
        if (t.branch[i].x == t.branch[j].x && t.branch[i].y == t.branch[j].y) continue;
        // any more duplicate?
        shared_ptr<salt::TreeNode> n1 = FindOrCreate(t.branch[i].x, t.branch[i].y);
        shared_ptr<salt::TreeNode> n2 = FindOrCreate(t.branch[j].x, t.branch[j].y);
        printlog("s", "%d - %d\n", n1->pin ? n1->pin->id : -1, n2->pin ? n2->pin->id : -1);
        n1->children.push_back(n2);
        n2->children.push_back(n1);
    }

    // Reverse parent-child orders
    saltTree.source = key2node[{net.source()->loc.x, net.source()->loc.y}];

    for (auto sink : vio_pin) {
        auto node = FindOrCreate(sink.get()->loc.x, sink.get()->loc.y);
        node->children.clear();
        node->parent = nullptr;
        node->children.push_back(saltTree.source);
        saltTree.source->children.push_back(node);
    }
    saltTree.SetParentFromUndirectedAdjList();
    saltTree.net = &net;

    auto length = fluteTree.length;

    free(fluteTree.branch);
    return length;
}

double DelayTreeBuilder::run_flute(const salt::Net& net, salt::Tree& saltTree) {
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
    for (size_t i = 0; i < d; ++i) {
        x[i] = net.pins[i]->loc.x;
        y[i] = net.pins[i]->loc.y;
    }
    if (fluteTree.branch) free(fluteTree.branch);
    fluteTree = flute::flute(d, x, y, ACCURACY);

    // Build adjacency list
    unordered_map<pair<DTYPE, DTYPE>, shared_ptr<salt::TreeNode>, boost::hash<pair<DTYPE, DTYPE>>> key2node;
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

    for (int i = 0; i < 2 * t.deg - 2; i++) {
        int j = t.branch[i].n;
        if (t.branch[i].x == t.branch[j].x && t.branch[i].y == t.branch[j].y) continue;
        // any more duplicate?
        shared_ptr<salt::TreeNode> n1 = FindOrCreate(t.branch[i].x, t.branch[i].y);
        shared_ptr<salt::TreeNode> n2 = FindOrCreate(t.branch[j].x, t.branch[j].y);
        n1->children.push_back(n2);
        n2->children.push_back(n1);
    }

    // Reverse parent-child orders
    saltTree.source = key2node[{net.source()->loc.x, net.source()->loc.y}];
    saltTree.SetParentFromUndirectedAdjList();
    saltTree.net = &net;

    auto length = fluteTree.length;

    free(fluteTree.branch);
    return length;
}

void DelayTreeBuilder::Run(const Net& net, Tree& tree, double eps, int refineLevel) {
    _eps = eps;

    auto best_WL = run_flute(net, tree);

    _best_WL = best_WL;

    int node_num = tree.UpdateId();
    _make_steiner_id = node_num;

    _source = tree.source;

    function<void(const shared_ptr<TreeNode>&, DTYPE)> traverse = [&](const shared_ptr<TreeNode>& node, DTYPE curDist) {
        _length_to_source[node->id] = curDist;
        auto node_to_source_manh_dis = utils::Dist(_source->loc, node->loc);
        _length_Mah_source[node->id] = node_to_source_manh_dis;
        for (auto c : node->children) traverse(c, curDist + c->WireToParent());
    };
    traverse(tree.source, 0);

    Break_points_1.clear();
    Init(tree, net.source());
    DFS(tree.source, _slSrc, eps - 1);

    int max_iter = 2 * node_num;
    int iter = 0;
    _current_WL = _best_WL;

#ifdef CRITICAL_SINK
    vector<std::pair<std::shared_ptr<salt::TreeNode>, double>> all_stretch_node;
    all_stretch_node = find_node_be_fix_PL(_source, tree);
    while (!all_stretch_node.empty()) {
        bool improve = PathLengthOpt(tree, all_stretch_node.back().first, eps, _current_WL);
        if (improve) {
            node_num = tree.UpdateId();

            function<void(const shared_ptr<TreeNode>&, DTYPE)> traverse = [&](const shared_ptr<TreeNode>& node,
                                                                              DTYPE curDist) {
                _length_to_source[node->id] = curDist;
                auto node_to_source_manh_dis = utils::Dist(_source->loc, node->loc);
                _length_Mah_source[node->id] = node_to_source_manh_dis;
                for (auto c : node->children) traverse(c, curDist + c->WireToParent());
            };
            traverse(tree.source, 0);
        }
        all_stretch_node.pop_back();
    }

#else
    vector<std::pair<std::shared_ptr<salt::TreeNode>, int>> all_stretch_node;
    if (!use_salt_class) {
        all_stretch_node = find_all_stretch_node(_source, tree);
    } else {
        all_stretch_node = Break_points_1;
        sort(all_stretch_node.begin(),
             all_stretch_node.end(),
             [](pair<shared_ptr<TreeNode>, int> p1, pair<shared_ptr<TreeNode>, int> p2) {
                 return p1.second > p2.second;
                 // return p1.second < p2.second;
             });
    }

    int stretch_node_num = all_stretch_node.size();
    int idx = 0;
    if (!all_stretch_node.empty()) {
        while (all_stretch_node[idx].first && iter < max_iter) {
            bool improve = PathLengthOpt(tree, all_stretch_node[idx].first, eps, _current_WL);
            if (improve) {
                node_num = tree.UpdateId();

                function<void(const shared_ptr<TreeNode>&, DTYPE)> traverse = [&](const shared_ptr<TreeNode>& node,
                                                                                  DTYPE curDist) {
                    _length_to_source[node->id] = curDist;
                    auto node_to_source_manh_dis = utils::Dist(_source->loc, node->loc);
                    _length_Mah_source[node->id] = node_to_source_manh_dis;
                    for (auto c : node->children) traverse(c, curDist + c->WireToParent());
                };
                traverse(tree.source, 0);

                if (!use_salt_class) {
                    all_stretch_node = find_all_stretch_node(_source, tree);
                } else {
                    Break_points_1.clear();
                    Init(tree, net.source());
                    DFS(tree.source, _slSrc, eps - 1);
                    sort(all_stretch_node.begin(),
                         all_stretch_node.end(),
                         [](pair<shared_ptr<TreeNode>, int> p1, pair<shared_ptr<TreeNode>, int> p2) {
                             return p1.second > p2.second;
                         });
                    all_stretch_node = Break_points_1;
                }

                idx = 0;
                stretch_node_num = all_stretch_node.size();
                iter++;

                // tree.Write("SELF");
                salt::WireLengthEval eval(tree);
                _current_WL = eval.wireLength;
                if (all_stretch_node.empty()) {
                    break;
                }

            } else {
                idx++;
                if (idx == stretch_node_num) {
                    break;
                }
            }
        }
    }
#endif

    _make_steiner_id = tree.UpdateId();

#ifdef CRITICAL_SINK
    Refine::DAES_C(tree, eps - 1, 25.35, true);
#else
    if (refineLevel == 1) {
        Refine::DAES_S(tree, eps - 1, 25.35, true);  // DAES-S
        // Refine::DAES_S_C(tree, eps - 1, 25.35, true);  // Used when there is timing information on the pin
    } else if (refineLevel == 2) {
        Refine::DAES(tree, eps - 1, 25.35, true);
        // Refine::DAES_C(tree, eps - 1, 25.35, true);  // Used when there is timing information on the pin
    }
#endif
}

bool DelayTreeBuilder::PathLengthOpt(Tree& tree, shared_ptr<salt::TreeNode> node, double eps, DTYPE current_WL) {
    bool improve = false;

    // 记录没优化前的长度信息
    vector<int> length_to_source;
    vector<int> length_Mah_source;
    length_to_source = _length_to_source;
    length_Mah_source = _length_Mah_source;

    int best_PL = INT_MAX;
    auto parent = node->parent;  // 父节点
    if (parent) {
        std::shared_ptr<salt::TreeNode> next_node = nullptr;

        next_node = find_ancestor_closest_source(node, tree);
        // next_node = find_father_closest_source_across_path(node, tree);
        assert(next_node);
        // 记录优化前的PL
        int pre_PL = std::accumulate(_length_to_source.begin(), _length_to_source.begin() + _pin_num, 0);

        shared_ptr<salt::TreeNode> best_remove_edge_1 = nullptr;
        shared_ptr<salt::TreeNode> best_remove_edge_2 = nullptr;
        DTYPE changed_best_WL = INT_MAX;  // 边删除过程中最好的WL

        set<int> unviolate_node_before;
        for (int i = 0; i < _pin_num; i++) {
            if (_length_to_source[i] <= _eps * _length_Mah_source[i]) {
                unviolate_node_before.insert(i);
            }
        }

        auto path = find_path(node, next_node);
        for (int edge_id = 0; edge_id < static_cast<int>(path.size() - 1); edge_id++) {  // 找到需要移除哪条边
            auto node_1 = path[edge_id];
            auto node_2 = path[edge_id + 1];
            if (!node_2) {
                continue;
            }

            DTYPE add_edge_length = utils::Dist(node->loc, next_node->loc);
            DTYPE remove_edge_length = utils::Dist(node_1->loc, node_2->loc);
            DTYPE changed_WL = current_WL + (add_edge_length - remove_edge_length);

            // 更新每个点到source的长度,只需要更新部分
            _length_to_source[node->id] = _length_to_source[next_node->id] + add_edge_length;

            auto updat_path = find_path(node, node_2);

            for (int i = 1; i < static_cast<int>(updat_path.size() - 1); i++) {
                auto n1 = path[i];
                auto n2 = path[i - 1];
                _length_to_source[n1->id] = _length_to_source[n2->id] + utils::Dist(n1->loc, n2->loc);
                // n1的子节点也要更新
                auto childs = n1->children;
                for (auto child : childs) {
                    if (child->id == n2->id) {
                        continue;
                    }
                    _length_to_source[child->id] = _length_to_source[n1->id] + utils::Dist(child->loc, n1->loc);
                    update_length(child);
                }
            }

            set<int> violate_node_after;
            for (int i = 0; i < _pin_num; i++) {
                if (_length_to_source[i] > _eps * _length_Mah_source[i]) {
                    violate_node_after.insert(i);
                }
            }
            // violate_node_after中不能包含unviolate_node_before中的点
            set<int> violate_node;
            set_intersection(violate_node_after.begin(),
                             violate_node_after.end(),
                             unviolate_node_before.begin(),
                             unviolate_node_before.end(),
                             inserter(violate_node, violate_node.begin()));
            if (!violate_node.empty()) {
                _length_to_source = length_to_source;
                continue;
            }

            auto post_PL = std::accumulate(_length_to_source.begin(), _length_to_source.begin() + _pin_num, 0);
            if ((post_PL < pre_PL && post_PL < best_PL) || (post_PL == best_PL && changed_WL < changed_best_WL) ||
                (post_PL < pre_PL && changed_WL < changed_best_WL) ||
                (post_PL - best_PL) < (changed_best_WL - changed_WL)) {
                best_PL = post_PL;
                best_remove_edge_1 = node_1;
                best_remove_edge_2 = node_2;
                changed_best_WL = changed_WL;
            }
            _length_to_source = length_to_source;
        }  // 找到需要移除哪条边

        if (best_remove_edge_1 && best_remove_edge_2) {
            if (best_remove_edge_1->parent == best_remove_edge_2) {
                flip_order(best_remove_edge_1, node, next_node);
            } else {
                flip_order(best_remove_edge_2, next_node, node);
            }
            improve = true;
            _current_WL = changed_best_WL;
        }
    }
    return improve;
}

vector<shared_ptr<TreeNode>> DelayTreeBuilder::find_path(shared_ptr<TreeNode> s_node, shared_ptr<TreeNode> e_node) {
    // 找到路径上的所有点
    std::vector<shared_ptr<TreeNode>> path;
    auto current_node = s_node;

    vector<bool> visited(2 * _pin_num, false);

    do {
        path.push_back(current_node);
        visited[current_node->id] = true;
        if (current_node->id == e_node->id) {
            break;
        }
        current_node = current_node->parent;
    } while (current_node);

    if (!current_node) {
        std::vector<shared_ptr<TreeNode>> path1;
        shared_ptr<TreeNode> common_ancestor = e_node;
        path1.push_back(e_node);
        while (!visited[common_ancestor->id]) {
            common_ancestor = common_ancestor->parent;
            path1.push_back(common_ancestor);
        }

        while (1) {
            auto n = path.back();
            path.pop_back();
            if (n == common_ancestor) {
                break;
            }
        }
        std::reverse(path1.begin(), path1.end());
        path.insert(path.end(), path1.begin(), path1.end());
    }
    return path;
}

void DelayTreeBuilder::flip_order(shared_ptr<TreeNode> p_node,
                                  shared_ptr<TreeNode> c_node,
                                  shared_ptr<TreeNode> c_node_father) {
    auto path = find_path(c_node, p_node);
    std::reverse(path.begin(), path.end());

    auto path_length = path.size();

    for (size_t i = 0; i < path_length - 1; ++i) {
        TreeNode::ResetParent(path[i]);
        TreeNode::SetParent(path[i], path[i + 1]);
    }
    TreeNode::ResetParent(path[path_length - 1]);
    TreeNode::SetParent(path[path_length - 1], c_node_father);
}

void DelayTreeBuilder::update_length(shared_ptr<TreeNode> node) {
    auto father_id = node->id;
    auto childs = node->children;
    for (auto child : childs) {
        auto child_id = child->id;
        _length_to_source[child_id] = _length_to_source[father_id] + child->WireToParent();
        update_length(child);
    }
}

vector<pair<shared_ptr<TreeNode>, int>> DelayTreeBuilder::find_all_stretch_node(shared_ptr<salt::TreeNode> source,
                                                                                Tree& tree) {
    vector<pair<shared_ptr<TreeNode>, int>> stretch_node;

    function<void(const shared_ptr<TreeNode>&, DTYPE)> traverse = [&](const shared_ptr<TreeNode>& node, DTYPE curDist) {
        _length_to_source[node->id] = curDist;
        auto node_to_source_manh_dis = _length_Mah_source[node->id];
        for (auto c : node->children) traverse(c, curDist + c->WireToParent());
        if (node->id < _pin_num) {
            if (curDist > _eps * node_to_source_manh_dis) {
                stretch_node.push_back(make_pair(node, curDist - node_to_source_manh_dis));
            }
        }
    };
    traverse(source, 0);
    sort(stretch_node.begin(),
         stretch_node.end(),
         [](pair<shared_ptr<TreeNode>, int> p1, pair<shared_ptr<TreeNode>, int> p2) { return p1.second > p2.second; });
    return stretch_node;
}

vector<pair<shared_ptr<TreeNode>, double>> DelayTreeBuilder::find_node_be_fix_PL(shared_ptr<salt::TreeNode> source,
                                                                                 Tree& tree) {
    vector<pair<shared_ptr<TreeNode>, double>> stretch_node;

    tree.PreOrder([&](const shared_ptr<TreeNode>& node) {
        if (node->pin) {
            cout << "pin id " << node->id << " slack " << node->pin->slack << endl;
            if (node->pin->slack < 0.0) {
                stretch_node.push_back(make_pair(node, node->pin->slack));
            }
        }
    });

    sort(stretch_node.begin(),
         stretch_node.end(),
         [](pair<shared_ptr<TreeNode>, double> p1, pair<shared_ptr<TreeNode>, double> p2) {
             return p1.second > p2.second;
             // return p1.second < p2.second;
         });

    return stretch_node;
}

shared_ptr<salt::TreeNode> DelayTreeBuilder::find_ancestor_closest_source(shared_ptr<salt::TreeNode> node, Tree& tree) {
    shared_ptr<salt::TreeNode> next_node = nullptr;
    shared_ptr<salt::TreeNode> next_node5 = nullptr;
    int extra_wl = INT_MAX;
    tree.PostOrderCopy([&](shared_ptr<TreeNode> nd) {
        // 两个点不能有连接
        if (nd->id != node->id && node->parent != nd && nd->parent != node) {
            auto dist = utils::Dist(nd->loc, node->loc);
            auto pl = dist + _length_to_source[nd->id];
            if (pl <= _eps * _length_Mah_source[node->id] && dist < extra_wl) {
                if (node->parent->id != nd->id) {
                    next_node = nd;
                    extra_wl = dist;
                }
            }
        }
    });

    return next_node;
}

shared_ptr<salt::TreeNode> DelayTreeBuilder::find_father_closest_source_across_path(shared_ptr<salt::TreeNode> node,
                                                                                    Tree& tree) {
    shared_ptr<salt::TreeNode> next_node = nullptr;
    shared_ptr<salt::TreeNode> next_node5 = nullptr;
    int extra_wl = INT_MAX;
    auto parent = node->parent;
    while (parent != nullptr) {
        auto nd = parent;
        // 两个点不能有连接
        if (nd->id != node->id) {
            auto dist = utils::Dist(nd->loc, node->loc);
            auto pl = dist + _length_to_source[nd->id];
            if (pl <= _eps * _length_Mah_source[node->id] && dist < extra_wl) {
                next_node = nd;
                extra_wl = dist;
            }
        }
        parent = nd->parent;
    }
    return next_node;
}

}  // namespace salt
