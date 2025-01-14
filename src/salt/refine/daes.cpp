#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "refine.h"
#include "salt/base/eval.h"
#include "salt/base/mst.h"

namespace salt {

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using BPoint = bg::model::point<DTYPE, 2, bg::cs::cartesian>;
using BSegment = bg::model::segment<BPoint>;
using BBox = bg::model::box<BPoint>;
using BPolygon = bg::model::polygon<BPoint>;
using RNode = pair<BBox, shared_ptr<TreeNode>>;  // R-Tree node
struct RNodeComp {
    bool operator()(const RNode& l, const RNode& r) const {
        return bg::equals(l.first, r.first) && l.second == r.second;
    }
};

void Refine::DAES_S(Tree& tree, double eps, double rd, bool useRTree) {
    bgi::rtree<RNode, bgi::rstar<8>, bgi::indexable<RNode>, RNodeComp> rtree;
    if (useRTree) {
        tree.PostOrder([&](const shared_ptr<TreeNode>& n) {
            if (n->parent) {
                BBox s;
                bg::envelope(BSegment(BPoint(n->loc.x, n->loc.y), BPoint(n->parent->loc.x, n->parent->loc.y)), s);
                rtree.insert({s, n});
            }
        });
    }
    auto Disconnect = [&](const shared_ptr<TreeNode>& n) {
        if (useRTree) {
            BBox s;
            bg::envelope(BSegment(BPoint(n->loc.x, n->loc.y), BPoint(n->parent->loc.x, n->parent->loc.y)), s);
            rtree.remove({s, n});
        }
        TreeNode::ResetParent(n);
    };
    auto Connect = [&](const shared_ptr<TreeNode>& n, const shared_ptr<TreeNode>& parent) {
        TreeNode::SetParent(n, parent);
        if (useRTree) {
            BBox s;
            bg::envelope(BSegment(BPoint(n->loc.x, n->loc.y), BPoint(parent->loc.x, parent->loc.y)), s);
            rtree.insert({s, n});
        }
    };

    // 预备
    salt::ElmoreDelayEval before_delayEval(rd, tree);
    auto before_sumDelay = before_delayEval.sumNorDelay;

    auto for_one_net_maxLB = before_delayEval._maxLb;

    int iter_num = 0;
    int max_iter = 5;

    while (true) {
        if (iter_num == max_iter) break;
        iter_num++;
        // Get nearest neighbors
        int num = tree.UpdateId();
        vector<shared_ptr<TreeNode>> nodes = tree.ObtainNodes(),
                                     orderedNodes(nodes.size());  // note: all pins should be covered

        vector<Point> points(nodes.size());
        for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
            orderedNodes[nodes[i]->id] = nodes[i];
            points[nodes[i]->id] = nodes[i]->loc;
        }
        nodes = orderedNodes;
        vector<vector<int>> nearestNeighbors;
        if (!useRTree) {
            MstBuilder mstB;
            mstB.GetAllNearestNeighbors(points, nearestNeighbors);
        } else {
            nearestNeighbors.resize(nodes.size());
            for (auto n : nodes) {
                if (n->parent) {
                    Point c = n->loc;  // center
                    DTYPE radius = n->WireToParent();
                    int tora = 1;
                    BBox queryBox{{c.x - tora * radius, c.y - tora * radius},
                                  {c.x + tora * radius, c.y + tora * radius}};
                    vector<RNode> cands;
                    rtree.query(bgi::intersects(queryBox), back_inserter(cands));
                    for (const auto& cand : cands) {
                        if (cand.second->id >= num) {
                            continue;
                        }
                        if (n->parent->id == cand.second->id) {
                            continue;
                        }
                        nearestNeighbors[n->id].push_back(cand.second->id);
                    }
                }
            }
        }

        // Prune descendants in nearest neighbors
        vector<int> preOrderIdxes(nodes.size(), -1);
        int globalPreOrderIdx = 0;
        function<void(const shared_ptr<TreeNode>&)> removeDescendants = [&](const shared_ptr<TreeNode>& node) {
            preOrderIdxes[node->id] = globalPreOrderIdx++;
            for (auto child : node->children) {
                removeDescendants(child);
            }
            for (auto& neighIdx : nearestNeighbors[node->id]) {
                int neighPreOrderIdx = preOrderIdxes[neighIdx];
                if (neighPreOrderIdx != -1 && neighPreOrderIdx >= preOrderIdxes[node->id]) {
                    neighIdx = -1;  // -1 stands for "descendant"
                }
            }
        };
        removeDescendants(tree.source);

        // Init path lengths and subtree slacks
        vector<DTYPE> pathLengths(nodes.size());
        vector<DTYPE> slacks(nodes.size());
        auto UpdatePathLengths = [&](const shared_ptr<TreeNode>& node) {
            if (node->parent) {
                pathLengths[node->id] = pathLengths[node->parent->id] + node->WireToParent();
            } else {
                pathLengths[node->id] = 0;
            }
        };
        auto UpdateSlacks = [&](const shared_ptr<TreeNode>& node) {
            if (node->children.empty()) {
                slacks[node->id] =
                    Dist(node->loc, tree.source->loc) * (1 + eps) - pathLengths[node->id];  // floor here...
            } else {
                DTYPE minSlack = Dist(node->loc, tree.source->loc) * (1 + eps) - pathLengths[node->id];
                for (auto child : node->children) {
                    minSlack = min(minSlack, slacks[child->id]);
                }
                slacks[node->id] = minSlack;
            }
        };
        tree.PreOrder(UpdatePathLengths);
        tree.PostOrder(UpdateSlacks);

        // Find legal candidate moves
        using MoveT = tuple<double, shared_ptr<TreeNode>, shared_ptr<TreeNode>>;
        vector<MoveT> candidateMoves;  // <wireLengthDelta, node, newParent>
        auto GetNearestPoint = [](const shared_ptr<TreeNode>& target, const shared_ptr<TreeNode>& neigh) {
            Box box(neigh->loc, neigh->parent->loc);
            box.Legalize();
            return box.GetNearestPointTo(target->loc);
        };

        vector<double> cap(num, 0);
        tree.PostOrder([&](const shared_ptr<TreeNode>& node) {
            if (node->pin && node != tree.source) cap[node->id] = node->pin->cap;
            for (auto c : node->children) {
                cap[node->id] += cap[c->id];
                cap[node->id] += c->WireToParent() * 8e-20;
            }
        });

        for (auto node : nodes) {
            std::shared_ptr<salt::TreeNode> tmp_steinerNode = nullptr;
            if (!(node->parent)) {
                continue;
            }

            auto originalParent = node->parent;

            double best_delay_delta = 0;

            shared_ptr<TreeNode> bestNewParent;
            for (int neighIdx : nearestNeighbors[node->id]) {
                if (neighIdx >= static_cast<int>(nodes.size())) continue;
                if (neighIdx == -1 || !nodes[neighIdx]->parent) continue;
                auto neigh = nodes[neighIdx];
                auto neighParent = neigh->parent;

                if (node == neighParent) continue;
                if (node == neigh) continue;
                auto steinerPt = GetNearestPoint(node, neigh);

                // pl不能违规
                DTYPE pathLengthDelta =
                    pathLengths[neighParent->id] + Dist(node->loc, neighParent->loc) - pathLengths[node->id];
                if (pathLengthDelta > slacks[node->id]) continue;

                tmp_steinerNode = nullptr;
                if (steinerPt == neigh->loc) {
                    TreeNode::ResetParent(node);
                    TreeNode::SetParent(node, neigh);
                } else if (steinerPt == neighParent->loc) {
                    TreeNode::ResetParent(node);
                    TreeNode::SetParent(node, neighParent);
                } else {
                    tmp_steinerNode = make_shared<TreeNode>(steinerPt);
                    TreeNode::SetParent(tmp_steinerNode, neighParent);
                    TreeNode::ResetParent(neigh);
                    TreeNode::SetParent(neigh, tmp_steinerNode);
                    TreeNode::ResetParent(node);
                    TreeNode::SetParent(node, tmp_steinerNode);
                    // for later moves
                    tmp_steinerNode->id = nodes.size();
                }

                salt::ElmoreDelayEval delayEval(rd, tree, for_one_net_maxLB);
                auto sumDelay = delayEval.sumNorDelay;

                if (sumDelay + 1e-12 < before_sumDelay) {
                    auto delay_delta = sumDelay - before_sumDelay;
                    if (delay_delta < best_delay_delta) {
                        best_delay_delta = delay_delta;
                        bestNewParent = neigh;
                    }
                }

                if (steinerPt == neigh->loc) {
                    TreeNode::ResetParent(node);
                    TreeNode::SetParent(node, originalParent);
                } else if (steinerPt == neighParent->loc) {
                    TreeNode::ResetParent(node);
                    TreeNode::SetParent(node, originalParent);
                } else {
                    TreeNode::ResetParent(tmp_steinerNode);
                    TreeNode::ResetParent(node);
                    TreeNode::ResetParent(neigh);
                    TreeNode::SetParent(node, originalParent);
                    TreeNode::SetParent(neigh, neighParent);
                }
            }
            if (bestNewParent) {
                candidateMoves.emplace_back(best_delay_delta, node, bestNewParent);
            }
            tmp_steinerNode = nullptr;
            tmp_steinerNode.reset();
        }
        if (candidateMoves.empty()) {
            break;
        }

        // Try candidate moves in the order of descending wire length savings
        // Note that earlier moves may influence the legality of later one
        sort(candidateMoves.begin(), candidateMoves.end(), [](const MoveT& lhs, const MoveT& rhs) {
            return get<0>(lhs) < get<0>(rhs);
        });

        for (const auto& move : candidateMoves) {
            auto node = get<1>(move), neigh = get<2>(move);
            auto neighParent = neigh->parent;
            // check due to earlier moves
            if (TreeNode::IsAncestor(node, neighParent)) continue;
            DTYPE pathLengthDelta =
                pathLengths[neighParent->id] + Dist(node->loc, neighParent->loc) - pathLengths[node->id];
            if (pathLengthDelta > slacks[node->id]) continue;
            auto steinerPt = GetNearestPoint(node, neigh);

            auto originalParent = node->parent;
            std::shared_ptr<salt::TreeNode> tmp_steinerNode = nullptr;
            if (steinerPt == neigh->loc) {
                TreeNode::ResetParent(node);
                TreeNode::SetParent(node, neigh);
            } else if (steinerPt == neighParent->loc) {
                TreeNode::ResetParent(node);
                TreeNode::SetParent(node, neighParent);
            } else {
                tmp_steinerNode = make_shared<TreeNode>(steinerPt);
                TreeNode::SetParent(tmp_steinerNode, neighParent);
                TreeNode::ResetParent(neigh);
                TreeNode::SetParent(neigh, tmp_steinerNode);
                TreeNode::ResetParent(node);
                TreeNode::SetParent(node, tmp_steinerNode);
                // for later moves
                tmp_steinerNode->id = nodes.size();
            }

            salt::ElmoreDelayEval delayEval(rd, tree, for_one_net_maxLB);
            auto sumDelay = delayEval.sumNorDelay;

            if (steinerPt == neigh->loc) {
                TreeNode::ResetParent(node);
                TreeNode::SetParent(node, originalParent);
            } else if (steinerPt == neighParent->loc) {
                TreeNode::ResetParent(node);
                TreeNode::SetParent(node, originalParent);
            } else {
                TreeNode::ResetParent(tmp_steinerNode);
                TreeNode::ResetParent(node);
                TreeNode::ResetParent(neigh);
                TreeNode::SetParent(node, originalParent);
                TreeNode::SetParent(neigh, neighParent);
            }

            tmp_steinerNode.reset();

            if (sumDelay + 1e-12 < before_sumDelay) {
                before_sumDelay = sumDelay;
            } else {
                continue;
            }

            // break
            Disconnect(node);
            // reroot
            if (steinerPt == neigh->loc) {
                Connect(node, neigh);
            } else if (steinerPt == neighParent->loc) {
                Connect(node, neighParent);
            } else {
                auto steinerNode = make_shared<TreeNode>(steinerPt);
                Connect(steinerNode, neighParent);
                Disconnect(neigh);
                Connect(neigh, steinerNode);
                Connect(node, steinerNode);
                // for later moves
                steinerNode->id = nodes.size();
                nodes.push_back(steinerNode);
                pathLengths.push_back(pathLengths[neighParent->id] + steinerNode->WireToParent());
                slacks.push_back(Dist(steinerNode->loc, tree.source->loc) * (1 + eps) - pathLengths.back());
            }
            // update slack for later moves: first subtree, then path to source
            TreeNode::PreOrder(neighParent, UpdatePathLengths);
            TreeNode::PostOrder(neighParent, UpdateSlacks);

            tree.UpdateId();

            auto tmp = neighParent;
            while (tmp->parent) {
                slacks[tmp->parent->id] = min(slacks[tmp->parent->id], slacks[tmp->id]);
                tmp = tmp->parent;
            }
        }

        // Finalize
        tree.RemoveTopoRedundantSteiner();
        tree.PostOrderCopy([&](const shared_ptr<TreeNode>& node) {
            // degree may change after post-order traversal of its children
            if (node->pin) return;
            if (node->children.empty()) {
                Disconnect(node);
            } else if (node->children.size() == 1) {
                auto oldParent = node->parent, oldChild = node->children[0];
                Disconnect(node);
                Disconnect(oldChild);
                Connect(oldChild, oldParent);
            }
        });
    }
}

// DAES
void Refine::DAES(Tree& tree, double eps, double rd, bool useRTree) {
    bgi::rtree<RNode, bgi::rstar<8>, bgi::indexable<RNode>, RNodeComp> rtree;
    if (useRTree) {
        tree.PostOrder([&](const shared_ptr<TreeNode>& n) {
            if (n->parent) {
                BBox s;
                bg::envelope(BSegment(BPoint(n->loc.x, n->loc.y), BPoint(n->parent->loc.x, n->parent->loc.y)), s);
                rtree.insert({s, n});
            }
        });
    }
    auto Disconnect = [&](const shared_ptr<TreeNode>& n) {
        if (useRTree) {
            BBox s;
            bg::envelope(BSegment(BPoint(n->loc.x, n->loc.y), BPoint(n->parent->loc.x, n->parent->loc.y)), s);
            rtree.remove({s, n});
        }
        TreeNode::ResetParent(n);
    };
    auto Connect = [&](const shared_ptr<TreeNode>& n, const shared_ptr<TreeNode>& parent) {
        TreeNode::SetParent(n, parent);
        if (useRTree) {
            BBox s;
            bg::envelope(BSegment(BPoint(n->loc.x, n->loc.y), BPoint(parent->loc.x, parent->loc.y)), s);
            rtree.insert({s, n});
        }
    };

    // 预备
    salt::ElmoreDelayEval before_delayEval(rd, tree);
    auto before_sumDelay = before_delayEval.sumNorDelay;

    auto for_one_net_maxLB = before_delayEval._maxLb;

    int iter_num = 0;
    int max_iter = 2;

    while (true) {
        if (iter_num == max_iter) break;
        iter_num++;
        // Get nearest neighbors
        int num = tree.UpdateId();
        vector<shared_ptr<TreeNode>> nodes = tree.ObtainNodes(),
                                     orderedNodes(nodes.size());  // note: all pins should be covered

        vector<Point> points(nodes.size());
        for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
            orderedNodes[nodes[i]->id] = nodes[i];
            points[nodes[i]->id] = nodes[i]->loc;
        }
        nodes = orderedNodes;
        vector<vector<int>> nearestNeighbors;
        if (!useRTree) {
            MstBuilder mstB;
            mstB.GetAllNearestNeighbors(points, nearestNeighbors);
        } else {
            nearestNeighbors.resize(nodes.size());
            for (auto n : nodes) {
                if (n->parent) {
                    Point c = n->loc;  // center
                    DTYPE radius = n->WireToParent();
                    int tora = 1;
                    BBox queryBox{{c.x - tora * radius, c.y - tora * radius},
                                  {c.x + tora * radius, c.y + tora * radius}};
                    vector<RNode> cands;
                    rtree.query(bgi::intersects(queryBox), back_inserter(cands));
                    for (const auto& cand : cands) {
                        if (cand.second->id >= num) {
                            continue;
                        }
                        if (n->parent->id == cand.second->id) {
                            continue;
                        }
                        nearestNeighbors[n->id].push_back(cand.second->id);
                    }
                }
            }
        }

        // Prune descendants in nearest neighbors
        vector<int> preOrderIdxes(nodes.size(), -1);
        int globalPreOrderIdx = 0;
        function<void(const shared_ptr<TreeNode>&)> removeDescendants = [&](const shared_ptr<TreeNode>& node) {
            preOrderIdxes[node->id] = globalPreOrderIdx++;
            for (auto child : node->children) {
                removeDescendants(child);
            }
            for (auto& neighIdx : nearestNeighbors[node->id]) {
                int neighPreOrderIdx = preOrderIdxes[neighIdx];
                if (neighPreOrderIdx != -1 && neighPreOrderIdx >= preOrderIdxes[node->id]) {
                    neighIdx = -1;  // -1 stands for "descendant"
                }
            }
        };
        removeDescendants(tree.source);

        // Init path lengths and subtree slacks
        vector<DTYPE> pathLengths(nodes.size());
        vector<DTYPE> slacks(nodes.size());
        auto UpdatePathLengths = [&](const shared_ptr<TreeNode>& node) {
            if (node->parent) {
                pathLengths[node->id] = pathLengths[node->parent->id] + node->WireToParent();
            } else {
                pathLengths[node->id] = 0;
            }
        };
        auto UpdateSlacks = [&](const shared_ptr<TreeNode>& node) {
            if (node->children.empty()) {
                slacks[node->id] =
                    Dist(node->loc, tree.source->loc) * (1 + eps) - pathLengths[node->id];  // floor here...
            } else {
                DTYPE minSlack = Dist(node->loc, tree.source->loc) * (1 + eps) - pathLengths[node->id];
                for (auto child : node->children) {
                    minSlack = min(minSlack, slacks[child->id]);
                }
                slacks[node->id] = minSlack;
            }
        };
        tree.PreOrder(UpdatePathLengths);
        tree.PostOrder(UpdateSlacks);

        // Find legal candidate moves
        using MoveT = tuple<double,
                            shared_ptr<TreeNode>,
                            shared_ptr<TreeNode>,
                            int>;  // 最后一个的int，1表示选neigh，2表示选neighParent，3表示选steinerPt
        int which_node = 0;
        vector<MoveT> candidateMoves;  // <wireLengthDelta, node, newParent>
        auto GetNearestPoint = [](const shared_ptr<TreeNode>& target, const shared_ptr<TreeNode>& neigh) {
            Box box(neigh->loc, neigh->parent->loc);
            box.Legalize();
            return box.GetNearestPointTo(target->loc);
        };

        vector<double> cap(num, 0);
        tree.PostOrder([&](const shared_ptr<TreeNode>& node) {
            if (node->pin && node != tree.source) cap[node->id] = node->pin->cap;
            for (auto c : node->children) {
                cap[node->id] += cap[c->id];
                cap[node->id] += c->WireToParent() * 8e-20;
            }
        });

        for (auto node : nodes) {
            std::shared_ptr<salt::TreeNode> tmp_steinerNode = nullptr;
            if (!(node->parent)) {
                continue;
            }

            auto originalParent = node->parent;

            double best_delay_delta = 0;

            shared_ptr<TreeNode> bestNewParent;
            for (int neighIdx : nearestNeighbors[node->id]) {
                if (neighIdx >= static_cast<int>(nodes.size())) continue;
                if (neighIdx == -1 || !nodes[neighIdx]->parent) continue;
                auto neigh = nodes[neighIdx];
                auto neighParent = neigh->parent;

                if (node == neighParent) continue;
                if (node == neigh) continue;

                // pl不能违规
                DTYPE pathLengthDelta =
                    pathLengths[neighParent->id] + Dist(node->loc, neighParent->loc) - pathLengths[node->id];
                if (pathLengthDelta > slacks[node->id]) continue;

                auto steinerPt = GetNearestPoint(node, neigh);

                bool exct_twice = (steinerPt == neigh->loc) || (steinerPt == neighParent->loc) ? 1 : 0;
                vector<shared_ptr<TreeNode>> tmp_nodes;
                tmp_nodes.emplace_back(neigh);
                tmp_nodes.emplace_back(neighParent);
                for (int i = 0; i < 2; i++) {
                    TreeNode::ResetParent(node);
                    TreeNode::SetParent(node, tmp_nodes[i]);

                    if (i == 0) {
                        // pl不能违规
                        DTYPE pathLengthDelta =
                            pathLengths[tmp_nodes[i]->id] + Dist(node->loc, tmp_nodes[i]->loc) - pathLengths[node->id];
                        if (pathLengthDelta > slacks[node->id]) continue;
                    }

                    salt::ElmoreDelayEval delayEval(rd, tree, for_one_net_maxLB);
                    auto sumDelay = delayEval.sumNorDelay;

                    if (sumDelay + 1e-12 < before_sumDelay) {
                        auto delay_delta = sumDelay - before_sumDelay;
                        if (delay_delta < best_delay_delta) {
                            best_delay_delta = delay_delta;
                            bestNewParent = neigh;
                            which_node = i + 1;
                        }
                    }

                    TreeNode::ResetParent(node);
                    TreeNode::SetParent(node, originalParent);
                }

                tmp_steinerNode = nullptr;
                if (!exct_twice) {
                    tmp_steinerNode = make_shared<TreeNode>(steinerPt);
                    TreeNode::SetParent(tmp_steinerNode, neighParent);
                    TreeNode::ResetParent(neigh);
                    TreeNode::SetParent(neigh, tmp_steinerNode);
                    TreeNode::ResetParent(node);
                    TreeNode::SetParent(node, tmp_steinerNode);
                    tmp_steinerNode->id = nodes.size();

                    salt::ElmoreDelayEval delayEval(rd, tree, for_one_net_maxLB);
                    auto sumDelay = delayEval.sumNorDelay;

                    if (sumDelay + 1e-12 < before_sumDelay) {
                        auto delay_delta = sumDelay - before_sumDelay;
                        if (delay_delta < best_delay_delta) {
                            best_delay_delta = delay_delta;
                            bestNewParent = neigh;
                            which_node = 3;
                        }
                    }

                    TreeNode::ResetParent(tmp_steinerNode);
                    TreeNode::ResetParent(node);
                    TreeNode::ResetParent(neigh);
                    TreeNode::SetParent(node, originalParent);
                    TreeNode::SetParent(neigh, neighParent);
                }
            }
            if (bestNewParent) {
                assert(which_node != 0);
                candidateMoves.emplace_back(best_delay_delta, node, bestNewParent, which_node);
            }
            tmp_steinerNode = nullptr;
            tmp_steinerNode.reset();
        }
        if (candidateMoves.empty()) {
            break;
        }

        // Try candidate moves in the order of descending wire length savings
        // Note that earlier moves may influence the legality of later one
        sort(candidateMoves.begin(), candidateMoves.end(), [](const MoveT& lhs, const MoveT& rhs) {
            return get<0>(lhs) < get<0>(rhs);
        });
        for (const auto& move : candidateMoves) {
            auto node = get<1>(move), neigh = get<2>(move);
            auto which_node = get<3>(move);
            auto neighParent = neigh->parent;
            // check due to earlier moves
            if (TreeNode::IsAncestor(node, neighParent)) continue;

            DTYPE pathLengthDelta =
                pathLengths[neighParent->id] + Dist(node->loc, neighParent->loc) - pathLengths[node->id];
            if (pathLengthDelta > slacks[node->id]) continue;

            auto steinerPt = GetNearestPoint(node, neigh);

            auto originalParent = node->parent;
            std::shared_ptr<salt::TreeNode> tmp_steinerNode = nullptr;
            if (which_node == 1) {
                TreeNode::ResetParent(node);
                TreeNode::SetParent(node, neigh);

                // pl不能违规
                DTYPE pathLengthDelta = pathLengths[neigh->id] + Dist(node->loc, neigh->loc) - pathLengths[node->id];
                if (pathLengthDelta > slacks[node->id]) continue;
            } else if (which_node == 2) {
                TreeNode::ResetParent(node);
                TreeNode::SetParent(node, neighParent);
            } else {
                tmp_steinerNode = make_shared<TreeNode>(steinerPt);
                TreeNode::SetParent(tmp_steinerNode, neighParent);
                TreeNode::ResetParent(neigh);
                TreeNode::SetParent(neigh, tmp_steinerNode);
                TreeNode::ResetParent(node);
                TreeNode::SetParent(node, tmp_steinerNode);
                // for later moves
                tmp_steinerNode->id = nodes.size();
            }

            salt::ElmoreDelayEval delayEval(rd, tree, for_one_net_maxLB);
            auto sumDelay = delayEval.sumNorDelay;

            if (which_node == 1) {
                TreeNode::ResetParent(node);
                TreeNode::SetParent(node, originalParent);
            } else if (which_node == 2) {
                TreeNode::ResetParent(node);
                TreeNode::SetParent(node, originalParent);
            } else {
                TreeNode::ResetParent(tmp_steinerNode);
                TreeNode::ResetParent(node);
                TreeNode::ResetParent(neigh);
                TreeNode::SetParent(node, originalParent);
                TreeNode::SetParent(neigh, neighParent);
            }

            tmp_steinerNode.reset();

            if (sumDelay + 1e-12 < before_sumDelay) {
                before_sumDelay = sumDelay;
            } else {
                continue;
            }

            // break
            Disconnect(node);
            // reroot
            if (which_node == 1) {
                Connect(node, neigh);
            } else if (which_node == 2) {
                Connect(node, neighParent);
            } else {
                auto steinerNode = make_shared<TreeNode>(steinerPt);
                Connect(steinerNode, neighParent);
                Disconnect(neigh);
                Connect(neigh, steinerNode);
                Connect(node, steinerNode);
                // for later moves
                steinerNode->id = nodes.size();
                nodes.push_back(steinerNode);
                pathLengths.push_back(pathLengths[neighParent->id] + steinerNode->WireToParent());
                slacks.push_back(Dist(steinerNode->loc, tree.source->loc) * (1 + eps) - pathLengths.back());
            }
            // update slack for later moves: first subtree, then path to source
            TreeNode::PreOrder(neighParent, UpdatePathLengths);
            TreeNode::PostOrder(neighParent, UpdateSlacks);

            tree.UpdateId();

            auto tmp = neighParent;
            while (tmp->parent) {
                slacks[tmp->parent->id] = min(slacks[tmp->parent->id], slacks[tmp->id]);
                tmp = tmp->parent;
            }
        }

        // Finalize
        tree.RemoveTopoRedundantSteiner();
        tree.PostOrderCopy([&](const shared_ptr<TreeNode>& node) {
            // degree may change after post-order traversal of its children
            if (node->pin) return;
            if (node->children.empty()) {
                Disconnect(node);
            } else if (node->children.size() == 1) {
                auto oldParent = node->parent, oldChild = node->children[0];
                Disconnect(node);
                Disconnect(oldChild);
                Connect(oldChild, oldParent);
            }
        });
    }
}

void Refine::DAES_C(Tree& tree, double eps, double rd, bool useRTree) {
    bgi::rtree<RNode, bgi::rstar<8>, bgi::indexable<RNode>, RNodeComp> rtree;
    if (useRTree) {
        tree.PostOrder([&](const shared_ptr<TreeNode>& n) {
            if (n->parent) {
                BBox s;
                bg::envelope(BSegment(BPoint(n->loc.x, n->loc.y), BPoint(n->parent->loc.x, n->parent->loc.y)), s);
                rtree.insert({s, n});
            }
        });
    }
    auto Disconnect = [&](const shared_ptr<TreeNode>& n) {
        if (useRTree) {
            BBox s;
            bg::envelope(BSegment(BPoint(n->loc.x, n->loc.y), BPoint(n->parent->loc.x, n->parent->loc.y)), s);
            rtree.remove({s, n});
        }
        TreeNode::ResetParent(n);
    };
    auto Connect = [&](const shared_ptr<TreeNode>& n, const shared_ptr<TreeNode>& parent) {
        TreeNode::SetParent(n, parent);
        if (useRTree) {
            BBox s;
            bg::envelope(BSegment(BPoint(n->loc.x, n->loc.y), BPoint(parent->loc.x, parent->loc.y)), s);
            rtree.insert({s, n});
        }
    };

    // 预备
    salt::ElmoreDelayEval before_delayEval(rd, tree);
    auto before_sumDelay = before_delayEval.sumNorDelay;

    auto for_one_net_maxLB = before_delayEval._maxLb;  // 用来正则化，不需要重复计算，这里直接保存一份

    int iter_num = 0;
    int max_iter = 2;

    while (true) {
        if (iter_num == max_iter) break;
        iter_num++;
        // Get nearest neighbors
        int num = tree.UpdateId();
        vector<shared_ptr<TreeNode>> nodes = tree.ObtainNodes(),
                                     orderedNodes(nodes.size());  // note: all pins should be covered

        vector<Point> points(nodes.size());
        for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
            orderedNodes[nodes[i]->id] = nodes[i];
            points[nodes[i]->id] = nodes[i]->loc;
        }
        nodes = orderedNodes;
        vector<vector<int>> nearestNeighbors;
        if (!useRTree) {
            MstBuilder mstB;
            mstB.GetAllNearestNeighbors(points, nearestNeighbors);
        } else {
            nearestNeighbors.resize(nodes.size());
            for (auto n : nodes) {
                if (n->parent) {
                    Point c = n->loc;  // center
                    DTYPE radius = n->WireToParent();
                    int tora = 1;
                    BBox queryBox{{c.x - tora * radius, c.y - tora * radius},
                                  {c.x + tora * radius, c.y + tora * radius}};
                    vector<RNode> cands;
                    rtree.query(bgi::intersects(queryBox), back_inserter(cands));
                    for (const auto& cand : cands) {
                        if (cand.second->id >= num) {
                            continue;
                        }
                        if (n->parent->id == cand.second->id) {
                            continue;
                        }
                        nearestNeighbors[n->id].push_back(cand.second->id);
                    }
                }
            }
        }

        // Prune descendants in nearest neighbors
        vector<int> preOrderIdxes(nodes.size(), -1);
        int globalPreOrderIdx = 0;
        function<void(const shared_ptr<TreeNode>&)> removeDescendants = [&](const shared_ptr<TreeNode>& node) {
            preOrderIdxes[node->id] = globalPreOrderIdx++;
            for (auto child : node->children) {
                removeDescendants(child);
            }
            for (auto& neighIdx : nearestNeighbors[node->id]) {
                int neighPreOrderIdx = preOrderIdxes[neighIdx];
                if (neighPreOrderIdx != -1 && neighPreOrderIdx >= preOrderIdxes[node->id]) {
                    neighIdx = -1;  // -1 stands for "descendant"
                }
            }
        };
        removeDescendants(tree.source);

        // Init path lengths and subtree slacks
        vector<DTYPE> pathLengths(nodes.size());
        vector<DTYPE> slacks(nodes.size());
        // 所有slack违例的点的id
        vector<int> slack_violation_nodes_id;  //////////////
        auto UpdatePathLengths = [&](const shared_ptr<TreeNode>& node) {
            if (node->pin) {
                if (node->pin->slack < 0.0) {
                    slack_violation_nodes_id.push_back(node->id);
                }
            }
            if (node->parent) {
                pathLengths[node->id] = pathLengths[node->parent->id] + node->WireToParent();
            } else {
                pathLengths[node->id] = 0;
            }
        };
        auto UpdateSlacks = [&](const shared_ptr<TreeNode>& node) {
            if (node->children.empty()) {
                slacks[node->id] =
                    Dist(node->loc, tree.source->loc) * (1 + eps) - pathLengths[node->id];  // floor here...
            } else {
                DTYPE minSlack = Dist(node->loc, tree.source->loc) * (1 + eps) - pathLengths[node->id];
                for (auto child : node->children) {
                    minSlack = min(minSlack, slacks[child->id]);
                }
                slacks[node->id] = minSlack;
            }
        };
        tree.PreOrder(UpdatePathLengths);
        tree.PostOrder(UpdateSlacks);

        // Find legal candidate moves
        using MoveT = tuple<double,
                            shared_ptr<TreeNode>,
                            shared_ptr<TreeNode>,
                            int>;  // 最后一个的int，1表示选neigh，2表示选neighParent，3表示选steinerPt
        int which_node = 0;
        vector<MoveT> candidateMoves;  // <wireLengthDelta, node, newParent>
        auto GetNearestPoint = [](const shared_ptr<TreeNode>& target, const shared_ptr<TreeNode>& neigh) {
            Box box(neigh->loc, neigh->parent->loc);
            box.Legalize();
            return box.GetNearestPointTo(target->loc);
        };

        for (auto node : nodes) {
            std::shared_ptr<salt::TreeNode> tmp_steinerNode = nullptr;
            if (!(node->parent)) {
                continue;
            }

            auto originalParent = node->parent;

            double best_delay_delta = 0;

            shared_ptr<TreeNode> bestNewParent;
            for (int neighIdx : nearestNeighbors[node->id]) {
                if (neighIdx >= static_cast<int>(nodes.size())) continue;
                if (neighIdx == -1 || !nodes[neighIdx]->parent) continue;
                auto neigh = nodes[neighIdx];
                auto neighParent = neigh->parent;

                if (node == neighParent) continue;
                if (node == neigh) continue;

                DTYPE pathLengthDelta =
                    pathLengths[neighParent->id] + Dist(node->loc, neighParent->loc) - pathLengths[node->id];
                if (pathLengthDelta > slacks[node->id]) continue;

                auto steinerPt = GetNearestPoint(node, neigh);

                bool exct_twice = (steinerPt == neigh->loc) || (steinerPt == neighParent->loc) ? 1 : 0;
                vector<shared_ptr<TreeNode>> tmp_nodes;
                tmp_nodes.emplace_back(neigh);
                tmp_nodes.emplace_back(neighParent);
                for (int i = 0; i < 2; i++) {
                    TreeNode::ResetParent(node);
                    TreeNode::SetParent(node, tmp_nodes[i]);

                    if (i == 0) {
                        // pl不能违规
                        DTYPE pathLengthDelta =
                            pathLengths[tmp_nodes[i]->id] + Dist(node->loc, tmp_nodes[i]->loc) - pathLengths[node->id];
                        if (pathLengthDelta > slacks[node->id]) continue;
                    }

                    salt::ElmoreDelayEval delayEval(rd, tree, for_one_net_maxLB);

                    double object = 0.0;
                    for (auto id : slack_violation_nodes_id) {
                        auto node = nodes[id];
                        auto node_delay = delayEval.nodeNorDelay[node->id];
                        auto node_slack = node->pin->slack;
                        auto len = Dist(node->loc, tree.net->source()->loc);
                        object += -len * node_slack * (node_slack / node_delay);
                    }
                    if (object < before_sumDelay) {
                        auto delay_delta = object - before_sumDelay;
                        if (delay_delta < best_delay_delta) {
                            best_delay_delta = delay_delta;
                            bestNewParent = neigh;
                            which_node = i + 1;
                        }
                    }

                    TreeNode::ResetParent(node);
                    TreeNode::SetParent(node, originalParent);
                }

                tmp_steinerNode = nullptr;
                if (!exct_twice) {
                    tmp_steinerNode = make_shared<TreeNode>(steinerPt);
                    TreeNode::SetParent(tmp_steinerNode, neighParent);
                    TreeNode::ResetParent(neigh);
                    TreeNode::SetParent(neigh, tmp_steinerNode);
                    TreeNode::ResetParent(node);
                    TreeNode::SetParent(node, tmp_steinerNode);
                    tmp_steinerNode->id = nodes.size();

                    salt::ElmoreDelayEval delayEval(rd, tree, for_one_net_maxLB);

                    double object = 0.0;
                    for (auto id : slack_violation_nodes_id) {
                        auto node = nodes[id];
                        auto node_delay = delayEval.nodeNorDelay[node->id];
                        auto node_slack = node->pin->slack;
                        auto len = Dist(node->loc, tree.net->source()->loc);
                        object += -len * node_slack * (node_slack / node_delay);
                        // object += len*(node_slack / node_delay);
                    }
                    if (object < before_sumDelay) {
                        auto delay_delta = object - before_sumDelay;
                        if (delay_delta < best_delay_delta) {
                            best_delay_delta = delay_delta;
                            bestNewParent = neigh;
                            which_node = 3;
                        }
                    }

                    TreeNode::ResetParent(tmp_steinerNode);
                    TreeNode::ResetParent(node);
                    TreeNode::ResetParent(neigh);
                    TreeNode::SetParent(node, originalParent);
                    TreeNode::SetParent(neigh, neighParent);
                }
            }
            if (bestNewParent) {
                assert(which_node != 0);
                candidateMoves.emplace_back(best_delay_delta, node, bestNewParent, which_node);
            }
            tmp_steinerNode = nullptr;
            tmp_steinerNode.reset();
        }
        if (candidateMoves.empty()) {
            break;
        }

        // Try candidate moves in the order of descending wire length savings
        // Note that earlier moves may influence the legality of later one
        sort(candidateMoves.begin(), candidateMoves.end(), [](const MoveT& lhs, const MoveT& rhs) {
            return get<0>(lhs) < get<0>(rhs);
        });
        for (const auto& move : candidateMoves) {
            auto node = get<1>(move), neigh = get<2>(move);
            auto which_node = get<3>(move);
            auto neighParent = neigh->parent;
            // check due to earlier moves
            if (TreeNode::IsAncestor(node, neighParent)) continue;

            DTYPE pathLengthDelta =
                pathLengths[neighParent->id] + Dist(node->loc, neighParent->loc) - pathLengths[node->id];
            if (pathLengthDelta > slacks[node->id]) continue;

            auto steinerPt = GetNearestPoint(node, neigh);

            auto originalParent = node->parent;
            std::shared_ptr<salt::TreeNode> tmp_steinerNode = nullptr;
            if (which_node == 1) {
                TreeNode::ResetParent(node);
                TreeNode::SetParent(node, neigh);

                // pl不能违规
                DTYPE pathLengthDelta = pathLengths[neigh->id] + Dist(node->loc, neigh->loc) - pathLengths[node->id];
                if (pathLengthDelta > slacks[node->id]) continue;
            } else if (which_node == 2) {
                TreeNode::ResetParent(node);
                TreeNode::SetParent(node, neighParent);
            } else {
                tmp_steinerNode = make_shared<TreeNode>(steinerPt);
                TreeNode::SetParent(tmp_steinerNode, neighParent);
                TreeNode::ResetParent(neigh);
                TreeNode::SetParent(neigh, tmp_steinerNode);
                TreeNode::ResetParent(node);
                TreeNode::SetParent(node, tmp_steinerNode);
                // for later moves
                tmp_steinerNode->id = nodes.size();
            }

            salt::ElmoreDelayEval delayEval(rd, tree, for_one_net_maxLB);
            auto sumDelay = delayEval.sumNorDelay;

            if (which_node == 1) {
                TreeNode::ResetParent(node);
                TreeNode::SetParent(node, originalParent);
            } else if (which_node == 2) {
                TreeNode::ResetParent(node);
                TreeNode::SetParent(node, originalParent);
            } else {
                TreeNode::ResetParent(tmp_steinerNode);
                TreeNode::ResetParent(node);
                TreeNode::ResetParent(neigh);
                TreeNode::SetParent(node, originalParent);
                TreeNode::SetParent(neigh, neighParent);
            }

            tmp_steinerNode.reset();

            if (sumDelay + 1e-12 < before_sumDelay) {
                before_sumDelay = sumDelay;
            } else {
                continue;
            }

            // break
            Disconnect(node);
            // reroot
            if (which_node == 1) {
                Connect(node, neigh);
            } else if (which_node == 2) {
                Connect(node, neighParent);
            } else {
                auto steinerNode = make_shared<TreeNode>(steinerPt);
                Connect(steinerNode, neighParent);
                Disconnect(neigh);
                Connect(neigh, steinerNode);
                Connect(node, steinerNode);
                // for later moves
                steinerNode->id = nodes.size();
                nodes.push_back(steinerNode);
                pathLengths.push_back(pathLengths[neighParent->id] + steinerNode->WireToParent());
                slacks.push_back(Dist(steinerNode->loc, tree.source->loc) * (1 + eps) - pathLengths.back());
            }
            // update slack for later moves: first subtree, then path to source
            TreeNode::PreOrder(neighParent, UpdatePathLengths);
            TreeNode::PostOrder(neighParent, UpdateSlacks);

            tree.UpdateId();

            auto tmp = neighParent;
            while (tmp->parent) {
                slacks[tmp->parent->id] = min(slacks[tmp->parent->id], slacks[tmp->id]);
                tmp = tmp->parent;
            }
        }

        // Finalize
        tree.RemoveTopoRedundantSteiner();
        // // cout << tree;
        tree.PostOrderCopy([&](const shared_ptr<TreeNode>& node) {
            // degree may change after post-order traversal of its children
            if (node->pin) return;
            if (node->children.empty()) {
                Disconnect(node);
            } else if (node->children.size() == 1) {
                auto oldParent = node->parent, oldChild = node->children[0];
                Disconnect(node);
                Disconnect(oldChild);
                Connect(oldChild, oldParent);
            }
        });
    }
}

void Refine::DAES_S_C(Tree& tree, double eps, double rd, bool useRTree) {
    bgi::rtree<RNode, bgi::rstar<8>, bgi::indexable<RNode>, RNodeComp> rtree;
    if (useRTree) {
        tree.PostOrder([&](const shared_ptr<TreeNode>& n) {
            if (n->parent) {
                BBox s;
                bg::envelope(BSegment(BPoint(n->loc.x, n->loc.y), BPoint(n->parent->loc.x, n->parent->loc.y)), s);
                rtree.insert({s, n});
            }
        });
    }
    auto Disconnect = [&](const shared_ptr<TreeNode>& n) {
        if (useRTree) {
            BBox s;
            bg::envelope(BSegment(BPoint(n->loc.x, n->loc.y), BPoint(n->parent->loc.x, n->parent->loc.y)), s);
            rtree.remove({s, n});
        }
        TreeNode::ResetParent(n);
    };
    auto Connect = [&](const shared_ptr<TreeNode>& n, const shared_ptr<TreeNode>& parent) {
        TreeNode::SetParent(n, parent);
        if (useRTree) {
            BBox s;
            bg::envelope(BSegment(BPoint(n->loc.x, n->loc.y), BPoint(parent->loc.x, parent->loc.y)), s);
            rtree.insert({s, n});
        }
    };

    // 预备
    salt::ElmoreDelayEval before_delayEval(rd, tree);
    auto before_sumDelay = before_delayEval.sumNorDelay;

    auto for_one_net_maxLB = before_delayEval._maxLb;

    int iter_num = 0;
    int max_iter = 5;

    while (true) {
        if (iter_num == max_iter) break;
        iter_num++;
        // Get nearest neighbors
        int num = tree.UpdateId();
        vector<shared_ptr<TreeNode>> nodes = tree.ObtainNodes(),
                                     orderedNodes(nodes.size());  // note: all pins should be covered
        vector<Point> points(nodes.size());
        for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
            orderedNodes[nodes[i]->id] = nodes[i];
            points[nodes[i]->id] = nodes[i]->loc;
        }
        nodes = orderedNodes;
        vector<vector<int>> nearestNeighbors;
        if (!useRTree) {
            MstBuilder mstB;
            mstB.GetAllNearestNeighbors(points, nearestNeighbors);
        } else {
            nearestNeighbors.resize(nodes.size());
            for (auto n : nodes) {
                if (n->parent) {
                    Point c = n->loc;  // center
                    DTYPE radius = n->WireToParent();
                    int tora = 1;
                    BBox queryBox{{c.x - tora * radius, c.y - tora * radius},
                                  {c.x + tora * radius, c.y + tora * radius}};
                    vector<RNode> cands;
                    rtree.query(bgi::intersects(queryBox), back_inserter(cands));
                    for (const auto& cand : cands) {
                        if (cand.second->id >= num) {
                            continue;
                        }
                        if (n->parent->id == cand.second->id) {
                            continue;
                        }
                        nearestNeighbors[n->id].push_back(cand.second->id);
                    }
                }
            }
        }

        // Prune descendants in nearest neighbors
        vector<int> preOrderIdxes(nodes.size(), -1);
        int globalPreOrderIdx = 0;
        function<void(const shared_ptr<TreeNode>&)> removeDescendants = [&](const shared_ptr<TreeNode>& node) {
            preOrderIdxes[node->id] = globalPreOrderIdx++;
            for (auto child : node->children) {
                removeDescendants(child);
            }
            for (auto& neighIdx : nearestNeighbors[node->id]) {
                int neighPreOrderIdx = preOrderIdxes[neighIdx];
                if (neighPreOrderIdx != -1 && neighPreOrderIdx >= preOrderIdxes[node->id]) {
                    neighIdx = -1;  // -1 stands for "descendant"
                }
            }
        };
        removeDescendants(tree.source);

        // Init path lengths and subtree slacks
        vector<DTYPE> pathLengths(nodes.size());
        vector<DTYPE> slacks(nodes.size());
        // 所有slack违例的点的id
        vector<int> slack_violation_nodes_id;
        auto UpdatePathLengths = [&](const shared_ptr<TreeNode>& node) {
            if (node->pin) {
                if (node->pin->slack < 0.0) {
                    slack_violation_nodes_id.push_back(node->id);
                }
            }
            if (node->parent) {
                pathLengths[node->id] = pathLengths[node->parent->id] + node->WireToParent();
            } else {
                pathLengths[node->id] = 0;
            }
        };
        auto UpdateSlacks = [&](const shared_ptr<TreeNode>& node) {
            if (node->children.empty()) {
                slacks[node->id] =
                    Dist(node->loc, tree.source->loc) * (1 + eps) - pathLengths[node->id];  // floor here...
            } else {
                DTYPE minSlack = Dist(node->loc, tree.source->loc) * (1 + eps) - pathLengths[node->id];
                for (auto child : node->children) {
                    minSlack = min(minSlack, slacks[child->id]);
                }
                slacks[node->id] = minSlack;
            }
        };
        tree.PreOrder(UpdatePathLengths);
        tree.PostOrder(UpdateSlacks);

        // Find legal candidate moves
        // using MoveT = tuple<DTYPE, shared_ptr<TreeNode>, shared_ptr<TreeNode>>;
        using MoveT = tuple<double, shared_ptr<TreeNode>, shared_ptr<TreeNode>>;
        vector<MoveT> candidateMoves;  // <wireLengthDelta, node, newParent>
        auto GetNearestPoint = [](const shared_ptr<TreeNode>& target, const shared_ptr<TreeNode>& neigh) {
            Box box(neigh->loc, neigh->parent->loc);
            box.Legalize();
            return box.GetNearestPointTo(target->loc);
        };

        for (auto node : nodes) {
            std::shared_ptr<salt::TreeNode> tmp_steinerNode = nullptr;
            if (!(node->parent)) {
                continue;
            }

            auto originalParent = node->parent;

            double best_delay_delta = 0;

            shared_ptr<TreeNode> bestNewParent;
            for (int neighIdx : nearestNeighbors[node->id]) {
                if (neighIdx >= static_cast<int>(nodes.size())) continue;
                if (neighIdx == -1 || !nodes[neighIdx]->parent) continue;
                auto neigh = nodes[neighIdx];
                auto neighParent = neigh->parent;

                if (node == neighParent) continue;
                if (node == neigh) continue;

                auto steinerPt = GetNearestPoint(node, neigh);

                // pl不能违规
                DTYPE pathLengthDelta =
                    pathLengths[neighParent->id] + Dist(node->loc, neighParent->loc) - pathLengths[node->id];
                if (pathLengthDelta > slacks[node->id]) continue;

                tmp_steinerNode = nullptr;
                if (steinerPt == neigh->loc) {
                    TreeNode::ResetParent(node);
                    TreeNode::SetParent(node, neigh);
                } else if (steinerPt == neighParent->loc) {
                    TreeNode::ResetParent(node);
                    TreeNode::SetParent(node, neighParent);
                } else {
                    tmp_steinerNode = make_shared<TreeNode>(steinerPt);
                    TreeNode::SetParent(tmp_steinerNode, neighParent);
                    TreeNode::ResetParent(neigh);
                    TreeNode::SetParent(neigh, tmp_steinerNode);
                    TreeNode::ResetParent(node);
                    TreeNode::SetParent(node, tmp_steinerNode);
                    // for later moves
                    tmp_steinerNode->id = nodes.size();
                }

                salt::ElmoreDelayEval delayEval(rd, tree, for_one_net_maxLB);

                // if (sumDelay + 1e-12 < before_sumDelay) {
                //     auto delay_delta = sumDelay - before_sumDelay;
                //     if (delay_delta < best_delay_delta) {
                //         best_delay_delta = delay_delta;
                //         bestNewParent = neigh;
                //     }
                // }

                double object = 0.0;
                for (auto id : slack_violation_nodes_id) {
                    auto node = nodes[id];
                    auto node_delay = delayEval.nodeNorDelay[node->id];
                    auto node_slack = node->pin->slack;
                    auto len = Dist(node->loc, tree.net->source()->loc);
                    object += -len * node_slack * (node_slack / node_delay);
                    // object += len*(node_slack / node_delay);
                }
                if (object < before_sumDelay) {
                    auto delay_delta = object - before_sumDelay;
                    if (delay_delta < best_delay_delta) {
                        best_delay_delta = delay_delta;
                        bestNewParent = neigh;
                    }
                }

                if (steinerPt == neigh->loc) {
                    TreeNode::ResetParent(node);
                    TreeNode::SetParent(node, originalParent);
                } else if (steinerPt == neighParent->loc) {
                    TreeNode::ResetParent(node);
                    TreeNode::SetParent(node, originalParent);
                } else {
                    TreeNode::ResetParent(tmp_steinerNode);
                    TreeNode::ResetParent(node);
                    TreeNode::ResetParent(neigh);
                    TreeNode::SetParent(node, originalParent);
                    TreeNode::SetParent(neigh, neighParent);
                }
            }
            if (bestNewParent) {
                candidateMoves.emplace_back(best_delay_delta, node, bestNewParent);
            }
            tmp_steinerNode = nullptr;
            tmp_steinerNode.reset();
        }
        if (candidateMoves.empty()) {
            break;
        }

        // Try candidate moves in the order of descending wire length savings
        // Note that earlier moves may influence the legality of later one
        sort(candidateMoves.begin(), candidateMoves.end(), [](const MoveT& lhs, const MoveT& rhs) {
            return get<0>(lhs) < get<0>(rhs);
        });
        for (const auto& move : candidateMoves) {
            auto node = get<1>(move), neigh = get<2>(move);
            auto neighParent = neigh->parent;
            // check due to earlier moves
            if (TreeNode::IsAncestor(node, neighParent)) continue;
            DTYPE pathLengthDelta =
                pathLengths[neighParent->id] + Dist(node->loc, neighParent->loc) - pathLengths[node->id];
            if (pathLengthDelta > slacks[node->id]) continue;
            auto steinerPt = GetNearestPoint(node, neigh);

            auto originalParent = node->parent;
            std::shared_ptr<salt::TreeNode> tmp_steinerNode = nullptr;
            if (steinerPt == neigh->loc) {
                TreeNode::ResetParent(node);
                TreeNode::SetParent(node, neigh);
            } else if (steinerPt == neighParent->loc) {
                TreeNode::ResetParent(node);
                TreeNode::SetParent(node, neighParent);
            } else {
                tmp_steinerNode = make_shared<TreeNode>(steinerPt);
                TreeNode::SetParent(tmp_steinerNode, neighParent);
                TreeNode::ResetParent(neigh);
                TreeNode::SetParent(neigh, tmp_steinerNode);
                TreeNode::ResetParent(node);
                TreeNode::SetParent(node, tmp_steinerNode);
                // for later moves
                tmp_steinerNode->id = nodes.size();
            }

            salt::ElmoreDelayEval delayEval(rd, tree, for_one_net_maxLB);
            auto sumDelay = delayEval.sumNorDelay;

            if (steinerPt == neigh->loc) {
                TreeNode::ResetParent(node);
                TreeNode::SetParent(node, originalParent);
            } else if (steinerPt == neighParent->loc) {
                TreeNode::ResetParent(node);
                TreeNode::SetParent(node, originalParent);
            } else {
                TreeNode::ResetParent(tmp_steinerNode);
                TreeNode::ResetParent(node);
                TreeNode::ResetParent(neigh);
                TreeNode::SetParent(node, originalParent);
                TreeNode::SetParent(neigh, neighParent);
            }

            tmp_steinerNode.reset();

            if (sumDelay + 1e-12 < before_sumDelay) {
                before_sumDelay = sumDelay;
            } else {
                continue;
            }

            // break
            Disconnect(node);
            // reroot
            if (steinerPt == neigh->loc) {
                Connect(node, neigh);
            } else if (steinerPt == neighParent->loc) {
                Connect(node, neighParent);
            } else {
                auto steinerNode = make_shared<TreeNode>(steinerPt);
                Connect(steinerNode, neighParent);
                Disconnect(neigh);
                Connect(neigh, steinerNode);
                Connect(node, steinerNode);
                // for later moves
                steinerNode->id = nodes.size();
                nodes.push_back(steinerNode);
                pathLengths.push_back(pathLengths[neighParent->id] + steinerNode->WireToParent());
                slacks.push_back(Dist(steinerNode->loc, tree.source->loc) * (1 + eps) - pathLengths.back());
            }
            // update slack for later moves: first subtree, then path to source
            TreeNode::PreOrder(neighParent, UpdatePathLengths);
            TreeNode::PostOrder(neighParent, UpdateSlacks);

            tree.UpdateId();

            auto tmp = neighParent;
            while (tmp->parent) {
                slacks[tmp->parent->id] = min(slacks[tmp->parent->id], slacks[tmp->id]);
                tmp = tmp->parent;
            }
        }

        // Finalize
        tree.RemoveTopoRedundantSteiner();
        // // cout << tree;
        tree.PostOrderCopy([&](const shared_ptr<TreeNode>& node) {
            // degree may change after post-order traversal of its children
            if (node->pin) return;
            if (node->children.empty()) {
                Disconnect(node);
            } else if (node->children.size() == 1) {
                auto oldParent = node->parent, oldChild = node->children[0];
                Disconnect(node);
                Disconnect(oldChild);
                Connect(oldChild, oldParent);
            }
        });
    }
}
}  // namespace salt