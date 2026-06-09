// NonDuplicationRecombiner implementation.
//
// Direct C++ port of benchmark/grg_recombination.py's NonDuplicationRecombination.
// Method semantics are preserved exactly so parity tests against the Python
// reference can drive correctness verification. See benchmark/native/DESIGN.md
// for the full porting decisions and benchmark/recombination_time_breakdown.md
// for the perf motivation.

#include "grgl/recombination.h"

#include "grgl/common.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <deque>
#include <limits>
#include <stdexcept>
#include <utility>

namespace grgl {

namespace {

// Sentinel used during span/anc-cov accumulation. Real BpPosition values are
// uint64_t, so we use the max as "no value yet" and accumulate via min/max.
constexpr BpPosition NO_POS_SENTINEL = std::numeric_limits<BpPosition>::max();

inline double elapsedSeconds(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

} // namespace

NonDuplicationRecombiner::NonDuplicationRecombiner(MutableGRGPtr grg, bool instrument)
    : m_grg(std::move(grg)),
      m_genomeLength(0),
      m_genVisited(0),
      m_genConnected(0),
      m_audit(),
      m_stats(),
      m_instrument(instrument),
      m_debug(false),
      m_deferSampleUpdates(false) {
    if (!m_grg) {
        throw std::invalid_argument("NonDuplicationRecombiner: null MutableGRGPtr");
    }
    m_genomeLength = m_grg->getBPRange().second;

    const size_t n = m_grg->numNodes();
    m_spanCache.resize(n);   // default-constructed -> Uninit
    m_ancCovCache.resize(n); // default-constructed -> Uninit
    m_visitedGen.resize(n, 0);
    m_connectedGen.resize(n, 0);

    buildAncestralCaches();
}

// One-pass O(V + E) Kahn's-sort initialization of every per-node cache.
// Direct port of Python _build_ancestral_caches; see its docstring for the
// algorithm and correctness argument.
//
// After this returns, every node in the input graph has:
//   m_spanCache[node]   = Set(min_pos, max_pos) over node's own + ancestors' mutations,
//                          or Empty if neither has any mutations
//   m_ancCovCache[node] = Set(anc_min, anc_max + 1) over parents' span unions,
//                          or Empty if node has no parents (root)
//   m_upEdgesCache[node], m_mutationCache[node], m_posCache[node] populated
//
// Lazy fallbacks in getNodeAndAncestorSpan / getAncestralCoverage handle
// bubbles + offspring created mid-recombination.
void NonDuplicationRecombiner::buildAncestralCaches() {
    std::chrono::steady_clock::time_point tStart;
    if (m_instrument) {
        tStart = std::chrono::steady_clock::now();
    }

    const size_t n = m_grg->numNodes();
    if (n == 0) {
        if (m_instrument) {
            m_stats.initCachesTime += elapsedSeconds(tStart);
        }
        return;
    }

    // Pass 1: populate _up_edges_cache and compute in-degrees.
    std::vector<uint32_t> inDegree(n, 0);
    for (NodeID node = 0; node < n; node++) {
        NodeIDList parents = m_grg->getUpEdges(node);
        const size_t fanin = parents.size();
        m_upEdgesCache[node] = std::vector<NodeID>(parents.begin(), parents.end());
        inDegree[node] = static_cast<uint32_t>(fanin);
    }

    // Seed: roots (no parents) are immediately processable.
    std::deque<NodeID> queue;
    for (NodeID node = 0; node < n; node++) {
        if (inDegree[node] == 0) {
            queue.push_back(node);
        }
    }

    // Pass 2: topological processing.
    while (!queue.empty()) {
        const NodeID node = queue.front();
        queue.pop_front();

        // ---- Populate mutation caches (inlined from getNodeMutations) ----
        std::vector<MutationId> mutIds = m_grg->getMutationsForNode<MutationId>(node, /*allowSort=*/false);
        std::vector<std::pair<MutationId, BpPosition>> mutations;
        mutations.reserve(mutIds.size());
        for (MutationId mid : mutIds) {
            const Mutation& mut = m_grg->getMutationById(mid);
            mutations.emplace_back(mid, mut.getPosition());
        }
        if (mutations.size() > 1) {
            std::sort(mutations.begin(),
                      mutations.end(),
                      [](const std::pair<MutationId, BpPosition>& a, const std::pair<MutationId, BpPosition>& b) {
                          return a.second < b.second;
                      });
        }
        std::vector<BpPosition> positions;
        positions.reserve(mutations.size());
        for (const auto& m : mutations) {
            positions.push_back(m.second);
        }

        // ---- Compute span_cache[node] = own muts (min,max) U each parent's span ----
        BpPosition minP = NO_POS_SENTINEL;
        BpPosition maxP = 0;
        bool anySet = false;
        if (!mutations.empty()) {
            minP = mutations.front().second;
            maxP = mutations.back().second;
            anySet = true;
        }

        const std::vector<NodeID>& parents = m_upEdgesCache[node];
        for (NodeID parent : parents) {
            const IntervalOpt& pSpan = m_spanCache[parent];
            if (pSpan.isSet()) {
                if (!anySet) {
                    minP = pSpan.lo;
                    maxP = pSpan.hi;
                    anySet = true;
                } else {
                    if (pSpan.lo < minP) {
                        minP = pSpan.lo;
                    }
                    if (pSpan.hi > maxP) {
                        maxP = pSpan.hi;
                    }
                }
            }
        }
        m_spanCache[node] = anySet ? IntervalOpt::makeSet(minP, maxP) : IntervalOpt::makeEmpty();

        // ---- Compute anc_cov_cache[node] = union of parents' span (with +1 on max) ----
        // Empty (Python None) if node has no parents (root).
        if (parents.empty()) {
            m_ancCovCache[node] = IntervalOpt::makeEmpty();
        } else {
            BpPosition ancMin = NO_POS_SENTINEL;
            BpPosition ancMax = 0;
            bool anyAnc = false;
            for (NodeID parent : parents) {
                const IntervalOpt& pSpan = m_spanCache[parent];
                if (pSpan.isSet()) {
                    if (!anyAnc) {
                        ancMin = pSpan.lo;
                        ancMax = pSpan.hi;
                        anyAnc = true;
                    } else {
                        if (pSpan.lo < ancMin) {
                            ancMin = pSpan.lo;
                        }
                        if (pSpan.hi > ancMax) {
                            ancMax = pSpan.hi;
                        }
                    }
                }
            }
            m_ancCovCache[node] = anyAnc ? IntervalOpt::makeSet(ancMin, ancMax + 1) : IntervalOpt::makeEmpty();
        }

        m_mutationCache[node] = std::move(mutations);
        m_posCache[node] = std::move(positions);

        // ---- Release children whose final parent just finished ----
        NodeIDList children = m_grg->getDownEdges(node);
        for (NodeID child : children) {
            inDegree[child]--;
            if (inDegree[child] == 0) {
                queue.push_back(child);
            }
        }
    }

    if (m_instrument) {
        m_stats.initCachesTime += elapsedSeconds(tStart);
    }
}

// ------------------------------------------------------------------
// Per-node array growth (mirrors Python _grow_node_arrays / _sync_to_grg)
// ------------------------------------------------------------------

void NonDuplicationRecombiner::growNodeArrays(NodeID nodeId) {
    const size_t target = static_cast<size_t>(nodeId) + 1;
    if (m_spanCache.size() < target) {
        m_spanCache.resize(target);   // default-constructed -> Uninit
        m_ancCovCache.resize(target); // default-constructed -> Uninit
    }
    if (m_visitedGen.size() < target) {
        m_visitedGen.resize(target, 0);
        m_connectedGen.resize(target, 0);
    }
}

void NonDuplicationRecombiner::syncToGrg() {
    const size_t n = m_grg->numNodes();
    if (n > 0) {
        growNodeArrays(static_cast<NodeID>(n - 1));
    }
}

// ------------------------------------------------------------------
// Cached graph access (lazy fallbacks for nodes created mid-recombination)
// ------------------------------------------------------------------

const std::vector<NodeID>& NonDuplicationRecombiner::getUpEdgesCached(NodeID nodeId) {
    auto it = m_upEdgesCache.find(nodeId);
    if (it != m_upEdgesCache.end()) {
        return it->second;
    }
    NodeIDList parents = m_grg->getUpEdges(nodeId);
    auto inserted = m_upEdgesCache.emplace(nodeId, std::vector<NodeID>(parents.begin(), parents.end()));
    return inserted.first->second;
}

const std::vector<std::pair<MutationId, BpPosition>>& NonDuplicationRecombiner::getNodeMutations(NodeID nodeId) {
    auto it = m_mutationCache.find(nodeId);
    if (it != m_mutationCache.end()) {
        return it->second;
    }

    std::vector<MutationId> mutIds = m_grg->getMutationsForNode<MutationId>(nodeId, /*allowSort=*/false);
    std::vector<std::pair<MutationId, BpPosition>> mutations;
    mutations.reserve(mutIds.size());
    for (MutationId mid : mutIds) {
        const Mutation& mut = m_grg->getMutationById(mid);
        mutations.emplace_back(mid, mut.getPosition());
    }
    // p50 mutations/node is 0 across all benchmarked files; skip the sort on
    // the 0/1-item common case to match the Python optimization.
    if (mutations.size() > 1) {
        std::sort(mutations.begin(),
                  mutations.end(),
                  [](const std::pair<MutationId, BpPosition>& a, const std::pair<MutationId, BpPosition>& b) {
                      return a.second < b.second;
                  });
    }
    std::vector<BpPosition> positions;
    positions.reserve(mutations.size());
    for (const auto& m : mutations) {
        positions.push_back(m.second);
    }

    m_posCache[nodeId] = std::move(positions);
    auto inserted = m_mutationCache.emplace(nodeId, std::move(mutations));
    return inserted.first->second;
}

// ------------------------------------------------------------------
// Iterative span / ancestral coverage (lazy DFS)
// ------------------------------------------------------------------

IntervalOpt NonDuplicationRecombiner::getNodeAndAncestorSpan(NodeID nodeId) {
    if (!m_spanCache[nodeId].isUninit()) {
        return m_spanCache[nodeId];
    }

    // Iterative post-order DFS — mirrors Python's two-phase stack. Each entry
    // is (nid, processed). processed=false means "push children, then re-push
    // self with processed=true"; processed=true means "all parents resolved,
    // compute and store own span_cache".
    std::vector<std::pair<NodeID, bool>> stack;
    std::unordered_set<NodeID> scheduled;
    stack.emplace_back(nodeId, false);
    scheduled.insert(nodeId);

    while (!stack.empty()) {
        const std::pair<NodeID, bool> entry = stack.back();
        stack.pop_back();
        const NodeID nid = entry.first;
        const bool processed = entry.second;

        if (processed) {
            BpPosition minP = NO_POS_SENTINEL;
            BpPosition maxP = 0;
            bool anySet = false;
            const auto& muts = getNodeMutations(nid);
            if (!muts.empty()) {
                minP = muts.front().second;
                maxP = muts.back().second;
                anySet = true;
            }
            for (NodeID parent : getUpEdgesCached(nid)) {
                const IntervalOpt& pSpan = m_spanCache[parent];
                if (pSpan.isSet()) {
                    if (!anySet) {
                        minP = pSpan.lo;
                        maxP = pSpan.hi;
                        anySet = true;
                    } else {
                        if (pSpan.lo < minP) {
                            minP = pSpan.lo;
                        }
                        if (pSpan.hi > maxP) {
                            maxP = pSpan.hi;
                        }
                    }
                }
            }
            m_spanCache[nid] = anySet ? IntervalOpt::makeSet(minP, maxP) : IntervalOpt::makeEmpty();
            scheduled.erase(nid);
        } else {
            stack.emplace_back(nid, true);
            for (NodeID parent : getUpEdgesCached(nid)) {
                if (m_spanCache[parent].isUninit() && scheduled.find(parent) == scheduled.end()) {
                    stack.emplace_back(parent, false);
                    scheduled.insert(parent);
                }
            }
        }
    }

    return m_spanCache[nodeId];
}

IntervalOpt NonDuplicationRecombiner::getAncestralCoverage(NodeID nodeId) {
    if (!m_ancCovCache[nodeId].isUninit()) {
        return m_ancCovCache[nodeId];
    }

    const std::vector<NodeID>& parents = getUpEdgesCached(nodeId);
    if (parents.empty()) {
        m_ancCovCache[nodeId] = IntervalOpt::makeEmpty();
        return m_ancCovCache[nodeId];
    }

    BpPosition minPos = NO_POS_SENTINEL;
    BpPosition maxPos = 0;
    bool anySet = false;
    for (NodeID parent : parents) {
        IntervalOpt pSpan = getNodeAndAncestorSpan(parent);
        if (pSpan.isSet()) {
            if (!anySet) {
                minPos = pSpan.lo;
                maxPos = pSpan.hi;
                anySet = true;
            } else {
                if (pSpan.lo < minPos) {
                    minPos = pSpan.lo;
                }
                if (pSpan.hi > maxPos) {
                    maxPos = pSpan.hi;
                }
            }
        }
    }
    m_ancCovCache[nodeId] = anySet ? IntervalOpt::makeSet(minPos, maxPos + 1) : IntervalOpt::makeEmpty();
    return m_ancCovCache[nodeId];
}

// ------------------------------------------------------------------
// Bubble extraction
// ------------------------------------------------------------------

NodeID
NonDuplicationRecombiner::extractBubble(NodeID nodeId, const std::vector<MutationId>& relMutIds, NodeID offspringId) {
    std::chrono::steady_clock::time_point t;
    if (m_instrument) {
        t = std::chrono::steady_clock::now();
    }
    const NodeID bubbleId = m_grg->makeNode();
    if (m_instrument) {
        m_stats.makeNodeTime += elapsedSeconds(t);
        m_stats.makeNodeCalls += 1;
    }
    m_audit.makeNodeCalls += 1;

    growNodeArrays(bubbleId);

    if (m_instrument) {
        t = std::chrono::steady_clock::now();
    }
    m_grg->connect(static_cast<SignedNodeID>(bubbleId), static_cast<SignedNodeID>(nodeId));
    m_grg->connect(static_cast<SignedNodeID>(bubbleId), -static_cast<SignedNodeID>(offspringId));
    if (m_instrument) {
        m_stats.connectTime += elapsedSeconds(t);
        m_stats.connectCalls += 2;
    }
    m_audit.connectCallsInExtract += 2;
    m_audit.extractBubbleCalls += 1;

    // node_id just gained a new up-edge; drop its cached up-edges list.
    m_upEdgesCache.erase(nodeId);

    BubbleOp op;
    op.nodeId = nodeId;
    op.bubbleId = bubbleId;
    op.relMutIds = relMutIds;
    m_pendingBubbles.push_back(std::move(op));

    // Track only nodes whose own caches actually become stale:
    // node_id loses mutations (mutation/pos/anc_cov are stale) and the
    // new bubble_id needs a clean slot. node_id's parents are NOT
    // invalidated -- their span and anc_cov depend only on themselves
    // and their own ancestors, neither of which changes here.
    m_modifiedNodes.insert(nodeId);
    m_modifiedNodes.insert(bubbleId);

    return bubbleId;
}

// ------------------------------------------------------------------
// Apply deferred work + cache eviction
// ------------------------------------------------------------------

void NonDuplicationRecombiner::applyPendingBubbles() {
    std::chrono::steady_clock::time_point tApplyStart;
    double tGet = 0.0;
    double tAdd = 0.0;
    double tRem = 0.0;
    size_t nMutsTotal = 0;
    const size_t nBubbles = m_pendingBubbles.size();
    if (m_instrument) {
        tApplyStart = std::chrono::steady_clock::now();
    }

    for (const BubbleOp& op : m_pendingBubbles) {
        if (m_instrument) {
            nMutsTotal += op.relMutIds.size();
        }
        for (MutationId mid : op.relMutIds) {
            std::chrono::steady_clock::time_point t;
            if (m_instrument) {
                t = std::chrono::steady_clock::now();
            }
            const Mutation& mut = m_grg->getMutationById(mid);
            if (m_instrument) {
                tGet += elapsedSeconds(t);
                t = std::chrono::steady_clock::now();
            }
            m_grg->addMutation(mut, op.bubbleId);
            if (m_instrument) {
                tAdd += elapsedSeconds(t);
                t = std::chrono::steady_clock::now();
            }
            m_grg->removeMutation(mid, op.nodeId);
            if (m_instrument) {
                tRem += elapsedSeconds(t);
            }
        }
    }
    m_pendingBubbles.clear();

    if (!m_deferSampleUpdates) {
        flushSampleUpdates();
    }

    if (m_instrument) {
        const double elapsed = elapsedSeconds(tApplyStart);
        m_stats.applyBubblesTime += elapsed;
        m_stats.getMutationByIdTime += tGet;
        m_stats.addMutationTime += tAdd;
        m_stats.removeMutationTime += tRem;
        m_stats.getMutationByIdCalls += nMutsTotal;
        m_stats.addMutationCalls += nMutsTotal;
        m_stats.removeMutationCalls += nMutsTotal;
        m_stats.bubblesCreated += nBubbles;
        m_stats.mutationsMoved += nMutsTotal;
    }
}

void NonDuplicationRecombiner::flushSampleUpdates() {
    if (m_pendingSampleRemovals.empty()) {
        return;
    }
    std::chrono::steady_clock::time_point t;
    if (m_instrument) {
        t = std::chrono::steady_clock::now();
    }
    NodeIDList currentList = m_grg->getSampleNodes();
    std::unordered_set<NodeID> current(currentList.begin(), currentList.end());
    for (NodeID rm : m_pendingSampleRemovals) {
        current.erase(rm);
    }
    NodeIDList nextSamples(current.begin(), current.end());
    m_grg->setSamples(nextSamples);
    m_pendingSampleRemovals.clear();
    if (m_instrument) {
        m_stats.flushSamplesTime += elapsedSeconds(t);
    }
}

void NonDuplicationRecombiner::clearModifiedCaches() {
    std::chrono::steady_clock::time_point t;
    if (m_instrument) {
        t = std::chrono::steady_clock::now();
    }
    for (NodeID nodeId : m_modifiedNodes) {
        m_mutationCache.erase(nodeId);
        m_posCache.erase(nodeId);
        if (static_cast<size_t>(nodeId) < m_spanCache.size()) {
            m_spanCache[nodeId] = IntervalOpt();   // back to Uninit
            m_ancCovCache[nodeId] = IntervalOpt(); // back to Uninit
        }
    }
    m_modifiedNodes.clear();
    if (m_instrument) {
        m_stats.clearCachesTime += elapsedSeconds(t);
    }
}

SignedNodeID NonDuplicationRecombiner::registerOffspring(NodeID offspringId) {
    auto it = m_negativeNodeIndex.find(offspringId);
    size_t idx;
    if (it == m_negativeNodeIndex.end()) {
        idx = m_negativeNodeIds.size();
        m_negativeNodeIndex.emplace(offspringId, idx);
        m_negativeNodeIds.push_back(offspringId);
    } else {
        idx = it->second;
    }
    return -static_cast<SignedNodeID>(idx + 1);
}

// ------------------------------------------------------------------
// Public entry points (recurseAttach + recombine* in subsequent tasks)
// ------------------------------------------------------------------

SignedNodeID NonDuplicationRecombiner::recombine(NodeID hapA, NodeID hapB, BpPosition breakpoint) {
    m_audit.recombineCalls++;
    m_pendingBubbles.clear();

    std::chrono::steady_clock::time_point t;
    const uint64_t genVBefore = m_genVisited;
    if (m_instrument) {
        t = std::chrono::steady_clock::now();
    }
    syncToGrg();
    if (m_instrument) {
        m_stats.syncToGrgTime += elapsedSeconds(t);
        t = std::chrono::steady_clock::now();
    }
    const NodeID offspringId = m_grg->makeNode(1, /*negative=*/true);
    if (m_instrument) {
        m_stats.makeNodeTime += elapsedSeconds(t);
        m_stats.makeNodeCalls += 1;
    }
    m_audit.makeNodeCalls++;
    growNodeArrays(static_cast<NodeID>(m_grg->numNodes() - 1));
    m_genConnected++;

    double recurseT = 0.0;
    std::chrono::steady_clock::time_point ts;
    if (m_instrument) {
        ts = std::chrono::steady_clock::now();
    }
    recurseAttach(hapA, offspringId, 0, breakpoint);
    if (m_instrument) {
        recurseT += elapsedSeconds(ts);
        m_stats.recurseAttachCalls += 1;
        m_stats.segmentsProcessed += 1;
        ts = std::chrono::steady_clock::now();
    }
    recurseAttach(hapB, offspringId, breakpoint, m_genomeLength);
    if (m_instrument) {
        recurseT += elapsedSeconds(ts);
        m_stats.recurseAttachCalls += 1;
        m_stats.segmentsProcessed += 1;
        m_stats.recurseAttachTime += recurseT;
    }

    applyPendingBubbles();
    clearModifiedCaches();

    if (m_instrument) {
        const uint64_t genVAfter = m_genVisited;
        uint64_t v = 0;
        for (uint64_t gen : m_visitedGen) {
            if (gen > genVBefore && gen <= genVAfter) {
                v++;
            }
        }
        m_stats.visitsTotal += v;
        m_stats.offspringCount += 1;
    }

    return registerOffspring(offspringId);
}

SignedNodeID NonDuplicationRecombiner::recombineMulti(const std::vector<std::pair<NodeID, BpPosition>>& segments) {
    m_audit.recombineCalls++;
    m_pendingBubbles.clear();

    std::chrono::steady_clock::time_point t;
    const uint64_t genVBefore = m_genVisited;
    if (m_instrument) {
        t = std::chrono::steady_clock::now();
    }
    syncToGrg();
    if (m_instrument) {
        m_stats.syncToGrgTime += elapsedSeconds(t);
        t = std::chrono::steady_clock::now();
    }
    const NodeID offspringId = m_grg->makeNode(1, /*negative=*/true);
    if (m_instrument) {
        m_stats.makeNodeTime += elapsedSeconds(t);
        m_stats.makeNodeCalls += 1;
    }
    m_audit.makeNodeCalls++;
    growNodeArrays(static_cast<NodeID>(m_grg->numNodes() - 1));
    m_genConnected++;

    double recurseT = 0.0;
    BpPosition start = 0;
    for (const auto& seg : segments) {
        const NodeID parentId = seg.first;
        const BpPosition end = seg.second;
        if (end > start) {
            if (m_debug) {
                std::fputs("BREAK\n", stderr);
            }
            std::chrono::steady_clock::time_point ts;
            if (m_instrument) {
                ts = std::chrono::steady_clock::now();
            }
            recurseAttach(parentId, offspringId, start, end);
            if (m_instrument) {
                recurseT += elapsedSeconds(ts);
                m_stats.recurseAttachCalls += 1;
                m_stats.segmentsProcessed += 1;
            }
        }
        start = end;
    }
    if (m_instrument) {
        m_stats.recurseAttachTime += recurseT;
    }

    applyPendingBubbles();
    clearModifiedCaches();

    if (m_instrument) {
        const uint64_t genVAfter = m_genVisited;
        uint64_t v = 0;
        for (uint64_t gen : m_visitedGen) {
            if (gen > genVBefore && gen <= genVAfter) {
                v++;
            }
        }
        m_stats.visitsTotal += v;
        m_stats.offspringCount += 1;
    }

    return registerOffspring(offspringId);
}

void NonDuplicationRecombiner::endGeneration() {
    std::chrono::steady_clock::time_point t;
    if (m_instrument) {
        t = std::chrono::steady_clock::now();
    }
    m_grg->sortMutations();
    if (m_instrument) {
        m_stats.sortMutationsTime += elapsedSeconds(t);
    }
    // sortMutations renumbers MutationIds (by (position, allele)), so any
    // cached (mut_id, position) pairs from this generation are now stale.
    // Position-derived caches (span, anc_cov, up_edges) survive intentionally.
    m_mutationCache.clear();
    m_posCache.clear();
}

void NonDuplicationRecombiner::clearPendingSampleRemovals() { m_pendingSampleRemovals.clear(); }

// ------------------------------------------------------------------
// Iterative recurseAttach (hot path) -- direct port of Python _recurse_attach.
//
// Every visit past both skip guards lands in exactly one of 13 audit buckets
// (9 standard + 4 root variants). The audit_check() invariants (see
// AuditCounters in the header) verify the implementation matches the
// algorithm spec by reconciling the per-case counts against connect /
// extract_bubble call totals.
// ------------------------------------------------------------------

void NonDuplicationRecombiner::recurseAttach(NodeID rootId,
                                             NodeID offspringId,
                                             BpPosition leftBound,
                                             BpPosition rightBound) {
    m_genVisited++;
    const uint64_t genV = m_genVisited;
    const uint64_t genC = m_genConnected;

    m_audit.recurseAttachCalls++;

    const SignedNodeID negOffspring = -static_cast<SignedNodeID>(offspringId);

    // Stack entries: (nodeId, L, R). Stored as separate vectors would be
    // marginally faster but a single vector-of-struct keeps the port direct.
    struct Frame {
        NodeID nodeId;
        BpPosition L;
        BpPosition R;
    };
    std::vector<Frame> stack;
    stack.push_back({rootId, leftBound, rightBound});

    while (!stack.empty()) {
        const Frame frame = stack.back();
        stack.pop_back();
        const NodeID nodeId = frame.nodeId;
        const BpPosition L = frame.L;
        const BpPosition R = frame.R;

        if (L >= R) {
            m_audit.skipEmptyInterval++;
            continue;
        }
        if (m_visitedGen[nodeId] == genV) {
            m_audit.skipAlreadyVisited++;
            continue;
        }
        m_visitedGen[nodeId] = genV;
        m_audit.visits++;

        // Cache-hit fast path: mutation/pos cache populated by buildAncestralCaches
        // for every input-graph node; lazy fallback handles bubbles/offspring.
        auto posIt = m_posCache.find(nodeId);
        if (posIt == m_posCache.end()) {
            getNodeMutations(nodeId);
            posIt = m_posCache.find(nodeId);
        }
        const std::vector<BpPosition>& positions = posIt->second;

        size_t left = 0;
        size_t right = 0;
        if (!positions.empty()) {
            left = static_cast<size_t>(std::lower_bound(positions.begin(), positions.end(), L) - positions.begin());
            right = static_cast<size_t>(std::lower_bound(positions.begin(), positions.end(), R) - positions.begin());
        }
        const size_t numRel = right - left;
        const size_t numAll = positions.size();

        const bool hasAllRelevant = numAll > 0 && numRel == numAll;
        const bool hasNoRelevant = numRel == 0;
        const bool hasPartialRelevant = numRel > 0 && numRel < numAll;

        // Resolve ancestral coverage (anc_cov_cache). Uninit -> lazy DFS.
        // Empty -> root (Python None). Set -> (Iu0, Iu1).
        IntervalOpt Iu = m_ancCovCache[nodeId];
        if (Iu.isUninit()) {
            Iu = getAncestralCoverage(nodeId);
        }

        if (Iu.isEmpty()) {
            // Root: only own muts matter. Four mutually-exclusive sub-cases.
            if (hasAllRelevant) {
                if (m_connectedGen[nodeId] != genC) {
                    m_grg->connect(static_cast<SignedNodeID>(nodeId), negOffspring);
                    m_connectedGen[nodeId] = genC;
                    m_pendingSampleRemovals.insert(nodeId);
                    m_audit.directAttachRoot++;
                    m_audit.connectCallsInAttach++;
                } else {
                    m_audit.directAttachDup++;
                }
            } else if (hasPartialRelevant) {
                const auto& mutsAtNode = m_mutationCache[nodeId];
                std::vector<MutationId> relMutIds;
                relMutIds.reserve(numRel);
                for (size_t i = left; i < right; i++) {
                    relMutIds.push_back(mutsAtNode[i].first);
                }
                const NodeID bubbleId = extractBubble(nodeId, relMutIds, offspringId);
                m_connectedGen[bubbleId] = genC;
                m_audit.bubbleStripPartialRt++;
            } else {
                m_audit.pruningRoot++;
            }
            continue;
        }

        // Non-root: classify by relevant-mutation status x interval overlap.
        const BpPosition Iu0 = Iu.lo;
        const BpPosition Iu1 = Iu.hi;
        const bool ancestralDisjoint = R <= Iu0 || L >= Iu1;
        const bool fullCoverage = Iu0 >= L && Iu1 <= R;

        if (hasAllRelevant && fullCoverage) {
            if (m_connectedGen[nodeId] != genC) {
                m_grg->connect(static_cast<SignedNodeID>(nodeId), negOffspring);
                m_connectedGen[nodeId] = genC;
                m_pendingSampleRemovals.insert(nodeId);
                m_audit.directAttach++;
                m_audit.connectCallsInAttach++;
            } else {
                m_audit.directAttachDup++;
            }
            continue;
        }

        if (hasNoRelevant && ancestralDisjoint) {
            m_audit.pruning++;
            continue;
        }

        if (hasNoRelevant) {
            const std::vector<NodeID>& parents = getUpEdgesCached(nodeId);
            // Early-attach optimization: multi-parent empty intermediary
            // fully covering [L, R) attaches directly instead of walking up.
            // Gated on num_all==0 + full_coverage + fanout>1; see Python
            // comments above the corresponding block for the correctness
            // argument (multitree preservation, etc.).
            if (numAll == 0 && fullCoverage && parents.size() > 1) {
                if (m_connectedGen[nodeId] != genC) {
                    m_grg->connect(static_cast<SignedNodeID>(nodeId), negOffspring);
                    m_connectedGen[nodeId] = genC;
                    m_pendingSampleRemovals.insert(nodeId);
                    m_audit.pathCompressionAttach++;
                    m_audit.connectCallsInAttach++;
                } else {
                    m_audit.directAttachDup++;
                }
                continue;
            }
            const BpPosition newL = (L > Iu0) ? L : Iu0;
            const BpPosition newR = (R < Iu1) ? R : Iu1;
            if (fullCoverage) {
                m_audit.pathCompression++;
            } else {
                m_audit.decomposition++;
            }
            if (newL >= newR) {
                m_audit.skipEmptyTrim++;
                continue;
            }
            // Skip parents already attached to this offspring (typically
            // bubbles from earlier segments). Iterate in reverse so the
            // first parent is processed first off the stack.
            for (auto it = parents.rbegin(); it != parents.rend(); ++it) {
                const NodeID parent = *it;
                if (m_connectedGen[parent] != genC) {
                    stack.push_back({parent, newL, newR});
                }
            }
            continue;
        }

        // Partial / all relevant, not full coverage -> bubble + maybe recurse.
        const auto& mutsAtNode = m_mutationCache[nodeId];
        std::vector<MutationId> relMutIds;
        relMutIds.reserve(numRel);
        for (size_t i = left; i < right; i++) {
            relMutIds.push_back(mutsAtNode[i].first);
        }
        const NodeID bubbleId = extractBubble(nodeId, relMutIds, offspringId);
        m_connectedGen[bubbleId] = genC;

        // Classify exactly one of the 5 non-root bubble cells.
        if (hasAllRelevant) {
            if (ancestralDisjoint) {
                m_audit.bubbleStrip++;
            } else {
                m_audit.bubbleSplit++;
            }
        } else {
            // hasPartialRelevant
            if (fullCoverage) {
                m_audit.bubbleFill++;
            } else if (ancestralDisjoint) {
                m_audit.bubbleStripPartial++;
            } else {
                m_audit.bubbleSplitPartial++;
            }
        }

        if (!ancestralDisjoint) {
            const BpPosition newL = (L > Iu0) ? L : Iu0;
            const BpPosition newR = (R < Iu1) ? R : Iu1;
            if (newL >= newR) {
                m_audit.skipEmptyTrim++;
                continue;
            }
            // extractBubble invalidated m_upEdgesCache[nodeId]; re-fetch.
            const std::vector<NodeID>& parents = getUpEdgesCached(nodeId);
            for (auto it = parents.rbegin(); it != parents.rend(); ++it) {
                const NodeID parent = *it;
                if (m_connectedGen[parent] != genC) {
                    stack.push_back({parent, newL, newR});
                }
            }
        }
    }
}

} // namespace grgl
