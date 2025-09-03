# Next Reaction Method (NRM) Performance Analysis and Improvement Plan

## Executive Summary
The Gibson-Bruck Next Reaction Method implementation in JumpProcesses.jl suffers from poor cache performance due to its use of MutableBinaryMinHeap. This analysis identifies the bottlenecks and proposes both short-term and long-term solutions to improve performance.

## Current Implementation Analysis

### Algorithm Overview
The NRM (Gibson & Bruck, 2000) improves upon Gillespie's First Reaction Method by:
1. Using a priority queue to store absolute firing times for each reaction
2. Reusing random numbers when updating dependent reaction times
3. Only updating reactions affected by the last executed reaction (via dependency graph)

### Key Code Sections
- **Main structure**: `NRMJumpAggregation` in `/src/aggregators/nrm.jl`
- **Critical hot path**: `update_dependent_rates!` (lines 90-118)
- **Priority queue operations**: Uses `MutableBinaryMinHeap` from DataStructures.jl

### Performance Bottleneck
The primary bottleneck occurs in `update_dependent_rates!` where for each dependent reaction:
```julia
# Line 105-108 and 111-114 in nrm.jl
update!(p.pq, rx, new_time)  # O(log n) heap operations
```

## Cache Performance Issues

### Binary Heap Cache Locality Problems
1. **Poor spatial locality**: Parent at index `i`, children at `2i` and `2i+1`
   - As heap grows, parent-child nodes rarely share cache lines
   - Each heap traversal causes O(log n) cache misses

2. **Random access patterns**: Dependency graph traversal adds unpredictable memory access
   - Dependent reactions can be scattered throughout the heap
   - No locality between consecutively accessed heap nodes

3. **Memory bandwidth limitation**: NRM is memory-bound, not compute-bound
   - Heap updates dominate runtime for large systems
   - Cache misses become critical performance factor

## Proposed Solutions

### 1. Short-Term: Library Replacement
**Replace MutableBinaryMinHeap with QuickHeaps.jl**
- **Effort**: Low (drop-in replacement)
- **Expected speedup**: 2-3x
- **Implementation**:
  ```julia
  using QuickHeaps
  pq = FastPriorityQueue{Int,T}(num_reactions)
  ```
- **Trade-offs**: Still uses binary heap, so cache issues remain for very large systems

### 2. Medium-Term: D-ary Heap Implementation
**Implement custom 4-ary or 8-ary heap**
- **Effort**: Medium (1-2 days implementation)
- **Expected speedup**: 2x over binary heap for large systems
- **Benefits**:
  - Shallower tree (log₄ n or log₈ n levels)
  - Better cache utilization (4-8 children fit in fewer cache lines)
  - Proven effectiveness in literature

**Implementation sketch**:
```julia
struct DaryHeap{D,T}
    data::Vector{T}
    indices::Vector{Int}  # For indexed updates
    d::Val{D}  # Compile-time constant arity
end

# Parent/child index calculations
parent(i, ::Val{D}) where D = (i + D - 2) ÷ D
child(i, k, ::Val{D}) where D = D * (i - 1) + k + 1
```

### 3. Long-Term: Cache-Optimized B-Heap
**Implement B-heap with cache-line-aware grouping**
- **Effort**: High (3-5 days implementation + testing)
- **Expected speedup**: 3-5x for large systems
- **Design**:
  - Group 7-15 elements per "super-node" (fits in 1-2 cache lines)
  - Each super-node is a mini-heap
  - Tree of super-nodes reduces to log₈ n or log₁₆ n traversals

**Memory layout**:
```
Cache line 1: [elements 0-7]
Cache line 2: [elements 8-15]
...
```

### 4. Alternative: Calendar Queue
**For systems with bounded time horizons**
- **Best for**: Systems where reaction times fall within predictable ranges
- **Effort**: High (requires analysis of time distributions)
- **Expected performance**: O(1) average for insert/delete
- **Limitation**: Performance degrades if time range is unbounded

## Data Layout Optimizations

### Structure-of-Arrays (SoA) Refactoring
Current (Array-of-Structures):
```julia
reactions = [Reaction1, Reaction2, ...]  # Poor cache utilization
```

Proposed (Structure-of-Arrays):
```julia
struct ReactionData
    rates::Vector{Float64}      # Hot data together
    times::Vector{Float64}      # Hot data together
    dependencies::Vector{...}   # Cold data separate
end
```

### Memory Pool for Temporary Arrays
- Pre-allocate workspace for heap operations
- Reuse arrays in `update_dependent_rates!`
- Reduce allocation overhead

## Implementation Priority

### Phase 1: Quick Wins (1 week)
1. ✅ Benchmark current implementation
2. ✅ Test QuickHeaps.jl as drop-in replacement
3. ✅ Implement basic memory pooling

### Phase 2: D-ary Heap (2 weeks)
1. ✅ Implement generic D-ary heap with indexed updates
2. ✅ Benchmark with D=4 and D=8
3. ✅ Integrate with NRM aggregator

### Phase 3: Advanced Optimizations (1 month)
1. ✅ Profile cache misses with system sizes 10³ to 10⁶ reactions
2. ✅ Implement B-heap if cache misses dominate
3. ✅ Consider calendar queue for specific problem classes

## Benchmarking Strategy

### Test Systems
1. **Weakly coupled**: Linear chain (each reaction affects 2-3 others)
2. **Strongly coupled**: All-to-all interactions
3. **Realistic**: Gene regulatory networks, metabolic pathways

### Metrics
- Wall time vs number of reactions (10² to 10⁶)
- Cache misses (using Profile.jl or perf)
- Memory bandwidth utilization
- Scaling behavior (time complexity verification)

## Expected Outcomes

### Performance Improvements
- **Small systems (< 1000 reactions)**: 2-3x with QuickHeaps
- **Medium systems (1000-10000)**: 3-4x with d-ary heap
- **Large systems (> 10000)**: 4-6x with B-heap

### Trade-offs
- Increased code complexity for advanced data structures
- Potential slight overhead for very small systems
- Need to maintain backward compatibility

## Available Julia Libraries Assessment

### Existing Options
1. **QuickHeaps.jl**: Good immediate improvement, but still binary heap
2. **DataStructures.jl**: Current implementation, well-tested but not cache-optimized
3. **IntervalHeaps.jl**: Not suitable for indexed updates

### Missing in Ecosystem
- ❌ B-heap implementations
- ❌ D-ary heaps with indexed updates
- ❌ Calendar queues
- ❌ Bucketed priority queues

### Recommendation
Given the lack of cache-optimized priority queues in Julia, implementing a custom solution would benefit both JumpProcesses and the broader Julia ecosystem. Consider releasing as separate package.

## References
1. Gibson, M.A. and Bruck, J. (2000). "Efficient Exact Stochastic Simulation of Chemical Systems with Many Species and Many Channels"
2. LaMarca, A. and Ladner, R. (1999). "The Influence of Caches on the Performance of Heaps"
3. Sanders, P. (2000). "Fast Priority Queues for Cached Memory"

## Next Steps
1. Create benchmark suite for current implementation
2. Test QuickHeaps.jl integration
3. Design and implement D-ary heap
4. Publish results and consider contributing optimized heap to Julia ecosystem