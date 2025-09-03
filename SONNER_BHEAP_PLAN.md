# Updated Cache-Friendly B-Heap Implementation Plan

## Interface Compatibility Analysis
Based on NRM usage, our B-heap needs to support:
- `MutableBinaryMinHeap(data)` constructor
- `top_with_handle(pq)` - get min value and handle
- `update!(pq, handle, new_value)` - update indexed element
- `pq[handle]` - access element by handle (indexing)
- **New requirement**: Efficient bulk reset from `fill_rates_and_get_times!`

## B-Heap Design Decisions

### 1. Fanout Selection: 8-way (B=8)
- **Rationale**: With 64-byte cache lines and 8-byte Float64 elements, exactly 8 elements fit per cache line
- **Benefits**: Minimizes cache line fetches, reduces tree height from log₂(n) to log₈(n)
- **Trade-off**: Slightly more comparisons per level, but cache gains dominate

### 2. Memory Layout: Segmented Super-Nodes with Reset Optimization
```
SuperNode = [8 values | 8 handles | metadata]
Cache Line 1: [val1, val2, ..., val8]          # 64 bytes
Cache Line 2: [hdl1, hdl2, ..., hdl8]          # 64 bytes  
```

### 3. Handle Management: Persistent Indirection Table
- Maintain `handle_to_position::Vector{Tuple{Int,Int}}` mapping handles to (super_node, slot)
- **Reset optimization**: Preserve handle mappings, only update values
- Enables O(1) handle lookups while allowing efficient bulk resets

## Implementation Structure

```julia
mutable struct CacheFriendlyBHeap{T} <: AbstractMutableHeap{T}
    super_nodes::Vector{SuperNode{T}}
    handle_to_pos::Vector{Tuple{Int,Int}}  # (super_node_idx, slot) - persistent
    pos_to_handle::Vector{Vector{Int}}     # reverse mapping - persistent
    size::Int
    next_handle::Int
    
    # Reset optimization fields
    initial_capacity::Int                  # Track original size for reset
    reset_buffer::Vector{T}               # Workspace for bulk operations
end

struct SuperNode{T}
    values::MVector{8,T}                  # Mutable for efficient reset
    children_start::Int                   # Index of first child super-node
    is_leaf::Bool
end
```

## Efficient Reset Implementation

### 1. **Bulk Reset Strategy**: `reset_heap!(heap, new_data::Vector{T})`
```julia
function reset_heap!(heap::CacheFriendlyBHeap{T}, new_data::Vector{T}) where T
    @assert length(new_data) == heap.initial_capacity "Size mismatch in reset"
    
    # Option A: In-place heapify (preserve structure, update values)
    copy_to_super_nodes!(heap, new_data)  # O(n) copy
    heapify_in_place!(heap)                # O(n) bottom-up heapify
    
    # Option B: Bulk rebuild (faster for heavily changed data)  
    if should_rebuild(heap, new_data)      # Heuristic based on changes
        rebuild_from_scratch!(heap, new_data)
    else
        incremental_heapify!(heap, new_data)
    end
end
```

### 2. **Memory-Efficient Copy Operations**
```julia
function copy_to_super_nodes!(heap::CacheFriendlyBHeap{T}, data::Vector{T}) where T
    # Use SIMD-friendly chunked copy
    @inbounds for (i, super_node) in enumerate(heap.super_nodes)
        start_idx = (i-1) * 8 + 1
        end_idx = min(start_idx + 7, length(data))
        
        # Copy 8 elements at once (cache-line optimized)
        for j in 1:8
            if start_idx + j - 1 <= length(data)
                super_node.values[j] = data[start_idx + j - 1]
            else
                super_node.values[j] = typemax(T)  # Sentinel for unused slots
            end
        end
    end
end
```

### 3. **Preserve Handle Mappings During Reset**
- **Key insight**: NRM uses reaction indices as handles (1, 2, 3, ..., n)
- Handle-to-position mapping remains stable across resets
- Only need to update values, not restructure handle system

### 4. **Smart Heapify Strategy**
```julia
function heapify_in_place!(heap::CacheFriendlyBHeap{T}) where T
    # Bottom-up heapify optimized for super-nodes
    for level in reverse(1:height(heap))
        @inbounds for super_node_idx in level_range(heap, level)
            super_node = heap.super_nodes[super_node_idx]
            
            # Heapify within super-node (8-way)
            heapify_super_node!(super_node)
            
            # Bubble down to children if needed
            bubble_down_super_node!(heap, super_node_idx)
        end
    end
end
```

## Key Implementation Choices for Reset

### 1. **Memory Pool Strategy**
```julia
mutable struct BHeapMemoryPool{T}
    reset_workspace::Vector{T}        # Pre-allocated for data copying
    comparison_workspace::Vector{Int} # For sorting within super-nodes  
    temp_super_node::SuperNode{T}     # Temporary storage for swaps
end
```

### 2. **Reset Performance Modes**
```julia
@enum ResetMode begin
    PRESERVE_STRUCTURE  # Minimal allocation, preserve super-node structure
    REBUILD_OPTIMAL     # Full rebuild for optimal cache layout
    ADAPTIVE           # Choose based on data change analysis
end
```

### 3. **Cache-Aligned Reset Operations**
- Process super-nodes sequentially for optimal cache utilization
- Use prefetch hints for next super-node during current processing
- Batch comparisons within super-nodes using SIMD when possible

### 4. **Handle Consistency Across Resets**
```julia
function verify_handle_consistency!(heap::CacheFriendlyBHeap)
    # Debug mode: verify handles still map correctly after reset
    @assert length(heap.handle_to_pos) == heap.initial_capacity
    # Ensure no handle points to invalid position
    for (handle, (sn_idx, slot)) in enumerate(heap.handle_to_pos)
        @assert 1 <= sn_idx <= length(heap.super_nodes)
        @assert 1 <= slot <= 8
    end
end
```

## Enhanced Interface Implementation

### Core Methods with Reset Support:
1. **Constructor**: `CacheFriendlyBHeap(data)` - Store `initial_capacity` for reset
2. **Bulk Reset**: `reset!(heap, new_data)` - Efficient O(n) reset without reallocation
3. **Top Access**: `top_with_handle()` - Return root's minimum and its handle  
4. **Updates**: `update!(heap, handle, value)` - Locate via handle table, bubble up/down
5. **Indexing**: `getindex(heap, handle)` - O(1) lookup through indirection

### Reset-Specific Optimizations:
- **Memory reuse**: All internal structures sized for `initial_capacity`
- **Handle preservation**: Reaction indices remain consistent across resets
- **Vectorized operations**: Use Julia's SIMD capabilities for bulk copying
- **Cache prefetch**: Hint next super-node access during processing

## Performance Expectations

### Reset Performance:
- **Traditional approach**: O(n log n) for rebuild from scratch
- **B-heap reset**: O(n) for in-place heapify with preserved structure
- **Memory benefit**: Zero allocation during reset (vs full heap reconstruction)

### Cache Benefits During Reset:
- **Sequential access**: Process super-nodes in order for optimal cache usage
- **Bulk operations**: Copy/compare 8 elements at once per cache line access
- **Prefetch optimization**: Pipeline cache loading for next super-node

### Expected Speedups for Reset Operations:
- **Reset time**: 5-10x faster than full reconstruction
- **Memory pressure**: Eliminated allocation spikes during reset
- **Cache efficiency**: 3-4x better cache utilization vs element-by-element rebuild

## Implementation Phases

### Phase 1: Basic B-Heap with Reset (1.5 weeks)
- Core super-node structure and heap property maintenance
- Efficient bulk reset with handle preservation
- Memory pool and workspace management

### Phase 2: Interface Compatibility (3-4 days)  
- Implement full MutableBinaryMinHeap API including reset
- Handle indirection and indexing with reset consistency
- Integration testing with NRM reset patterns

### Phase 3: Reset Optimizations (1 week)
- SIMD-optimized bulk operations
- Adaptive reset strategies and cache prefetching
- Performance benchmarking vs existing heap reset

The B-heap will provide substantial performance gains for both regular operations and the critical reset functionality used by NRM during simulation initialization and reinitialization.