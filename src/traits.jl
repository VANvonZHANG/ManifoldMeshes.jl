abstract type TopologyStyle end
struct IsGrid <: TopologyStyle end
struct IsMesh <: TopologyStyle end

TopologyStyle(::Type{IsGrid}) = IsGrid()
TopologyStyle(::Type{IsMesh}) = IsMesh()
TopologyStyle(::Type{T}) where {T} = IsMesh()
TopologyStyle(m::T) where {T} = TopologyStyle(T)

abstract type AbstractLocation end
struct NodeLoc <: AbstractLocation end
struct CellLoc <: AbstractLocation end
struct EdgeLoc <: AbstractLocation end
