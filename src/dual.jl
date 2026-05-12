# Dual mesh framework
#
# Every mesh type that supports dual construction stores a
# RefValue{Union{Nothing, AbstractManifoldMesh{M}}} field named _dual. The dual() function
# computes lazily on first call.

"""
    AbstractDualMesh{M} <: AbstractManifoldMesh{M}

Marker abstract type for meshes that are explicitly constructed as duals.
Concrete dual mesh types subtype `AbstractDualMesh{M}`, which itself extends
`AbstractManifoldMesh{M}`.
"""
abstract type AbstractDualMesh{M} <: AbstractManifoldMesh{M} end

# -- Default no-dual implementations --

has_dual(g::AbstractManifoldMesh) = false

# Mesh types that support duals must override these.
# The pattern is:
#
#   struct MyMesh{M} <: AbstractManifoldMesh{M}
#       # ... geometry fields ...
#       _dual::Base.RefValue{Union{Nothing, AbstractManifoldMesh{M}}}
#   end
#
#   # In constructor: M is the manifold instance, so use typeof(M) for the type parameter
#   _dual = Ref{Union{Nothing, AbstractManifoldMesh{typeof(M)}}}(nothing)
#
#   function dual(g::MyMesh{M}) where {M}
#       if g._dual[] === nothing
#           g._dual[] = _compute_dual(g)
#       end
#       return g._dual[]
#   end
#
#   has_dual(g::MyMesh) = g._dual[] !== nothing
