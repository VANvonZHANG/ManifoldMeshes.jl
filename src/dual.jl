# Dual mesh framework
#
# Every mesh type that supports dual construction stores a RefValue{Nothing}
# or RefValue{AbstractManifoldMesh} field named _dual. The dual() function
# computes lazily on first call.

"""
    AbstractDualMesh{M}

Marker abstract type for meshes that are explicitly constructed as duals.
Concrete dual mesh types subtype both this and AbstractManifoldMesh{M}.
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
#   function dual(g::MyMesh{M}) where {M}
#       if g._dual[] === nothing
#           g._dual[] = _compute_dual(g)
#       end
#       return g._dual[]
#   end
#
#   has_dual(g::MyMesh) = g._dual[] !== nothing
