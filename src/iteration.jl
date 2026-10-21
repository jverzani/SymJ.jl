# Expressions and polynomials behave as 0-dimensional, single-element collections, the
# way `Number` does in Base, so generic code that iterates or indexes a scalar works.

const _Scalar = Union{Sym,Poly}

Base.iterate(x::_Scalar) = (x, nothing)
Base.iterate(x::_Scalar, ::Any) = nothing
Base.isempty(x::_Scalar) = false
Base.length(x::_Scalar) = 1
Base.size(x::_Scalar) = ()
Base.size(x::_Scalar, d::Integer) = d < 1 ? throw(BoundsError()) : 1
Base.axes(x::_Scalar) = ()
Base.axes(x::_Scalar, d::Integer) = d < 1 ? throw(BoundsError()) : Base.OneTo(1)
Base.ndims(x::_Scalar) = 0
Base.ndims(::Type{<:_Scalar}) = 0
Base.eltype(::Type{T}) where {T<:_Scalar} = T
Base.IteratorSize(::Type{<:_Scalar}) = Base.HasShape{0}()
Base.keys(::_Scalar) = Base.OneTo(1)
Base.firstindex(x::_Scalar) = 1
Base.firstindex(x::_Scalar, d::Int) = d < 1 ? throw(BoundsError()) : 1
Base.lastindex(x::_Scalar) = 1
Base.lastindex(x::_Scalar, d::Int) = d < 1 ? throw(BoundsError()) : 1
Base.first(x::_Scalar) = x
Base.last(x::_Scalar) = x
Base.copy(x::_Scalar) = x

Base.getindex(x::_Scalar) = x
function Base.getindex(x::_Scalar, i::Integer)
    @inline
    @boundscheck i == 1 || throw(BoundsError(x, i))
    x
end
function Base.getindex(x::_Scalar, I::Integer...)
    @inline
    @boundscheck all(isone, I) || throw(BoundsError(x, I))
    x
end
