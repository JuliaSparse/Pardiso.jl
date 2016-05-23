# https://github.com/JuliaGPU/VulkanCore.jl/blob/master/src/CEnum.jl
# The VulkanCore.jl package is licensed under the MIT "Expat" License:
# Copyright (c) 2016: Valentin Churavy.

module CEnum

using Compat

abstract Cenum{T}
@compat Base.:|{T<:Cenum}(a::T, b::T) = T(Int(a) | Int(b))
@compat Base.:&{T<:Cenum}(a::T, b::T) = T(Int(a) & Int(b))
# typemin and typemax won't change for an enum, so we might as well inline them per type
function Base.typemax{T<:Cenum}(::Type{T})
    last(enum_values(T))
end
function Base.typemin{T<:Cenum}(::Type{T})
    first(enum_values(T))
end
Base.convert{T<:Integer}(::Type{T}, x::Cenum) = convert(T, Base.box(Int32, x))
Base.write(io::IO, x::Cenum) = write(io, Int32(x))
Base.read{T<:Cenum}(io::IO, ::Type{T}) = T(read(io, Int32))

enum_values{T<:Cenum}(::T) = enum_values(T)
enum_names{T<:Cenum}(::T) = enum_names(T)

function is_member{T<:Cenum}(::Type{T}, x::Integer)
    is_member(T, enum_values(T), x)
end
@inline is_member{T<:Cenum}(::Type{T}, r::UnitRange, x::Integer) = x in r
@inline function is_member{T<:Cenum}(::Type{T}, values::Tuple, x::Integer)
    lo, hi = typemin(T), typemax(T)
    x<lo || x>hi && return false
    for val in values
        val == x && return true
        val > x && return false # is sorted
    end
    return false
end

function enum_name{T<:Cenum}(x::T)
    index = findfirst(enum_values(T), Int(x))
    if index != 0
        return enum_names(T)[index]
    end
    error("Invalid enum: $x, name not found")
end
function Base.show(io::IO, x::Cenum)
    print(io, enum_name(x), "($(Int(x)))")
end

function islinear(array)
    isempty(array) && return false # false, really? it's kinda undefined?
    lastval = first(array)
    for val in rest(array, 2)
        val-lastval == 1 || return false
    end
    return true
end


macro cenum(name, args...)
    if !isa(name, Symbol)
        error("Name must be symbol or Name{Type}. Found: $name")
    end
    lastval = -1
    name_values = map([args...]) do arg
        if isa(arg, Symbol)
            lastval += 1
            val = lastval
            sym = arg
        elseif arg.head == :(=) || arg.head == :kw
            sym,val = arg.args
        else
            error("Expression of type $arg not supported. Try only symbol or name = value")
        end
        (sym, val)
    end
    sort!(name_values, by=last) # sort for values
    values = map(last, name_values)

    if islinear(values) # optimize for linear values
        values = :($(first(values)):$(last(values)))
    else
        values = :(tuple($(values...)))
    end
    value_block = Expr(:block)
    typename = esc(name)
    for (ename, value) in name_values
        push!(value_block.args, :(const $(esc(ename)) = $typename($value)))
    end

    expr = quote
        bitstype 32 $typename <: CEnum.Cenum{UInt32}
        function Base.convert(::Type{$typename}, x::Integer)
            is_member($typename, x) || Base.Enums.enum_argument_error($(Expr(:quote, name)), x)
            Base.box($typename, convert(Int32, x))
        end
        CEnum.enum_names(::Type{$typename}) = tuple($(map(x-> Expr(:quote, first(x)), name_values)...))
        CEnum.enum_values(::Type{$typename}) = $values
        $value_block
    end
    expr
end
export @cenum

end # module