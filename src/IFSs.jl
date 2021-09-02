
abstract type AbstractIFS <: Function end

abstract type AbstractIFSComplex <: AbstractIFS end


"""
An Iterated Function System formed by general Complex Contractions.
"""
struct IFSComplex{F <: Function} <: AbstractIFSComplex
  contractions::Vector{F}

  function IFSComplex{F}(funs::Vector{F}) where F <: Function
    new(deepcopy(funs))
  end
end

IFSComplex(funs::Vector{F}) where F <: Function = IFSComplex{F}(funs)


function (ifs::IFSComplex)(z::Number)
  n = rand(1:size(ifs))
  ifs.contractions[n](z), n
end


"""
An Iterated Function System formed by Complex Affine Contractions.
"""
struct IFSComplexAffine{F <: AbstractAffineTransformation} <: AbstractIFSComplex
  contractions::Vector{F}

  function IFSComplexAffine{F}(funs::Vector{F}) where F <: AbstractAffineTransformation
    new(deepcopy(funs))
  end
end

IFSComplexAffine(funs::Vector{F}) where F <: AbstractAffineTransformation =
IFSComplexAffine{F}(funs)


function (ifs::IFSComplexAffine)(z::Number)
  n = rand(1:size(ifs))
  ifs.contractions[n](z), n
end


contractions(ifs::AbstractIFS) = ifs.contractions
size(ifs::AbstractIFS) = length(ifs.contractions)
