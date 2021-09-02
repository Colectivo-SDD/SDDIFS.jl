
"""
Drawings and calculations useful for research in the theory of
Iterated Functions Systems (IFSs) in the euclidean and complex planes.
"""
module SDDIFS

using Reexport

@reexport using SDDCore, SDDGeometry, SDDGraphics


include("IFSs.jl")

export
    AbstractIFS,
    AbstractIFSComplex,
    IFSComplex,
    IFSComplexAffine,
    #WIFSComplex,
    #WIFSComplexAffine,
    #AbstractIFSR2,
    #IFSR2,
    #IFSR2Affine,
    #WIFSR2,
    #WIFSR2Affine,
    #AbstractIFSR3,
    #IFSR3,
    #IFSR3Affine,
    #WIFSR3,
    #WIFSR3Affine,
    contractions,
    size


include("attractors.jl")

export
    drawattractorC
    #drawattractorR2
    #drawattractorR3


#include("densistymaps.jl")
#include("trappedpoints.jl")
#include("fractaldimensions.jl")

end # module
