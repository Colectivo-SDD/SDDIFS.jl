
"""
Drawings and calculations useful for research in the theory of
Iterated Functions Systems (IFSs) in the euclidean spaces and complex plane.
"""
module SDDIFS

using Reexport

@reexport using SDDCore, SDDGeometry, StaticArrays, Colors, ColorSchemes, Images, Makie
#, SDDGraphics


#include("IFSs.jl") # Deprecated!!!!

#export
    #AbstractIFS,
    #AbstractIFSComplex,
    #IFSComplex,
    #IFSComplexAffine,
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
    #contractions,
    #size


include("randindex.jl")

include("attractors.jl")

export
    imgattractor,
    #PlotAttractor,
    plotattractor,
    plotattractor!
    #plotattractor3D, # Trabajo de Edgar
    #plotattractor3D! # Trabajo de Edgar


include("densitymaps.jl")

export
    plotdensitymap3d,
    plotdensitymap3d!
#    imgdensitymap, # Trabajo de Linda
#    plotdensitymap, # Trabajo de Linda
#    imgstructuraldensitymap, # Trabajo de Linda
#    plotstructuraldensitymap, # Trabajo de Linda
#    plotstructuraldensitymap3D # Es posible???    


#include("distancemaps.jl")

#export
#    imgdistancemap,
#    plotdistancemap,
#    plotdistancemap3D # Trabajo de Edgar

#include("fractaldimensions.jl")

#export
#    boxdimension


# Extendend Makie plot functions
export
    plot,
    plot!


end # module
