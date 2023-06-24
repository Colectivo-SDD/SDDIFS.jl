
"""
Drawings and calculations useful for research in the theory of
Iterated Functions Systems (IFSs) in the euclidean spaces and complex plane.
"""
module SDDIFS

using Reexport

@reexport using SDDCore, SDDGeometry, StaticArrays, Colors, ColorSchemes, Images, Makie
#, SDDGraphics # Deprecated!!!


#include("IFSs.jl") # Deprecated!!!

#export # Deprecated!!!
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
    plotimgattractor,
    plotimgattractor!,
    plotattractor,
    plotattractor!,
    plotscatterattractor,
    plotscatterattractor!,
    plotattractor3d, # Trabajo de Edgar
    plotattractor3d!, # Trabajo de Edgar
    plotscatterattractor3d, # Trabajo de Edgar
    plotscatterattractor3d! # Trabajo de Edgar


include("densitymaps.jl")

export
    imgdensitymap, # Trabajo de Linda
    plotdensitymap, # Trabajo de Linda
    plotdensitymap!, # Trabajo de Linda
    plotimgdensitymap, # Trabajo de Linda
    plotimgdensitymap!, # Trabajo de Linda
    imgstructuraldensitymap, # Trabajo de Linda
    plotimgstructuraldensitymap, # Trabajo de Linda
#    plotstructuraldensitymap, # Imposible con heatmap!!!
    imgstructuraldensitymapff, # Fractal Flame!
    plotimgstructuraldensitymapff, # Fractal Flame!
#    plotstructuraldensitymapff, # Imposible con heatmap!!!
    plotdensitymap3d, 
    plotdensitymap3d!,
    plotstructuraldensitymap3d, # ???
    plotstructuraldensitymap3d! # ???


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
