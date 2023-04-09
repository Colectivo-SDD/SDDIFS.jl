#
# Aquí va la implementación de los mapas de densidad 2D
#


"""
Plot the density map of a 3D IFS.
"""
@recipe(PlotDensityMap3D) do scene
    Attributes(
        weights = Float64[],
        seed = nothing,
        iterations = 10000,
        preiterations = 0,
        kind = :volume, # Makie's plot volume by default
        densityscale = log #log-density by default
    )
end

function Makie.plot!(
    plt::PlotDensityMap3D{<:Tuple{AbstractVector{<:Function},
    AbstractVector{<:Real}, AbstractVector{<:Real}, AbstractVector{<:Real}}})

    # Recipe attributes
    ifs = plt[1][] # IFS, array of functions
    obs_xs = plt[2] 
    xs = obs_xs[]
    obs_ys = plt[3]
    ys = obs_ys[]
    obs_zs = plt[4]
    zs = obs_zs[]

    # Plot keyword arguments
    ws = plt.weights[]
    # Cheking weights (probabilites) array
    if length(ws) > 0
        if length(ws) != length(ifs) && !checkindex(ws)
          @error "Non-suitable weights vector."
        end
    end

    sd = plt.seed[]
    nits = plt.iterations[]
    npits = plt.preiterations[]
    knd = plt.kind[]
    denscl = plt.densityscale[]
    colormap = plt.colormap[]

    ### Density map algorithm

    # Rectangular 3D region and dimensions
    w, h, d = length(xs), length(ys), length(zs)
    rr = RectRegion3D(xs, ys, zs)

    # Cube to contain final densities 
    #densitiescube = StaticArray{Tuple{w,h,d}}(fill(1f0, w, h, d)) # Aparently an array to large...
    densitiescube = fill(1f0, w, h, d)

    ## "Chaos game" algorithm, accumulating densities

    # Initial data for random orbit
    #pk = isnothing(sd) ? SVector(rand(xs), rand(ys), rand(zs)) : sd
    pk = isnothing(sd) ? [rand(xs), rand(ys), rand(zs)] : sd
    nfs = length(ifs) # Number of functions in the system
    randindex = createrandindex(nfs, ws)
    minval, maxval = 1f0, 1f0
    
    # Random orbit without "drawing point"
    for n in 1:npits
        k = randindex()
        pk = ifs[k](pk)
    end

    # Random orbit with "drawing point"
    for n in 1:nits
        k = randindex()
        pk = ifs[k](pk)

        if isinsideclosed(pk, rr)
            cindex = tocubeindex(pk,w,h,d,rr)
            densitiescube[cindex...] += 1f0
            if densitiescube[cindex...] > maxval
                maxval = densitiescube[cindex...]
            end
        end
    end

    ## end of "chaos game" algorithm

    # Scaling density (log is default)
    if denscl != nothing
        densitiescube = Float32.(denscl.(densitiescube))
        minval = Float32(denscl(minval))
        maxval = Float32(denscl(maxval))
    end

    # Normalizing data
    if minval > 0f0
        densitiescube = (densitiescube .- minval) ./ (maxval-minval)
    else
        densitiescube = densitiescube ./ maxval
    end
    
    ### End of density map algorithm

    # Force first color to totally transparent
    cm = to_colormap(colormap)
    cm[1] = RGBAf(red(cm[1]), green(cm[1]), blue(cm[1]), 0)

    if knd == :contour
        contour!(plt, obs_xs, obs_ys, obs_zs, densitiescube; plt.attributes.attributes..., colormap = cm)
    elseif knd == :volumeslices
        volumeslices!(plt, obs_xs, obs_ys, obs_zs, densitiescube; plt.attributes.attributes..., colormap = cm)
    else # knd == :volume
        volume!(plt, obs_xs, obs_ys, obs_zs, densitiescube; plt.attributes.attributes..., colormap = cm)
    end

    # Freeing memory
    densitiescube = 0
    GC.gc()

    plt
end


#=
"""
Plot the structural density map of a 3D IFS. Structural refers to the asignation of different
colors to each function in the system.
"""
@recipe(PlotStructuralDensityMap3D) do scene
Attributes(
    weights = Float64[],
    seed = nothing,
    iterations = 10000,
    preiterations = 0,
    kind = :volume, # Makie's plot volume by default
    densityscale = log #log-density by default
)
end

function Makie.plot!(
    plt::PlotStructuralDensityMap3D{<:Tuple{AbstractVector{<:Function},
    AbstractVector{<:Real}, AbstractVector{<:Real}, AbstractVector{<:Real}}})

    # Recipe attributes
    ifs = plt[1][] # IFS, array of functions
    obs_xs = plt[2] 
    xs = obs_xs[]
    obs_ys = plt[3]
    ys = obs_ys[]
    obs_zs = plt[4]
    zs = obs_zs[]

    # Plot keyword arguments
    ws = plt.weights[]
    sd = plt.seed[]
    nits = plt.iterations[]
    npits = plt.preiterations[]
    knd = plt.kind[]
    dens = plt.densityscale[]

    # Structural density-map algorithm

    # To do!!!

    densitycube = [ 0 for x in xs, y in ys, z in zs] # Cube to contain final
    #densitycube = dens.(densitycube)

    if knd == :contour
        contour!(plt, obs_xs, obs_ys, obs_zs, densitycube; plt.attributes.attributes...)
    elseif knd == :volumeslices
        volumeslices!(plt, obs_xs, obs_ys, obs_zs, densitycube; plt.attributes.attributes...)
    else # knd == :volume
        volume!(plt, obs_xs, obs_ys, obs_zs, densitycube; plt.attributes.attributes...)
    end

    plt
end

=#