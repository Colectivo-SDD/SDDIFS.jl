#
# 2D IFS Attractors with Chas Game, Density Map and without Structural Coloring
#

"""
Implementation of the "Chaos Game" algorithm to "draw" an IFS attractor in a matrix,
given a rectangular region, using the density map technique.
"""
function matrixdensitymap(ifs::AbstractVector{<:Function},
  xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real};
  weights::AbstractVector{<:Real} = Float64[], seed = nothing,
  preiterations::Int = 100, iterations::Int = 10000, value = nothing, #::Function
  rot90::Bool = false, scale::Function = log)

  # Cheking weights (probabilites) array
  if length(weights) > 0
    if length(weights) != length(ifs) && !checkindex(weights)
      @error "Non-suitable weights vector."
    end
  end

  # Checking for a suitable 2D function and creating a "topoint" function
  topoint = createtopoint2D(functionkind2D(ifs[1]))

  # Matrix and rectangular region
  w, h = length(xs), length(ys)
  rr = RectRegion(xs, ys)
  #mtrx = rot90 ? fill(value(1), h, w) : fill(value(1), w, h)
  mtrx = rot90 ? fill(1.0, h, w) : fill(1.0, w, h)
  mtrxindx = rot90 ? torot90matrixindex : tomatrixindex

  # Initial data for random orbit
  pk = isnothing(seed) ? topoint(rand(xs), rand(ys)) : seed
  nfs = length(ifs) # Number of functions in the system
  randindex = createrandindex(nfs, weights)

  # The "Chaos Game" Algorithm

  # Random orbit without "drawing point"
  for n in 1:preiterations
    k = randindex()
    pk = ifs[k](pk)
  end

  nmax = 1

  # Random orbit with "drawing point"
  for n in 1:iterations
    k = randindex()
    pk = ifs[k](pk)

    if isinsideclosed(pk, rr)
      indx = mtrxindx(pk, w, h, rr)
      mtrx[indx...] += 1
      if nmax < mtrx[indx...]
        nmax = mtrx[indx...]
      end
    end
  end
  
  # Return matrix
  if isnothing(value)
    return scale.(mtrx)./scale(nmax)
  end

  value.(scale.(mtrx)./scale(nmax))
end


"""
Draw an IFS attractor using the chaos game and density map in an image.
"""
function imgdensitymap(ifs::AbstractVector{<:Function}, xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real};
  weights::AbstractVector{<:Real} = Float64[], seed = nothing,
  preiterations::Int = 100, iterations::Int = 10000, scale::Function = log,
  colormap::Union{Symbol, Vector{<:Colorant}} = :viridis)

  cm = typeof(colormap) == Symbol ? colorschemes[colormap] : ColorScheme(colormap)
  #nfs = length(ifs)
  #colors = [ cm[k/nfs] for k in 0:nfs ]

  matrixdensitymap(ifs, xs, ys, weights = weights, seed = seed,
    iterations = iterations, preiterations = preiterations, scale = scale,
    rot90 = true, value = t::Real -> cm[t])
end


"""
    plotdensitymap(ifs, xs, ys; kwargs)

Plot an IFS attractor using the "Chaos Game" algorithm with the density map technique,
extending the Makie's `heatmap`.

## Attributes
  weights
  seed
  iterations
  preiterations
  scale
"""
@recipe(PlotDensityMap) do scene
  Attributes(
    weights = Float64[],
    seed = nothing,
    iterations = 10000,
    preiterations = 100,
    scale = log
  )
end

function Makie.plot!(
  plt::PlotDensityMap{<:Tuple{AbstractVector{<:Function},
  AbstractVector{<:Real}, AbstractVector{<:Real}}})

  # Recipe attributes
  ifs = plt[1][] # IFS, array of functions
  obs_xs = plt[2] 
  xs = obs_xs[]
  obs_ys = plt[3]
  ys = obs_ys[]

  # Plot keyword arguments
  ws = plt.weights[]
  sd = plt.seed[]
  nits = plt.iterations[]
  npits = plt.preiterations[]
  scl = plt.scale[]
   
  heatmap!(plt, obs_xs, obs_ys,
      matrixdensitymap(ifs, xs, ys, weights = ws, seed = sd,
        iterations = nits, preiterations = npits, scale=scl);
      plt.attributes.attributes...)
  
  plt
end


"""
    plotimgdensitymap(ifs, xs, ys; kwargs)

Plot an IFS attractor using the "Chaos Game" algorithm with the density map technique,
extending the Makie's `image`.

## Attributes
  weights
  seed
  iterations
  preiterations
"""
@recipe(PlotImgDensityMap) do scene
  Attributes(
    weights = Float64[],
    seed = nothing,
    iterations = 10000,
    preiterations = 100,
    scale = log
  )
end

function Makie.plot!(
  plt::PlotImgDensityMap{<:Tuple{AbstractVector{<:Function},
  AbstractVector{<:Real}, AbstractVector{<:Real}}})

  # Recipe attributes
  ifs = plt[1][] # IFS, array of functions
  obs_xs = plt[2] 
  xs = obs_xs[]
  obs_ys = plt[3]
  ys = obs_ys[]

  # Plot keyword arguments
  ws = plt.weights[]
  sd = plt.seed[]
  nits = plt.iterations[]
  npits = plt.preiterations[]
  scl = plt.scale[]
 
  nfs = length(ifs)  
  pltcm = plt.colormap[]
  cm = typeof(pltcm) == Symbol ? colorschemes[pltcm] : ColorScheme(pltcm)
  #colors = [ cm[k/nfs] for k in 0:nfs ]

  image!(plt, obs_xs, obs_ys,
      matrixdensitymap(ifs, xs, ys, weights = ws, seed = sd,
        iterations = nits, preiterations = npits, scale = scl,
        value = t::Real -> cm[t]);
      plt.attributes.attributes...)
  
  plt
end


#
# 2D IFS Attractors with Chas Game, Density Map and Structural Coloring
#

"""
Implementation of the "Chaos Game" algorithm to "draw" an IFS attractor in a matrix,
given a rectangular region, using the density map technique and structural coloring.
"""
function matrixstructuraldensitymap(ifs::AbstractVector{<:Function},
  xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real};
  weights::AbstractVector{<:Real} = Float64[], seed = nothing,
  preiterations::Int = 100, iterations::Int = 10000, #value = nothing, #::Function
  rot90::Bool = false, scale::Function = log)

  # Cheking weights (probabilites) array
  if length(weights) > 0
    if length(weights) != length(ifs) && !checkindex(weights)
      @error "Non-suitable weights vector."
    end
  end

  # Checking for a suitable 2D function and creating a "topoint" function
  topoint = createtopoint2D(functionkind2D(ifs[1]))

  # Matrix and rectangular region
  w, h = length(xs), length(ys)
  rr = RectRegion(xs, ys)
  nfs = length(ifs) # Number of functions in the system
  mtrxs = AbstractArray{<:Real}[]
  for k in 1:nfs
    push!(mtrxs, rot90 ? fill(1.0, h, w) : fill(1.0, w, h))
  end
  mtrxindx = rot90 ? torot90matrixindex : tomatrixindex

  # Initial data for random orbit
  pk = isnothing(seed) ? topoint(rand(xs), rand(ys)) : seed
  randindex = createrandindex(nfs, weights)

  # The "Chaos Game" Algorithm

  # Random orbit without "drawing point"
  for n in 1:preiterations
    k = randindex()
    pk = ifs[k](pk)
  end

  #nmaxs = fill(1, 1, nfs)
  nmax = 0

  # Random orbit with "drawing point"
  for n in 1:iterations
    k = randindex()
    pk = ifs[k](pk)

    if isinsideclosed(pk, rr)
      indx = mtrxindx(pk, w, h, rr)
      mtrxs[k][indx...] += 1
      if nmax < mtrxs[k][indx...]
        nmax = mtrxs[k][indx...]
      end
    end
  end
  
  for k in 1:length(ifs)
    mtrxs[k] = scale.(mtrxs[k])./(length(ifs)*scale(nmax))
  end

  mtrxs
end


"""
Draw an IFS attractor using the chaos game and density map in an image.
"""
function imgstructuraldensitymap(ifs::AbstractVector{<:Function}, xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real};
  weights::AbstractVector{<:Real} = Float64[], seed = nothing,
  preiterations::Int = 100, iterations::Int = 10000, scale::Function = log,
  colormap::Union{Symbol, Vector{<:Colorant}} = :viridis)

  cm = typeof(colormap) == Symbol ? colorschemes[colormap] : ColorScheme(colormap)
  nfs = length(ifs)
  colors = [ cm[k/nfs] for k in 0:nfs ]

  mtrxs = matrixstructuraldensitymap(ifs, xs, ys, weights = weights, seed = seed,
    iterations = iterations, preiterations = preiterations, scale = scale,
    rot90 = true) #, value = val)

  H, W = length(ys), length(xs)
  img = fill(RGB(0,0,0), H, W)

  for h in 1:H
    for w in 1:W        
        s = 0.0
        for k in 1:nfs 
            s += mtrxs[k][h,w]
        end

        if s < 0.00001
            img[h,w] = colors[1]
        else
            for k in 1:nfs 
                img[h,w] += mtrxs[k][h,w]*colors[k+1]
            end
        end
    end
  end

  mtrxs = 0
  GC.gc()

  img
end

#=

# ¿Cómo hacer mapas de densidad con coloreado estructural y mapa de densidad?

"""
    plotstructuraldensitymap(ifs, xs, ys; kwargs)

Plot an IFS attractor using the "Chaos Game" algorithm with the density map technique,
extending the Makie's `heatmap`.

## Attributes
  weights
  seed
  iterations
  preiterations
  scale
"""
@recipe(PlotStructuralDensityMap) do scene
  Attributes(
    weights = Float64[],
    seed = nothing,
    iterations = 10000,
    preiterations = 100,
    scale = log
  )
end

function Makie.plot!(
  plt::PlotStructuralDensityMap{<:Tuple{AbstractVector{<:Function},
  AbstractVector{<:Real}, AbstractVector{<:Real}}})

  # Recipe attributes
  ifs = plt[1][] # IFS, array of functions
  obs_xs = plt[2] 
  xs = obs_xs[]
  obs_ys = plt[3]
  ys = obs_ys[]

  # Plot keyword arguments
  ws = plt.weights[]
  sd = plt.seed[]
  nits = plt.iterations[]
  npits = plt.preiterations[]
  scl = plt.scale[]
   
  heatmap!(plt, obs_xs, obs_ys,
      matrixstructuraldensitymap(ifs, xs, ys, weights = ws, seed = sd,
        iterations = nits, preiterations = npits, scale=scl);
      plt.attributes.attributes...)
  
  plt
end


"""
    plotimgstructuraldensitymap(ifs, xs, ys; kwargs)

Plot an IFS attractor using the "Chaos Game" algorithm with the density map technique,
extending the Makie's `image`.

## Attributes
  weights
  seed
  iterations
  preiterations
"""
@recipe(PlotImgStructuralDensityMap) do scene
  Attributes(
    weights = Float64[],
    seed = nothing,
    iterations = 10000,
    preiterations = 100,
    scale = log
  )
end

function Makie.plot!(
  plt::PlotImgStructuralDensityMap{<:Tuple{AbstractVector{<:Function},
  AbstractVector{<:Real}, AbstractVector{<:Real}}})

  # Recipe attributes
  ifs = plt[1][] # IFS, array of functions
  obs_xs = plt[2] 
  xs = obs_xs[]
  obs_ys = plt[3]
  ys = obs_ys[]

  # Plot keyword arguments
  ws = plt.weights[]
  sd = plt.seed[]
  nits = plt.iterations[]
  npits = plt.preiterations[]
  scl = plt.scale[]
 
  nfs = length(ifs)  
  pltcm = plt.colormap[]
  cm = typeof(pltcm) == Symbol ? colorschemes[pltcm] : ColorScheme(pltcm)
  #colors = [ cm[k/nfs] for k in 0:nfs ]

  image!(plt, obs_xs, obs_ys,
      matrixdensitymap(ifs, xs, ys, weights = ws, seed = sd,
        iterations = nits, preiterations = npits, scale = scl,
        value = t::Real -> cm[t]);
      plt.attributes.attributes...)
  
  plt
end

=#


#
# 3D IFS Attractors with chaos Game and  Density Map
#

"""
Plot the density map of a 3D IFS.
"""
@recipe(PlotImgStructuralDensityMap3D) do scene
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
    plt::PlotImgStructuralDensityMap3D{<:Tuple{AbstractVector{<:Function},
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