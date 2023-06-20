
"""
Implementation of the "Chaos Game" algorithm to "draw" an IFS attractor in a matrix,
given a rectangular region.
"""
function matrixattractor(ifs::AbstractVector{<:Function}, xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real};
  weights::AbstractVector{<:Real} = Float64[], seed = nothing,
  preiterations::Int = 100, iterations::Int = 10000, 
  rot90::Bool = false, value::Function = k::Int -> k)

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
  mtrx = rot90 ? fill(value(0), h, w) : fill(value(0), w, h)
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

  # Random orbit with "drawing point"
  for n in 1:iterations
    k = randindex()
    pk = ifs[k](pk)

    if isinsideclosed(pk, rr)
      mtrx[mtrxindx(pk, w, h, rr)...] = value(k)
    end
  end
  
  # Return matrix
  mtrx
end


"""
Draw an IFS attractor in an image.
"""
function imgattractor(ifs::AbstractVector{<:Function}, xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real};
  weights::AbstractVector{<:Real} = Float64[], seed = nothing,
  preiterations::Int = 100, iterations::Int = 10000,
  colormap::Union{Symbol, Vector{<:Colorant}} = :viridis)

  cm = typeof(colormap) == Symbol ? colorschemes[colormap] : ColorScheme(colormap)
  nfs = length(ifs)
  colors = [ cm[k/nfs] for k in 0:nfs ]

  matrixattractor(ifs, xs, ys, weights = weights, seed = seed,
    iterations = iterations, preiterations = preiterations,
    #value = k -> cm[k/numcolors], rot90 = true)
    value = k::Int -> colors[k+1], rot90 = true)
end


"""
    plotattractor(ifs, xs, ys; kwargs)

Plot an IFS attractor using the "Chaos Game" algorithm, extending the Makie's `heatmap`.

## Attributes
  weights
  seed
  iterations
  preiterations
"""
@recipe(PlotAttractor) do scene
  Attributes(
    weights = Float64[],
    seed = nothing,
    iterations = 10000,
    preiterations = 100
  )
end

function Makie.plot!(
  plt::PlotAttractor{<:Tuple{AbstractVector{<:Function},
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
   
  heatmap!(plt, obs_xs, obs_ys,
      matrixattractor(ifs, xs, ys, weights = ws, seed = sd,
        iterations = nits, preiterations = npits);
      plt.attributes.attributes...)
  
  plt
end


"""
    plotimgattractor(ifs, xs, ys; kwargs)

Plot an IFS attractor using the "Chaos Game" algorithm, extending the Makie's `image`.

## Attributes
  weights
  seed
  iterations
  preiterations
"""
@recipe(PlotImgAttractor) do scene
  Attributes(
    weights = Float64[],
    seed = nothing,
    iterations = 10000,
    preiterations = 100
  )
end

function Makie.plot!(
  plt::PlotImgAttractor{<:Tuple{AbstractVector{<:Function},
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
 
  nfs = length(ifs)  
  pltcm = plt.colormap[]
  cm = typeof(pltcm) == Symbol ? colorschemes[pltcm] : ColorScheme(pltcm)
  colors = [ cm[k/nfs] for k in 0:nfs ]

  image!(plt, obs_xs, obs_ys,
      matrixattractor(ifs, xs, ys, weights = ws, seed = sd,
        iterations = nits, preiterations = npits,
        value = k::Int -> colors[k+1]);
      plt.attributes.attributes...)
  
  plt
end


#
# Aquí va la implementación de plotattractor3D
#


#=
# SDDGraphics DEPRECATED!!!
function drawattractorC(ifs::AbstractIFSComplex, z0::Number=0;
  preiterations::Int=200, iterations::Int=1000)

  SDDGraphics.newdrawing()

  SDDGraphics.updatecolorarray(size(ifs))

  zn = z0
  lastused = 0

  for n in 1:preiterations
    zn, lastused = ifs(zn)
  end

  for n in 1:iterations
    zn, lastused = ifs(zn)
    if SDDGraphics.insiderectregion(zn)
      SDDGraphics.color(lastused)
      SDDGraphics.drawpoint(zn)
    end
  end

  SDDGraphics.drawing()
end

drawattractorC(funs::Vector{F}, z0::Number=0; preiterations::Int=200,
  iterations::Int=1000) where F <: Function =
drawattractorC(IFSComplex(funs),z0,preiterations=preiterations, iterations=iterations)

drawattractorC(funs::Vector{F}, z0::Number=0; preiterations::Int=200,
  iterations::Int=1000) where F <: AbstractAffineTransformation =
drawattractorC(IFSComplexAffine(funs),z0,preiterations=preiterations, iterations=iterations)
=#