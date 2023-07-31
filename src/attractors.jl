
#####
# 2D
#####

"""
Implementation of the "Chaos Game" algorithm to "draw" an IFS attractor in a matrix,
given a rectangular region.
"""
function matrixattractor(ifs::AbstractVector{<:Function}, xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real};
  weights::AbstractVector{<:Real} = Float64[], seed = nothing,
  preiterations::Int = 100, iterations::Int = 10000, 
  rot90::Bool = false, value::Function = k::Int -> k, fillvalue=0)

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
  fv = fillvalue #convert(typeof(value(1)), fillvalue)
  mtrx = rot90 ? fill(fv, h, w) : fill(fv, w, h)
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
    imgattractor(ifs, xs, ys[; kwargs])

Given an IFS

\$\\{f_k:X \\rightarrow X\\}_{k=1}^{K},\$

where \$X=\\mathbb{R}^2\$ or \$X=\\mathbb{C}\$, draw their attractor in an image,
using the an implementation of the "Chaos Game" algorithm with structural coloring.

#### Arguments
- `ifs::Vector{Function}`: An Iterated Functions System \$\\{f_k:X \\rightarrow X\\}\$.
- `xs::AbstractVector{Real}`: Base \$x\$ coordinates.
- `ys::AbstractVector{Real}`: Base \$y\$ coordinates.

#### Keyword Arguments
- `weights::Vector{Real} = []`: Probabilities associated to the IFS
- `seed = nothing`. Initial point for random orbit.
- `preiterations::Int = 100`: Iterations without drawing.
- `iterations::Int = 10000`: Iterations for drawing.
- `colormap::Union{Symbol, Vector{Colorant}} = :viridis`: Colormap.
- `bgcolor = RGB(1,1,1)`: Background color.
"""
function imgattractor(ifs::AbstractVector{<:Function},
  xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real};
  weights::AbstractVector{<:Real} = Float64[], seed = nothing,
  preiterations::Int = 100, iterations::Int = 10000,
  colormap::Union{Symbol, Vector{<:Colorant}} = :viridis,
  bgcolor::Colorant = RGBf(1,1,1))

  cm = typeof(colormap) == Symbol ? colorschemes[colormap] : ColorScheme(colormap)
  nfs = length(ifs)
  colors = [ cm[k/(nfs-1)] for k in 0:(nfs-1) ]

  matrixattractor(ifs, xs, ys, weights = weights, seed = seed,
    iterations = iterations, preiterations = preiterations,
    rot90 = true, value = k::Int -> colors[k], fillvalue = bgcolor)
end


"""
    plotimgattractor(ifs, xs, ys[; kwargs])

Given an IFS

\$\\{f_k:X \\rightarrow X\\}_{k=1}^{K},\$

where \$X=\\mathbb{R}^2\$ or \$X=\\mathbb{C}\$,
plot their attractor with the **Makie**'s `image`,
using the an implementation of the "Chaos Game" algorithm with structural coloring.

#### Arguments
- `ifs::Vector{Function}`: An Iterated Functions System \$\\{f_k:X \\rightarrow X\\}\$.
- `xs::AbstractVector{Real}`: Base \$x\$ coordinates.
- `ys::AbstractVector{Real}`: Base \$y\$ coordinates.

#### Keyword Arguments
- `weights::Vector{Real} = []`: Probabilities associated to the IFS
- `seed = nothing`. Initial point for random orbit.
- `preiterations::Int = 100`: Iterations without drawing.
- `iterations::Int = 10000`: Iterations for drawing.
- `colormap::Union{Symbol, Vector{Colorant}} = :viridis`: Colormap.
- `bgcolor = RGBA(1,1,1,0)`: Background color.    
"""
@recipe(PlotImgAttractor) do scene
  Attributes(
    weights = Float64[],
    seed = nothing,
    iterations = 10000,
    preiterations = 100,
    bgcolor = RGBA(1,1,1,0)
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
  bgc = plt.bgcolor[]
 
  nfs = length(ifs)  
  pltcm = plt.colormap[] # From Makie
  cm = typeof(pltcm) == Symbol ? colorschemes[pltcm] : ColorScheme(pltcm)
  colors = [ cm[k/(nfs-1)] for k in 0:(nfs-1) ]

  image!(plt, obs_xs, obs_ys,
      matrixattractor(ifs, xs, ys, weights = ws, seed = sd,
        iterations = nits, preiterations = npits,
        value = k::Int -> colors[k], fillvalue=bgc);
      plt.attributes.attributes...)
  
  plt
end


"""
    plotattractor(ifs, xs, ys[; kwargs])

Given an IFS

\$\\{f_k:X \\rightarrow X\\}_{k=1}^{K},\$

where \$X=\\mathbb{R}^2\$ or \$X=\\mathbb{C}\$,
plot their attractor with the **Makie**'s `heatmap`,
using the an implementation of the "Chaos Game" algorithm with structural coloring.

#### Arguments
- `ifs::Vector{Function}`: An Iterated Functions System \$\\{f_k:X \\rightarrow X\\}\$.
- `xs::AbstractVector{Real}`: Base \$x\$ coordinates.
- `ys::AbstractVector{Real}`: Base \$y\$ coordinates.

#### Keyword Arguments
- `weights::Vector{Real} = []`: Probabilities associated to the IFS
- `seed = nothing`. Initial point for random orbit.
- `preiterations::Int = 100`: Iterations without drawing.
- `iterations::Int = 10000`: Iterations for drawing.
- `colormap::Union{Symbol, Vector{Colorant}} = :viridis`: Colormap.

#### Notes
The background color is taken as from the first color of the colormap.
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
  #bgc = plt.bgcolor[]

  # Remove non Makie keyword arguments to avoid errors
  delete!(plt.attributes.attributes, :seed)
  delete!(plt.attributes.attributes, :weights)
  delete!(plt.attributes.attributes, :iterations)
  delete!(plt.attributes.attributes, :preiterations)
   
  heatmap!(plt, obs_xs, obs_ys,
      matrixattractor(ifs, xs, ys, weights = ws, seed = sd,
        iterations = nits, preiterations = npits);
      plt.attributes.attributes...)
  
  plt
end


"""
    plotscatterattractor(ifs[; kwargs])

Given an IFS

\$\\{f_k:X \\rightarrow X\\}_{k=1}^{K},\$

where \$X=\\mathbb{R}^2\$ or \$X=\\mathbb{C}\$,
plot their attractor with the **Makie**'s `scatter`,
using the an implementation of the "Chaos Game" algorithm with structural coloring.

#### Arguments
- `ifs::Vector{Function}`: An Iterated Functions System \$\\{f_k:X \\rightarrow X\\}\$.

#### Keyword Arguments
- `weights::Vector{Real} = []`: Probabilities associated to the IFS
- `seed = nothing`. Initial point for random orbit.
- `preiterations::Int = 100`: Iterations without drawing.
- `iterations::Int = 10000`: Iterations for drawing.
- `coloring::Symbol = :structural`: One of
  - `:structural`: Coloring by function index.
  - `:ordered`: Coloring by random orbit index.
  - `:random`: Random coloring.
"""
@recipe(PlotScatterAttractor) do scene
  Attributes(
    weights = Float64[],
    seed = nothing,
    iterations = 10000,
    preiterations = 100,
    coloring = :structural
  )
end

function Makie.plot!(plt::PlotScatterAttractor{<:Tuple{<:AbstractVector{<:Function}}})
  # Recipe attributes
  ifs = plt[1][] # IFS, array of functions

  topoint = createtopoint2D(functionkind2D(ifs[1]))

  # Plot keyword arguments
  sd = plt.seed[]
  nits = plt.iterations[]
  npits = plt.preiterations[]
  clr = plt.coloring[]
  ws = plt.weights[]

  # Cheking weights (probabilites) array
  if length(ws) > 0
      if length(ws) != length(ifs) && !checkindex(ws)
        @error "Non-suitable weights vector."
      end
  end

  ## "Chaos game" algorithm

  # Initial data for random orbit
  pk = isnothing(sd) ? topoint(rand(), rand()) : sd
  nfs = length(ifs) # Number of functions in the system
  randindex = createrandindex(nfs, ws)
    
  # Random orbit without "drawing point"
  for n in 1:npits
    k = randindex()
    pk = ifs[k](pk)
  end

  # Random orbit
  pts = typeof(pk)[]

  # Structural coloring
  colors = Float64[]

  # Random orbit with "drawing point"
  if clr == :structural
    for n in 1:nits
      k = randindex()
      pk = ifs[k](pk)
      push!(pts, pk)
      push!(colors, k)
    end
  else
    if clr == :ordered
      colors = 1:nits
    else
      colors = rand(nits)
    end
    
    for n in 1:nits
      k = randindex()
      pk = ifs[k](pk)
      push!(pts, pk)
      #push!(colors, n/nits)
    end
  end    
  ## end of "chaos game" algorithm

  # Remove non Makie keyword arguments to avoid errors
  delete!(plt.attributes.attributes, :seed)
  delete!(plt.attributes.attributes, :weights)
  delete!(plt.attributes.attributes, :iterations)
  delete!(plt.attributes.attributes, :preiterations)
  delete!(plt.attributes.attributes, :coloring)

  # Drawing the random orbit
  if typeof(pk) <: Number # Complex numbers points
    scatter!(plt, real.(pts), imag.(pts); plt.attributes.attributes..., color = colors)
  else # Bidimensinal cartesian plane points
    scatter!(plt, Point2f.(pts); plt.attributes.attributes..., color = colors)
  end
  
  # Freeing memory
  pts = 0
  colors = 0
  GC.gc()

  plt
end


#####
# 3D
#####

"""
    plotattractor3d(ifs,xs, ys, zs [; kwargs])

Given an IFS

\$\\{f_k:\\mathbb{R}^3 \\rightarrow \\mathbb{R}^3\\}_{k=1}^{K},\$

plot their attractor with the **Makie**'s volume,
using the an implementation of the "Chaos Game" algorithm with structural coloring.

#### Arguments
- `ifs::Vector{Function}`: An Iterated Functions System \$\\{f_k:\\mathbb{R}^3 \\rightarrow \\mathbb{R}^3\\}\$.
- `xs::AbstractVector{Real}`: Base \$x\$ coordinates.
- `ys::AbstractVector{Real}`: Base \$y\$ coordinates.
- `zs::AbstractVector{Real}`: Base \$z\$ coordinates.

#### Keyword Arguments
- `weights::Vector{Real} = []`: Probabilities associated to the IFS
- `seed = nothing`. Initial point for random orbit.
- `preiterations::Int = 100`: Iterations without drawing.
- `iterations::Int = 10000`: Iterations for drawing.
- `coloring::Symbol = :structural`: One of
  - `:structural`: Coloring by function index.
  - `:ordered`: Coloring by random orbit index.
  - `:random`: Random coloring.
- `kind::Symbol = :volume`: **Makie**'s plot kind
  - `:volume`: `volume` for *volume render*.
  - `:contour`: `contour` for multiple isosurfaces.
  - `:volumeslices`: `volumeslices` for volume slices.

#### Notes
Looks better with plot kind `:volume` and algorithm `:mip`.
"""
@recipe(PlotAttractor3D) do scene
    Attributes(
        weights = Float64[],
        seed = nothing,
        iterations = 10000,
        preiterations = 0,
        coloring = :structural,
        kind = :volume, # Makie's plot volume by default
        #scale = log #log-density-scale by default
    )
end

function Makie.plot!(
    plt::PlotAttractor3D{<:Tuple{AbstractVector{<:Function},
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
    clr = plt.coloring[]
    knd = plt.kind[]
    colormap = plt.colormap[]

    # Rectangular 3D region and dimensions
    w, h, d = length(xs), length(ys), length(zs)
    rr = RectRegion3D(xs, ys, zs)

    # Cube to contain values
    cube = fill(0f0, w, h, d)

    #
    nfs = length(ifs) # Number of functions in the system
    randindex = createrandindex(nfs, ws)

    # Coloring function selection
    funclrdict = Dict(
      :structural => (k::Int, n::Int) -> Float32(k),#,/nfs),
      :ordered => (k::Int, n::Int) -> Float32(n/nits),
      :random => (k::Int, n::Int) -> rand()
    )
    funclr = funclrdict[ clr == :structural ? :structural : ( clr == :orderd ? :ordered : :random ) ]

    ## "Chaos game" algorithm
    
    # Initial data for random orbit
    pk = isnothing(sd) ? [rand(xs), rand(ys), rand(zs)] : sd

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
            cube[cindex...] = funclr(k,n)
        end
    end

    ## end of "chaos game" algorithm

    # Force first color to totally transparent
    cm = to_colormap(colormap)
    cm[1] = RGBAf(red(cm[1]), green(cm[1]), blue(cm[1]), 0)

  # Remove non Makie keyword arguments to avoid errors
  delete!(plt.attributes.attributes, :seed)
  delete!(plt.attributes.attributes, :weights)
  delete!(plt.attributes.attributes, :iterations)
  delete!(plt.attributes.attributes, :preiterations)
  delete!(plt.attributes.attributes, :coloring)
  delete!(plt.attributes.attributes, :kind)

    if knd == :contour
        contour!(plt, obs_xs, obs_ys, obs_zs, cube; plt.attributes.attributes..., colormap = cm)
    elseif knd == :volumeslices
        volumeslices!(plt, obs_xs, obs_ys, obs_zs, cube; plt.attributes.attributes..., colormap = cm)
    else # knd == :volume
        volume!(plt, obs_xs, obs_ys, obs_zs, cube; plt.attributes.attributes..., colormap = cm)
    end

    # Freeing memory
    cube = 0
    GC.gc()

    plt
end


"""
    plotscatterattractor3d(ifs[; kwargs])

Given an IFS

\$\\{f_k:\\mathbb{R}^3 \\rightarrow \\mathbb{R}^3\\}_{k=1}^{K},\$

plot their attractor with the **Makie**'s `scatter`,
using the an implementation of the "Chaos Game" algorithm with structural coloring.

#### Arguments
- `ifs::Vector{Function}`: An Iterated Functions System \$\\{f_k:\\mathbb{R}^3 \\rightarrow \\mathbb{R}^3\\}\$.

#### Keyword Arguments
- `weights::Vector{Real} = []`: Probabilities associated to the IFS
- `seed = nothing`. Initial point for random orbit.
- `preiterations::Int = 100`: Iterations without drawing.
- `iterations::Int = 10000`: Iterations for drawing.
- `coloring::Symbol = :structural`: One of
  - `:structural`: Coloring by function index.
  - `:ordered`: Coloring by random orbit index.
  - `:random`: Random coloring.
- `usemesh::Bool = false`: Use `meshscatter` instead of `scatter`.

#### Notes
There are some strange bug with **(W)GLMakie** and `@recipe`... Works perfectly with **CairoMakie**.
"""
@recipe(PlotScatterAttractor3D) do scene
  Attributes(
    weights = Float64[],
    seed = nothing,
    iterations = 10000,
    preiterations = 100,
    coloring = :structural,
    usemesh = false
  )
end

function Makie.plot!(plt::PlotScatterAttractor3D{<:Tuple{AbstractVector{<:Function}}})
  # Recipe attributes
  ifs = plt[1][] # IFS, array of functions

  # Plot keyword arguments
  sd = plt.seed[]
  nits = plt.iterations[]
  npits = plt.preiterations[]
  clr = plt.coloring[]
  umsh = plt.usemesh[]
  ws = plt.weights[]

  # Cheking weights (probabilites) array
  if length(ws) > 0
      if length(ws) != length(ifs) && !checkindex(ws)
        @error "Non-suitable weights vector."
      end
  end

  ## "Chaos game" algorithm

  # Initial data for random orbit
  pk = isnothing(sd) ? [rand(), rand(), rand()] : sd
  nfs = length(ifs) # Number of functions in the system
  randindex = createrandindex(nfs, ws)
    
  # Random orbit without "drawing point"
  for n in 1:npits
    k = randindex()
    pk = ifs[k](pk)
  end

  # Separated coordinates for scatter
  pxs = Float64[]
  pys = Float64[]
  pzs = Float64[]

  # Coloring
  colors = Float64[]

  # Random orbit with "drawing point"
  if clr == :structural
    for n in 1:nits
      k = randindex()
      pk = ifs[k](pk)
      push!(pxs, pk[1])
      push!(pys, pk[2])
      push!(pzs, pk[3])
      push!(colors, k)
    end
  else
    colors = (clr == :ordered) ? (1:nits) : rand(nits)

    for n in 1:nits
      k = randindex()
      pk = ifs[k](pk)
      push!(pxs, pk[1])
      push!(pys, pk[2])
      push!(pzs, pk[3])
    end
  end

  ## end of "chaos game" algorithm

  # Remove non Makie keyword arguments to avoid errors
  delete!(plt.attributes.attributes, :seed)
  delete!(plt.attributes.attributes, :weights)
  delete!(plt.attributes.attributes, :iterations)
  delete!(plt.attributes.attributes, :preiterations)
  delete!(plt.attributes.attributes, :coloring)
  delete!(plt.attributes.attributes, :usemesh)

  # Drawing the random orbit
  if umsh == true
    meshscatter!(plt, pxs, pys, pzs; plt.attributes.attributes..., color = colors)
  else
    scatter!(plt, pxs, pys, pzs; plt.attributes.attributes..., color = colors)
  end

  # Freeing memory
  pxs = 0
  pys = 0
  pzs = 0
  colors = 0
  GC.gc()

  plt
end


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