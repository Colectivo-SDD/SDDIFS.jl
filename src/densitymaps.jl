
#####
# 2D
#####

#
# 2D IFS Attractors with Chaos Game, Density Map and without Structural Coloring
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
  histogram = rot90 ? fill(1, h, w) : fill(1, w, h)
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
      histogram[indx...] += 1
      if nmax < histogram[indx...]
        nmax = histogram[indx...]
      end
    end
  end

  val0 = isnothing(value) ? 0.0 : value(0.0)
  mtrx = rot90 ? fill(val0, h, w) : fill(val0, w, h)
  snmax = scale(nmax)

  if isnothing(value)
    for j in 1:h
      for i in 1:w
        if histogram[j,i] > 1
          mtrx[j,i] = scale(histogram[j,i])/snmax
        end
      end
    end
  else
    for j in 1:h
      for i in 1:w
        if densmtrx[j,i] > 1
          mtrx[j,i] = value(scale(histogram[j,i])/snmax)
        end
      end
    end
  end

  histogram = 0
  GC.gc()

  # Return values matrix
  mtrx
end


"""
    imgdensitymap(ifs, xs, ys[; kwargs])

Given an IFS

\$\\{f_k:X \\rightarrow X\\}_{k=1}^{K},\$

where \$X=\\mathbb{R}^2\$ or \$X=\\mathbb{C}\$, draw their attractor in an image,
using the an implementation of the "Chaos Game" algorithm with the density map technique.

#### Arguments
- `ifs::Vector{Function}`: An Iterated Functions System \$\\{f_k:X \\rightarrow X\\}\$.
- `xs::AbstractVector{Real}`: Base \$x\$ coordinates.
- `ys::AbstractVector{Real}`: Base \$y\$ coordinates.

#### Keyword Arguments
- `weights::Vector{Real} = []`: Probabilities associated to the IFS
- `seed = nothing`. Initial point for random orbit.
- `preiterations::Int = 100`: Iterations without drawing.
- `iterations::Int = 10000`: Iterations for drawing.
- `scale::Function = log`: Density scaling function \$(1,\\infty)\\rightarrow (0,\\infty)\$.
- `colormap::Union{Symbol, Vector{Colorant}} = :viridis`: Colormap.
"""
function imgdensitymap(ifs::AbstractVector{<:Function}, xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real};
  weights::AbstractVector{<:Real} = Float64[], seed = nothing,
  preiterations::Int = 100, iterations::Int = 10000, scale::Function = log,
  colormap::Union{Symbol, Vector{<:Colorant}} = :viridis)

  cm = typeof(colormap) == Symbol ? colorschemes[colormap] : ColorScheme(colormap)

  matrixdensitymap(ifs, xs, ys, weights = weights, seed = seed,
    iterations = iterations, preiterations = preiterations, scale = scale,
    rot90 = true, value = t::Real -> cm[t])
end


"""
    plotimgdensitymap(ifs, xs, ys[; kwargs])

Given an IFS

\$\\{f_k:X \\rightarrow X\\}_{k=1}^{K},\$

where \$X=\\mathbb{R}^2\$ or \$X=\\mathbb{C}\$, draw their attractor with **Makie**'s `image`,
using the an implementation of the "Chaos Game" algorithm with the density map technique.

#### Arguments
- `ifs::Vector{Function}`: An Iterated Functions System \$\\{f_k:X \\rightarrow X\\}\$.
- `xs::AbstractVector{Real}`: Base \$x\$ coordinates.
- `ys::AbstractVector{Real}`: Base \$y\$ coordinates.

#### Keyword Arguments
- `weights::Vector{Real} = []`: Probabilities associated to the IFS
- `seed = nothing`. Initial point for random orbit.
- `preiterations::Int = 100`: Iterations without drawing.
- `iterations::Int = 10000`: Iterations for drawing.
- `scale::Function = log`: Density scaling function \$(1,\\infty)\\rightarrow (0,\\infty)\$.
- `colormap::Union{Symbol, Vector{Colorant}} = :viridis`: Colormap.
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

  image!(plt, obs_xs, obs_ys,
      matrixdensitymap(ifs, xs, ys, weights = ws, seed = sd,
        iterations = nits, preiterations = npits, scale = scl,
        value = t::Real -> cm[t]);
      plt.attributes.attributes...)
  
  plt
end


"""
    plotdensitymap(ifs, xs, ys[; kwargs])

Given an IFS

\$\\{f_k:X \\rightarrow X\\}_{k=1}^{K},\$

where \$X=\\mathbb{R}^2\$ or \$X=\\mathbb{C}\$, draw their attractor with **Makie**'s `heatmap`,
using the an implementation of the "Chaos Game" algorithm with the density map technique.

#### Arguments
- `ifs::Vector{Function}`: An Iterated Functions System \$\\{f_k:X \\rightarrow X\\}\$.
- `xs::AbstractVector{Real}`: Base \$x\$ coordinates.
- `ys::AbstractVector{Real}`: Base \$y\$ coordinates.

#### Keyword Arguments
- `weights::Vector{Real} = []`: Probabilities associated to the IFS
- `seed = nothing`. Initial point for random orbit.
- `preiterations::Int = 100`: Iterations without drawing.
- `iterations::Int = 10000`: Iterations for drawing.
- `scale::Function = log`: Density scaling function \$(1,\\infty)\\rightarrow (0,\\infty)\$.
- `colormap::Union{Symbol, Vector{Colorant}} = :viridis`: Colormap.
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


#
# 2D IFS Attractors with Chaos Game, Density Map and Structural Coloring.
#

"""
Implementation of the "Chaos Game" algorithm to "draw" an IFS attractor in a matrix,
given a rectangular region, using the density map technique and structural coloring.
"""
function matrixstructuraldensitymap(ifs::AbstractVector{<:Function},
  xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real};
  weights::AbstractVector{<:Real} = Float64[], seed = nothing,
  preiterations::Int = 100, iterations::Int = 10000, #value = nothing, #::Function
  rot90::Bool = false, scale::Function = log, rgbscale::Function = sqrt,
  colormap::Union{Symbol, Vector{<:Colorant}} = :viridis, bgcolor::Union{Symbol, Colorant} = RGB(0,0,0),
  usealpha::Bool = false, gamma::Real = 1.0)

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
  mtrxs = AbstractArray{<:Real}[] # Histogram matrices
  for k in 1:nfs
    push!(mtrxs, rot90 ? fill(0, h, w) : fill(0, w, h)) # "Colors"
  end
  push!(mtrxs, rot90 ? fill(1, h, w) : fill(1, w, h)) # "Alpha"
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

  nmax = 0

  # Random orbit with "drawing point"
  for n in 1:iterations
    k = randindex()
    pk = ifs[k](pk)

    if isinsideclosed(pk, rr)
      indx = mtrxindx(pk, w, h, rr)
      mtrxs[k][indx...] += 1
      mtrxs[end][indx...] += 1
      if nmax < mtrxs[end][indx...]
        nmax = mtrxs[end][indx...]
      end      
    end
  end
  
  cm = typeof(colormap) == Symbol ? colorschemes[colormap] : ColorScheme(colormap)
  colors = [ cm[k/(nfs-1)] for k in 0:(nfs-1) ]
  bgc = RGB{Float64}(bgcolor) # Needed for deal with floating point arithmetic errors...
  if usealpha
    bgc = RGBA{Float64}(bgcolor)
  end
  mtrx = rot90 ? fill(bgc, h, w) : fill(bgc, w, h)
  dalpha = scale(nmax)
  brfct = gamma > 1.0 ? 1.0/gamma : 1.0

  for j in 1:h
    for i in 1:w
      if mtrxs[end][j,i] > 1
        d = mtrxs[end][j,i] - 1
        mtrx[j,i] = (mtrxs[1][j,i]/d)*colors[1]
        for k in 2:nfs
          mtrx[j,i] += (mtrxs[k][j,i]/d)*colors[k]
        end
        mtrx[j,i] = RGB{Float64}(
          rgbscale(  red(mtrx[j,i])),
          rgbscale(green(mtrx[j,i])),
          rgbscale( blue(mtrx[j,i]))
        ) * (scale(mtrxs[end][j,i])/dalpha)^brfct
      end
    end
  end  

  if usealpha
    for j in 1:h
      for i in 1:w
        if mtrxs[end][j,i] > 1
          d = mtrxs[end][j,i] - 1
          mtrx[j,i] = (mtrxs[1][j,i]/d)*colors[1]
          for k in 2:nfs
            mtrx[j,i] += (mtrxs[k][j,i]/d)*colors[k]
          end
          d = (scale(mtrxs[end][j,i])/dalpha)^brfct
          mtrx[j,i] = RGBA{Float64}(
            rgbscale(  red(mtrx[j,i])) * d,
            rgbscale(green(mtrx[j,i])) * d,
            rgbscale( blue(mtrx[j,i])) * d,
            d
          ) 
        end
      end
    end  
  end

  mtrxs = 0
  GC.gc()

  mtrx
end


"""
    imgstructuraldensitymap(ifs, xs, ys[; kwargs])

Given an IFS

\$\\{f_k:X \\rightarrow X\\}_{k=1}^{K},\$

where \$X=\\mathbb{R}^2\$ or \$X=\\mathbb{C}\$, draw their attractor in an image,
using the an implementation of the "Chaos Game" algorithm with
structural coloring and the density map technique.

#### Arguments
- `ifs::Vector{Function}`: An Iterated Functions System \$\\{f_k:X \\rightarrow X\\}\$.
- `xs::AbstractVector{Real}`: Base \$x\$ coordinates.
- `ys::AbstractVector{Real}`: Base \$y\$ coordinates.

#### Keyword Arguments
- `weights::Vector{Real} = []`: Probabilities associated to the IFS
- `seed = nothing`. Initial point for random orbit.
- `preiterations::Int = 100`: Iterations without drawing.
- `iterations::Int = 10000`: Iterations for drawing.
- `scale::Function = log`: Density scaling function \$(1,\\infty)\\rightarrow (0,\\infty)\$.
- `rgbscale::Function = sqrt`: Density RGB scaling function \$(0,1)\\rightarrow (0,1)\$.
- `colormap::Union{Symbol, Vector{Colorant}} = :viridis`: Colormap.
- `bgcolor::Colorant = RGBf(0,0,0)`: Background color.
- `gamma::Real = 1.0`: Gamma correction parameter.
- `usealpha::Bool = false`: Use alpha component.
"""
function imgstructuraldensitymap(ifs::AbstractVector{<:Function}, xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real};
  weights::AbstractVector{<:Real} = Float64[], seed = nothing,
  preiterations::Int = 100, iterations::Int = 10000,
  scale::Function = log, rgbscale::Function = sqrt,
  colormap::Union{Symbol, Vector{<:Colorant}} = :viridis, bgcolor::Colorant = RGBf(0,0,0),
  gamma::Real = 1.0, usealpha::Bool = false)

  matrixstructuraldensitymap(ifs, xs, ys, weights = weights, seed = seed,
    iterations = iterations, preiterations = preiterations,
    scale = scale, rgbscale = rgbscale,
    colormap = colormap, bgcolor = bgcolor, gamma = gamma, usealpha = usealpha,
    rot90 = true)
end


"""
    plotimgstructuraldensitymap(ifs, xs, ys[; kwargs])

Given an IFS

\$\\{f_k:X \\rightarrow X\\}_{k=1}^{K},\$

where \$X=\\mathbb{R}^2\$ or \$X=\\mathbb{C}\$, draw their attractor with **Makie**'s `image`,
using the an implementation of the "Chaos Game" algorithm with
structural coloring and the density map technique.

#### Arguments
- `ifs::Vector{Function}`: An Iterated Functions System \$\\{f_k:X \\rightarrow X\\}\$.
- `xs::AbstractVector{Real}`: Base \$x\$ coordinates.
- `ys::AbstractVector{Real}`: Base \$y\$ coordinates.

#### Keyword Arguments
- `weights::Vector{Real} = []`: Probabilities associated to the IFS
- `seed = nothing`. Initial point for random orbit.
- `preiterations::Int = 100`: Iterations without drawing.
- `iterations::Int = 10000`: Iterations for drawing.
- `scale::Function = log`: Density scaling function \$(1,\\infty)\\rightarrow (0,\\infty)\$.
- `rgbscale::Function = sqrt`: Density RGB scaling function \$(0,1)\\rightarrow (0,1)\$.
- `colormap::Union{Symbol, Vector{Colorant}} = :viridis`: Colormap.
- `bgcolor::Colorant = RGBf(0,0,0)`: Background color.
- `gamma::Real = 1.0`: Gamma correction parameter.
- `usealpha::Bool = false`: Use alpha component.
"""
@recipe(PlotImgStructuralDensityMap) do scene
  Attributes(
    weights = Float64[],
    seed = nothing,
    iterations = 10000,
    preiterations = 100,
    scale = log,
    rgbscale = sqrt,
    bgcolor = RGBf(0,0,0),
    gamma = 1.0,
    usealpha = false
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
  rgbscl = plt.rgbscale[]
  bgclr = plt.bgcolor[]
  gam = plt.gamma[]
  usea = plt.usealpha[]
  pltcm = plt.colormap[]

  image!(plt, obs_xs, obs_ys,
      matrixstructuraldensitymap(ifs, xs, ys, weights = ws, seed = sd,
        iterations = nits, preiterations = npits,
        scale = scl, rgbscale = rgbscl, bgcolor = bgclr, gamma = gam, usealpha = usea,
        colormap = pltcm);
      plt.attributes.attributes...)
  
  plt
end


#=
# ¿Cómo hacer mapas de densidad con coloreado estructural y mapa de calor (heatmap)?
# Me parece que es imposible...
# Se puede simular como en la implementación de plotstructuraldensitymap3d
=#


#
# 2D IFS Attractors with Chas Game, Density Map, Structural Coloring, Fractal Flame.
#

"""
Implementation of the "Chaos Game" algorithm to "draw" an IFS attractor in a matrix,
given a rectangular region, using the density map technique and structural coloring.
Fractal flame version.
"""
function matrixstructuraldensitymapff(ifs::AbstractVector{<:Function},
  xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real};
  weights::AbstractVector{<:Real} = Float64[], seed = nothing,
  preiterations::Int = 100, iterations::Int = 10000, #value = nothing, #::Function
  rot90::Bool = false, scale::Function = log,
  colormap::Union{Symbol, Vector{<:Colorant}} = :viridis, bgcolor::Union{Symbol, Colorant} = RGB(0,0,0),
  usealpha::Bool = false, gamma::Real = 1.0)

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
  histogram = rot90 ? fill(1, h, w) : fill(1, w, h)
  mtrxindx = rot90 ? torot90matrixindex : tomatrixindex

  cm = typeof(colormap) == Symbol ? colorschemes[colormap] : ColorScheme(colormap)
  colors = [ cm[k/(nfs-1)] for k in 0:(nfs-1) ]

  bgc = RGBf(bgcolor)
  if usealpha
    bgc = RGBAf(bgcolor)
  end
  mtrx = rot90 ? fill(bgc, h, w) : fill(bgc, w, h)

  # Initial data for random orbit
  pk = isnothing(seed) ? topoint(rand(xs), rand(ys)) : seed
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

      mtrx[indx...] = (mtrx[indx...] + colors[k])/2

      histogram[indx...] += 1
      if nmax < histogram[indx...]
        nmax = histogram[indx...]
      end
    end
  end

  d = scale(nmax)
  brfct = gamma > 1.0 ? 1.0/gamma : 1.0 # Brightness factor from gamma correction

  for j in 1:h
    for i in 1:w
      if histogram[j,i] > 1
        mtrx[j,i] *= (scale(histogram[j,i])/d)^brfct
      end
    end
  end

  mtrx
end


"""
    imgstructuraldensitymap(ifs, xs, ys[; kwargs])

Given an IFS

\$\\{f_k:X \\rightarrow X\\}_{k=1}^{K},\$

where \$X=\\mathbb{R}^2\$ or \$X=\\mathbb{C}\$, draw their attractor in an image,
using the an implementation of the "Chaos Game" algorithm with
structural coloring and the density map technique. Fractal flame version.

#### Arguments
- `ifs::Vector{Function}`: An Iterated Functions System \$\\{f_k:X \\rightarrow X\\}\$.
- `xs::AbstractVector{Real}`: Base \$x\$ coordinates.
- `ys::AbstractVector{Real}`: Base \$y\$ coordinates.

#### Keyword Arguments
- `weights::Vector{Real} = []`: Probabilities associated to the IFS
- `seed = nothing`. Initial point for random orbit.
- `preiterations::Int = 100`: Iterations without drawing.
- `iterations::Int = 10000`: Iterations for drawing.
- `scale::Function = log`: Density scaling function \$(1,\\infty)\\rightarrow (0,\\infty)\$.
- `colormap::Union{Symbol, Vector{Colorant}} = :viridis`: Colormap.
- `bgcolor::Colorant = RGBf(0,0,0)`: Background color.
- `gamma::Real = 1.0`: Gamma correction parameter.
- `usealpha::Bool = false`: Use alpha component.
"""
function imgstructuraldensitymapff(ifs::AbstractVector{<:Function}, xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real};
  weights::AbstractVector{<:Real} = Float64[], seed = nothing,
  preiterations::Int = 100, iterations::Int = 10000, scale::Function = log,
  colormap::Union{Symbol, Vector{<:Colorant}} = :viridis, bgcolor::Union{Symbol, Colorant} = RGB(0,0,0),
  usealpha::Bool = false, gamma::Real = 1.0)

  matrixstructuraldensitymapff(ifs, xs, ys, weights = weights, seed = seed,
    iterations = iterations, preiterations = preiterations, scale = scale,
    rot90 = true,
    colormap = colormap, bgcolor = bgcolor,
    usealpha = usealpha, gamma = gamma)
end


"""
    plotimgstructuraldensitymapff(ifs, xs, ys[; kwargs])

Given an IFS

\$\\{f_k:X \\rightarrow X\\}_{k=1}^{K},\$

where \$X=\\mathbb{R}^2\$ or \$X=\\mathbb{C}\$, draw their attractor with **Makie**'s `image`,
using the an implementation of the "Chaos Game" algorithm with
structural coloring and the density map technique. Fractal flame version.

#### Arguments
- `ifs::Vector{Function}`: An Iterated Functions System \$\\{f_k:X \\rightarrow X\\}\$.
- `xs::AbstractVector{Real}`: Base \$x\$ coordinates.
- `ys::AbstractVector{Real}`: Base \$y\$ coordinates.

#### Keyword Arguments
- `weights::Vector{Real} = []`: Probabilities associated to the IFS
- `seed = nothing`. Initial point for random orbit.
- `preiterations::Int = 100`: Iterations without drawing.
- `iterations::Int = 10000`: Iterations for drawing.
- `scale::Function = log`: Density scaling function \$(1,\\infty)\\rightarrow (0,\\infty)\$.
- `colormap::Union{Symbol, Vector{Colorant}} = :viridis`: Colormap.
- `bgcolor::Colorant = RGBf(0,0,0)`: Background color.
- `gamma::Real = 1.0`: Gamma correction parameter.
- `usealpha::Bool = false`: Use alpha component.
"""
@recipe(PlotImgStructuralDensityMap) do scene
  Attributes(
    weights = Float64[],
    seed = nothing,
    iterations = 10000,
    preiterations = 100,
    scale = log,
    bgcolor = RGBf(0,0,0),
    gamma = 1.0,
    usealpha = false
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
  bgclr = plt.bgcolor[]
  gam = plt.gamma[]
  usea = plt.usealpha[]
  pltcm = plt.colormap[]

  image!(plt, obs_xs, obs_ys,
      matrixstructuraldensitymapff(ifs, xs, ys, weights = ws, seed = sd,
        iterations = nits, preiterations = npits,
        scale = scl, bgcolor = bgclr, gamma = gam, usealpha = usea,
        colormap = pltcm);
      plt.attributes.attributes...)
  
  plt
end


#####
# 3D 
#####

"""
    plotdensitymap3d(ifs,xs, ys, zs [; kwargs])

Given an IFS

\$\\{f_k:\\mathbb{R}^3 \\rightarrow \\mathbb{R}^3\\}_{k=1}^{K},\$

plot their attractor with the **Makie**'s volume,
using the an implementation of the "Chaos Game" algorithm with density map.

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
- `scale::Function = log`: Density scaling function.
- `kind::Symbol = :volume`: **Makie**'s plot kind
  - `:volume`: `volume` for *volume render*.
  - `:contour`: `contour` for multiple *isosurfaces*.
  - `:volumeslices`: `volumeslices` for *volume slices*.

#### Notes
Looks better with plot kind `:volume` and algorithm `:absorption` with `absorption=8f0`.
"""
@recipe(PlotDensityMap3D) do scene
    Attributes(
        weights = Float64[],
        seed = nothing,
        iterations = 10000,
        preiterations = 0,
        kind = :volume, # Makie's plot volume by default
        scale = log #log-density-scale by default
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
    scl = plt.scale[]
    colormap = plt.colormap[]

    ### Density map algorithm

    # Rectangular 3D region and dimensions
    w, h, d = length(xs), length(ys), length(zs)
    rr = RectRegion3D(xs, ys, zs)

    # Cube to contain final densities 
    denscube = fill(1f0, w, h, d)

    ## "Chaos game" algorithm, accumulating densities

    # Initial data for random orbit
    pk = isnothing(sd) ? [rand(xs), rand(ys), rand(zs)] : sd
    nfs = length(ifs) # Number of functions in the system
    randindex = createrandindex(nfs, ws)
    #minval = 0f0
    maxval = 1f0
    
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
            denscube[cindex...] += 1f0
            if denscube[cindex...] > maxval
                maxval = denscube[cindex...]
            end
        end
    end

    ## end of "chaos game" algorithm

    # Scaling density (log is default)
    d = scl(maxval)
    for k in 1:length(zs)
      for j in 1:length(ys)
        for i in 1:length(xs)
          denscube[k,j,i] = denscube[k,j,i] > 1f0 ? Float32(scl(denscube[k,j,i])/d) : 0f0
        end
      end
    end
    
    ### End of density map algorithm

    # Force first color to totally transparent
    cm = to_colormap(colormap)
    cm[1] = RGBAf(red(cm[1]), green(cm[1]), blue(cm[1]), 0)

    if knd == :contour
        contour!(plt, obs_xs, obs_ys, obs_zs, denscube; plt.attributes.attributes..., colormap = cm)
    elseif knd == :volumeslices
        volumeslices!(plt, obs_xs, obs_ys, obs_zs, denscube; plt.attributes.attributes..., colormap = cm)
    else # knd == :volume
        volume!(plt, obs_xs, obs_ys, obs_zs, denscube; plt.attributes.attributes..., colormap = cm)
    end

    # Freeing memory
    denscube = 0
    GC.gc()

    plt
end


"""
    plotstructuraldensitymap3d(ifs, xs, ys, zx [; kwargs])

Given an IFS

\$\\{f_k:\\mathbb{R}^3 \\rightarrow \\mathbb{R}^3,\\}_{k=1}^{K},\$

plot their attractor with the **Makie**'s volume,
using the an implementation of the "Chaos Game" algorithm
with density map and simulated structural coloring.
The structural coloring is obtained with an adaptation from the fractal flame algorithm,
taht is, a simulated structural coloring.

#### Arguments
- `ifs::Vector{Function}`: An Iterated Functions System \$\\{f_k:X \\rightarrow X\\}\$.
- `xs::AbstractVector{Real}`: Base \$x\$ coordinates.
- `ys::AbstractVector{Real}`: Base \$y\$ coordinates.
- `zs::AbstractVector{Real}`: Base \$z\$ coordinates.

#### Keyword Arguments
- `weights::Vector{Real} = []`: Probabilities associated to the IFS
- `seed = nothing`. Initial point for random orbit.
- `preiterations::Int = 100`: Iterations without drawing.
- `iterations::Int = 10000`: Iterations for drawing.
- `kind::Symbol = :volume`: **Makie**'s plot kind
  - `:volume`: `volume` for *volume render*.
  - `:contour`: `contour` for multiple isosurfaces.
  - `:volumeslices`: `volumeslices` for volume slices.
- `scale::Function = log`: Scale function.

#### Notes
Looks better with plot kind `:volume` and algorithm `:mip`.
"""
@recipe(PlotStructuralDensityMap3D) do scene
    Attributes(
        weights = Float64[],
        seed = nothing,
        iterations = 10000,
        preiterations = 0,
        kind = :volume, # Makie's plot volume by default
        scale = log #log-density-scale by default
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
    scl = plt.scale[]
    colormap = plt.colormap[]

    # Rectangular 3D region and dimensions
    w, h, d = length(xs), length(ys), length(zs)
    rr = RectRegion3D(xs, ys, zs)

    # Cube to contain values and histogram
    cube = fill(0f0, w, h, d)
    histogram = fill(1, w, h, d)

    ## "Chaos game" algorithm

    # Initial data for random orbit
    pk = isnothing(sd) ? [rand(xs), rand(ys), rand(zs)] : sd
    nfs = length(ifs) # Number of functions in the system
    randindex = createrandindex(nfs, ws)
    maxval = 1
    
    # Random orbit without "drawing point"
    for n in 1:npits
        k = randindex()
        pk = ifs[k](pk)
    end

    v = Float32[ k/nfs for k in 1:nfs ] #sqrt
    v[1] = 1.5f0*v[1] # Trick to reduce transparece fo first color

    # Random orbit with "drawing point"
    for n in 1:nits
        k = randindex()
        pk = ifs[k](pk)

        if isinsideclosed(pk, rr)
            cindex = tocubeindex(pk,w,h,d,rr)
            cube[cindex...] = (cube[cindex...] + v[k])/2f0
            if histogram[cindex...] > maxval
                maxval = histogram[cindex...]
            end
        end
    end

    ## end of "chaos game" algorithm

    # Scaling density (log is default)
    smax = scl(maxval)
    for k in 1:d
      for j in 1:h
        for i in 1:w
          if histogram[k,j,i] > 1
            cube[k,j,i] *= Float32(scl(histogram[k,j,i])/smax)
          end
        end
      end
    end
    
    ### End of density map algorithm

    # Force first color to totally transparent
    cm = to_colormap(colormap)
    cm[1] = RGBAf(red(cm[1]), green(cm[1]), blue(cm[1]), 0)

    if knd == :contour
        contour!(plt, obs_xs, obs_ys, obs_zs, cube; plt.attributes.attributes..., colormap = cm)
    elseif knd == :volumeslices
        volumeslices!(plt, obs_xs, obs_ys, obs_zs, cube; plt.attributes.attributes..., colormap = cm)
    else # knd == :volume
        volume!(plt, obs_xs, obs_ys, obs_zs, cube; plt.attributes.attributes..., colormap = cm)
    end

    # Freeing memory
    denscube = 0
    histogram = 0
    GC.gc()

    plt
end
