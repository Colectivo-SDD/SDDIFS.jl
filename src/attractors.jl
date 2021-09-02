
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
