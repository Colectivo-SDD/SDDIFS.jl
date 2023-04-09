
"""
Verify that the sum of propabilities is 1.
"""
checkweights(weights::Vector{<:Real}) = ( 1.0-eps() < sum(weights) < 1.0+eps() )


"""
Create a function to create random indexes of functions in the IFS.
"""
function createrandindex(nfs::Int, weights::Vector{<:Real} = Float64[])
    if length(weights) == 0
        return () -> rand(1:nfs)
    end
    function idxrnd()
        prob = rand()
        wsum = 0
        for n in 1:length(weights)
            wsum += weights[n]
            if prob <= wsum
                return n
            end
        end
        1
    end
end
