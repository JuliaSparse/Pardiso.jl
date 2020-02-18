# remove deps.jl if it exists, in case build.jl fails
isfile("deps.jl") && rm("deps.jl")

#################################################

println("\nMKL Pardiso")
println("=============")
function find_mklparadiso()
    if haskey(ENV, "MKLROOT")
        println("found MKLROOT environment varaible, enabling local MKL")
        return true
    end
    println("did not find MKLROOT environment variable, using provided MKL")
    return false
end

found_mklpardiso = find_mklparadiso()

open("deps.jl", "w") do f
    print(f,
"""
const LOCAL_MKL_FOUND = $found_mklpardiso
"""
)

end
