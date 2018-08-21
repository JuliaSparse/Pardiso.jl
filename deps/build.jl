# remove deps.jl if it exists, in case build.jl fails
isfile("deps.jl") && rm("deps.jl")

using Libdl

const LIBPARDISONAMES = [
    "libpardiso600-WIN-X86-64.dll",
    "libpardiso600-MACOS-X86-64.dylib",
    "libpardiso600-GNU720-X86-64",
    "libpardiso500-WIN-X86-64.dll",
    "libpardiso500-MACOS-X86-64.dylib",
    "libpardiso500-GNU461-X86-64",
    "libpardiso500-GNU472-X86-64",
    "libpardiso500-GNU481-X86-64",
]

const PATH_PREFIXES = [
   @__DIR__,
]

# print to stderr, since that is where Pkg prints its messages
eprintln(x...) = println(stderr, x...)

pardiso_version = 0
function find_paradisolib()
    found_lib = false
    for prefix in PATH_PREFIXES
        for libname in LIBPARDISONAMES
            local path
            try
                path = joinpath(prefix, libname)
                Libdl.dlopen(path, Libdl.RTLD_GLOBAL)
                global PARDISO_LIB_FOUND = true
                eprintln("found libpardiso at $(abspath(path)), using it")
                if occursin("600", libname)
                    global pardiso_version = 6
                else
                    global pardiso_version = 5
                end
                return path, true
            catch e
                if isfile(path)
                    eprintln("found library but it failed to load due to:")
                    Base.showerror(stderr, e)
                end
            end
        end
    end
    eprintln("did not find libpardiso, assuming PARDISO 5/6 is not installed")
    return "", false
end

function find_mklparadiso()
    if haskey(ENV, "MKLROOT")
        eprintln("found MKLROOT key, using it")
        return ENV["MKLROOT"], true
    end
    eprintln("did not find MKLROOT key, assuming MKL is not installed")
    return "", false
end

pardisopath, found_pardisolib = find_paradisolib()
mklroot, found_mklpardiso = find_mklparadiso()

if !(found_mklpardiso || found_pardisolib)
    @warn("no Pardiso library managed to load")
end

open("deps.jl", "w") do f
    print(f,
"""
const MKL_PARDISO_LIB_FOUND = $found_mklpardiso
const PARDISO_LIB_FOUND = $found_pardisolib
const PARDISO_VERSION = $pardiso_version
const MKLROOT = $(repr(mklroot))
const PARDISO_PATH = raw"$pardisopath"
"""
)

end
