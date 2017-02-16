# remove deps.jl if it exists, in case build.jl fails
using Compat

isfile("deps.jl") && rm("deps.jl")

@unix_only const LIBPARDISONAMES = [
    "libpardiso500-GNU461-X86-64",
    "libpardiso500-GNU472-X86-64",
    "libpardiso500-GNU481-X86-64",
    "libpardiso"
]

@windows_only const LIBPARDISONAMES = ["libpardiso500-WIN-X86-64.dll", "libpardiso"]

const PATH_PREFIXES = [
    dirname(@__FILE__),
    ""
]

# print to stderr, since that is where Pkg prints its messages
eprintln(x...) = println(STDERR, x...)

function find_paradisolib()
    found_lib = false
    for prefix in PATH_PREFIXES
        for libname in LIBPARDISONAMES
            try
                path = joinpath(prefix, libname)
                Libdl.dlopen(path, Libdl.RTLD_GLOBAL)
                global PARDISO_LIB_FOUND = true
                eprintln("found libpardiso at $path, using it")
                return path, true
            end
        end
    end
    eprintln("did not find libpardiso, assuming PARDISO 5.0 is not installed")
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
    warn("no Pardiso library managed to load")
end

open("deps.jl", "w") do f
    print(f,
"""
const MKL_PARDISO_LIB_FOUND = $found_mklpardiso
const PARDISO_LIB_FOUND = $found_pardisolib
const MKLROOT = R"$mklroot"
const PARDISO_PATH = R"$pardisopath"
"""
)

end
