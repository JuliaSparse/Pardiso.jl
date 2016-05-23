# remove deps.jl if it exists, in case build.jl fails
isfile("deps.jl") && rm("deps.jl")

const LIBPARDISONAMES = [
    "libpardiso500-GNU461-X86-64",
    "libpardiso500-GNU472-X86-64",
    "libpardiso500-GNU481-X86-64",
    "libpardiso"
]

const PATH_PREFIXES = [
    "",
    joinpath(Pkg.dir(), "Pardiso", "deps")
]

# print to stderr, since that is where Pkg prints its messages
eprintln(x...) = println(STDERR, x...)

PARDISO_LIB_FOUND = false
MKL_PARDISO_LIB_FOUND = false

function find_paradisolib()
    found_lib = false
    for prefix in PATH_PREFIXES
        for libname in LIBPARDISONAMES
            try
                path = joinpath(prefix, libname)
                Libdl.dlopen(path, Libdl.RTLD_GLOBAL)
                global PARDISO_LIB_FOUND = true
                eprintln("found libpardiso at $path")
                return path, true
            end
        end
    end
    return "", false
end

pardisopath, loaded = find_paradisolib()


if !(MKL_PARDISO_LIB_FOUND || PARDISO_LIB_FOUND)
    warn("no Pardiso library managed to load")
end

open("deps.jl", "w") do f
    print(f,
"""
const MKL_PARDISO_LIB_FOUND = $MKL_PARDISO_LIB_FOUND
const PARDISO_LIB_FOUND = $PARDISO_LIB_FOUND
const PARDISO_PATH = "$pardisopath"
"""
)

end
