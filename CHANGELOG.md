# Changelog

## [1.0.0] - 2025-01-24

### Features

- Allow for MKL v2025
- Bump version to 1.0
- add Changelog created with git-cliff

## [0.5.7] - 2024-04-22

### Features

- Remove unused convert in pardiso (#100)
- Fix CI badge + bump patch version

## [0.5.6] - 2024-03-04

### Features

- Adaptations for Panua Pardiso (#75)

Adaptations for PanuaPardiso, the successor of "project pardiso"

* rewrote corresponding parts of README

* Advise to set JULIA_PARDISO to `panua-pardiso-yyyymmdd-os/lib`, created by unzipping  the PanuaPardiso distribution file

* Additional lookup for libpardiso.so/dylib/dll in the JULIA_PARDISO directory

* export `mkl_is_available` and `panua_is_available` API call

Co-authored-by: Kristoffer Carlsson <kcarlsson89@gmail.com>
- Bump patch version to 0.5.6

## [0.5.5] - 2024-02-24

### Features

- Test on 1.6 as well (#92)
- Update TagBot.yml
- Bump MKL compat
- Update Project.toml
- Merge branch 'master' into kc/mkl_compat
- Merge pull request #91 from JuliaSparse/kc/mkl_compat

bump MKL compat
- Fix loading Pardiso on mac when MKL does not support it (#95)

* allow Pardiso.jl to load when MKL is not available

* add a small note about MKL versions on mac
- Removing hard-coded MKL in tests
- Using Pardiso.MKL_jll.is_available()
- Update tests for macOS

Trying to pin MKL_jll to 2023 for macOS.

Furthermore making the examples functions to take in the available solvers.
- Fixing mistakes in the test examples

Fixing mistakes in the test examples
- Update runtests.jl
- Updating readme to refer to MKL_jll 2023
- Merge pull request #97 from mipals/fixing-macos-tests

Removing hard-coded MKL in tests
- Fixes for self-installed MKL from intel oneAPI (#99)

* Fixes for self-installed MKL from intel oneAPI

* introduces Pardiso.mkl_is_available() for use in runtests.jl
* test load mkl in __init__ (like pardiso)
* some link updates, hint at LD_LIBRARY_PATH
* use mkl_is_available() everywhere
* make error msg less specific if it is not available.
* introduce MKL_LOAD_FAILED flag
* set to true in __init__ if libmkl_rt cannot be loaded
* checked by mkl_is_available(), so MKLPardisoSolver()
  cannot be instantiated if load failed

---------

Co-authored-by: Kristoffer Carlsson <kcarlsson89@gmail.com>

## [0.5.4] - 2022-03-01

### Features

- Update to allow MKL 2022 (#87)
- Bump version

## [0.5.3] - 2021-07-21

### Features

- Bump version

## [0.5.2] - 2021-07-09

### Features

- Allow StridedVecOrMat for RHS (#43)
- Update supported versions
- Adding support for gcc-6 to the Pardiso MKL solver (#46)
- Modify gfortran lib check (#47)
- Add Project.toml (#49)

* add Project.toml

* fix uuid

* address comments
- Check for libgfortran v9 (#54)
- Use MKL_jll if MKLROOT is not set
- Drop support for pardiso version 5
- Modernize travis
- Only test threads when running on a machine that supports it
- This finalizers are buggy
- Only try load gfortran on UNIX
- Merge pull request #56 from JuliaSparse/kc/mkl_jll

Use MKL_jll unless MKLROOT is set and drop support for version 5
- Better default if OMP_NUM_THREADS is not set
- Some typo fixes
- Merge pull request #57 from JuliaSparse/kc/threads

better default if OMP_NUM_THREADS is not set
- Warn about using MKL pardiso when Pardiso 6.0 is loaded
- Added schur_complement, pardisogetschur, wrapper for pardiso_get_schur (#39)

* added schur_complement, pardisogetschur, wrapper for pardiso_get_schur

* added docs for schur, tests (complex failing)

* docs

* removed schur_complement(...rows...) which never existed

* conjugate transpose

* version specification for schur

* low-level api readme

* testing for schur comp

* test for Schur comp

* only test if Pardiso 6.0 is loaded

Co-authored-by: Kristoffer Carlsson <kristoffer.carlsson@chalmers.se>
- Remove mentions of 5.0
- Always pardisoinit after changing matrixtype (#59)
- Add Travis badge
- Better ignoring of cov and mem files
- Remove workaround that is fixed in new MKL
- Improve and fix problems with examples
- Merge pull request #62 from JuliaSparse/improvements

Various Improvements
- Add Julia Computing to LICENSE
- Make Pardiso work with MKL in sysimage using 64 bit integers (#65)

* make Pardiso work with MKL in sysimage using 64 bit integers

* fix
- Bump version
- Test using travis_wait (#66)
- Fix nprocs acting weird in MKL when libpardiso is loaded (#67)
- Move to GHA (#79)
- Tweak computing the libmkl_rt path (#78)
- Update badge
- Update version

## [0.4.2] - 2019-03-26

### Features

- Added support for Pardiso compiled with gcc 8.0 on linux (#40)

## [0.4.1] - 2018-12-11

### Features

- Fix MKL Pardiso support on Linux. Code was trying to load libgomp twice, once in the way that works on 1.0 and once the old way. (#38)

## [0.4.0] - 2018-09-03

### Features

- Wip upgrade to 0.7 (#36)

* wip upgrade to 0.7

* fixes

* fixes

* fix missing variable

* fix Int32

* updates

* add .so extension

* updates

## [0.3.2] - 2018-01-30

### Features

- Fixes and readme updates for macOS (#30)

## [0.3.1] - 2017-08-25

### Features

- Add back loading libgomp (#26)

## [0.3.0] - 2017-08-01

### Features

- Remove mentions of BaseTestNext (#18)

no longer in test/REQUIRE
- Update README.md
- Use version number instead of release in .travis.yml (#19)

release will change over time, but your REQUIRE file says this
package supports julia 0.5 so it should continue to be tested
- Simplify loading (#22)

* simplify loading

* make mkl work on macos

* use showerror

* Update README.md

* fixes

* fixes

* precomp

* fix deprecation in example
- Add back windows dll (#23)

## [0.2.0] - 2017-02-24

### Features

- Fix deprecations on 0.6 (#16)
- Fix REQUIRE
- Loosen Compat requirement to 0.17 (#17)

parametric type aliases not used here currently,
only new abstract type syntax

## [0.1.3] - 2017-02-17

### Features

- Update README.md
- Update README.md
- Update README.md
- Fixed deprecation warnings thrown in Julia 0.5-rc2. Warnings were from @unix_only and @windows_only macros. (#12)
- Typo in LICENSE.md (#14)
- Fixes (#15)

* fixes

* fix

## [0.1.2] - 2016-05-27

### Features

- Fix typos
- Remove debug print
- Fixes and examples
- Dont load libgfortran on windows
- Ignore windows libs

## [0.1.1] - 2016-05-26

### Features

- Update README.md
- Update README.md
- Add forgotton BaseTestNext REQUIRE
- Avoid overwriting variables in loading

## [0.1.0] - 2016-05-24

### Features

- Rename pardiso.jl and update deprecations
- Fix deprecation
- Merge pull request #6 from JuliaSparse/kc/deprecs

fix deprecations and file nameing
- Document OMP_NUM_THREADS requirement

fixes #7
- Fix some julia highlighting
- Fix integer to int
- Use enums + build script + precompilation (#10)

use enums, support windows, add build script, enable precompilation
- Enums seems to not support precompilation
- Use @__FILE__ over Pkg.dir
- Fix working when MKL not loaded
- Update README.md
- Update README.md
- Change to matrixtype
- Update README.md
- Update README.md
- Add back opening gfortran because it is needed for prints

## [0.0.2] - 2015-05-25

### Features

- Fix some missed @compats
- Warn on no loaded pardiso library
- Update tests
- Export set_nprocs
- Update README.md
- Merge branch 'master' of https://github.com/KristofferC/Pardiso.jl
- Fix not loading threaded mkl

## [0.0.1] - 2015-05-21

### Features

- Initial empty commit
- Pardiso.jl generated files.

    license:  MIT
    authors:  Kristoffer Carlsson
    years:    2015
    user:

Julia Version 0.4.0-dev+4179 [c8c967c*]
- General updates
- Update README
- Small sentence fix
- Fix spelling of pardiso
- Fix reading of OMP_GET_THREADS variable
- Some clarifications
- Update to new API
- Remove unused variable
- Fix typos
- Removed unneeded commands from example
- Correct getter and setter for IPARM and DPARM
- Fix typo
- General updates
- Update README
- General fixes

- Release memory after solve
- Add missing phase to phase -> msg
- Remove B from checkmatrix arguments
- Add different number of factorizations to PardisoSolver type
- Add user permutation to PardisoSolver type
- Fix unclear variable name
- Update library list
- Update with MNUM, MAXFCT, PERM
- Do not return X from padiso
- Merge branch 'master' of https://github.com/KristofferC/Pardiso.jl

Conflicts:
	README.md
	src/Pardiso.jl
- Started on MKL
- Started on MKL
- Moved some things around
- Continued with MKL, add automatic mtype setting
- Started on MKL
- Moved some things around
- Continued with MKL, add automatic mtype setting
- Type stabilize error checks
- Fix hardcoded library path
- Improved matrix type detection etc
- Add LAPACK
- Fix typo
- Merge branch 'kc/mkl' of https://github.com/KristofferC/Pardiso.jl into kc/mkl

Conflicts:
	README.md
	src/Pardiso.jl
	src/mkl_pardiso.jl
	src/pardiso.jl
	test/runtests.jl
- Fixed diff between MKL and 5.0
- Merge branch 'kc/mkl'

Conflicts:
	src/Pardiso.jl
- Simplified logic
- Correct test description
- Update README.md
- Update README.md
- Update README.md
- Update README.md
- Update README.md
- Add info when pardiso fails to load
- Fix typo

<!-- generated by git-cliff -->
