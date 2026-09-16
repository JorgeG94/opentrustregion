![License](https://img.shields.io/github/license/eriksen-lab/opentrustregion)
![CI](https://github.com/eriksen-lab/opentrustregion/actions/workflows/main.yml/badge.svg)
[![codecov](https://codecov.io/github/eriksen-lab/opentrustregion/graph/badge.svg?token=NJIM6FDADD)](https://codecov.io/github/eriksen-lab/opentrustregion)
[![DOI](https://zenodo.org/badge/855894860.svg)](https://doi.org/10.5281/zenodo.17142554)

# OpenTrustRegion: A Reusable Library for Second-Order Trust Region Orbital Optimization

This library provides a robust and flexible implementation for second-order trust region orbital optimization, with extensive customization options to suit various use cases.

The following paper documents the theory and implementation of the methodology in OpenTrustRegion, and should be cited in any work using OpenTrustRegion:  

- Greiner, J.; Høyvik, I.-M.; Lehtola, S.; Eriksen, J. J. 
  A Reusable Library for Second-Order Orbital Optimization Using the Trust Region Method. Journal of Chemical Theory and Computation 2026, 22(2), 881–895. 
  DOI: [10.1021/acs.jctc.5c01576](https://doi.org/10.1021/acs.jctc.5c01576). 
  arXiv: [2509.13931](https://arxiv.org/abs/2509.13931).

## Installation

### Fortran or C Usage
To build the library for Fortran or C:

```sh
mkdir build
cd build
cmake ..
cmake --build .
```

The installation can be tested by running the ```testsuite.py``` file in the ```pyopentrustregion``` directory.

### Python Usage
To install the library for use with Python:

```sh
pip install .
```

The installation can be tested by running

```sh
python3 -m pyopentrustregion.testsuite
```

### CMake Configuration Options

The build process can be customized using the following CMake options:

| Option | Type | Default | Description |
|--------|------|----------|------------|
| **BUILD_SHARED_LIBS** | `BOOL` | `OFF` | Build shared libraries (`.so`, `.dylib`) instead of static ones. |
| **OpenTrustRegion_BUILD_TESTING** | `BOOL` | `ON` | Build the project’s testsuite. |
| **OpenTrustRegion_INSTALL_CMAKEDIR** | `STRING` | (auto) | Project install directory. |
| **CMAKE_BUILD_TYPE** | `STRING` | `Release` | Choose the build type (`Debug`, `Release`, etc.). |
| **INTEGER_SIZE** | `STRING` | *(auto)* | Set the integer precision to `4` (32-bit) or `8` (64-bit). Required when providing custom BLAS/LAPACK libraries. Otherwise defaults to 32-bit integers and tries to locate compatible BLAS and LAPACK libraries. Falls back to 64-bit integers if 32-bit libraries cannot be found. The resulting library name reflects the chosen integer precision (`libopentrustregion_32.*` or `libopentrustregion_64.*`) |
| **BLAS_LIBRARIES** | `PATH` | *(auto)* | Path(s) to BLAS libraries. If not provided, CMake attempts to locate a suitable BLAS automatically. |
| **LAPACK_LIBRARIES** | `PATH` | *(auto)* | Path(s) to LAPACK libraries. If not provided, CMake attempts to locate a suitable LAPACK automatically. |
| **OpenTrustRegion_HOST_PROVIDES_BLAS** | `BOOL` | `OFF` | When enabled, OpenTrustRegion will not attempt to detect or link BLAS/LAPACK and the testsuite is automatically disabled. The calling program must provide BLAS/LAPACK routines that expose the unsuffixed symbol names (for example `ddot`, `dsyev`) with an integer width matching `INTEGER_SIZE`; otherwise linking will fail loudly. |
| **OpenTrustRegion_ENABLE_XHOST** | `BOOL` | `ON` | Optimize the Release build for the current machine's instruction set (`-march=native` for GNU, `-xHost` for Intel). Automatically disabled when cross-compiling, regardless of this setting. Turn off when building for a different machine than the one compiling (e.g. packaging/conda builds). |

## Program Interfaces

The PySCF interface is available as an extension hosted at https://github.com/eriksen-lab/pyscf_opentrustregion. To install it, simply add its path to the **`PYSCF_EXT_PATH`** environment variable:
```sh
export PYSCF_EXT_PATH=path/to/pyscf_opentrustregion
```
Usage examples can be found in the **`examples`** directory of the PySCF interface repository.

The interface supports Hartree–Fock and DFT calculations via the **`mf_to_otr`** function, which wraps PySCF **`HF`** and **`KS`** objects into their OpenTrustRegion counterparts. Similarly, localization methods are available through the **`BoysOTR`**, **`PipekMezeyOTR`**, and **`EdmistonRuedenbergOTR`** classes, and state-specific CASSCF calculations are supported via the **`casscf_to_otr`** function applied to a PySCF **`CASSCF`** object. All returned objects are fully compatible with the original PySCF classes and can be used interchangeably.

Optional settings can be adjusted by modifying object attributes directly. Orbital optimization and internal stability analysis are performed using the **`kernel`** and **`stability_check`** member functions, respectively.

## Usage

The optimization process is initiated by calling a `solver` subroutine. This routine requires the following input arguments:

### Required Arguments

- **`update_orbs`** (subroutine):  
  Accepts and applies a variable update (e.g., orbital rotation), updates the internal state, and provides:
  - Objective function value (real)
  - Gradient (real array, written in-place)
  - Hessian diagonal (real array, written in-place)
  - A **`hess_x`** subroutine that performs Hessian-vector products:
    - Accepts a trial vector and writes the result of the Hessian transformation into an output array (real array, written in-place)
    - Returns an integer error code (0 for success, positive integers < 100 for errors)
  - Returns an integer error code (0 for success, positive integers < 100 for errors)
- **`obj_func`** (function):  
  Accepts and applies a variable update (e.g., orbital rotation) and returns:
  - Objective function value (real)
  - An integer error code (0 for success, positive integers < 100 for errors)
- **`n_param`** (integer): Specifies the number of parameters to be optimized.
- **`error`** (integer): An integer code indicating the success or failure of the solver. The error code structure is explained below.
- **`settings`** (settings_type): Settings object which controls optional arguments as described below.

---

The following Fortran snippet demonstrates how to use the `solver` interface:

```fortran
use opentrustregion, only: ip, rp, update_orbs_type, obj_func_type, solver_settings_type, solver

procedure(update_orbs_type), pointer :: update_orbs_funptr
procedure(obj_func_type), pointer :: obj_func_funptr
integer(ip) :: n_param, error
type(solver_settings_type) :: settings

! set callback function pointers to existing implementations
update_orbs_funptr => update_orbs
obj_func_funptr => obj_func

! initialize settings
call settings%init(error)

! override default settings
settings%conv_tol = 1e-6_rp
settings%n_macro = 100
settings%subsystem_solver = "tcg"

! run solver
call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
```

- Callback function pointers (`update_orbs_funptr`, `obj_func_funptr`) point to existing implementations elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- Solver settings are initialized using the `init()` method of the derived type, and default settings can be overridden (here, `conv_tol` and `n_macro`).
- Finally, the `solver` is called with the initialized settings and callback functions.

---

The callback functions above receive no host data, so a host program has to reach its own data through module-level variables. The `solver_ctx` entry point therefore accepts an opaque context which is handed to every callback function as its first argument. The context is an unlimited polymorphic argument, so it can be any type the host program defines, and the callback functions recover it with `select type`:

```fortran
use opentrustregion, only: ip, rp, update_orbs_ctx_type, obj_func_ctx_type, &
                           hess_x_ctx_type, solver_settings_type, solver_ctx

! host derived type holding everything the callback functions need
type :: host_type
    real(rp), allocatable :: mo_coeff(:, :), ints(:, :)
end type

type(host_type), target :: host
procedure(update_orbs_ctx_type), pointer :: update_orbs_funptr
procedure(obj_func_ctx_type), pointer :: obj_func_funptr
integer(ip) :: n_param, error
type(solver_settings_type) :: settings

! set callback function pointers to existing implementations
update_orbs_funptr => update_orbs
obj_func_funptr => obj_func

! initialize settings
call settings%init(error)

! run solver with the host data supplied through the context
call solver_ctx(update_orbs_funptr, obj_func_funptr, host, n_param, error, settings)
```

The callback functions have the same arguments as before, preceded by the context:

```fortran
subroutine update_orbs(context, kappa, func, grad, h_diag, hess_x_funptr, error)
    class(*), intent(inout), target :: context
    real(rp), intent(in), target :: kappa(:)
    real(rp), intent(out) :: func
    real(rp), intent(out), target :: grad(:), h_diag(:)
    procedure(hess_x_ctx_type), intent(out), pointer :: hess_x_funptr
    integer(ip), intent(out) :: error

    ! initialize error flag
    error = 0

    ! recover the host data from the context
    select type (host => context)
    type is (host_type)
        ! apply the variable update to host%mo_coeff and evaluate func, grad and
        ! h_diag from host%mo_coeff and host%ints
    class default
        error = 1
        return
    end select

    ! the Hessian linear transformation is itself context-carrying
    hess_x_funptr => hess_x

end subroutine update_orbs
```

- The context is passed to `solver_ctx` as an argument and is never stored by the library, so this feature adds no shared state of its own.
- The context dummy arguments carry the `target` attribute even though nothing in the library points at the context. This is an affordance for a host program which wants to keep a pointer to its own context inside a callback function, which is only valid if the actual argument passed by the host program has the `target` attribute as well.
- The abstract interfaces `update_orbs_ctx_type`, `obj_func_ctx_type`, `hess_x_ctx_type`, `precond_ctx_type`, `project_ctx_type`, and `conv_check_ctx_type` describe the context-carrying callback functions. The logging function needs no host data and therefore has no context-carrying counterpart.
- The optional callback functions are supplied as the `precond_ctx`, `project_ctx`, and `conv_check_ctx` settings instead of `precond`, `project`, and `conv_check`. Setting both flavours of the same callback function is refused with an error rather than silently resolved.
- `solver` is a thin adapter which bundles the callback functions of the plain interfaces into a context and calls `solver_ctx`, so both entry points run the same algorithm.
- Since that bundle is the context which `solver` passes to the callback functions, the context-carrying flavours of the optional callback functions require a context-carrying entry point: `solver` refuses to run when `precond_ctx`, `project_ctx`, `conv_check_ctx`, `stability_settings%precond_ctx`, or `stability_settings%project_ctx` is set, and `stability_check` refuses to run when `precond_ctx` or `project_ctx` is set.
- The reverse combination is supported: `solver_ctx` and `stability_check_ctx` accept the plain `precond`, `project`, and `conv_check` callback functions, which are called without a context.

---

The following C snippet demonstrates the equivalent usage through the C interface:

```c
#include <string.h>
#include "opentrustregion.h"

c_int n_param;

// set callback function pointers to existing implementations
update_orbs_fp update_orbs_funptr = update_orbs;
obj_func_fp obj_func_funptr = obj_func;

// initialize settings
solver_settings_type settings = solver_settings_init();

// override default settings
settings.conv_tol = 1e-6;
settings.n_macro = 100;
strcpy(settings.subsystem_solver, "tcg");

// run solver
c_int error = solver(update_orbs_funptr, obj_func_funptr, n_param, settings);
```

- Callback function pointers (`update_orbs_funptr`, `obj_func_funptr`) point to existing implementations elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- Solver settings are initialized via a small helper function `solver_settings_init()`, which returns a struct with default values. Individual settings (here, `conv_tol` and `n_macro`) can then be overridden.
- Finally, the `solver` is called with the initialized settings and callback functions and directly returns an error code in typical C fashion.

---

The following Python snippet demonstrates the equivalent usage through the Python interface:

```python
from pyopentrustregion import SolverSettings, solver

# initialize settings
settings = SolverSettings()

# override default settings
settings.conv_tol = 1e-6
settings.n_macro = 100
settings.subsystem_solver = "tcg"

# run solver
solver(update_orbs, obj_func, n_param, settings)
```

- Callback functions (`update_orbs`, `obj_func`) are defined elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- Solver settings are initialized via the `SolverSettings` class, which returns an object with default values; individual settings (here, `conv_tol`, and `n_macro`) can then be overridden.
- Finally, the `solver` is called with the initialized settings and callback functions and errors can be caught in pythonic fashion in the form of a `RuntimeException`.

### Optional Settings
The optimization process can be fine-tuned using the following settings:

- **`precond`** (subroutine): Applies a preconditioner to a residual vector. Writes the result in-place into a provided array and returns an integer error code (0 for success, positive integers < 100 for errors).
- **`project`** (subroutine): Applies a projection in-place to a provided vector and returns an integer error code (0 for success, positive integers < 100 for errors). Required for optimization using non-redundant parameters. When this is used, all other passed routines (`update_orbs`, `hess_x`, and `precond`) must be self-projecting.
- **`conv_check`** (function): Returns whether the optimization has converged due to some supplied convergence criterion. Additionally, outputs an integer code indicating the success or failure of the function, positive integers less than 100 represent error conditions.
- **`precond_ctx`**, **`project_ctx`**, **`conv_check_ctx`** (subroutine/function): Context-carrying counterparts of `precond`, `project`, and `conv_check`, which receive the context passed to `solver_ctx` as their first argument. Only one flavour of each callback function can be set.
- **`stability`** (boolean): Determines whether a stability check is performed upon convergence.
- **`line_search`** (boolean): Determines whether a line search is performed after every macro iteration.
- **`subsystem_solver`** (string): Specifies which subsystem solver to use. Options include:
  - `"davidson"`: standard Davidson method,
  - `"jacobi-davidson"`: Davidson method with fallback to Jacobi-Davidson if convergence is difficult, or automatically after `jacobi_davidson_start` micro iterations,
  - `"tcg"`: truncated conjugate gradient method.
- **`conv_tol`** (real): Specifies the convergence criterion for the RMS gradient.
- **`n_random_trial_vectors`** (integer): Number of random trial vectors used to initialize the micro iterations.
- **`start_trust_radius`** (real): Initial trust radius.
- **`n_macro`** (integer): Maximum number of macro iterations.
- **`n_micro`** (integer): Maximum number of micro iterations.
- **`jacobi_davidson_start`** (integer): Number of micro iterations after which the subsystem solver switches to the Jacobi-Davidson method.
- **`global_red_factor`** (real): Reduction factor for the residual during micro iterations in the global region.
- **`local_red_factor`** (real): Reduction factor for the residual during micro iterations in the local region.
- **`verbose`** (integer): Controls the verbosity of output during optimization.
- **`seed`** (integer): Seed value for generating random trial vectors.
- **`logger`** (subroutine): Accepts a log message. Logging is otherwise routed to stdout.
- **`stability_settings`** (stability_settings_type): Settings object controlling the internal stability check that is automatically performed upon convergence when `stability` is `True` or when starting at a stationary point (see the Stability Check section below). If `stability_settings%precond`, `stability_settings%project`, `stability_settings%precond_ctx`, `stability_settings%project_ctx`, or `stability_settings%logger` are left unset, they default to the corresponding `precond`, `project`, `precond_ctx`, `project_ctx`, and `logger` supplied to `solver`. `stability_settings%verbose` is raised to at least the solver's own `verbose` level.

## Stability Check
A separate `stability_check` subroutine is available to verify whether the current solution corresponds to a minimum. If not, it returns a boolean indicating instability and optionally, writes the eigenvector corresponding to the negative eigenvalue in-place to the provided memory.

### Required Arguments

- **`h_diag`** (real array): Represents the Hessian diagonal at the current point.
- **`hess_x`** (subroutine): Performs Hessian-vector products at the current point:
  - Accepts a trial vector and writes the result of the Hessian transformation into an output array (real array, written in-place)
  - Returns an integer error code (0 for success, positive integers < 100 for errors)
- **`stable`** (boolean): Returns whether the current point is stable.
- **`error`** (integer): An integer code indicating the success or failure of the solver. The error code structure is explained below.
- **`kappa`** (real array): If the memory is provided and the current point is not stable (as can be checked from return code of `stable`), the descent direction is written in-place in this array.
- **`settings`** (settings_type): Settings object which controls optional arguments as described below.

---

The following Fortran snippet demonstrates how to use the stability check interface:

```fortran
use opentrustregion, only: ip, rp, stability_settings_type, hess_x_type, stability_check

real(rp), allocatable :: h_diag(:), kappa(:)
procedure(hess_x_type), pointer :: hess_x_funptr
integer(ip) :: n_param, error
logical :: stable
type(stability_settings_type) :: settings

! set callback function pointer to existing implementation
hess_x_funptr => hess_x

! initialize settings
call settings%init(error)

! override default settings
settings%conv_tol = 1e-6_rp
settings%n_iter = 100
settings%diag_solver = "jacobi-davidson"

! run stability check
call stability_check(h_diag, hess_x_funptr, n_param, stable, error, settings, kappa=kappa)
```

- `hess_x_funptr` points to an existing Hessian-vector product implementation elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- Stability settings are initialized via the `init()` method of the derived type and can be overridden (here, `conv_tol` and `n_iter`).
- The `stable` logical output receives the result of the stability check.
- The descent direction `kappa` is optional and is only returned if provided.

---

As for the solver, a `stability_check_ctx` entry point accepts an opaque context which is handed to the Hessian linear transformation as its first argument:

```fortran
use opentrustregion, only: ip, rp, stability_settings_type, hess_x_ctx_type, &
                           stability_check_ctx

type(host_type), target :: host
real(rp), allocatable :: h_diag(:), kappa(:)
procedure(hess_x_ctx_type), pointer :: hess_x_funptr
integer(ip) :: error
logical :: stable
type(stability_settings_type) :: settings

! set callback function pointer to existing implementation
hess_x_funptr => hess_x

! initialize settings
call settings%init(error)

! run stability check with the host data supplied through the context
call stability_check_ctx(h_diag, hess_x_funptr, host, stable, error, settings, &
                         kappa=kappa)
```

- The callback function has the same arguments as the one of the plain `hess_x_type` interface, preceded by `class(*), intent(inout), target :: context`, and recovers the host data with `select type`.
- The optional `precond` and `project` callback functions are supplied as the `precond_ctx` and `project_ctx` settings, as described for the solver above.

---

The following C snippet demonstrates the equivalent usage through the C interface:

```c
#include <string.h>
#include "opentrustregion.h"

c_int n_param;
c_bool stable;

// set callback function pointer to existing implementation
hess_x_fp hess_x_funptr = hess_x;

// initialize settings
stability_settings_type settings = stability_settings_init();

// override default settings
settings.conv_tol = 1e-6;
settings.n_iter = 100;
strcpy(settings.diag_solver, "jacobi-davidson");

// pointers to Hessian diagonal and descent direction
double* h_diag;
double* kappa;

// run stability check
c_int error = stability_check(h_diag, hess_x_funptr, n_param, &stable, settings, kappa);
```

- `hess_x_funptr` points to an existing Hessian-vector product implementation elsewhere in the program.
- `n_param` and `h_diag` are assumed to be defined elsewhere.
- Stability settings are initialized via a small helper function `stability_settings_init()`, which returns a struct with default values; individual settings (here, `conv_tol` and `n_iter`) can then be overridden.
- The `stable` output receives the result of the stability check which directly returns an error code in typical C fashion.
- The descent direction `kappa` can be defined elsewhere if needed; otherwise, it can be set to `nullptr`.

---

The following Python snippet demonstrates the equivalent usage through the Python interface:

```python
from pyopentrustregion import StabilitySettings, stability_check

# initialize settings
settings = StabilitySettings()

# override default settings
settings.conv_tol = 1e-6
settings.n_iter = 100
settings.diag_solver = "jacobi-davidson"

# Hessian diagonal and descent direction arrays
h_diag = np.asarray(h_diag, dtype=np.float64)
kappa = np.empty(n_param, dtype=np.float64)

# run stability check
stable = stability_check(h_diag, hess_x, n_param, settings, kappa=kappa)
```

- `hess_x` is an existing Hessian-vector product implementation elsewhere in the program.
- `n_param` and `h_diag` are assumed to be defined elsewhere.
- Stability settings are initialized via the `StabilitySettings` class, which returns an object with default values; individual settings (here, `conv_tol`) can then be overridden.
- The `stable` output receives the result of the stability check and errors can be caught in pythonic fashion in the form of a `RuntimeException`.
- The descent direction `kappa` is optional and is only returned if provided.

### Optional Settings
The stability check can be fine-tuned using the following settings:

- **`precond`** (subroutine): Applies a preconditioner to a residual vector. Writes the result in-place into a provided array and returns an integer error code (0 for success, positive integers < 100 for errors).
- **`project`** (subroutine): Applies a projection in-place to a provided vector and returns an integer error code (0 for success, positive integers < 100 for errors). Required for stability check using non-redundant parameters. When this is used, all other passed routines (`hess_x` and `precond`) must be self-projecting.
- **`precond_ctx`**, **`project_ctx`** (subroutine): Context-carrying counterparts of `precond` and `project`, which receive the context passed to `stability_check_ctx` as their first argument. Only one flavour of each callback function can be set.
- **`diag_solver`** (string): Specifies which diagonalization solver to use. Options include:
  - `"davidson"`: standard Davidson method,
  - `"jacobi-davidson"`: Davidson method with fallback to Jacobi-Davidson if convergence is difficult, or automatically after `jacobi_davidson_start` micro iterations.
- **`conv_tol`** (real): Convergence criterion for the residual norm.
- **`n_random_trial_vectors`** (integer): Number of random trial vectors used to start the Davidson iterations.
- **`n_iter`** (integer): Maximum number of Davidson iterations.
- **`jacobi_davidson_start`** (integer): Number of micro iterations after which the subsystem solver switches to the Jacobi-Davidson method.
- **`verbose`** (integer): Controls the verbosity of output during the stability check.
- **`seed`** (integer): Seed value for generating random trial vectors.
- **`logger`** (function): Accepts a log message. Logging is otherwise routed to stdout.

## Error Code Structure

The library uses structured integer return codes to indicate whether a function has encountered an error. These codes follow the format **`OOEE`**, where:

- **`OO`** = Origin of the error (which component/function reported the error)
- **`EE`** = Specific error code

### General Rules

- A return code of `0` means success.
- Return codes between `1` and `99` are currently unused.
- All current error codes start from `100` and follow the `OOEE` structure.

### Origins (`OO`)

| Code Prefix (`OO`) | Component           |
|--------------------|---------------------|
| `01`               | `solver`            |
| `02`               | `stability_check`   |
| `11`               | `obj_func`          |
| `12`               | `update_orbs`       |
| `13`               | `hess_x`            |
| `14`               | `precond`           |
| `15`               | `conv_check`        |
| `16`               | `project`           |

### Error Codes (`EE`)

The error field (`EE`) is `01` for a general, unspecified error. Some origins define additional, more specific codes where the host program can plausibly react differently to them (e.g. by adjusting a setting and retrying):

| Error Code | Meaning                                                                |
|------------|------------------------------------------------------------------------|
| `0101`     | General error in `solver` |
| `0102`     | Orbital optimization did not converge within the maximum number of macro iterations (`n_macro`) |
| `0201`     | General error in `stability_check` |
| `0202`     | Stability check did not converge within the maximum number of iterations (`n_iter`) |

Future versions may define more specific codes for other actionable failure modes.

### Example Error Codes

| Error Code | Meaning                   |
|------------|---------------------------|
| `0101`     | General error in `solver` |
| `1201`     | Error in `update_orbs`    |

