# ifx error #8682: `declare target` in both interface and body of a separate module procedure

## Summary

`ifx` rejects a separate module procedure that has `!$omp declare target` both in
its interface (in the parent module) and in its body (in a submodule):

```fortran
module parent
  interface
    elemental module function f(x) result(a)
      real, intent(in) :: x
      real :: a
      !$omp declare target
    end function f
  end interface
end module parent

submodule (parent) child
contains
  module procedure f
    !$omp declare target          ! error #8682
    a = 2.0 * x
  end procedure f
end submodule child
```

```
reproducer.F90(36): error #8682: Multiple !$OMP DECLARE TARGET declarations.
```

This is a front-end error, not an ICE. It happens with or without
`-fopenmp-targets`, and at both `-O0` and `-O3`.

## Why this looks like a compiler bug

OpenMP's Fortran restrictions for `declare target` say that a directive specified
in an interface block for a procedure must match a `declare target` directive in
the definition of that procedure. For a separate module procedure, the parent
module's interface body is that interface block and the `module procedure` body
is the definition. So writing the directive in both places should be the
conforming form, not an error.

Please check this against the OpenMP version you report under (4.5 through 6.0
contain this restriction) before filing.

## Environment

- Compiler: `ifx (IFX) 2026.1.0 20260617`
- System: Aurora (Intel GPU, `spir64` offload target)

## Reproduce

```bash
./run.sh
```

| Build           | Flags                                          | Result                     |
|-----------------|------------------------------------------------|----------------------------|
| bug             | `-qopenmp -fopenmp-targets=spir64`             | **error #8682**            |
| workaround      | `-qopenmp -fopenmp-targets=spir64 -DWORKAROUND`| compiles (run on a GPU node) |
| workaround_host | `-qopenmp -DWORKAROUND`                        | compiles, prints `PASS`    |

| Directive placement     | Result          |
|-------------------------|-----------------|
| interface and body      | **error #8682** |
| interface only          | OK              |
| body only               | OK              |

## Impact on MOM6

`exp_repro` (`src/framework/MOM_intrinsic_functions.F90` interface,
`src/framework/intrinsics/MOM_exp.F90` submodule body) uses this pattern, so
`MOM_exp.F90` fails to compile.

## Workaround

Under `-DIFX_WORKAROUNDS_MOM6`, the directive in the submodule body
(`MOM_exp.F90`) is omitted. The interface still declares `exp_repro` as a device
procedure.
