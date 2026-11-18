# Compiler bug reports

Standalone reproducers for compiler bugs hit while porting MOM6 to GPUs. Each
subdirectory has a `README.md`, a runnable `reproducer.F90` and a `run.sh`.

## ifx workarounds in MOM6: `IFX_WORKAROUNDS_MOM6`

Building with `-DIFX_WORKAROUNDS_MOM6` enables source workarounds for ifx bugs,
so compilation can get further and expose the next bug. Each workaround is
guarded like this:

```fortran
#ifdef IFX_WORKAROUNDS_MOM6
  ! Avoid an ifx ICE, see bug_reports/<report_dir>
  <workaround>
#else
  <original code>
#endif
```

Without the macro the code is unchanged. Search for a report's directory name
to find every place its workaround is applied.

## Reports

| Directory | Compiler | Symptom | Workaround sites |
|---|---|---|---|
| [ifx_dc_dtype_bounds_ice](ifx_dc_dtype_bounds_ice/) | ifx 2026.1.0 | ICE #5633 when an offloaded `do concurrent` takes bounds from 2+ components of one derived-type object | 34 loops in 11 files |
| [ifx_map_extended_type_array_ice](ifx_map_extended_type_array_ice/) | ifx 2026.1.0 | ICE #5623 when OpenMP maps an array of an extended type whose allocatable component is in the extension (also host-only) | `MOM_diag_buffers.F90` (2 types) |
