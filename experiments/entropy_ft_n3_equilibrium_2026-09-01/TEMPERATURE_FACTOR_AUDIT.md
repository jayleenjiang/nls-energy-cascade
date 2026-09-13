# Source-level temperature-factor audit

Audited source: `flux/NLS_entropy_ft.cpp` at commit
`3d89659432fac0e512a1cc86fea2b63f8f849762`, SHA-256
`9ae5835ed708c8794c8b00ba799b23761482953aaf0ed47cd0b4ba3966d4eaf2`.

The boundary update in lines 361--377 uses

```cpp
const double noise_scale = std::sqrt(2.0 * GAMMA * temperature);
state.x[site] +=
    -GAMMA * state.force_real[site] * dt +
    noise_scale * sqrt_dt * normal_x;
state.y[site] +=
    -GAMMA * state.force_imag[site] * dt +
    noise_scale * sqrt_dt * normal_y;
```

The source convention in lines 14--19 identifies
`force_real=partial_x E` and `force_imag=partial_y E`.  Thus the continuous
bath SDE for either Cartesian boundary component is

`dX = -gamma partial_X(E) dt + sqrt(2 gamma T) dW`.

Its forward dissipative operator is

`L_bath^* p = gamma div(p grad E + T grad p)`.

For `p_eq proportional to exp(-E/T)`, `grad p_eq=-(p_eq/T) grad E`, so
`L_bath^* p_eq=0` exactly.  The diffusion coefficient in the forward equation
is `(1/2)*(sqrt(2 gamma T))^2=gamma T`; therefore the same input `T` appears in
the Gibbs exponent.  There is no hidden factor two.  Equivalently,
`T_effective=noise_scale^2/(2 gamma)=T`, giving exactly `6` and `10` for the
two controls.

This is an identity for the target continuous SDE, not an empirical validation
of the finite-`dt` split integrator.  The equal-temperature production runs
provide that independent numerical known-answer test.
