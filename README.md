# maximum-principle_CAC_solver

MATLAB implementation of a maximum principle-preserving explicit finite difference
method for the conservative Allen--Cahn (CAC) equation with a space-time-dependent
Lagrange multiplier.

## Requirements

* MATLAB (no additional toolboxes required)

## Files

* `Fig1_Stability_test_1D.m` — Maximum-principle test in one dimension.
* `Fig2_Stability_test_1D_growth_rate.m` — Linear stability analysis, 1D growth rate.
* `Fig3_Stability_test_1D_growth_rate_profile.m` — 1D growth rate profile.
* `Fig4_Stability_test_2D.m` — Maximum-principle test in two dimensions.
* `Fig5a_Energy_dissipation_test_1D.m`, `Fig5b_Energy_dissipation_test_2D.m`, `Fig5c_Energy_dissipation_test_3D.m` — Evolution of the discrete Ginzburg--Landau energy.
* `Fig6_Tab1_Convergence_test.m` — Temporal convergence test for the explicit Euler method (1D/2D/3D).
* `Fig7_MBM_2D.m` — Mass-conserving motion by mean curvature (two circles, 2D).
* `Fig8_MBMCM_2D.m` — Mass-conserving motion by mean curvature (star-shaped interface, 2D).
* `Fig8_plot_MBMCM_2D.m` — Plots the saved 2D star-shaped interface snapshots.
* `Fig9_MBMCM_3D.m` — Mass-conserving motion by mean curvature (3D).
* `Fig9_plot_MBMCM_3D.m` — Plots the saved 3D interface snapshots.

## Numerical method

The conservative Allen--Cahn equation is discretized using a finite difference method
in space and the explicit Euler method in time. The space-time-dependent Lagrange
multiplier preserves the total mass.

The time-step condition used to preserve the discrete maximum principle is

$$
\Delta t \leq \frac{h^2\epsilon^2}{4h^2+2d\epsilon^2},
$$

where $h$ is the spatial grid size, $\epsilon$ is the interfacial-width parameter, and
$d$ is the spatial dimension.

## Usage

Open MATLAB, move to the repository directory, and run the desired script. For example,

```
Fig1_Stability_test_1D
```

The parameters, initial conditions, and output settings can be changed near the
beginning of each script.

## Reference

Youngjin Hwang and Junseok Kim, "Maximum principle-preserving explicit method for the
conservative Allen--Cahn equation," manuscript.

## Author

Youngjin Hwang and Junseok Kim
