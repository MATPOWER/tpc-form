# TPC-Form
## Approximate formulations of AC power flows in Transformed Polar Coordinates (TPC)

[TPC-Form](https://github.com/MATPOWER/tpc-form) is a package of MATLAB M-files and classes for using a novel quadratic approximation of the AC power flow equations in [MATPOWER](https://github.com/MATPOWER/matpower). TPC-Form is based on code written by Wilson González Vanegas (`wgv`) using the new Object-Oriented Programming architecture of MATPOWER for developers. The current implementation allows for using a quadratic model where the power flow equations are approximated as a set of quadratic forms in transformed polar coordinates consisting of voltage angles and logarithm of the voltage magnitudes.

## System requirements
- [MATPOWER](https://github.com/MATPOWER/matpower) 8.1 or later.

## Installation
Installation and use of TPC-Form requires familiarity with the basic operation of MATLAB, including setting up your MATLAB path.

1. Clone the repository or download and extract the zip file of the TPC-Form distribution from the [TPC-Form project page](https://github.com/MATPOWER/tpc-form) to the location of your choice. We will use `<TPC>` to denote the path to this directory.
   
2. Add the following directory to your MATLAB path:
   
    - `<TPC>/lib/`

## Sample usage
We use `case118` as a test system to compare the power flow results obtained using the AC power flow equations and those calculated with the Quadratic Approximation in Transformed Polar Coordinates (QTPC). In the following code, notice that the `tpc` formulation is used via a [MATPOWER extension](https://matpower.org/doc/dev-manual/customizing.html#matpower-extensions) called `xt_tpc`:

```matlab
%% Load case
mpc = loadcase('case118');

%% Create an options struct with QTPC formulation enabled
mpopt_qtpc = mpoption('model','TPC');
mpopt_qtpc.pf.tpc.form = 'QUAD';

%% Run the exact AC power flow with default options
res_ac = run_pf(mpc);

%% Run the approximated power flow with QTPC formulation using the tpc extension
res_qtpc = run_pf(mpc, mpopt_qtpc, 'mpx', wgv.xt_tpc);
```
You should see something like the following printing in the Command Window:
```
MATPOWER Version 8.1, 12-Jul-2025
Power Flow -- AC-polar-power formulation

Newton's method converged in 3 iterations.
PF successful
.
.
.

MATPOWER Version 8.1, 12-Jul-2025
Power Flow -- TPC-QUAD formulation

Newton's method converged in 3 iterations.
PF successful
.
.
.
```
Finally, some metrics to compare the results with both formulations:
```matlab
%% Extract voltage angles and magnitudes with AC formulation
va_ac = res_ac.dm.elements.bus.tab.va;
vm_ac = res_ac.dm.elements.bus.tab.vm;

%% Extract voltage angles and magnitudes with QTPC formulation
va_qtpc = res_qtpc.dm.elements.bus.tab.va;
vm_qtpc = res_qtpc.dm.elements.bus.tab.vm;

%% Root Mean Square Error
nb = length(va_ac);
rmse_va = 1/sqrt(nb) * norm(va_qtpc - va_ac)
rmse_vm = 1/sqrt(nb) * norm(vm_qtpc - vm_ac)

%% Maximum percentage realtive error
mpre = max(abs( (vm_qtpc - vm_ac) ./ vm_ac)) * 100
```
Which should output the following results:
```
rmse_va =

    0.0527


rmse_vm =

   8.5975e-05


mpre =

    0.0617
```

## Documentation
Some Live Scripts and Live Functions included in `<TPC>/examples` and `<TPC>/other` can serve as illustrations of the theoretical background behind the QTPC formulation. 

## Citing
The mathematical formulation of QTPC is included in the draft paper ***A Quadratic Approximation of AC Power Flows in Transformed Polar Coordinates*** by Wilson González-Vanegas and Carlos E. Murillo-Sánchez, which is being considered for publication in the journal [ e-Prime - Advances in Electrical Engineering, Electronics and Energy](https://www.sciencedirect.com/journal/e-prime-advances-in-electrical-engineering-electronics-and-energy). The exact reference will be updated in the future in case of journal acceptance.
