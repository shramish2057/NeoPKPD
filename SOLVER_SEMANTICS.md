# Solver semantics

Solver behavior in OpenPKPD is treated as versioned numerical semantics.

## Version

Current solver semantics version: 1.0.0

## Supported solvers

The following solver identifiers are supported:

- :Tsit5 → OrdinaryDiffEq Tsitouras 5th order explicit Runge-Kutta
- :Rosenbrock23 → OrdinaryDiffEq Rosenbrock23 linearly implicit method

The mapping between symbol and algorithm is fixed within a solver semantics version.

## Grid handling

- saveat is mandatory
- Solver output is interpolated only at saveat points
- Returned solution time vector equals saveat exactly

## Tolerances

- reltol and abstol are passed directly to the solver
- No hidden defaults are applied by the engine

## Iteration limits

- maxiters is enforced explicitly
- Solver termination due to maxiters is treated as an error

## Rationale

Explicit solver semantics prevent silent numerical changes across versions
and allow long-term reproducibility of published simulations.
