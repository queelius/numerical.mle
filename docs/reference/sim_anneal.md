# sim_anneal

This function implements the simulated annealing algorithm, which is a
global optimization algorithm that is useful for finding a good starting
point for a local optimization algorithm. We do not return this as an
MLE object because, to be a good estimate of the MLE, the gradient of
\`f\` evaluated at its solution should be close to zero, assuming the
MLE is interior to the domain of \`f\`. However, since this algorithm is
not guided by gradient information, it is not sensitive to the gradient
of \`f\` and instead only seeks to maximize \`f\`.

## Usage

``` r
sim_anneal(x0 = NULL, obj_fn = NULL, options = list(), ...)
```

## Arguments

- x0:

  Initial guess, default is NULL (must be specified in options)

- obj_fn:

  Objective function to maximize, default is NULL (must be specified in
  options)

- options:

  List of optional arguments

- ...:

  Additional arguments that may be passed to \`options\$neigh\`

## Value

list with best solution (argmax) and its corresponding objective
function value (max), and optionally path

## Functions

- `sim_anneal()`: options

## Fields

- `t_init`:

  Initial temperature

- `t_end`:

  Final temperature

- `alpha`:

  Cooling factor

- `iter_per_temp`:

  Number of iterations per temperature

- `max_iter`:

  Maximum number of iterations, used instead of t_end if not NULL,
  defaults to NULL

- `debug`:

  If TRUE, print debugging information to the console

- `trace`:

  If TRUE, track the history of positions and values

- `sup`:

  Support function, returns TRUE if x is in the domain of f

- `neigh`:

  Neighborhood function, returns a random neighbor of x

- `debug_freq`:

  Frequency of debug output, defaults to 10

- `obj_fn`:

  Objective function to maximize, if not specified in formal parameter
  \`obj_fn\`

- `x0`:

  Initial guess, if not specified in formal parameter \`x0\`
