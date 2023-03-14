# [Frequently Asked Questions](@id faq)

This page is a compilation of frequently asked questions and answers.

## [Stability and Divergence of ODE Solves](@id faq_stability)

For guidelines on debugging ODE solve issues, see
[PSA: How to help yourself debug differential equation solving issues](https://discourse.julialang.org/t/psa-how-to-help-yourself-debug-differential-equation-solving-issues/62489).

#### My model is reporting unstable results. What can I do?

First of all, don't panic. You may have experienced one of the following warnings:

> dt <= dtmin. Aborting. There is either an error in your model specification or the true solution is unstable.
> 
> NaN dt detected. Likely a NaN value in the state, parameters, or derivative value caused this outcome.
> 
> Instability detected. Aborting

These are all pointing to a similar behavior: for some reason or another, the
ODE solve is diverging to infinity. As it diverges to infinity, the `dt` of the
integrator will drop (trying to control the speed and error), so it will either
hit the minimum `dt`, hit `dt=NaN`, or have a value in the ODE hit `Inf`. Whichever
one occurs first will throw the respective warning.

How to handle this? 99.99% of the time this has been debugged, it has turned out
to be an error in the user's model! A missing minus sign, an incorrect term,
etc. There are many other behaviors to watch out for. In some ODEs, increasing
a parameter can cause a bifurcation so that the solution diverges. With
`u'=a*u`, if `a` is negative then, it nicely falls to zero, but if `a` is positive
the solution quickly diverges to infinity! This means, double-check your parameters
are indexed correctly!

**Note: if you see these warnings during a parameter estimation process, this is
likely the underlying problem. Simply check `sol.retcode != :Success` and throw
an `Inf` cost. Most optimizers will then reject steps in those parameter regimes!**

There are a few other things to check as well. Often, the stability of
an ODE solve improves as you decrease the tolerance, so you may want to try a
smaller `abstol` and `reltol`. One behavior to watch out for is that if your
model is a differential-algebraic equation and your DAE is of high index (say
index>1), this can impact the numerical solution. In this case, you may want to
use the [ModelingToolkit.jl index reduction tools](https://mtk.sciml.ai/dev/mtkitize_tutorials/modelingtoolkitize_index_reduction/)
to improve the numerical stability of a solve. In addition, if it's a highly
stiff ODE/DAE that is large, and you're using a matrix-free solver (such as GMRES),
make sure the tolerance of the GMRES is well-tuned and an appropriate preconditioner
is applied. Finally, try other solvers. They all have different stability, so try
`Tsit5()`, `Vern7()`, `QNDF()`, `Rodas5()`, `TRBDF2()`, `KenCarp4()`, `Sundials.CVODE_BDF()`,
etc. and see what works.

If none of this works out, double-check that your ODE truly has the behavior
that you believe it should. This is one of the most common issues: your intuition
may be deceiving. For example, `u' = -sqrt(u)` with `u(0)=1` cannot hit zero
because its derivative shrinks to zero, right? Wrong! [It will hit zero in a finite
time, after which the solution is undefined and does not have a purely real solution](https://www.wolframalpha.com/input/?i=u%27%3D-sqrt%28u%29).
`u' = u^2 - 100u` will “usually” go to zero, but if `u(0)>10` then it will go to
infinity. Plot out your diverging solution and see whether the asymptotics are
correct: if `u[i]` gets big, do your equations make `u'[i]` positive and growing?
That would be a problem!

Let's say you don't believe you made an error at all, and you want to file a bug
report. To do so, you'll first want to prove that it's isolated to a solver.
If it's a solver issue, then you shouldn't see it happen with every single solver.
Do you think it's an issue with the Julia solvers? Well fortunately,
DifferentialEquations.jl offers direct unmodified wrappers to almost all previously
built solvers, so if you think it's a Julia issue, try running your ODE through:

  - Sundials.jl, a wrapper for the C++ SUNDIALS library though `CVODE_Adams`,
    `CVODE_BDF`, `IDA`, and `ARKODE`.
  - ODEInterfaceDiffEq.jl, a wrapper for the classic Hairer Fortran codes like
    `dorpi5`, `dop853`, `radau`, `rodas`, etc.
  - LSODA.jl, a wrapper for the classic `lsoda` algorithm.
  - MATLABDiffEq.jl, a wrapper for the MATLAB ODE solvers `ode45`, `ode15s`, etc.
  - SciPyDiffEq.jl, a wrapper for SciPy's `odeint` (LSODA) and other methods (LSODE, etc.).
  - deSolveDiffEq.jl, a wrapper for the commonly used R library.

And many more. Testing this is as simple as changing `solve(prob,Tsit5())` to
`solve(prob,lsoda())`, so please give this a try. If you translated your code
from another language, like Python or MATLAB, use the direct wrapper to double
check the steps are the same. If they are not, then your ODE is not the same,
because it's using a direct call to the solvers of those packages!

If your ODE diverges to infinity with every ODE solver ever made, the problem is
most likely not the ODE solvers. Or rather, to put it in meme form:

![](https://user-images.githubusercontent.com/1814174/120933617-eb65ac80-c6c8-11eb-85f7-ef98688d054c.jpg)

Don't be like Patrick. If after trying these ideas, your ODE solve still seems
to have issues, and you haven't narrowed it down, feel free to ask on the
[Julia Discourse](https://discourse.julialang.org/) to get some help diagnosing
it. If you did find a solver issue, please open an issue on the GitHub repository.

#### A larger maxiters seems to be needed, but it's already high?

If you see:

> Interrupted. Larger maxiters is needed.

Note that it could quite possibly arise just from having a very long timespan.
If you check `sol.t` from the returned object, and it looks like it's stepping
at reasonable lengths, feel free to just pass `maxiters=...` into solve to
bump it up from the default of `Int(1e5)`.

But if your `maxiters` is already high, then the problem is likely that your model
is stiff. A stiff ODE requires very small timesteps from many explicit solvers,
such as `Tsit5()`, `Vern7()`, etc., and thus those methods are not appropriate
for this kind of problem. You will want to change to a different method, like
`Rodas5()`, `Rosenbrock23()`, `TRBDF2()`, `KenCarp4()`, or `QNDF()`.

#### My ODE goes negative but should stay positive, what tools can help?

There are many tools to help! However, let's first focus on one piece first:
when you say “should” be positive, what do you mean by “should”? If you mean
“mathematically you can prove that the ODE with these values and these initial
conditions will have a solution that is positive for all time”, then yes, you're
looking in the right place. If by “should” you mean “it's a model of biochemical
reactions, so the concentration should always be positive”, then ask yourself
first, did you write down a model where it will always be positive?

The following set of tools are designed to accuracy enforce positivity in ODE
models, which mathematically should be positive in the true solution. If they
encounter a model that is actually going negative, they will work really hard
to get a positive but correct solution, which is impossible, so they will simply
error out. This can be more subtle than you think. Solving `u'=-sqrt(u)` is
not guaranteed to stay positive, even though the derivative goes to zero as `u`
goes to zero (check the analytical solution if you're curious). Similarly,
analyzing nonlinear models can showcase all sorts of behavior. A common cause
for accidental negativity is Hill functions in systems biology models: just
because derivatives go to zero doesn't mean they are going to zero fast enough
to keep things positive!

With that in mind, let's see the options.

The simplest trick is to change the solver tolerance. Reduce `abstol` (and maybe
`reltol`) a bit. That can help reduce the error and thus keep the solution positive.
For some more difficult equations, changing to a stiff ODE solver like `Rosenbrock23()`
`QNDF`, or `TRBDF2()` can be helpful.

If those don't work, call out the big guns. One of them is `isoutofdomain`, where
you can define a boolean function which will cause step rejections whenever it
is not satisfied. For example, `isoutofdomain = (u,p,t)->any(x->x<0,u)` will make
the solver reject any step which cases any variable `u` to go negative. Now, using
any pure-Julia solver with this option, it's impossible to get a negative in the
result! One thing you may see though is:

> dt <= dtmin. Aborting. There is either an error in your model specification or the true solution is unstable.

or

> Interrupted. Larger maxiters is needed.

What this means is that enforcing positivity is not possible. It keeps rejecting
steps that go negative, reducing `dt`, taking another step, rejecting, reducing,
repeat until `dt` hits `dtmin` or it hits maxiters. This means that even when
trying to solve the problem with the most accurate infinitesimal `dt`, the
solution still goes negative. Are you sure the true solution is supposed to be
positive? If you see this, check for issues like a missing minus sign in your
equations.

If that works but is a little slow, the domain handling callbacks in
[the callback library](@ref callback_library) are designed to function similarly
but in a way that gets better performance. Instead of repeating lots of steps
through rejections, it interpolates back to still take a smaller step, always
progressing forwards. However, this can be a bit less stable, so its applicability
depends on the equation, and once again this requires that the solution is
truly positive. If the true solution goes negative, it will repeatedly try
interpolating backwards until it can no longer and end with a `dtmin` issue.

Finally, note that ODE solvers will not be more correct than tolerance, and so
one should expect that if the solution is supposed to be positive but `abstol=1e-12`,
you may end up with `u[i]=-1e-12`. That is okay,
[that is expected behavior of numerical solvers](https://www.radford.edu/%7Ethompson/RP/nonnegative.pdf),
the ODE solver is still doing its job. If this is a major issue for your application,
you may want to write your model to be robust to this behavior, such as changing
`sqrt(u[i])` to `sqrt(max(0,u[i]))`. You should also consider transforming your
values, like solving for `u^2` or `exp(u)` instead of `u`, which mathematically
can only be positive. Look into using a tool like [ModelingToolkit.jl](https://mtk.sciml.ai/dev/)
for automatically transforming your equations.

## [Performance](@id faq_performance)

#### GPUs, multithreading and distributed computation support

Yes. The `*`DiffEq.jl libraries (OrdinaryDiffEq.jl, StochasticDiffEq.jl, and
DelayDiffEq.jl) are all written to be generic to the array and number types.
This means they will adopt the implementation that is given by the array type.
The in-place algorithms internally utilize Julia's broadcast (with some exceptions
due to a Julia bug for now, see [this issue](https://github.com/SciML/OrdinaryDiffEq.jl/issues/106))
and Julia's `mul!` in-place matrix multiplication function. The out-of-place
algorithms utilize standard arithmetical functions. Both additionally utilize
the user's norm specified via the common interface options and, if a stiff
solver, ForwardDiff/DiffEqDiffTools for the Jacobian calculation, and Base linear
factorizations for the linear solve. For your type, you may likely need to give
a [better form of the norm](https://docs.sciml.ai/DiffEqDocs/stable/extras/timestepping/#Common-Setup),
[Jacobian](@ref performance_overloads),
or [linear solve calculations](@ref linear_nonlinear)
to fully utilize parallelism.

GPUArrays.jl (CuArrays.jl), ArrayFire.jl, DistributedArrays.jl have been tested and work in
various forms, where the last one is still not recommended for common use
yet.

The next question is whether it matters. Generally, your system has to be large
for parallelism to matter. Using a multithreaded array for broadcast we find
helpful around `N>1000`, though the Sundials manual says `N>100,000`. For high
order Runge-Kutta methods it's likely lower than the Sundials estimate because
of more operations packed into each internal step, but as always, that will need
more benchmarks to be precise and will depend on the problem being solved. GPUs
generally require some intensive parallel operation in the user's `f` function
to be viable, for example, a matrix multiplication for a stencil computation
in a PDE. If you're simply solving some ODE element-wise on a big array, it likely
won't do much, or it will slow things down just due to how GPUs work.
DistributedArrays require parallel linear solves to really matter, and thus are
only recommended when you have a problem that cannot fit into memory or are using
a stiff solver with a Krylov method for the linear solves.

### Note About Setting Up Your Julia Installation for Speed: BLAS Choices

Julia uses an underlying BLAS implementation for its matrix multiplications
and factorizations. This library is automatically multithreaded and accelerates
the internal linear algebra of DifferentialEquations.jl. However, for optimality,
you should make sure that the number of BLAS threads that you are using matches
the number of physical cores and not the number of logical cores. See
[this issue for more details](https://github.com/JuliaLang/julia/issues/33409).

To check the number of BLAS threads, use:

```julia
ccall((:openblas_get_num_threads64_, Base.libblas_name), Cint, ())
```

If I want to set this directly to 4 threads, I would use:

```julia
using LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(4)
```

Additionally, sometimes Intel's MKL might be a faster BLAS than the standard
BLAS that ships with Julia (OpenBLAS). To switch your BLAS implementation, you
can use [MKL.jl](https://github.com/JuliaLinearAlgebra/MKL.jl), which will accelerate
the linear algebra routines. This is done via:

```julia
using MKL
```

#### My ODE is solving really slow

First, check for bugs. These solvers go through a ton of convergence tests and
so if there's a solver issue, it's either just something to do with how numerical
methods work or it's a user-error (generally the latter, though check the later
part of the FAQ on normal numerical errors). User-errors in the `f` function
causing a divergence of the solution is the most common reason for reported
slow codes.

If you have no bugs, great! The standard tricks for optimizing Julia code then
apply. Take a look at the [Optimizing DiffEq Code tutorial](https://tutorials.sciml.ai/html/introduction/03-optimizing_diffeq_code.html)
for some tips and pointers.

What you want to do first is make sure your function does not allocate.
If your system is small (`<=100` ODEs/SDEs/DDEs/DAEs?), then you should set your
system up to use [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl).
This is demonstrated
[in the ODE tutorial](@ref ode_other_types)
with static matrices. Static vectors/arrays are stack-allocated, and thus creating
new arrays is free, and the compiler doesn't have to heap-allocate any of the
temporaries (that's the expensive part!). These have specialized superfast
dispatches for arithmetic operations and extra things like LU-factorizations,
and thus they are preferred when possible. However, they lose efficiency if they
grow too large.

For anything larger, you should use the `in-place` syntax `f(du,u,p,t)` and make
sure that your function doesn't allocate. Assuming you know of a `u0`, you
should be able to do:

```julia
du = similar(u0)
@time f(du, u0, p, t)
```

and see close to zero allocations and close to zero memory allocated. If you see
more, then you might have a type-instability or have temporary arrays. To find
type-instabilities, you should do:

```julia
@code_warntype f(du, u, p, t)
```

and read the printout to see if there are any types that aren't inferred by the
compiler, and fix them. If you have any global variables, you should make them
`const`. As for allocations, some common things that allocate
are:

  - Array slicing, like `u[1:5]`. Instead, use `@view u[1:5]`
  - Matrix multiplication with `*`. Instead of `A*b`, use `mul!(c,A,b)` for some
    pre-allocated cache vector `c`.
  - Non-broadcasted expressions. Every expression on arrays should `.=` into another
    array, or it should be re-written to loop and do computations with scalar (or
    static array) values.

For an example of optimizing a function resulting from a PDE discretization, see
[this blog post](http://www.stochasticlifestyle.com/solving-systems-stochastic-pdes-using-gpus-julia/).

#### The stiff solver takes forever to take steps for my PDE discretization

The solvers for stiff solvers require solving a nonlinear equation each step.
To do so, they have to do a few Newton steps. By default, these methods
assume that the Jacobian is dense, automatically calculate the Jacobian for you,
and do a dense factorization. However, in many cases you may want to use alternatives
that are more tuned for your problem.

First of all, when available, it's recommended that you pass a function for computing
your Jacobian. This is discussed in the [performance overloads](@ref performance_overloads)
section. Jacobians are especially helpful for Rosenbrock methods.

Secondly, if your Jacobian isn't dense, you shouldn't use a dense Jacobian!
Instead, if you're using  a `*DiffEq` library you should
[specify a linear solver](@ref linear_nonlinear) and/or a `jac_prototype` for the
matrix form, and for Sundials.jl, you should change the `linear_solver` option. See
[the ODE solve Sundials portion](@ref ode_solve_sundials)
for details on that.

Right now, `QNDF` is the recommended method for stiff
problems with large sparse Jacobians. You should specify `jac_prototype` as a
special matrix, such as a banded or tridiagonal matrix, if it satisfies a special
structure. If you only know the Jacobian is sparse, using automated sparsity
detection can help with identifying the sparsity pattern. See the [stiff ODE tutorial](@ref stiff)
for more details. Lastly, using `LinSolveGMRES()` can help if a sparsity pattern
cannot be obtained, but the matrix is large, or if the sparsity cannot fit into
memory. Once again, a good reference for how to handle PDE discretizations can be found
[at this blog post](http://www.stochasticlifestyle.com/solving-systems-stochastic-pdes-using-gpus-julia/).

#### My Problem Has Discontinuities and is Unstable / Slow

[This Discourse post](https://discourse.julialang.org/t/handling-instability-when-solving-ode-problems/9019/5)
goes into detail for how to handle discontinuities in your ODE function and how
to use that extra information to speed up the solver.

## Complicated Models

#### Switching ODE functions in the middle of integration

There are a few ways to do this. The simplest way is to just have a parameter to
switch between the two. For example:

```julia
function f(du, u, p, t)
    if p == 0
        du[1] = 2u[1]
    else
        du[1] = -2u[1]
    end
    du[2] = -u[2]
end
```

Then in a callback, you can make the `affect!` function modify `integrator.prob.p`.
For example, we can make it change when `u[2]<0.5` via:

```julia
condition(t, u, integrator) = u[2] - 0.5
affect!(integrator) = integrator.p = 1
```

Then it will change between the two ODE choices for `du1` at that moment.
Another way to do this is to make the ODE functions all be the same type
via FunctionWrappers.jl, but that is unnecessary. With the way that modern
processors work, there exists branch prediction and thus execution of a conditional
is free if it's predictable which branch will be taken. In this case, almost every
call to `f` takes the `p==0` route until the callback, at which point it is
almost always the `else` route. Therefore, the processor will effectively get
rid of the computational cost associated with this, so you're likely over-optimizing
if you're going further (unless this change happens every step, but even then,
this is probably the cheapest part of the computation…).

## Numerical Error

#### What does tolerance mean and how much error should I expect

The most useful options are the tolerances `abstol` and `reltol`. These tell the
internal adaptive time stepping engine how precise of a solution you want.
Generally, `reltol` is the relative accuracy while `abstol` is the accuracy when
`u` is near zero. *These tolerances are local tolerances and thus are not global
guarantees*. However, a good rule of thumb is that the total solution accuracy
is 1-2 digits less than the relative tolerances. Thus, for the defaults
`abstol=1e-6` and `reltol=1e-3`, you can expect a global accuracy of about 1-2
digits. This is standard across the board and applies to the native Julia methods,
the wrapped Fortran and C++ methods, the calls to MATLAB/Python/R, etc.

#### The solver doesn't obey physical law X (e.g. conservation of energy)

Yes, this is because the numerical solution of the ODE is not the exact solution.
There are a few ways that you can handle this problem. One way is to get a more
exact solution. Thus instead of

```julia
sol = solve(prob, alg)
```

use

```julia
sol = solve(prob, alg, abstol = 1e-10, reltol = 1e-10)
```

Of course, there's always a tradeoff between accuracy and efficiency, so play
around to find out what's right for your problem.

Another thing you can do is use a callback. There are some
[premade callbacks in the callback library](@ref callback_library) which
handle these sorts of things like projecting to manifolds and preserving positivity.

##### Symplectic integrators don't conserve energy

Yes, symplectic integrators do not exactly conserve energy. It is a common
misconception that they do. What symplectic integrators actually do is solve
for a trajectory which rests on a symplectic manifold that is perturbed from
the true solution's manifold by the truncation error. This means that symplectic
integrators do not experience (very much) longtime drift, but their orbit is
not exactly the same as the true solution in phase space, and thus you will
see differences in energy that tend to look periodic. There is a small drift
which grows linearly and is related to floating-point error, but this drift
is much less than standard methods. This is why symplectic methods are recommended
for longtime integration.

For conserving energy, there are a few things you can do. First of all, the energy
error is related to the integration error, so simply solving with higher accuracy
will reduce the error. The results in the
[SciMLBenchmarks](https://github.com/SciML/SciMLBenchmarks.jl) show
that using a `DPRKN` method with low tolerance can be a great choice. Another
thing you can do is use
[the ManifoldProjection callback from the callback library](@ref callback_library).

#### How to get to zero error

You can't. For floating-point numbers, you shouldn't use below `abstol=1e-14`
and `reltol=1e-14`. If you need lower than that, use arbitrary precision numbers
like BigFloats or [ArbFloats.jl](https://github.com/JuliaArbTypes/ArbFloats.jl).

## Autodifferentiation and Dual Numbers

#### Native Julia solvers compatibility with autodifferentiation

Yes, they are compatible with automatic differentiation! Take a look at the
[sensitivity analysis](https://docs.sciml.ai/SciMLSensitivity/stable/) page for more details.

If the algorithm does not have differentiation of parameter-dependent events,
then you simply need to make the initial condition have elements of Dual numbers.
If the algorithm uses Dual numbers, you need to make sure that time is also
given by Dual numbers.

To show this in action, let's say we want to find the Jacobian of solution
of the Lotka-Volterra equation at `t=10` with respect to the parameters.

```@example faq1
using DifferentialEquations
function func(du, u, p, t)
    du[1] = p[1] * u[1] - p[2] * u[1] * u[2]
    du[2] = -3 * u[2] + u[1] * u[2]
end
function f(p)
    prob = ODEProblem(func, eltype(p).([1.0, 1.0]), (0.0, 10.0), p)
    # Lower tolerances to show the methods converge to the same value
    solve(prob, Tsit5(), save_everystep = false, abstol = 1e-12, reltol = 1e-12)[end]
end
```

This function takes in new parameters and spits out the solution at the end.
We make the initial condition `eltype(p).([1.0,1.0])` so that way it's typed to
be Dual numbers whenever `p` is an array of `Dual` numbers, and we do the same
for the timespan just to show what you'd do if there were parameters-dependent events.
Then we can take the Jacobian via ForwardDiff.jl:

```@example faq1
using ForwardDiff
ForwardDiff.jacobian(f, [1.5, 1.0])
```

and compare it to FiniteDiff.jl:

```@example faq1
using FiniteDiff
FiniteDiff.finite_difference_jacobian(f, [1.5, 1.0])
```

#### I get Dual number errors when I solve my ODE with Rosenbrock or SDIRK methods

This is because you're using a cache which is incompatible with autodifferentiation
via ForwardDiff.jl. For example, if we use the ODE function:

```julia
using LinearAlgebra, OrdinaryDiffEq
function foo(du, u, (A, tmp), t)
    mul!(tmp, A, u)
    @. du = u + tmp
    nothing
end
prob = ODEProblem(foo, ones(5, 5), (0.0, 1.0), (ones(5, 5), zeros(5, 5)))
solve(prob, Rosenbrock23())
```

Here we use a cached temporary array to avoid the allocations of matrix
multiplication. When autodifferentiation occurs, the element type of `u` is
`Dual` numbers, so `A*u` produces `Dual` numbers, so the error arises when it
tries to write into `tmp`. There are two ways to avoid this. The first way,
the easy way, is to just turn off autodifferentiation with the `autodiff=false`
option in the solver. Every solver which uses autodifferentiation has this option.
Thus, we'd solve this with:

```julia
prob = ODEProblem(f, ones(5, 5), (0.0, 1.0))
sol = solve(prob, Rosenbrock23(autodiff = false))
```

and it will use a numerical differentiation fallback (DiffEqDiffTools.jl) to
calculate Jacobians.

We could use `get_tmp` and `dualcache` functions from
[PreallocationTools.jl](https://github.com/SciML/PreallocationTools.jl)
to solve this issue, e.g.,

```julia
using LinearAlgebra, OrdinaryDiffEq, PreallocationTools
function foo(du, u, (A, tmp), t)
    tmp = get_tmp(tmp, first(u) * t)
    mul!(tmp, A, u)
    @. du = u + tmp
    nothing
end
prob = ODEProblem(foo, ones(5, 5), (0.0, 1.0),
                  (ones(5, 5), PreallocationTools.dualcache(zeros(5, 5))))
solve(prob, TRBDF2())
```

## Sparse Jacobians

#### I get errors when I try to solve my problem using sparse Jacobians

This is likely because you're using a Jacobian matrix with a sparsity structure that changes, which is incompatible with the default linear solver for sparse matrices.  If the linear solver catches the issue, you'll see the error message

```
ERROR: ArgumentError: The pattern of the original matrix must match the pattern of the refactor.
```

or

```
ERROR: ArgumentError: pattern of the matrix changed
```

though, an `Error: SingularException` is also possible if the linear solver fails to detect that the sparsity structure changed. To address this issue, you'll need to disable caching the symbolic factorization, e.g.,

```julia
solve(prob, Rodas4(linsolve = KLUFactorization(; reuse_symbolic = false)))
```

For more details about possible linear solvers, consult the [LinearSolve.jl documentation](https://docs.sciml.ai/LinearSolve/stable/)

## Odd Error Messages

#### “Error Exception: `llvmcall` must be compiled to be called” when running the debugger?

The debugger is incompatible with `llvmcall` which is used in the `AutoSpecialize` form
that is used to reduce the compile times.
In order to make use of the debugger, make use of the `FullSpecialize` form.
I.e., change `prob = ODEProblem(lorenz!,u0,tspan)`
to `prob = ODEProblem{true, SciMLBase.FullSpecialize}(lorenz!,u0,tspan)`.
We plan to have a fix for this, but for now,
the workaround should be sufficient for all cases.
