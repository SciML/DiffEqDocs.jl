var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#DifferentialEquations.jl-Documentation-1",
    "page": "Home",
    "title": "DifferentialEquations.jl Documentation",
    "category": "section",
    "text": "This is a suite for numerically solving differential equations in Julia. The purpose of this package is to supply efficient Julia implementations of solvers for various differential equations. Equations within the realm of this package include:Discrete equations (function maps, discrete stochastic (Gillespie/Markov) simulations)\nOrdinary differential equations (ODEs)\nSplit and Partitioned ODEs (Symplectic integrators, IMEX Methods)\nStochastic ordinary differential equations (SODEs or SDEs)\nRandom differential equations (RODEs or RDEs)\nDifferential algebraic equations (DAEs)\nDelay differential equations (DDEs)\nMixed discrete and continuous equations (Hybrid Equations, Jump Diffusions)\n(Stochastic) partial differential equations ((S)PDEs) (with both finite difference and finite element methods)The well-optimized DifferentialEquations solvers benchmark as the some of the fastest implementations, using classic algorithms and ones from recent research which routinely outperform the \"standard\" C/Fortran methods, and include algorithms optimized for high-precision and HPC applications. At the same time, it wraps the classic C/Fortran methods, making it easy to switch over to them whenever necessary. Solving differential equations with different methods from  different languages and packages can be done by changing one line of code,  allowing for easy benchmarking to ensure you are using the fastest method possible.DifferentialEquations.jl integrates with the Julia package sphere with: GPU accleration through CUDAnative.jl and CuArrays.jl\nAutomated sparsity detection with SparsityDetection.jl\nAutomatic Jacobian coloring with SparseDiffTools.jl, allowing for fast solutions to problems with sparse or structured (Tridiagonal, Banded, BlockBanded, etc.) Jacobians\nAllowing the specification of linear solvers for maximal efficiency\nProgress meter integration with the Juno IDE for estimated time to solution \nAutomatic plotting of time series and phase plots \nBuilt-in interpolations\nWraps for common C/Fortran methods like Sundials and Hairer\'s radau\nArbitrary precision with BigFloats and Arbfloats\nArbitrary array types, allowing the definition of differential equations on  matrices and distributed arrays\nUnit checked arithmetic with UnitfulAdditionally, DifferentialEquations.jl comes with built-in analysis features, including:Forward and adjoint local sensitivity analysis for fast gradient computations\nOptimization-based and Bayesian parameter estimation\nNeural differential equations with DiffEqFlux.jl for efficient scientific machine learning (scientific ML) and scientific AI.\nAutomatic distributed, multithreaded, and GPU parallelism of ensemble trajectories\nGlobal sensitivity analysis\nUncertainty quantificationIf you have any questions, or just want to chat about solvers/using the package, please feel free to use the Gitter channel. For bug reports, feature requests, etc., please submit an issue. If you\'re interested in contributing, please see the Developer Documentation."
},

{
    "location": "index.html#Supporting-and-Citing-1",
    "page": "Home",
    "title": "Supporting and Citing",
    "category": "section",
    "text": "The software in this ecosystem was developed as part of academic research. If you would like to help support it, please star the repository as such metrics may help us secure funding in the future. If you use JuliaDiffEq software as part of your research, teaching, or other activities, we would be grateful if you could cite our work. Please see our citation page for guidelines."
},

{
    "location": "index.html#Getting-Started:-Installation-And-First-Steps-1",
    "page": "Home",
    "title": "Getting Started: Installation And First Steps",
    "category": "section",
    "text": "To install the package, use the following command inside the Julia REPL:using Pkg\nPkg.add(\"DifferentialEquations\")To load the package, use the command:using DifferentialEquationsThis will add solvers and dependencies for all kinds of Differential Equations (e.g. ODEs or SDEs etc., see the Supported Equations section below). If you are interested in only one type of equation solvers of DifferentialEquations.jl or simply want a more lightweight version, see the Low Dependency Usage page.To understand the package in more detail, check out the following tutorials in this manual. It is highly recommended that new users start with the ODE tutorial. Example IJulia notebooks can also be found in DiffEqTutorials.jl. If you find any example where there seems to be an error, please open an issue.For the most up to date information on using the package, please join the Gitter channel.Using the bleeding edge for the latest features and development is only recommended for power users. Information on how to get to the bleeding edge is found in the developer documentation."
},

{
    "location": "index.html#IJulia-Notebook-Tutorials-1",
    "page": "Home",
    "title": "IJulia Notebook Tutorials",
    "category": "section",
    "text": "You can access extra tutorials supplied in the DiffEqTutorials.jl repository via the commands:using Pkg\npkg\"add https://github.com/JuliaDiffEq/DiffEqTutorials.jl\"\nusing DiffEqTutorials\nDiffEqTutorials.open_notebooks()Or you can view the webpages for the rendered tutorials at the links found in the repository."
},

{
    "location": "index.html#Video-Tutorial-1",
    "page": "Home",
    "title": "Video Tutorial",
    "category": "section",
    "text": "(Image: Video Tutorial)"
},

{
    "location": "index.html#Tutorials-1",
    "page": "Home",
    "title": "Tutorials",
    "category": "section",
    "text": "The following tutorials will introduce you to the functionality of DifferentialEquations.jl. More examples can be found by checking out the IJulia notebooks in the examples folder.Pages = [\n    \"tutorials/ode_example.md\",\n    \"tutorials/sde_example.md\",\n    \"tutorials/dde_example.md\",\n    \"tutorials/dae_example.md\",\n    \"tutorials/discrete_stochastic_example.md\",\n    \"tutorials/jump_diffusion.md\",\n    \"tutorials/bvp_example.md\",\n    \"tutorials/additional.md\"\n    ]\nDepth = 2"
},

{
    "location": "index.html#Basics-1",
    "page": "Home",
    "title": "Basics",
    "category": "section",
    "text": "These pages introduce you to the core of DifferentialEquations.jl and the common interface. It explains the general workflow, options which are generally available, and the general tools for analysis.Pages = [\n    \"basics/overview.md\",\n    \"basics/common_solver_opts.md\",\n    \"basics/solution.md\",\n    \"basics/plot.md\",\n    \"basics/integrator.md\",\n    \"basics/problem.md\",\n    \"basics/faq.md\",\n    \"basics/compatibility_chart.md\"\n    ]\nDepth = 2"
},

{
    "location": "index.html#Problem-Types-1",
    "page": "Home",
    "title": "Problem Types",
    "category": "section",
    "text": "These pages describe building the problem types to define differential equations for the solvers, and the special features of the different solution types.Pages = [\n  \"types/discrete_types.md\",\n  \"types/ode_types.md\",\n  \"types/dynamical_types.md\",\n  \"types/split_ode_types.md\",\n  \"types/steady_state_types.md\",\n  \"types/bvp_types.md\",\n  \"types/sde_types.md\",\n  \"types/rode_types.md\",\n  \"types/dde_types.md\",\n  \"types/dae_types.md\",\n  \"types/jump_types.md\",\n]\nDepth = 2"
},

{
    "location": "index.html#Solver-Algorithms-1",
    "page": "Home",
    "title": "Solver Algorithms",
    "category": "section",
    "text": "These pages describe the solvers and available algorithms in detail.Pages = [\n  \"solvers/discrete_solve.md\",\n  \"solvers/ode_solve.md\",\n  \"solvers/dynamical_solve.md\",\n  \"solvers/split_ode_solve.md\",\n  \"solvers/steady_state_solve.md\",\n  \"solvers/bvp_solve.md\",\n  \"solvers/jump_solve.md\",\n  \"solvers/sde_solve.md\",\n  \"solvers/rode_solve.md\",\n  \"solvers/dde_solve.md\",\n  \"solvers/dae_solve.md\",\n  \"solvers/benchmarks.md\"\n]\nDepth = 2"
},

{
    "location": "index.html#Additional-Features-1",
    "page": "Home",
    "title": "Additional Features",
    "category": "section",
    "text": "These sections discuss extra performance enhancements, event handling, and other in-depth features.Pages = [\n    \"features/performance_overloads.md\",\n    \"features/diffeq_arrays.md\",\n    \"features/diffeq_operator.md\",\n    \"features/noise_process.md\",\n    \"features/linear_nonlinear.md\",\n    \"features/callback_functions.md\",\n    \"features/callback_library.md\",\n    \"features/ensemble.md\",\n    \"features/io.md\",\n    \"features/low_dep.md\",\n    \"features/progress_bar.md\"\n]\nDepth = 2"
},

{
    "location": "index.html#Analysis-Tools-1",
    "page": "Home",
    "title": "Analysis Tools",
    "category": "section",
    "text": "Because DifferentialEquations.jl has a common interface on the solutions, it is easy to add functionality to the entire DiffEq ecosystem by developing it to the solution interface. These pages describe the add-on analysis tools which are available.Pages = [\n    \"analysis/parameterized_functions.md\",\n    \"analysis/parameter_estimation.md\",\n    \"analysis/bifurcation.md\",\n    \"analysis/sensitivity.md\",\n    \"analysis/global_sensitivity.md\",\n    \"analysis/uncertainty_quantification.md\",\n    \"analysis/neural_networks.md\",\n    \"analysis/dev_and_test.md\"\n]\nDepth = 2"
},

{
    "location": "index.html#Modeling-Tools-1",
    "page": "Home",
    "title": "Modeling Tools",
    "category": "section",
    "text": "While DifferentialEquations.jl can be used to directly build any differential or difference equation (/ discrete stochastic) model, in many cases it can be helpful to have a tailored-built API for making certain types of common models easier. This is provided by the modeling functionality.Pages = [\n    \"models/multiscale.md\",\n    \"models/physical.md\",\n    \"models/financial.md\",\n    \"models/biological.md\",\n    \"models/external_modeling.md\"\n]\nDepth = 2"
},

{
    "location": "index.html#APIs-1",
    "page": "Home",
    "title": "APIs",
    "category": "section",
    "text": "Some DifferentialEquations.jl packages provide documented APIs, these include:Pages = [\n    \"apis/diffeqbio.md\"\n]\nDepth = 2"
},

{
    "location": "index.html#Extra-Details-1",
    "page": "Home",
    "title": "Extra Details",
    "category": "section",
    "text": "These are just assorted extra explanations for the curious.Pages = [\n    \"extras/timestepping.md\"\n]\nDepth = 2"
},

{
    "location": "index.html#Acknowledgements-1",
    "page": "Home",
    "title": "Acknowledgements",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Core-Contributors-1",
    "page": "Home",
    "title": "Core Contributors",
    "category": "section",
    "text": "JuliaDiffEq and DifferentialEquations.jl has been a collaborative effort by many individuals. Significant contributions have been made by the following individuals:Chris Rackauckas (@ChrisRackauckas) (lead developer)\nYingbo Ma (@YingboMa)\nDavid Widmann (@devmotion)\nHendrik Ranocha (@ranocha)\nEthan Levien (@elevien)\nTom Short (@tshort)\n@dextorious\nSamuel Isaacson (@isaacsas)"
},

{
    "location": "index.html#Google-Summer-of-Code-Alumni-1",
    "page": "Home",
    "title": "Google Summer of Code Alumni",
    "category": "section",
    "text": "Yingbo Ma (@YingboMa)\nShivin Srivastava (@shivin9)\nAyush Pandey (@Ayush-iitkgp)\nXingjian Guo (@MSeeker1340)\nShubham Maddhashiya (@sipah00)\nVaibhav Kumar Dixit (@Vaibhavdixit02)"
},

{
    "location": "tutorials/ode_example.html#",
    "page": "Ordinary Differential Equations",
    "title": "Ordinary Differential Equations",
    "category": "page",
    "text": ""
},

{
    "location": "tutorials/ode_example.html#Ordinary-Differential-Equations-1",
    "page": "Ordinary Differential Equations",
    "title": "Ordinary Differential Equations",
    "category": "section",
    "text": "This tutorial will introduce you to the functionality for solving ODEs. Other introductions can be found by checking out DiffEqTutorials.jl. Additionally, a video tutorial walks through this material."
},

{
    "location": "tutorials/ode_example.html#Example-1-:-Solving-Scalar-Equations-1",
    "page": "Ordinary Differential Equations",
    "title": "Example 1 : Solving Scalar Equations",
    "category": "section",
    "text": "In this example we will solve the equationfracdudt = f(upt)on the time interval tin01 where f(upt)=αu. We know by Calculus that the solution to this equation is u(t)=u₀exp(αt).The general workflow is to define a problem, solve the problem, and then analyze the solution. The full code for solving this problem is:using DifferentialEquations\nf(u,p,t) = 1.01*u\nu0=1/2\ntspan = (0.0,1.0)\nprob = ODEProblem(f,u0,tspan)\nsol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)\nusing Plots\nplot(sol,linewidth=5,title=\"Solution to the linear ODE with a thick line\",\n     xaxis=\"Time (t)\",yaxis=\"u(t) (in μm)\",label=\"My Thick Line!\") # legend=false\nplot!(sol.t, t->0.5*exp(1.01t),lw=3,ls=:dash,label=\"True Solution!\")where the pieces are described below."
},

{
    "location": "tutorials/ode_example.html#Step-1:-Defining-a-Problem-1",
    "page": "Ordinary Differential Equations",
    "title": "Step 1: Defining a Problem",
    "category": "section",
    "text": "To solve this numerically, we define a problem type by giving it the equation, the initial condition, and the timespan to solve over:using DifferentialEquations\nf(u,p,t) = 1.01*u\nu0=1/2\ntspan = (0.0,1.0)\nprob = ODEProblem(f,u0,tspan)Note that DifferentialEquations.jl will choose the types for the problem based on the types used to define the problem type. For our example, notice that u0 is a Float64, and therefore this will solve with the dependent variables being Float64. Since tspan = (0.0,1.0) is a tuple of Float64\'s, the independent variables will be solved using Float64\'s (note that the start time and end time must match types). You can use this to choose to solve with arbitrary precision numbers, unitful numbers, etc. Please see the notebook tutorials for more examples.The problem types include many other features, including the ability to define mass matrices and hold callbacks for events. Each problem type has a page which details its constructor and the available fields. For ODEs, the appropriate page is here. In addition, a user can specify additional functions to be associated with the function in order to speed up the solvers. These are detailed at the performance overloads page."
},

{
    "location": "tutorials/ode_example.html#Step-2:-Solving-a-Problem-1",
    "page": "Ordinary Differential Equations",
    "title": "Step 2: Solving a Problem",
    "category": "section",
    "text": ""
},

{
    "location": "tutorials/ode_example.html#Controlling-the-Solvers-1",
    "page": "Ordinary Differential Equations",
    "title": "Controlling the Solvers",
    "category": "section",
    "text": "After defining a problem, you solve it using solve.sol = solve(prob)The solvers can be controlled using the available options are described on the Common Solver Options manual page. For example, we can lower the relative tolerance (in order to get a more correct result, at the cost of more timesteps) by using the command reltol:sol = solve(prob,reltol=1e-6)There are many controls for handling outputs. For example, we can choose to have the solver save every 0.1 time points by setting saveat=0.1. Chaining this with the tolerance choice looks like:sol = solve(prob,reltol=1e-6,saveat=0.1)More generally, saveat can be any collection of time points to save at. Note that this uses interpolations to keep the timestep unconstrained to speed up the solution. In addition, if we only care about the endpoint, we can turn off intermediate saving in general:sol = solve(prob,reltol=1e-6,save_everystep=false)which will only save the final time point."
},

{
    "location": "tutorials/ode_example.html#Choosing-a-Solver-Algorithm-1",
    "page": "Ordinary Differential Equations",
    "title": "Choosing a Solver Algorithm",
    "category": "section",
    "text": "DifferentialEquations.jl has a method for choosing the default solver algorithm which will find an efficient method to solve your problem. To help users receive the right algorithm, DifferentialEquations.jl offers a method for choosing algorithms through hints. This default chooser utilizes the precisions of the number types and the keyword arguments (such as the tolerances) to select an algorithm. Additionally one can provide alg_hints to help choose good defaults using properties of the problem and necessary features for the solution. For example, if we have a stiff problem where we need high accuracy, but don\'t know the best stiff algorithm for this problem, we can use:sol = solve(prob,alg_hints=[:stiff],reltol=1e-8,abstol=1e-8)You can also explicitly choose the algorithm to use. DifferentialEquations.jl offers a much wider variety of solver algorithms than traditional differential equations libraries. Many of these algorithms are from recent research and have been shown to be more efficient than the \"standard\" algorithms. For example, we can choose a 5th order Tsitouras method:sol = solve(prob,Tsit5())Note that the solver controls can be combined with the algorithm choice. Thus we can for example solve the problem using Tsit5() with a lower tolerance via:sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)In DifferentialEquations.jl, some good \"go-to\" choices for ODEs are:AutoTsit5(Rosenbrock23()) handles both stiff and non-stiff equations. This is a good algorithm to use if you know nothing about the equation.\nBS3() for fast low accuracy non-stiff.\nTsit5() for standard non-stiff. This is the first algorithm to try in most cases.\nVern7() for high accuracy non-stiff.\nRodas4() for stiff equations with Julia-defined types, events, etc.\nradau() for really high accuracy stiff equations (requires installing ODEInterfaceDiffEq.jl)For a comprehensive list of the available algorithms and detailed recommendations, Please see the solver documentation. Every problem type has an associated page detailing all of the solvers associated with the problem."
},

{
    "location": "tutorials/ode_example.html#Step-3:-Analyzing-the-Solution-1",
    "page": "Ordinary Differential Equations",
    "title": "Step 3: Analyzing the Solution",
    "category": "section",
    "text": ""
},

{
    "location": "tutorials/ode_example.html#Handling-the-Solution-Type-1",
    "page": "Ordinary Differential Equations",
    "title": "Handling the Solution Type",
    "category": "section",
    "text": "The result of solve is a solution object. We can access the 5th value of the solution with:sol[5] #.637or get the time of the 8th timestep by:sol.t[8] #.438Convenience features are also included. We can build an array using a comprehension over the solution tuples via:[t+u for (u,t) in tuples(sol)]or more generally[t+2u for (u,t) in zip(sol.u,sol.t)]allows one to use more parts of the solution type. The object that is returned by default acts as a continuous solution via an interpolation. We can access the interpolated values by treating sol as a function, for example:sol(0.45) # The value of the solution at t=0.45Note the difference between these: indexing with [i] is the value at the ith step, while (t) is an interpolation at time t!If in the solver dense=true (this is the default unless saveat is used), then this interpolation is a high order interpolation and thus usually matches the error of the solution time points. The interpolations associated with each solver is detailed at the solver algorithm page. If dense=false (unless specifically set, this only occurs when save_everystep=false or saveat is used) then this defaults to giving a linear interpolation.For details on more handling the output, see the solution handling page."
},

{
    "location": "tutorials/ode_example.html#Plotting-Solutions-1",
    "page": "Ordinary Differential Equations",
    "title": "Plotting Solutions",
    "category": "section",
    "text": "While one can directly plot solution time points using the tools given above, convenience commands are defined by recipes for Plots.jl. To plot the solution object, simply call plot:#]add Plots # You need to install Plots.jl before your first time using it!\nusing Plots\n#plotly() # You can optionally choose a plotting backend\nplot(sol)(Image: ode_tutorial_linear_plot)If you are in Juno, this will plot to the plot pane. To open an interactive GUI (dependent on the backend), use the gui command:gui()The plot function can be formatted using the attributes available in Plots.jl. Additional DiffEq-specific controls are documented at the plotting page.For example, from the Plots.jl attribute page we see that the line width can be set via the argument linewidth. Additionally, a title can be set with title. Thus we add these to our plot command to get the correct output, fix up some axis labels, and change the legend (note we can disable the legend with legend=false) to get a nice looking plot:plot(sol,linewidth=5,title=\"Solution to the linear ODE with a thick line\",\n     xaxis=\"Time (t)\",yaxis=\"u(t) (in μm)\",label=\"My Thick Line!\") # legend=falseWe can then add to the plot using the plot! command:plot!(sol.t,t->0.5*exp(1.01t),lw=3,ls=:dash,label=\"True Solution!\")(Image: ode_tutorial_thick_linear)"
},

{
    "location": "tutorials/ode_example.html#Example-2:-Solving-Systems-of-Equations-1",
    "page": "Ordinary Differential Equations",
    "title": "Example 2: Solving Systems of Equations",
    "category": "section",
    "text": "In this example we will solve the Lorenz equations:beginalign\nfracdxdt = σ(y-x) \nfracdydt = x(ρ-z) - y \nfracdzdt = xy - βz \nendalignDefining your ODE function to be in-place updating can have performance benefits. What this means is that, instead of writing a function which outputs its solution, you write a function which updates a vector that is designated to hold the solution. By doing this, DifferentialEquations.jl\'s solver packages are able to reduce the amount of array allocations and achieve better performance.The way we do this is we simply write the output to the 1st input of the function. For example, our Lorenz equation problem would be defined by the function:function lorenz(du,u,p,t)\n du[1] = 10.0*(u[2]-u[1])\n du[2] = u[1]*(28.0-u[3]) - u[2]\n du[3] = u[1]*u[2] - (8/3)*u[3]\nendand then we can use this function in a problem:u0 = [1.0;0.0;0.0]\ntspan = (0.0,100.0)\nprob = ODEProblem(lorenz,u0,tspan)\nsol = solve(prob)Using the plot recipe tools defined on the plotting page, we can choose to do a 3D phase space plot between the different variables:plot(sol,vars=(1,2,3))(Image: Lorenz System)Note that the default plot for multi-dimensional systems is an overlay of each timeseries. We can plot the timeseries of just the second component using the variable choices interface once more:plot(sol,vars=(0,2))(Image: Lorenz Timeseries)Note that here \"variable 0\" corresponds to the independent variable (\"time\")."
},

{
    "location": "tutorials/ode_example.html#Defining-Parameterized-Functions-1",
    "page": "Ordinary Differential Equations",
    "title": "Defining Parameterized Functions",
    "category": "section",
    "text": "In many cases you may want to explicitly have parameters associated with your differential equations. This can be used by things like parameter estimation routines. In this case, you use the p values via the syntax:function parameterized_lorenz(du,u,p,t)\n du[1] = p[1]*(u[2]-u[1])\n du[2] = u[1]*(p[2]-u[3]) - u[2]\n du[3] = u[1]*u[2] - p[3]*u[3]\nendand then we add the parameters to the ODEProblem:u0 = [1.0,0.0,0.0]\ntspan = (0.0,1.0)\np = [10.0,28.0,8/3]\nprob = ODEProblem(parameterized_lorenz,u0,tspan,p)We can make our functions look nicer by doing a few tricks. For example:function parameterized_lorenz(du,u,p,t)\n  x,y,z = u\n  σ,ρ,β = p\n  du[1] = dx = σ*(y-x)\n  du[2] = dy = x*(ρ-z) - y\n  du[3] = dz = x*y - β*z\nendNote that the type for the parameters p can be anything: you can use arrays, static arrays, named tuples, etc. to enclose your parameters in a way that is sensible for your problem.Additionally, there exists a @ode_def macro allows for \"defining your ODE in pseudocode\" and getting a function which is efficient and runnable. To use the macro, you write out your system of equations with the left-hand side being d_ and those variables will be parsed as the dependent variables. The independent variable is t, and the other variables are parameters which you pass at the end. For example, we can write the Lorenz system as:#]add ParameterizedFunctions\nusing ParameterizedFunctions\ng = @ode_def begin\n  dx = σ*(y-x)\n  dy = x*(ρ-z) - y\n  dz = x*y - β*z\nend σ ρ βDifferentialEquations.jl will automatically translate this to be exactly the same as f. The result is more legible code with no performance loss. For more information on the macro Domain Specific Language (DSL) and its limitations, please see the parameterized function page The result is that g is a function which you can now use to define the Lorenz problem.u0 = [1.0;0.0;0.0]\ntspan = (0.0,1.0)\np = [10.0,28.0,8/3]\nprob = ODEProblem(g,u0,tspan,p)The macro does \"behind-the-scenes\" symbolic calculations to pre-compute things like the Jacobian, inverse Jacobian, etc. in order to speed up calculations. Thus not only will this lead to legible ODE definitions, but \"unfairly fast\" code! We can turn off some of the calculations by using a more specific macro, like @ode_def_bare. See ParameterizedFunctions.jl for more details.Since the parameters exist within the function, functions defined in this manner can also be used for sensitivity analysis, parameter estimation routines, and bifurcation plotting. This makes DifferentialEquations.jl a full-stop solution for differential equation analysis which also achieves high performance."
},

{
    "location": "tutorials/ode_example.html#Example-3:-Using-Other-Types-for-Systems-of-Equations-1",
    "page": "Ordinary Differential Equations",
    "title": "Example 3: Using Other Types for Systems of Equations",
    "category": "section",
    "text": "DifferentialEquations.jl can handle many different dependent variable types (generally, anything with a linear index should work!). So instead of solving a vector equation, let\'s let u be a matrix! To do this, we simply need to have u0 be a matrix, and define f such that it takes in a matrix and outputs a matrix. We can define a matrix of linear ODEs as follows:A  = [1. 0  0 -5\n      4 -2  4 -3\n     -4  0  0  1\n      5 -2  2  3]\nu0 = rand(4,2)\ntspan = (0.0,1.0)\nf(u,p,t) = A*u\nprob = ODEProblem(f,u0,tspan)Here our ODE is on a 4x2 matrix, and the ODE is the linear system defined by multiplication by A. To solve the ODE, we do the same steps as before.sol = solve(prob)\nplot(sol)(Image: ODE System Solution)We can instead use the in-place form by using Julia\'s in-place matrix multiplication function mul!:using LinearAlgebra\nf(du,u,p,t) = mul!(du,A,u)Additionally, we can use non-traditional array types as well. For example, StaticArrays.jl offers immutable arrays which are stack-allocated, meaning that their usage does not require any (slow) heap-allocations that arrays normally have. This means that they can be used to solve the same problem as above, with the only change being the type for the initial condition and constants:using StaticArrays, DifferentialEquations\nA  = @SMatrix [ 1.0  0.0 0.0 -5.0\n                4.0 -2.0 4.0 -3.0\n               -4.0  0.0 0.0  1.0\n                5.0 -2.0 2.0  3.0]\nu0 = @SMatrix rand(4,2)\ntspan = (0.0,1.0)\nf(u,p,t) = A*u\nprob = ODEProblem(f,u0,tspan)\nsol = solve(prob)\nusing Plots; plot(sol)Note that the analysis tools generalize over to systems of equations as well.sol[4]still returns the solution at the fourth timestep. It also indexes into the array as well. The last value is the timestep, and the beginning values are for the component. This meanssol[5,3]is the value of the 5th component (by linear indexing) at the 3rd timepoint, orsol[2,1,:]is the timeseries for the component which is the 2nd row and 1 column."
},

{
    "location": "tutorials/ode_example.html#Going-Beyond-ODEs:-How-to-Use-the-Documentation-1",
    "page": "Ordinary Differential Equations",
    "title": "Going Beyond ODEs: How to Use the Documentation",
    "category": "section",
    "text": "Not everything can be covered in the tutorials. Instead, this tutorial will end by pointing you in the directions for the next steps."
},

{
    "location": "tutorials/ode_example.html#Common-API-for-Defining,-Solving,-and-Plotting-1",
    "page": "Ordinary Differential Equations",
    "title": "Common API for Defining, Solving, and Plotting",
    "category": "section",
    "text": "One feature of DifferentialEquations.jl is that this pattern for solving equations is conserved across the different types of differential equations. Every equation has a problem type, a solution type, and the same solution handling (+ plotting) setup. Thus the solver and plotting commands in the Basics section applies to all sorts of equations, like stochastic differential equations and delay differential equations. Each of these different problem types are defined in the Problem Types section of the docs. Every associated solver algorithm is detailed in the Solver Algorithms section, sorted by problem type. The same steps for ODEs can then be used for the analysis of the solution."
},

{
    "location": "tutorials/ode_example.html#Additional-Features-and-Analysis-Tools-1",
    "page": "Ordinary Differential Equations",
    "title": "Additional Features and Analysis Tools",
    "category": "section",
    "text": "In many cases, the common workflow only starts with solving the differential equation. Many common setups have built-in solutions in DifferentialEquations.jl. For example, check out the features for:Handling, parallelizing, and analyzing large Monte Carlo experiments\nSaving the output to tabular formats like DataFrames and CSVs\nEvent handling\nParameter estimation (inverse problems)\nQuantification of numerical uncertainty and errorMany more are defined in the relevant sections of the docs. Please explore the rest of the documentation, including tutorials for getting started with other types of equations. In addition, to get help, please either file an issue at the main repository or come have an informal discussion at our Gitter chatroom."
},

{
    "location": "tutorials/sde_example.html#",
    "page": "Stochastic Differential Equations",
    "title": "Stochastic Differential Equations",
    "category": "page",
    "text": ""
},

{
    "location": "tutorials/sde_example.html#Stochastic-Differential-Equations-1",
    "page": "Stochastic Differential Equations",
    "title": "Stochastic Differential Equations",
    "category": "section",
    "text": "This tutorial will introduce you to the functionality for solving SDEs. Other introductions can be found by checking out DiffEqTutorials.jl. This tutorial assumes you have read the Ordinary Differential Equations tutorial."
},

{
    "location": "tutorials/sde_example.html#Example-1:-Scalar-SDEs-1",
    "page": "Stochastic Differential Equations",
    "title": "Example 1: Scalar SDEs",
    "category": "section",
    "text": "In this example we will solve the equationdu = f(upt)dt + g(upt)dWwhere f(upt)=αu and g(upt)=βu. We know via Stochastic Calculus that the solution to this equation isu(tWₜ)=u₀exp((α-fracβ^22)t+βWₜ)To solve this numerically, we define a problem type by giving it the equation and the initial condition:using DifferentialEquations\nα=1\nβ=1\nu₀=1/2\nf(u,p,t) = α*u\ng(u,p,t) = β*u\ndt = 1//2^(4)\ntspan = (0.0,1.0)\nprob = SDEProblem(f,g,u₀,(0.0,1.0))The solve interface is then the same as with ODEs. Here we will use the classic Euler-Maruyama algorithm EM and plot the solution:sol = solve(prob,EM(),dt=dt)\nusing Plots; plotly() # Using the Plotly backend\nplot(sol)(Image: Basic Solution)"
},

{
    "location": "tutorials/sde_example.html#Using-Higher-Order-Methods-1",
    "page": "Stochastic Differential Equations",
    "title": "Using Higher Order Methods",
    "category": "section",
    "text": "One unique feature of DifferentialEquations.jl is that higher-order methods for stochastic differential equations are included. For reference, let\'s also give the SDEProblem the analytical solution. We can do this by making a test problem. This can be a good way to judge how accurate the algorithms are, or is used to test convergence of the algorithms for methods developers. Thus we define the problem object with:f_analytic(u₀,p,t,W) = u₀*exp((α-(β^2)/2)*t+β*W)\nff = SDEFunction(f,g,analytic=f_analytic)\nprob = SDEProblem(ff,g,u₀,(0.0,1.0))and then we pass this information to the solver and plot:#We can plot using the classic Euler-Maruyama algorithm as follows:\nsol = solve(prob,EM(),dt=dt)\nplot(sol,plot_analytic=true)(Image: SDE Solution)We can choose a higher-order solver for a more accurate result:sol = solve(prob,SRIW1(),dt=dt,adaptive=false)\nplot(sol,plot_analytic=true)(Image: Better SDE Solution)By default, the higher order methods have adaptivity. Thus one can usesol = solve(prob,SRIW1())\nplot(sol,plot_analytic=true)(Image: Better Automatic Solution)Here we allowed the solver to automatically determine a starting dt. This estimate at the beginning is conservative (small) to ensure accuracy. We can instead start the method with a larger dt by passing in a value for the starting dt:sol = solve(prob,SRIW1(),dt=dt)\nplot(sol,plot_analytic=true)(Image: Better Automatic Solution)"
},

{
    "location": "tutorials/sde_example.html#Ensemble-Simulations-1",
    "page": "Stochastic Differential Equations",
    "title": "Ensemble Simulations",
    "category": "section",
    "text": "Instead of solving single trajectories, we can turn our problem into a EnsembleProblem to solve many trajectories all at once. This is done by the EnsembleProblem constructor:ensembleprob = EnsembleProblem(prob)The solver commands are defined at the Parallel Ensemble Simulations page. For example we can choose to have 1000 trajectories via num_monte=1000. In addition, this will automatically parallelize using Julia native parallelism if extra processes are added via addprocs(), but we can change this to use multithreading via EnsembleThreads(). Together, this looks like:sol = solve(ensembleprob,EnsembleThreads(),trajectories=1000)Many more controls are defined at the Monte Carlo page, including analysis tools. A very simple analysis can be done with the MonteCarloSummary, which builds mean/var statistics and has an associated plot recipe. For example, we can get the statistics at every 0.01 timesteps and plot the average + error using:summ = EnsembleSummary(sol,0:0.01:1)\nplot(summ,labels=\"Middle 95%\")\nsumm = EnsembleSummary(sol,0:0.01:1;quantiles=[0.25,0.75])\nplot!(summ,labels=\"Middle 50%\",legend=true)(Image: sde_tutorial_monte)Additionally we can easily calculate the correlation between the values at t=0.2 and t=0.7 viatimepoint_meancor(sim,0.2,0.7) # Gives both means and then the correlation coefficient"
},

{
    "location": "tutorials/sde_example.html#Example-2:-Systems-of-SDEs-with-Diagonal-Noise-1",
    "page": "Stochastic Differential Equations",
    "title": "Example 2: Systems of SDEs with Diagonal Noise",
    "category": "section",
    "text": "More generally, an SDEdu = f(upt)dt + g(upt)dWgeneralizes to systems of equations is done in the same way as ODEs. Here, g is now a matrix of values. One common case, and the default for DifferentialEquations.jl, is diagonal noise where g is a diagonal matrix. This means that every function in the system gets a different random number. Instead of handling matrices in this case, we simply define both f and g as in-place functions. Thus f(du,u,p,t) gives a vector of du which is the deterministic change, and g(du2,u,p,t) gives a vector du2 for which du2.*W is the stochastic portion of the equation.For example, the Lorenz equation with additive noise has the same deterministic portion as the Lorenz equations, but adds an additive noise, which is simply 3*N(0,dt) where N is the normal distribution dt is the time step, to each step of the equation. This is done via:function lorenz(du,u,p,t)\n  du[1] = 10.0(u[2]-u[1])\n  du[2] = u[1]*(28.0-u[3]) - u[2]\n  du[3] = u[1]*u[2] - (8/3)*u[3]\nend\n\nfunction σ_lorenz(du,u,p,t)\n  du[1] = 3.0\n  du[2] = 3.0\n  du[3] = 3.0\nend\n\nprob_sde_lorenz = SDEProblem(lorenz,σ_lorenz,[1.0,0.0,0.0],(0.0,10.0))\nsol = solve(prob_sde_lorenz)\nplot(sol,vars=(1,2,3))(Image: stochastic_3d_lorenz)Note that it\'s okay for the noise function to mix terms. For examplefunction σ_lorenz(du,u,p,t)\n  du[1] = sin(u[3])*3.0\n  du[2] = u[2]*u[1]*3.0\n  du[3] = 3.0\nendis a valid noise function, which will once again give diagonal noise by du2.*W. Note also that in this format, it is fine to use ParameterizedFunctions. For example, the Lorenz equation could have been defined as:#]add ParameterizedFunctions\nusing ParameterizedFunctions\n\nf = @ode_def begin\n  dx = σ*(y-x)\n  dy = x*(ρ-z) - y\n  dz = x*y - β*z\nend σ ρ β\n\ng = @ode_def begin\n  dx = α\n  dy = α\n  dz = α\nend α"
},

{
    "location": "tutorials/sde_example.html#Example-3:-Systems-of-SDEs-with-Scalar-Noise-1",
    "page": "Stochastic Differential Equations",
    "title": "Example 3: Systems of SDEs with Scalar Noise",
    "category": "section",
    "text": "In this example we\'ll solve a system of SDEs with scalar noise. This means that the same noise process is applied to all SDEs. First we need to define a scalar noise process using the Noise Process interface. Since we want a WienerProcess that starts at 0.0 at time 0.0, we use the command W = WienerProcess(0.0,0.0,0.0) to define the Brownian motion we want, and then give this to the noise option in the SDEProblem. For a full example, let\'s solve a linear SDE with scalar noise using a high order algorithm:f(du,u,p,t) = (du .= u)\ng(du,u,p,t) = (du .= u)\nu0 = rand(4,2)\n\nW = WienerProcess(0.0,0.0,0.0)\nprob = SDEProblem(f,g,u0,(0.0,1.0),noise=W)\nsol = solve(prob,SRIW1())(Image: Scalar Noise)"
},

{
    "location": "tutorials/sde_example.html#Example-4:-Systems-of-SDEs-with-Non-Diagonal-Noise-1",
    "page": "Stochastic Differential Equations",
    "title": "Example 4: Systems of SDEs with Non-Diagonal Noise",
    "category": "section",
    "text": "In the previous examples we had diagonal noise, that is a vector of random numbers dW whose size matches the output of g where the noise is applied element-wise, and scalar noise where a single random variable is applied to all dependent variables. However, a more general type of noise allows for the terms to linearly mixed via g being a matrix.(Note that nonlinear mixings are not SDEs but fall under the more general class of random ordinary differential equations (RODEs) which have a separate set of solvers.Let\'s define a problem with four Wiener processes and two dependent random variables. In this case, we will want the output of g to be a 2x4 matrix, such that the solution is g(u,p,t)*dW, the matrix multiplication. For example, we can do the following:f(du,u,p,t) = du .= 1.01u\nfunction g(du,u,p,t)\n  du[1,1] = 0.3u[1]\n  du[1,2] = 0.6u[1]\n  du[1,3] = 0.9u[1]\n  du[1,4] = 0.12u[2]\n  du[2,1] = 1.2u[1]\n  du[2,2] = 0.2u[2]\n  du[2,3] = 0.3u[2]\n  du[2,4] = 1.8u[2]\nend\nprob = SDEProblem(f,g,ones(2),(0.0,1.0),noise_rate_prototype=zeros(2,4))In our g we define the functions for computing the values of the matrix. We can now think of the SDE that this solves as the system of equationsdu_1 = f_1(upt)dt + g_11(upt)dW_1 + g_12(upt)dW_2 + g_13(upt)dW_3 + g_14(upt)dW_4 \ndu_2 = f_2(upt)dt + g_21(upt)dW_1 + g_22(upt)dW_2 + g_23(upt)dW_3 + g_24(upt)dW_4meaning that for example du[1,1] and du[2,1] correspond to stochastic changes with the same random number in the first and second SDEs.Note that this problem can only be solved my SDE methods which are compatible with non-diagonal noise. This is discussed in the SDE solvers page.The matrix itself is determined by the keyword argument noise_rate_prototype in the SDEProblem constructor. This is a prototype for the type that du will be in g. This can be any AbstractMatrix type. Thus for example, we can define the problem as\n# Define a sparse matrix by making a dense matrix and setting some values as not zero\nA = zeros(2,4)\nA[1,1] = 1\nA[1,4] = 1\nA[2,4] = 1\nA=sparse(A)\n\n# Make `g` write the sparse matrix values\nfunction g(du,u,p,t)\n  du[1,1] = 0.3u[1]\n  du[1,4] = 0.12u[2]\n  du[2,4] = 1.8u[2]\nend\n\n# Make `g` use the sparse matrix\nprob = SDEProblem(f,g,ones(2),(0.0,1.0),noise_rate_prototype=A)and now g(u,p,t) writes into a sparse matrix, and g(u,p,t)*dW is sparse matrix multiplication."
},

{
    "location": "tutorials/sde_example.html#Example-4:-Colored-Noise-1",
    "page": "Stochastic Differential Equations",
    "title": "Example 4: Colored Noise",
    "category": "section",
    "text": "Colored noise can be defined using the Noise Process interface. In that portion of the docs, it is shown how to define your own noise process my_noise, which can be passed to the SDEProblemSDEProblem(f,g,u0,tspan,noise=my_noise)Note that general colored noise problems are only compatible with the EM and EulerHeun methods. This is discussed in the SDE solvers page."
},

{
    "location": "tutorials/sde_example.html#Example:-Spatially-Colored-Noise-in-the-Heston-Model-1",
    "page": "Stochastic Differential Equations",
    "title": "Example: Spatially-Colored Noise in the Heston Model",
    "category": "section",
    "text": "Let\'s define the Heston equation from financial mathematics:dS = μSdt + sqrtvSdW_1 \ndv = κ(Θ-v)dt + σsqrtvdW_2 \ndW_1 dW_2 = ρ dtIn this problem, we have a diagonal noise problem given by:function f(du,u,p,t)\n  du[1] = μ*u[1]\n  du[2] = κ*(Θ-u[2])\nend\nfunction g(du,u,p,t)\n  du[1] = √u[2]*u[1]\n  du[2] = Θ*√u[2]\nendHowever, our noise has a correlation matrix for some constant ρ. Choosing ρ=0.2:Γ = [1 ρ;ρ 1]To solve this, we can define a CorrelatedWienerProcess which starts at zero (W(0)=0) via:heston_noise = CorrelatedWienerProcess!(Γ,tspan[1],zeros(2),zeros(2))This is then used to build the SDE:SDEProblem(f,g,u0,tspan,noise=heston_noise)Of course, to fully define this problem we need to define our constants. Constructors for making common models like this easier to define can be found in the modeling toolkits. For example, the HestonProblem is pre-defined as part of the financial modeling tools."
},

{
    "location": "tutorials/rode_example.html#",
    "page": "Random Ordinary Differential Equations",
    "title": "Random Ordinary Differential Equations",
    "category": "page",
    "text": ""
},

{
    "location": "tutorials/rode_example.html#Random-Ordinary-Differential-Equations-1",
    "page": "Random Ordinary Differential Equations",
    "title": "Random Ordinary Differential Equations",
    "category": "section",
    "text": "This tutorial will introduce you to the functionality for solving RODEs. Other introductions can be found by checking out DiffEqTutorials.jl."
},

{
    "location": "tutorials/rode_example.html#Example-1:-Scalar-RODEs-1",
    "page": "Random Ordinary Differential Equations",
    "title": "Example 1: Scalar RODEs",
    "category": "section",
    "text": "In this example we will solve the equationdu = f(uptW)dtwhere f(uptW)=2usin(W) and W(t) is a Wiener process (Gaussian process).using DifferentialEquations\nfunction f(u,p,t,W)\n  2u*sin(W)\nend\nu0 = 1.00\ntspan = (0.0,5.0)\nprob = RODEProblem(f,u0,tspan)\nsol = solve(prob,RandomEM(),dt=1/100)(Image: intro_rode)The random process defaults to a Gaussian/Wiener process, so there is nothing else required here! See the documentation on NoiseProcesses for details on how to define other noise proceses."
},

{
    "location": "tutorials/rode_example.html#Example-2:-Systems-of-RODEs-1",
    "page": "Random Ordinary Differential Equations",
    "title": "Example 2: Systems of RODEs",
    "category": "section",
    "text": "As with the other problem types, there is an in-place version which is more efficient for systems. The signature is f(du,u,p,t,W). For example,using DifferentialEquations\nfunction f(du,u,p,t,W)\n  du[1] = 2u[1]*sin(W[1] - W[2])\n  du[2] = -2u[2]*cos(W[1] + W[2])\nend\nu0 = [1.00;1.00]\ntspan = (0.0,5.0)\nprob = RODEProblem(f,u0,tspan)\nsol = solve(prob,RandomEM(),dt=1/100)(Image: rode_system)By default, the size of the noise process matches the size of u0. However, you can use the rand_prototype keyword to explicitly set the size of the random process:function f(du,u,p,t,W)\n  du[1] = -2W[3]*u[1]*sin(W[1] - W[2])\n  du[2] = -2u[2]*cos(W[1] + W[2])\nend\nu0 = [1.00;1.00]\ntspan = (0.0,5.0)\nprob = RODEProblem(f,u0,tspan,rand_prototype=zeros(3))\nsol = solve(prob,RandomEM(),dt=1/100)(Image: noise_choice)"
},

{
    "location": "tutorials/dde_example.html#",
    "page": "Delay Differential Equations",
    "title": "Delay Differential Equations",
    "category": "page",
    "text": ""
},

{
    "location": "tutorials/dde_example.html#Delay-Differential-Equations-1",
    "page": "Delay Differential Equations",
    "title": "Delay Differential Equations",
    "category": "section",
    "text": "This tutorial will introduce you to the functionality for solving delay differential equations. This tutorial assumes you have read the Ordinary Differential Equations tutorial.Delay differential equations are equations which have a delayed argument. To allow for specifying the delayed argument, the function definition for a delay differential equation is expanded to include a history function h(p, t) which uses interpolations throughout the solution\'s history to form a continuous extension of the solver\'s past and depends on parameters p and time t. The function signature for a delay differential equation is f(u, h, p, t) for not in-place computations, and f(du, u, h, p, t) for in-place computations.In this example we will solve a model of breast cancer growth kinetics:beginalign\ndx_0 = fracv_01+beta_0left(x_2(t-tau)right)^2left(p_0-q_0right)x_0(t)-d_0x_0(t)\ndx_1 = fracv_01+beta_0left(x_2(t-tau)right)^2left(1-p_0+q_0right)x_0(t)\n       + fracv_11+beta_1left(x_2(t-tau)right)^2left(p_1-q_1right)x_1(t)-d_1x_1(t) notag \ndx_2 = fracv_11+beta_1left(x_2(t-tau)right)^2left(1-p_1+q_1right)x_1(t)-d_2x_2(t)\nendalignFor this problem we note that tau is constant, and thus we can use a method which exploits this behavior. We first write out the equation using the appropriate function signature. Most of the equation writing is the same, though we use the history function by first interpolating and then choosing the components. Thus the ith component at time t-tau is given by h(p, t-tau)[i]. Components with no delays are written as in the ODE.Thus, the function for this model is given by:const p0 = 0.2; const q0 = 0.3; const v0 = 1; const d0 = 5\nconst p1 = 0.2; const q1 = 0.3; const v1 = 1; const d1 = 1\nconst d2 = 1; const beta0 = 1; const beta1 = 1; const tau = 1\nfunction bc_model(du,u,h,p,t)\n  du[1] = (v0/(1+beta0*(h(p, t-tau)[3]^2))) * (p0 - q0)*u[1] - d0*u[1]\n  du[2] = (v0/(1+beta0*(h(p, t-tau)[3]^2))) * (1 - p0 + q0)*u[1] +\n          (v1/(1+beta1*(h(p, t-tau)[3]^2))) * (p1 - q1)*u[2] - d1*u[2]\n  du[3] = (v1/(1+beta1*(h(p, t-tau)[3]^2))) * (1 - p1 + q1)*u[2] - d2*u[3]\nendTo use the constant lag model, we have to declare the lags. Here we will use tau=1.lags = [tau]Now we build a DDEProblem. The signatureprob = DDEProblem(f, u0, h, tspan, p=nothing;\n                  constant_lags=[], dependent_lags=[], kwargs...)is very similar to ODEs, where we now have to give the lags and a function h. h is the history function that declares what the values were before the time the model starts. Here we will assume that for all time before t0 the values were 1 and define h as an out-of-place function:h(p, t) = ones(3)Next, we choose to solve on the timespan (0.0,10.0) and create the problem type:tspan = (0.0,10.0)\nu0 = [1.0,1.0,1.0]\nprob = DDEProblem(bc_model,u0,h,tspan; constant_lags=lags)An efficient way to solve this problem (given the constant lags) is with the MethodOfSteps solver. Through the magic that is Julia, it translates an OrdinaryDiffEq.jl ODE solver method into a method for delay differential equations which is highly efficient due to sweet compiler magic. A good choice is the order 5 method Tsit5():alg = MethodOfSteps(Tsit5())For lower tolerance solving, one can use the BS3() algorithm to good effect (this combination is similar to the MATLAB dde23, but more efficient tableau), and for high tolerances the Vern6() algorithm will give an 6th order solution.To solve the problem with this algorithm, we do the same thing we\'d do with other methods on the common interface:sol = solve(prob,alg)Note that everything available to OrdinaryDiffEq.jl can be used here, including event handling and other callbacks. The solution object has the same interface as for ODEs. For example, we can use the same plot recipes to view the results:using Plots; plot(sol)(Image: DDE Example Plot)"
},

{
    "location": "tutorials/dde_example.html#Speeding-Up-Interpolations-with-Idxs-1",
    "page": "Delay Differential Equations",
    "title": "Speeding Up Interpolations with Idxs",
    "category": "section",
    "text": "We can speed up the previous problem in two different ways. First of all, if we need to interpolate multiple values from a previous time, we can use the in-place form for the history function h(out, p, t) which writes the output to out. In this case, we must supply the history initial conditions as in-place as well. For the previous example, that\'s simplyh(out, p, t) = (out.=1.0)and then our DDE is:const out = zeros(3) # Define a cache variable\nfunction bc_model(du,u,h,p,t)\n  h(out, p, t-tau) # updates out to be the correct history function\n  du[1] = (v0/(1+beta0*(out[3]^2))) * (p0 - q0)*u[1] - d0*u[1]\n  du[2] = (v0/(1+beta0*(out[3]^2))) * (1 - p0 + q0)*u[1] +\n          (v1/(1+beta1*(out[3]^2))) * (p1 - q1)*u[2] - d1*u[2]\n  du[3] = (v1/(1+beta1*(out[3]^2))) * (1 - p1 + q1)*u[2] - d2*u[3]\nendHowever, we can do something even slicker in most cases. We only ever needed to interpolate past values at index 3. Instead of generating a bunch of arrays, we can instead ask specifically for that value by passing the keyword idxs = 3. The DDE function bc_model is now:function bc_model(du,u,h,p,t)\n  u3_past_sq = h(p, t-tau; idxs=3)^2\n  du[1] = (v0/(1+beta0*(u3_past_sq))) * (p0 - q0)*u[1] - d0*u[1]\n  du[2] = (v0/(1+beta0*(u3_past_sq))) * (1 - p0 + q0)*u[1] +\n          (v1/(1+beta1*(u3_past_sq))) * (p1 - q1)*u[2] - d1*u[2]\n  du[3] = (v1/(1+beta1*(u3_past_sq))) * (1 - p1 + q1)*u[2] - d2*u[3]\nendNote that this requires that we define the historical valuesh(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(3)where idxs can be an integer for which variable in the history to compute, and here for any number idxs we give back 1.0. Note that if we wanted to use past values of the ith derivative then we would call the history function h(p, t, Val{i}) in our DDE function and would have to define a dispatch likeh(p, t, ::Type{Val{1}}) = zeros(3)to say that derivatives before t0 are zero for any index. Again, we could use an in-place function instead or only compute specific indices by passing an idxs keyword.The functional forms for the history function are discussed also on the DDEProblem page."
},

{
    "location": "tutorials/dde_example.html#Undeclared-Delays-and-State-Dependent-Delays-via-Residual-Control-1",
    "page": "Delay Differential Equations",
    "title": "Undeclared Delays and State-Dependent Delays via Residual Control",
    "category": "section",
    "text": "You might have noticed DifferentialEquations.jl allows you to solve problems with undeclared delays since you can interpolate h at any value. This is a feature, but use it with caution. Undeclared delays can increase the error in the solution. It\'s recommended that you use a method with a residual control, such as MethodOfSteps(RK4()) whenever there are undeclared delays. With this you can use interpolated derivatives, solve functional differential equations by using quadrature on the interpolant, etc. However, note that residual control solves with a low level of accuracy, so the tolerances should be made very small and the solution should not be trusted for more than 2-3 decimal places.Note: MethodOfSteps(RK4()) with undeclared delays is similar to MATLAB\'s ddesd. Thus, for example, the following is similar to solving the example from above with residual control:prob = DDEProblem(bc_model,u0,h,tspan)\nalg = MethodOfSteps(RK4())\nsol = solve(prob,alg)Note that this method can solve problems with state-dependent delays."
},

{
    "location": "tutorials/dde_example.html#State-Dependent-Delay-Discontinuity-Tracking-1",
    "page": "Delay Differential Equations",
    "title": "State-Dependent Delay Discontinuity Tracking",
    "category": "section",
    "text": "State-dependent delays are problems where the delay is allowed to be a function of the current state. They can be more efficiently solved with discontinuity tracking. To do this in DifferentialEquations.jl, requires to pass lag functions g(u,p,t) as keyword dependent_lags to the DDEProblem definition. Other than that, everything else is the same, and one solves that problem using the common interface.We can solve the above problem with dependent delay tracking by declaring the dependent lags and solving with a MethodOfSteps algorithm:prob = DDEProblem(bc_model,u0,h,tspan; dependent_lags = ((u,p,t) -> tau,))\nalg = MethodOfSteps(Tsit5())\nsol = solve(prob,alg)Here we treated the single lag t-tau as a state-dependent delay. Of course, you can then replace that tuple of functions with whatever functions match your lags."
},

{
    "location": "tutorials/dae_example.html#",
    "page": "Differential Algebraic Equations",
    "title": "Differential Algebraic Equations",
    "category": "page",
    "text": ""
},

{
    "location": "tutorials/dae_example.html#Differential-Algebraic-Equations-1",
    "page": "Differential Algebraic Equations",
    "title": "Differential Algebraic Equations",
    "category": "section",
    "text": "This tutorial will introduce you to the functionality for solving DAEs. Other introductions can be found by checking out DiffEqTutorials.jl. This tutorial assumes you have read the Ordinary Differential Equations tutorial.In this example we will solve the implicit ODE equationf(duupt) = 0where f is the a variant of the Roberts equation. This equation is of the formbeginalign\ndu = f(upt) \n 0 = g(upt) \n endalignor is also known as a constrained differential equation where g is the constraint equation. The Roberts model can be written in the form:beginalign\ndy_1 = -004y₁ + 10^4 y_2 y_3 \ndy_2 = 004 y_1 - 10^4 y_2 y_3 - 3*10^7 y_2^2 \n1 =  y_1 + y_2 + y_3 \nendalignwith initial conditions y_1(0) = 1, y_2(0) = 0, y_3(0) = 0, dy_1 = - 004, dy_2 = 004, and dy_3 = 00.The workflow for DAEs is the same as for the other types of equations, where all you need to know is how to define the problem. A DAEProblem is specified by defining an in-place update f(out,du,u,p,t) which uses the values to mutate out as the output. To makes this into a DAE, we move all of the variables to one side. Thus we can define the function:function f(out,du,u,p,t)\n  out[1] = - 0.04u[1]              + 1e4*u[2]*u[3] - du[1]\n  out[2] = + 0.04u[1] - 3e7*u[2]^2 - 1e4*u[2]*u[3] - du[2]\n  out[3] = u[1] + u[2] + u[3] - 1.0\nendwith initial conditionsu₀ = [1.0, 0, 0]\ndu₀ = [-0.04, 0.04, 0.0]\ntspan = (0.0,100000.0)and make the DAEProblem:using DifferentialEquations\ndifferential_vars = [true,true,false]\nprob = DAEProblem(f,du₀,u₀,tspan,differential_vars=differential_vars)differential_vars is an option which states which of the variables are differential, i.e. not purely algebraic (which means that their derivative shows up in the residual equations). This is required for the algorithm to be able to find consistant initial conditions. Notice that the first two variables are determined by their changes, but the last is simply determined by the conservation equation. Thus we use differential_vars = [true,true,false].As with the other DifferentialEquations problems, the commands are then to solve and plot. Here we will use the IDA solver from Sundials:using Sundials\nsol = solve(prob,IDA())In order to clearly see all the features of this solution, it should be plotted on a logarithmic scale. We\'ll also plot each on a different subplot to allow scaling the y-axis appropriately.using Plots; plotly() # Using the Plotly backend\nplot(sol, xscale=:log10, tspan=(1e-6, 1e5), layout=(3,1))This gives the following plot:(Image: IntroDAEPlot)"
},

{
    "location": "tutorials/discrete_stochastic_example.html#",
    "page": "Discrete Stochastic (Gillespie) Equations",
    "title": "Discrete Stochastic (Gillespie) Equations",
    "category": "page",
    "text": ""
},

{
    "location": "tutorials/discrete_stochastic_example.html#Discrete-Stochastic-(Gillespie)-Equations-1",
    "page": "Discrete Stochastic (Gillespie) Equations",
    "title": "Discrete Stochastic (Gillespie) Equations",
    "category": "section",
    "text": "In this tutorial we will describe how to define and simulate continuous-time jump processes, also known in biological fields as Gillespie-type models. This tutorial assumes you have read the Ordinary Differential Equations tutorial. The discrete stochastic simulations we consider are a form of jump equation with a \"trivial\" (non-existent) differential equation. We will first demonstrate how to build these types of models using the biological modeling functionality, then describe how to build them directly and more generally using jumps, and finally show how to add discrete stochastic simulations to differential equation models."
},

{
    "location": "tutorials/discrete_stochastic_example.html#Defining-a-Model-using-Reactions-1",
    "page": "Discrete Stochastic (Gillespie) Equations",
    "title": "Defining a Model using Reactions",
    "category": "section",
    "text": "For our example, we will build an SIR model which matches the tutorial from Gillespie.jl. SIR stands for susceptible, infected, and recovered, and is a model of disease spread. When a susceptible person comes in contact with an infected person, the disease has a chance of infecting the susceptible person. This \"chance\" is determined by the number of susceptible persons and the number of infected persons, since in larger populations there is a greater chance that two people come into contact. Normally, the rate is modeled as the amountrate_constant*num_of_susceptible_people*num_of_infected_peopleThe rate_constant is determined by factors like the type of the disease. It can be interpreted as the probability per time one pair of susceptible and infected people encounter each other, with the susceptible person becoming sick. The overall rate (i.e. probability per time) that some susceptible person gets sick is then given by the rate constant multiplied by the number of possible pairs of susceptible and infected people. This formulation is known as the law of mass action.Let s be the number of susceptible persons, i be the number of infected persons, and r be the number of recovered persons. In this case, we can re-write our overall rate as:rate_constant*s*iThus we have that our \"reactants\" are components 1 and 2. When this \"reaction\" occurs, the result is that one susceptible person turns into an infected person. We can think of this as doing:s -= 1\ni += 1that is, we decrease the number of susceptible persons by 1 and increase the number of infected persons by 1.These are the facts the are encoded in the reaction:c1, s + i --> 2iThis \"reaction\" encodes that a susceptible person and an infected person can interact, resulting in two infected persons (i.e. the susceptible person was infected). Here, c1 is the reaction constant.To finish the model, we define one more reaction. Over time, infected people become less infected. The chance that any one person heals during some time unit depends on the number of people who are infected. Thus the rate at which infected people turn into recovered people israte_constant*iWhen this happens, we lose one infected person and gain a recovered person. This reaction is thus modeled as:c2, i --> rThus our full reaction network is:# ]add DiffEqBiological\nusing DiffEqBiological\nsir_model = @reaction_network SIR begin\n    c1, s + i --> 2i\n    c2, i --> r\nend c1 c2Notice that the order the variables are introduced in the model is s, then i, then r, and thus this is the canonical ordering of the variables."
},

{
    "location": "tutorials/discrete_stochastic_example.html#Building-and-Solving-the-Problem-1",
    "page": "Discrete Stochastic (Gillespie) Equations",
    "title": "Building and Solving the Problem",
    "category": "section",
    "text": "First, we have to define some kind of differential equation that we can \"solve\" to simulate the jump process. Since we do not want any continuous changes in the numbers of the different types of people, we will build a DiscreteProblem. We do this by giving the constructor u0, the initial condition, and tspan, the timespan. Here, we will start with 999 susceptible people, 1 infected person, and 0 recovered people, and solve the problem from t=0.0 to t=250.0. We use the parameters c1 = 0.1/1000 and c2 = 0.01. Thus we build the problem via:p = (0.1/1000,0.01)\nprob = DiscreteProblem([999,1,0],(0.0,250.0),p)The reaction network can be converted into various differential equations like JumpProblem, ODEProblem, or an SDEProblem. To turn it into a jump problem, we simply do:jump_prob = JumpProblem(prob,Direct(),sir_model)This is now a problem that can be solved using the differential equations solvers. Since our problem is discrete, we will use the FunctionMap() method.sol = solve(jump_prob,FunctionMap())This solve command takes the standard commands of the common interface, and the solution object acts just like any other differential equation solution. Thus there exists a plot recipe, which we can plot with:using Plots; plot(sol)(Image: gillespie_solution)"
},

{
    "location": "tutorials/discrete_stochastic_example.html#SSAStepper-1",
    "page": "Discrete Stochastic (Gillespie) Equations",
    "title": "SSAStepper",
    "category": "section",
    "text": "The previous example used FunctionMap() to perform the jump process simulation. FunctionMap is a DiscreteProblem algorithm in OrdinaryDiffEq.jl. This shows that any common interface algorithm can be used to perform the timestepping since this is implemented over the callback interface. In many cases we may have a pure jump system that only involves ConstantRateJumps and/or MassActionJumps. When that\'s the case, a substantial performance benefit may be gained by using SSAStepper()sol = solve(jump_prob,SSAStepper())SSAStepper is a barebones SSA method which doesn\'t allow defining events or integrating simultaneous ODEs, but is very efficient for pure jump/SSA problems."
},

{
    "location": "tutorials/discrete_stochastic_example.html#Controlling-Saving-Behavior-1",
    "page": "Discrete Stochastic (Gillespie) Equations",
    "title": "Controlling Saving Behavior",
    "category": "section",
    "text": "Note that jumps act via the callback interface which defaults to saving at each event. The reason is because this is required in order to accurately resolve every discontinuity exactly (and this is what allows for perfectly vertical lines!). However, in many cases when using jump problems you may wish to decrease the saving pressure given by large numbers of jumps. To do this, you set save_positions in the JumpProblem. Just like for other callbacks, this is a tuple (bool1,bool2) which saves whether to save before or after a jump. If we do not want to save at every jump, we would thus pass:jump_prob = JumpProblem(prob,Direct(),sir_model,save_positions=(false,false))Now the saving controls associated with the integrator are the only ones to note. For example, we can use saveat=0.5 to save at an evenly spaced grid:sol = solve(jump_prob,FunctionMap(),saveat=0.5)"
},

{
    "location": "tutorials/discrete_stochastic_example.html#Defining-the-Jumps-Directly:-ConstantRateJump-1",
    "page": "Discrete Stochastic (Gillespie) Equations",
    "title": "Defining the Jumps Directly: ConstantRateJump",
    "category": "section",
    "text": "Instead of using the biological modeling functionality of @reaction_network, we can directly define jumps. This allows for more general types of rates, at the cost of some modeling friendliness. The constructor for a ConstantRateJump is:jump = ConstantRateJump(rate,affect!)where rate is a function rate(u,p,t) and affect! is a function of the integrator affect!(integrator) (for details on the integrator, see the integrator interface docs). Thus, to define the jump equivalents to the above reactions, we can use:rate1(u,p,t) = (0.1/1000.0)*u[1]*u[2]\nfunction affect1!(integrator)\n  integrator.u[1] -= 1\n  integrator.u[2] += 1\nend\njump = ConstantRateJump(rate1,affect1!)\n\nrate2(u,p,t) = 0.01u[2]\nfunction affect2!(integrator)\n  integrator.u[2] -= 1\n  integrator.u[3] += 1\nend\njump2 = ConstantRateJump(rate2,affect2!)We can then use JumpProblem to augment a problem with jumps. To add the jumps to the DiscreteProblem and solve it, we would simply do:jump_prob = JumpProblem(prob,Direct(),jump,jump2)\nsol = solve(jump_prob,FunctionMap())Note, in systems with more than a few jumps (more than ~10), it can be advantageous to use a different internal representation for the jump collection. For such systems it is recommended to use DirectFW(), which should offer better performance than Direct()."
},

{
    "location": "tutorials/discrete_stochastic_example.html#Defining-the-Jumps-Directly:-MassActionJump-1",
    "page": "Discrete Stochastic (Gillespie) Equations",
    "title": "Defining the Jumps Directly: MassActionJump",
    "category": "section",
    "text": "For systems that can be represented as mass action reactions, a further specialization of the jump type is possible that offers improved computational performance; MasssActionJump. Suppose the system has N chemical species S_1dotsS_N. A general mass action reaction has the formR_1 S_1 + R_2 S_2 + dots + R_N S_N oversetkrightarrow P_1 S_1 + P_2 S_2 + dots + P_N S_Nwhere the non-negative integers (R_1dotsR_N) denote the reactant stoichiometry of the reaction, and the non-negative integers (P_1dotsP_N) the product stoichiometry. The net stoichiometry is the net change in each chemical species from the reaction occurring one time, given by (P_1-R_1dotsP_N-R_N).As an example, consider again the SIR model defined in the @reaction_network above. The species are then (s,i,r). The first reaction has rate c1, reactant stoichiometry (1,1,0), product stoichiometry (0,2,0), and net stoichiometry (-1,1,0). The second reaction has rate c2, reactant stoichiometry (0,1,0), product stoichiometry (0,0,1), and net stoichiometry (0,-1,1).We can encode this system as a mass action jump by specifying the rates, reactant stoichiometry, and the net stoichiometry as follows:rates = [0.1/1000, 0.01]    # i.e. [c1,c2]\nreactant_stoich =\n[\n  [1 => 1, 2 => 1],         # 1*s and 1*i\n  [2 => 1]                  # 1*i\n]\nnet_stoich =\n[\n  [1 => -1, 2 => 1],        # -1*s and 1*i\n  [2 => -1, 3 => 1]         # -1*i and 1*r\n]\nmass_act_jump = MassActionJump(rates, reactant_stoich, net_stoich)Just like for ConstantRateJumps, to then simulate the system we create a JumpProblem and call solve:jump_prob = JumpProblem(prob, Direct(), mass_act_jump)\nsol = solve(jump_prob, SSAStepper())"
},

{
    "location": "tutorials/discrete_stochastic_example.html#Defining-the-Jumps-Directly:-Mixing-ConstantRateJump-and-MassActionJump-1",
    "page": "Discrete Stochastic (Gillespie) Equations",
    "title": "Defining the Jumps Directly: Mixing ConstantRateJump and MassActionJump",
    "category": "section",
    "text": "Suppose we now want to add in to the SIR model another jump that can not be represented as a mass action reaction. We can create a new ConstantRateJump and simulate a hybrid system using both the MassActionJump for the two previous reactions, and the new ConstantRateJump. Let\'s suppose we want to let susceptible people be born with the following jump rate:birth_rate(u,p,t) = 10.*u[1]/(200. + u[1]) + 10.\nfunction birth_affect!(integrator)\n  integrator.u[1] += 1\nend\nbirth_jump = ConstantRateJump(birth_rate, birth_affect!)We can then simulate the hybrid system asjump_prob = JumpProblem(prob, Direct(), mass_act_jump, birth_jump)\nsol = solve(jump_prob, SSAStepper())(Image: gillespie_hybrid_jumps)"
},

{
    "location": "tutorials/discrete_stochastic_example.html#Adding-Jumps-to-a-Differential-Equation-1",
    "page": "Discrete Stochastic (Gillespie) Equations",
    "title": "Adding Jumps to a Differential Equation",
    "category": "section",
    "text": "Notice that if we instead used some form of differential equation instead of a DiscreteProblem, we would add the jumps/reactions to the differential equation. Let\'s define an ODE problem, where the continuous part only acts on some new 4th component:function f(du,u,p,t)\n  du[4] = u[2]*u[3]/100000 - u[1]*u[2]/100000\nend\n\nprob = ODEProblem(f,[999.0,1.0,0.0,100.0],(0.0,250.0))Notice we gave the 4th component a starting value of 100. The same steps as above will allow us to solve this hybrid equation when using ConstantRateJumps. For example, we can solve it using the Tsit5() method via:jump_prob = JumpProblem(prob,Direct(),jump,jump2)\nsol = solve(jump_prob,Tsit5())(Image: gillespie_ode)"
},

{
    "location": "tutorials/discrete_stochastic_example.html#Caution-about-Constant-Rate-Jumps-1",
    "page": "Discrete Stochastic (Gillespie) Equations",
    "title": "Caution about Constant Rate Jumps",
    "category": "section",
    "text": "Note that the assumption which is required for constant rate jumps is that their reaction rates must be constant on the interval between any constant rate jumps. Thus in the examples above,rate(u,p,t) = (0.1/1000.0)*u[1]*u[2]\nrate(u,p,t) = 0.01u[2]both must be constant other than changes due to some constant rate jump (the same applies to reactions). Since these rates only change when u[1] or u[2] is changed, and u[1] and u[2] only change when one of the jumps occur, this setup is valid. However, t*(0.1/1000.0)*u[1]*u[2] would not be valid because the rate would change during the interval, as would (0.1/1000.0)*u[1]*u[4]. Thus one must be careful about to follow this rule when choosing rates.(but note that it\'s okay for u[4] to depend on the other variables because its updated in a continuous manner!)If your problem must have the rates depend on a continuously changing quantity, you need to use the VariableRateJump or VariableRateReaction instead."
},

{
    "location": "tutorials/discrete_stochastic_example.html#Adding-a-VariableRateJump-1",
    "page": "Discrete Stochastic (Gillespie) Equations",
    "title": "Adding a VariableRateJump",
    "category": "section",
    "text": "Now let\'s consider adding a reaction whose rate changes continuously with the differential equation. To continue our example, let\'s let there be a new reaction which has the same effect as r2, but now is dependent on the amount of u[4].rate3(u,p,t) = 1e-2u[4]\nfunction affect3!(integrator)\n  integrator.u[2] -= 1\n  integrator.u[3] += 1\nend\njump3 = VariableRateJump(rate3,affect3!)We would expect this reaction to increase the amount of transitions from state 2 to 3. Solving the equation is exactly the same:prob = ODEProblem(f,[999.0,1.0,0.0,1.0],(0.0,250.0))\njump_prob = JumpProblem(prob,Direct(),jump,jump2,jump3)\nsol = solve(jump_prob,Tsit5())(Image: variable_rate_gillespie)Notice that this increases the amount of 3 at the end, reducing the falloff in the rate (though this model is kind of nonsensical).Note that even if the problem is a DiscreteProblem, VariableRateJumps and VariableRateReactions require a continuous solver, like an ODE/SDE/DDE/DAE solver.Lastly, we are not restricted to ODEs. For example, we can solve the same jump problem except with multiplicative noise on u[4] by using an SDEProblem instead:function g(du,u,p,t)\n  du[4] = 0.1u[4]\nend\n\nprob = SDEProblem(f,g,[999.0,1.0,0.0,1.0],(0.0,250.0))\njump_prob = JumpProblem(prob,Direct(),jump,jump2,jump3)\nsol = solve(jump_prob,SRIW1())(Image: sde_gillespie)"
},

{
    "location": "tutorials/discrete_stochastic_example.html#RegularJumps-and-Tau-Leaping-1",
    "page": "Discrete Stochastic (Gillespie) Equations",
    "title": "RegularJumps and Tau-Leaping",
    "category": "section",
    "text": "The previous parts described how to use ConstantRateJump and VariableRateJump to add jumps to differential equation algorithms over the callback interface. However, in many cases you do not need to step to every jump time. Instead, regular jumping allows you to pool together jumps and perform larger updates in a statistically-correct but more efficient manner.For RegularJumps, we pool together the jumps we wish to perform. Here our rate is a vector equation which computes the rates of each jump process together:function rate(out,u,p,t)\n    out[1] = (0.1/1000.0)*u[1]*u[2]\n    out[2] = 0.01u[2]\nendand then we compute the total change matrix cfunction c(dc,u,p,t,mark)\n    dc[1,1] = -1\n    dc[2,1] = 1\n    dc[2,2] = -1\n    dc[3,2] = 1\nendwhere each column is a different jump process. We then declare the form of dc and build a RegularJump:dc = zeros(3,2)\nrj = RegularJump(rate,c,dc;constant_c=true)From there we build a JumpProblem:prob = DiscreteProblem([999.0,1.0,0.0],(0.0,250.0))\njump_prob = JumpProblem(prob,Direct(),rj)Note that when a JumpProblem has a RegularJump, special algorithms are required. This is detailed on the jump solvers page. One such algorithm is SimpleTauLeaping, which we use as follows:sol = solve(jump_prob,SimpleTauLeaping();dt=1.0)"
},

{
    "location": "tutorials/jump_diffusion.html#",
    "page": "Jump Diffusion Equations",
    "title": "Jump Diffusion Equations",
    "category": "page",
    "text": ""
},

{
    "location": "tutorials/jump_diffusion.html#Jump-Diffusion-Equations-1",
    "page": "Jump Diffusion Equations",
    "title": "Jump Diffusion Equations",
    "category": "section",
    "text": "This tutorial assumes you have read the Ordinary Differential Equations tutorial.Jump Diffusion equations are stochastic diffeential equations with discontinuous jumps. These can be written as:fracdudt = f(upt) + Σgᵢ(ut)dWⁱ + Σ h_i(upt)N_i(t)where N_i is a Poisson-counter which denotes jumps of size h_i. In this tutorial we will show how to solve problems with even more general jumps."
},

{
    "location": "tutorials/jump_diffusion.html#Defining-a-ConstantRateJump-Problem-1",
    "page": "Jump Diffusion Equations",
    "title": "Defining a ConstantRateJump Problem",
    "category": "section",
    "text": "To start, let\'s solve an ODE with constant rate jumps. A jump is defined as being \"constant rate\" if the rate is only dependent on values from other constant rate jumps, meaning that its rate must not be coupled with time or the solution to the differential equation. However, these types of jumps are cheaper to compute.(Note: if your rate is only \"slightly\" dependent on the solution of the differential equation, then it may be okay to use a ConstantRateJump. Accuracy loss will be related to the percentage that the rate changes over the jump intervals.)Let\'s solve the following problem. We will have a linear ODE with a Poisson counter of rate 2 (which is the mean and variance), where at each jump the current solution will be halved. To solve this problem, we first define the ODEProblem:function f(du,u,p,t)\n  du[1] = u[1]\nend\n\nprob = ODEProblem(f,[0.2],(0.0,10.0))Notice that, even though our equation is on 1 number, we define it using the in-place array form. Variable rate jump equations will require this form. Note that for this tutorial we solve a one-dimensional problem, but the same syntax applies for solving a system of differential equations with multiple jumps.Now we define our rate equation for our jump. Since it\'s just the constant value 2, we do:rate(u,p,t) = 2Now we define the affect! of the jump. This is the same as an affect! from a DiscreteCallback, and thus acts directly on the integrator. Therefore, to make it halve the current value of u, we do:affect!(integrator) = (integrator.u[1] = integrator.u[1]/2)Then we build our jump:jump = ConstantRateJump(rate,affect!)Next, we extend our ODEProblem to a JumpProblem by attaching the jump:jump_prob = JumpProblem(prob,Direct(),jump)We can now solve this extended problem using any ODE solver:sol = solve(jump_prob,Tsit5())\nplot(sol)(Image: constant_rate_jump)"
},

{
    "location": "tutorials/jump_diffusion.html#Variable-Rate-Jumps-1",
    "page": "Jump Diffusion Equations",
    "title": "Variable Rate Jumps",
    "category": "section",
    "text": "Now let\'s define a jump which is coupled to the differential equation. Let\'s let the rate be the current value of the solution, that is:rate(u,p,t) = u[1]Using the same affect!affect!(integrator) = (integrator.u[1] = integrator.u[1]/2)we build a VariableRateJump:jump = VariableRateJump(rate,affect!)To make things interesting, let\'s copy this jump:jump2 = deepcopy(jump)so that way we have two independent jump processes. We now couple these jumps to the ODEProblem:jump_prob = JumpProblem(prob,Direct(),jump,jump2)which we once again solve using an ODE solver:sol = solve(jump_prob,Tsit5())\nplot(sol)(Image: variable_rate_jump)"
},

{
    "location": "tutorials/jump_diffusion.html#Jump-Diffusion-1",
    "page": "Jump Diffusion Equations",
    "title": "Jump Diffusion",
    "category": "section",
    "text": "Now we will finally solve the jump diffusion problem. The steps are the same as before, except we now start with a SDEProblem instead of an ODEProblem. Using the same drift function f as before, we add multiplicative noise via:function g(du,u,p,t)\n  du[1] = u[1]\nend\n\nprob = SDEProblem(f,g,[0.2],(0.0,10.0))and couple it to the jumps:jump_prob = JumpProblem(prob,Direct(),jump,jump2)We then solve it using an SDE algorithm:sol = solve(jump_prob,SRIW1())\nplot(sol)(Image: jump_diffusion)"
},

{
    "location": "tutorials/jump_diffusion.html#Coupling-Jump-Problems-1",
    "page": "Jump Diffusion Equations",
    "title": "Coupling Jump Problems",
    "category": "section",
    "text": "In many applications one is interested in coupling two stochastic processes. This has applications in Monte Carlo simulations and sensitivity analysis, for example. Currently, the coupling that is implemented for jump processes is known as the split coupling. The split coupling couples two jump processes by coupling the underlying Poisson processes driving the jump components.Suppose prob and prob_control are two problems we wish to couple. Then the coupled problem is obtained byprob_coupled =  SplitCoupledJumpProblem(jump_prob,jump_prob_control,Direct(),coupling_map)Here, coupling_map specifies which jumps to couple. If (j,i) is in coupling_map, then the ith jump in prob will be coupled to the jth jump in prob_control. Note that currently SplitCoupledJumpProblem is only implemented for constant rate jump problems.As an example, consider a doubly stochastic Poisson process, that is, a Poisson process whose rate is itself a stochastic process. In particular, we will take the rate to randomly switch between zero and 10 at unit rates:rate(u,p,t) = u[2]*10\naffect!(integrator) = integrator.u[1] += 1.\njump1 = ConstantRateJump(rate,affect!)\nrate(u,p,t) = u[2]\naffect!(integrator) = (integrator.u[2] -= 1.;integrator.u[3] += 1.)\njump2 = ConstantRateJump(rate,affect!)\n\nrate(u,p,t) = u[3]\naffect!(integrator) = (integrator.u[2] += 1.;integrator.u[3] -= 1.)\njump3 = ConstantRateJump(rate,affect!)\nprob = DiscreteProblem(u0,tspan)\njump_prob = JumpProblem(prob,Direct(),jump1,jump2,jump3)The doubly stochastic poisson process has two sources of randomness: one due to the Poisson process, and another due to random evolution of the rate. This is typical of many multiscale stochastic processes appearing in applications, and it is often useful to compare such a process to one obtained by removing one source of randomness. In present context, this means looking at an ODE with constant jump rates, where the deterministic evolution between jumps is given by the expected value of the Poisson process:function f(du,u,p,t)\n  du[1] = u[2]*10\n  du[2] = 0.\n  du[3] = 0.\nend\nprob_control = ODEProblem(f,u0,tspan)\njump_prob_control = JumpProblem(prob_control,Direct(),jump2,jump3)Let\'s couple the two problems by coupling the jumps corresponding the switching of the rate:coupling_map = [(2,1),(3,2)]\nprob_coupled =  SplitCoupledJumpProblem(jump_prob,jump_prob_control,Direct(),coupling_map)Now prob_coupled will be dealt with like any other JumpProblem:sol = solve(coupled_prob,Tsit5())(Image: jump_diffusion)"
},

{
    "location": "tutorials/bvp_example.html#",
    "page": "Boundary Value Problems",
    "title": "Boundary Value Problems",
    "category": "page",
    "text": ""
},

{
    "location": "tutorials/bvp_example.html#Boundary-Value-Problems-1",
    "page": "Boundary Value Problems",
    "title": "Boundary Value Problems",
    "category": "section",
    "text": "This tutorial will introduce you to the functionality for solving BVPs. Other introductions can be found by checking out DiffEqTutorials.jl. This tutorial assumes you have read the Ordinary Differential Equations tutorial.In this example we will solve the ODE that satisfies the boundary condition in the form ofbeginalign\nfracddt = f(t u) \ng(u) = vec0\nendalign"
},

{
    "location": "tutorials/bvp_example.html#Example-1:-Simple-Pendulum-1",
    "page": "Boundary Value Problems",
    "title": "Example 1: Simple Pendulum",
    "category": "section",
    "text": "The concrete example that we are solving is the simple pendulum ddotu+fracgLu=0 on the time interval tin0fracpi2. First, we need to define the ODEusing BoundaryValueDiffEq\nconst g = 9.81\nL = 1.0\ntspan = (0.0,pi/2)\nfunction simplependulum!(du,u,p,t)\n    θ  = u[1]\n    dθ = u[2]\n    du[1] = dθ\n    du[2] = -(g/L)*sin(θ)\nend"
},

{
    "location": "tutorials/bvp_example.html#Boundary-Condition-1",
    "page": "Boundary Value Problems",
    "title": "Boundary Condition",
    "category": "section",
    "text": "There are two problem types available:A problem type for general boundary conditions BVProblem ( including conditions that may be anywhere/ everywhere on the integration interval ).\nA problem type for boundaries that are specified at the beginning and the end of the integration interval TwoPointBVProblem"
},

{
    "location": "tutorials/bvp_example.html#BVProblem-1",
    "page": "Boundary Value Problems",
    "title": "BVProblem",
    "category": "section",
    "text": "The boundary conditions are specified by a function that calculates the residual in-place from the problem solution, such that the residual is vec0 when the boundary condition is satisfied.function bc1!(residual, u, p, t)\n    residual[1] = u[end÷2][1] + pi/2 # the solution at the middle of the time span should be -pi/2\n    residual[2] = u[end][1] - pi/2 # the solution at the end of the time span should be pi/2\nend\nbvp1 = BVProblem(simplependulum!, bc1!, [pi/2,pi/2], tspan)\nsol1 = solve(bvp1, GeneralMIRK4(), dt=0.05)\nplot(sol1)(Image: BVP Example Plot1)The third argument of BVProblem  is the initial guess of the solution, which is constant in this example. <!– add examples of more general initial conditions –> We need to use GeneralMIRK4 or Shooting methods to solve BVProblem. GeneralMIRK4 is a collocation method, whereas Shooting treats the problem as an IVP and varies the initial conditions until the boundary conditions are met. If you can have a good initial guess, Shooting method works very well.using OrdinaryDiffEq\nu₀_2 = [-1.6, -1.7] # the initial guess\nfunction bc3!(residual, sol, p, t)\n    residual[1] = sol(pi/4)[1] + pi/2 # use the interpolation here, since indexing will be wrong for adaptive methods\n    residual[2] = sol(pi/2)[1] - pi/2\nend\nbvp3 = BVProblem(simplependulum!, bc3!, u₀_2, tspan)\nsol3 = solve(bvp3, Shooting(Vern7()))The initial guess can also be supplied via a function of t or a previous solution type, this is espacially handy for parameter analysis. We changed u to sol to emphasize the fact that in this case the boundary condition can be written on the solution object. Thus all of the features on the solution type such as interpolations are available when using the Shooting method (i.e. you can have a boundary condition saying that the maximum over the interval is 1 using an optimization function on the continuous output). Note that user has to import the IVP solver before it can be used. Any common interface ODE solver is acceptable.plot(sol3)(Image: BVP Example Plot3)"
},

{
    "location": "tutorials/bvp_example.html#TwoPointBVProblem-1",
    "page": "Boundary Value Problems",
    "title": "TwoPointBVProblem",
    "category": "section",
    "text": "Defining a similar problem as TwoPointBVProblem is shown in the following example. At the moment MIRK4 is the only solver for TwoPointBVProblems.function bc2!(residual, u, p, t) # u[1] is the beginning of the time span, and u[end] is the ending\n    residual[1] = u[1][1] + pi/2 # the solution at the beginning of the time span should be -pi/2\n    residual[2] = u[end][1] - pi/2 # the solution at the end of the time span should be pi/2\nend\nbvp2 = TwoPointBVProblem(simplependulum!, bc2!, [pi/2,pi/2], tspan)\nsol2 = solve(bvp2, MIRK4(), dt=0.05) # we need to use the MIRK4 solver for TwoPointBVProblem\nplot(sol2)Note that u is a tuple of ( u[1], u[end] ) just like t is ( t[1], t[end] ) and p holds the parameters of the given problem.(Image: BVP Example Plot2)"
},

{
    "location": "tutorials/additional.html#",
    "page": "Additional Tutorials",
    "title": "Additional Tutorials",
    "category": "page",
    "text": ""
},

{
    "location": "tutorials/additional.html#Additional-Tutorials-1",
    "page": "Additional Tutorials",
    "title": "Additional Tutorials",
    "category": "section",
    "text": "Additional tutorials can be found at DiffEqTutorials.jl. These include interactive introductions, optimizing code, modeling examples, and deeper examples for extra features."
},

{
    "location": "basics/overview.html#",
    "page": "Overview of DifferentialEquations.jl",
    "title": "Overview of DifferentialEquations.jl",
    "category": "page",
    "text": ""
},

{
    "location": "basics/overview.html#Overview-of-DifferentialEquations.jl-1",
    "page": "Overview of DifferentialEquations.jl",
    "title": "Overview of DifferentialEquations.jl",
    "category": "section",
    "text": "The general workflow for using the package is as follows:Define a problem\nSolve the problem\nAnalyze the output"
},

{
    "location": "basics/overview.html#Defining-Problems-1",
    "page": "Overview of DifferentialEquations.jl",
    "title": "Defining Problems",
    "category": "section",
    "text": "Problems are specified via a type interface. The problem types are designed to contain the necessary information to fully define their associated differential equation. Each problem type has a page explaining their problem type and the special features associated with them. For example, an ordinary differential equation is defined byfracdudt = f(upt)over some time interval tspan with some initial condition u0, and therefore the ODEProblem is defined by those components:prob = ODEProblem(f,u0,tspan)\nprob = ODEProblem(f,u0,tspan,p)Note that the number types in the solution will match the types you designate in the problem. For example, if one uses Rational{BigInt} for specifying the timespan and BigFloat for specifying the initial condition, then the solution will solve using Rational{BigInt} for the timesteps and BigFloat for the independent variables. A wide variety of number types are compatible with the solvers such as complex numbers, unitful numbers (via Unitful.jl), decimals (via DecFP), dual numbers, and many more which may not have been tested yet (thanks to the power of multiple dispatch!). For information on type-compatibilty, please see the solver pages for the specific problems."
},

{
    "location": "basics/overview.html#Solving-the-Problems-1",
    "page": "Overview of DifferentialEquations.jl",
    "title": "Solving the Problems",
    "category": "section",
    "text": "Each type of differential equation has its own problem type which allow the solvers to dispatch to the right methods. The common interface for calling the solvers is:sol = solve(prob,alg;kwargs)Into the command, one passes the differential equation problem that they defined prob, optionally choose an algorithm alg (a default is given if not chosen), and change the properties of the solver using keyword arguments. The common arguments which are accepted by most methods is defined in the common solver options manual page. The solver returns a solution object sol which hold all of the details for the solution."
},

{
    "location": "basics/overview.html#Analyzing-the-Solution-1",
    "page": "Overview of DifferentialEquations.jl",
    "title": "Analyzing the Solution",
    "category": "section",
    "text": "With the solution object, you do the analysis as you please! The solution type has a common interface which makes handling the solution similar between the different types of differential equations. Tools such as interpolations are seamlessly built into the solution interface to make analysis easy. This interface is described in the solution handling manual page.Plotting functionality is provided by a recipe to Plots.jl. To use plot solutions, simply call the plot(sol) and the plotter will generate appropriate plots. If save_everystep was used, the plotters can generate animations of the solutions to evolution equations using the animate(sol) command. Plots can be customized using all of the keyword arguments provided by Plots.jl. Please see Plots.jl\'s documentation for more information."
},

{
    "location": "basics/overview.html#Add-on-Tools-1",
    "page": "Overview of DifferentialEquations.jl",
    "title": "Add-on Tools",
    "category": "section",
    "text": "One of the most compelling features of DifferentialEquations.jl is that the common solver interface allows one to build tools which are \"algorithm and problem agnostic\". For example, one of the provided tools allows for performing parameter estimation on ODEProblems. Since the solve interface is the same for the different algorithms, one can use any of the associated solving algorithms. This modular structure allows one to mix and match overarching analysis tools with specialized algorithms to one\'s problem, leading to high performance with a large feature base. Isn\'t that the promise of Julia just being fulfilled?"
},

{
    "location": "basics/overview.html#Development-and-Testing-Tools-1",
    "page": "Overview of DifferentialEquations.jl",
    "title": "Development and Testing Tools",
    "category": "section",
    "text": "Lastly, one unique feature of DifferentialEquations.jl is the existence of algorithm development and testing functionality. This suite was designed by researchers in the field of numerical differential equations to both try out new ideas and distribute finalized results to large audiences. The tools for algorithm development allow for easy convergence testing, benchmarking, and higher order analysis (stability plotting, etc.). This is one of the reasons why DifferentialEquations.jl contains many algorithms which are unique and the results of recent publications! Please check out the developer documentation for more information on using the development tools.Note that DifferentialEquations.jl allows for distributed development, meaning that algorithms which \"plug-into ecosystem\" don\'t have to be a part of the major packages. If you are interested in adding your work to the ecosystem, checkout the developer documentation for more information."
},

{
    "location": "basics/common_solver_opts.html#",
    "page": "Common Solver Options",
    "title": "Common Solver Options",
    "category": "page",
    "text": ""
},

{
    "location": "basics/common_solver_opts.html#Common-Solver-Options-1",
    "page": "Common Solver Options",
    "title": "Common Solver Options",
    "category": "section",
    "text": "The DifferentialEquations.jl universe has a large set of common arguments available for the solve function. These arguments apply to solve on any problem type and are only limited by limitations of the specific implementations.Many of the defaults depend on the algorithm or the package the algorithm derives from. Not all of the interface is provided by every algorithm. For more detailed information on the defaults and the available options for specific algorithms / packages, see the manual pages for the solvers of specific problems. To see whether a specific package is compaible with the use of a given option, see the Solver Compatibility Chart"
},

{
    "location": "basics/common_solver_opts.html#Default-Algorithm-Hinting-1",
    "page": "Common Solver Options",
    "title": "Default Algorithm Hinting",
    "category": "section",
    "text": "To help choose the default algorithm, the keyword argument alg_hints is provided to solve. alg_hints is a Vector{Symbol} which describe the problem at a high level to the solver. The options are::auto vs :nonstiff vs :stiff - Denotes the equation as nonstiff/stiff. :auto allow the default handling algorithm to choose stiffness detection algorithms. The default handling defaults to using :auto.Currently unused options include::interpolant - Denotes that a high-precision interpolation is important.\n:memorybound - Denotes that the solver will be memory bound.This functionality is derived via the benchmarks in DiffEqBenchmarks.jl"
},

{
    "location": "basics/common_solver_opts.html#SDE-Specific-Alghints-1",
    "page": "Common Solver Options",
    "title": "SDE Specific Alghints",
    "category": "section",
    "text": ":additive - Denotes that the underlying SDE has additive noise.\n:stratonovich - Denotes that the solution should adhere to the Stratonovich interpretation."
},

{
    "location": "basics/common_solver_opts.html#Output-Control-1",
    "page": "Common Solver Options",
    "title": "Output Control",
    "category": "section",
    "text": "These arguments control the output behavior of the solvers. It defaults to maximum output to give the best interactive user experience, but can be reduced all the way to only saving the solution at the final timepoint.The following options are all related to output control. See the \"Examples\" section at the end of this page for some example usage.dense: Denotes whether to save the extra pieces required for dense (continuous) output. Default is true for algorithms which have the ability to produce dense output. If dense is false, the solution still acts like a function, and sol(t) is a linear interpolation between the saved time points.\nsaveat: Denotes specific times to save the solution at, during the solving phase. The solver will save at each of the timepoints in this array in the most efficient manner available to the solver. Note that this can be used even if dense=false. If only saveat is given, then the arguments save_everystep and dense are false by default.    If saveat is given a number, then it will automatically expand to  tspan[1]:saveat:tspan[2]. For methods where interpolation is not possible,  saveat may be equivalent to tstops. The default value is [].\nsave_idxs: Denotes the indices for the components of the equation to save. Defaults to saving all indices. For example, if you are solving a 3-dimensional ODE, and given save_idxs = [1, 3], only the first and third components of the solution will be outputted. Notice that of course in this case the outputed solution will be two-dimensional.\ntstops: Denotes extra times that the timestepping algorithm must step to. This should be used to help the solver deal with discontinuities and singularities, since stepping exactly at the time of the discontinuity will improve accuracy. If a method cannot change timesteps (fixed timestep multistep methods), then tstops will use an interpolation, matching the behavior of saveat. If a method cannot change timesteps and also cannot interpolate, then tstops must be a multiple of dt or else an error will be thrown. Default is [].\nd_discontinuities: Denotes locations of discontinuities in low order derivatives. This will force FSAL algorithms which assume derivative continuity to re-evaluate the derivatives at the point of discontinuity. The default is [].\nsave_everystep: Saves the result at every step.     Default is true if isempty(saveat).\nsave_on: Denotes whether intermediate solutions are saved. This overrides the settings of dense, saveat and save_everystep and is used by some applicatioins to manually turn off saving temporarily. Everyday use of the solvers should leave this unchanged. Defaults to true.\nsave_start: Denotes whether the initial condition should be included in the solution type as the first timepoint. Defaults to true.\nsave_end: Denotes whether the final timepoint is forced to be saved, regardless of the other saving settings. Defaults to true.\ninitialize_save: Denotes whether to save after the callback initialization phase (when u_modified=true). Defaults to true.Note that dense requires save_everystep=true and saveat=false. If you need additional saving while keeping dense output, see the SavingCallback in the Callback Library."
},

{
    "location": "basics/common_solver_opts.html#Stepsize-Control-1",
    "page": "Common Solver Options",
    "title": "Stepsize Control",
    "category": "section",
    "text": "These arguments control the timestepping routines."
},

{
    "location": "basics/common_solver_opts.html#Basic-Stepsize-Control-1",
    "page": "Common Solver Options",
    "title": "Basic Stepsize Control",
    "category": "section",
    "text": "These are the standard options for controlling stepping behavior. Error estimates do the comparisonerr_scaled = err(abstol + max(uprevu)*reltol)The scaled error is guaranteed to be <1 for a given local error estimate (note: error estimates are local unless the method specifies otherwise). abstol controls the non-scaling error and thus can be though of as the error around zero. reltol scales with the size of the dependent variables and so one can interpret reltol=1e-3 as roughly being (locally) correct to 3 digits. Note tolerances can be specified element-wise by passing a vector whose size matches u0.adaptive: Turns on adaptive timestepping for appropriate methods. Default is true.\nabstol: Absolute tolerance in adaptive timestepping. This is the tolerance on local error estimates, not necessarily the global error (though these quantities are related). Defaults to 1e-6 on deterministic equations (ODEs/DDEs/DAEs) and 1e-2 on stochastic equations (SDEs/RODEs).\nreltol: Relative tolerance in adaptive timestepping.  This is the tolerance on local error estimatoes, not necessarily the global error (though these quantities are related). Defaults to 1e-3 on deterministic equations (ODEs/DDEs/DAEs) and 1e-2 on stochastic equations (SDEs/RODEs).\ndt: Sets the initial stepsize. This is also the stepsize for fixed timestep methods. Defaults to an automatic choice if the method is adaptive.\ndtmax: Maximum dt for adaptive timestepping. Defaults are package-dependent.\ndtmin: Minimum dt for adaptive timestepping. Defaults are package-dependent.\nforce_dtmin: Declares whether to continue, forcing the minimum dt usage. Default is false, which has the solver throw a warning and exit early when encountering the minimum dt. Setting this true allows the solver to continue, never letting dt go below dtmin (and ignoring error tolerances in those cases). Note that true is not compatible with most interop packages."
},

{
    "location": "basics/common_solver_opts.html#Fixed-Stepsize-Usage-1",
    "page": "Common Solver Options",
    "title": "Fixed Stepsize Usage",
    "category": "section",
    "text": "Note that if a method does not have adaptivity, the following rules apply:If dt is set, then the algorithm will step with size dt each iteration.\nIf tstops and dt are both set, then the algorithm will step with either a size dt, or use a smaller step to hit the tstops point.\nIf tstops is set without dt, then the algorithm will step directly to each value in tstops\nIf neither dt nor tstops are set, the solver will throw an error."
},

{
    "location": "basics/common_solver_opts.html#Advanced-Adaptive-Stepsize-Control-1",
    "page": "Common Solver Options",
    "title": "Advanced Adaptive Stepsize Control",
    "category": "section",
    "text": "These arguments control more advanced parts of the internals of adaptive timestepping and are mostly used to make it more efficient on specific problems. For detained explanations of the timestepping algorithms, see the timestepping descriptionsinternalnorm: The norm function internalnorm(u,t) which error estimates are calculated. Required are two dispatches: one dispatch for the state variable and the other on the elements of the state variable (scalar norm). Defaults are package-dependent.\ngamma: The risk-factor γ in the q equation for adaptive timestepping. Default is algorithm dependent.\nbeta1: The Lund stabilization α parameter. Defaults are algorithm-dependent.\nbeta2: The Lund stabilization β parameter. Defaults are algorithm-dependent.\nqmax: Defines the maximum value possible for the adaptive q. Defaults are algorithm-dependent.\nqmin: Defines the minimum value possible for the adaptive q. Defaults are algorithm-dependent.\nqsteady_min: Defines the minimum for the range around 1 where the timestep is held constant. Defaults are algorithm-dependent.\nqsteady_max: Defines the maximum for the range around 1 where the timestep is held constant. Defaults are algorithm-dependent.\nqoldinit: The initial qold in stabilization stepping. Defaults are algorithm-dependent.\nfailfactor: The amount to decrease the timestep by if the Newton iterations of an implicit method fail. Default is 2."
},

{
    "location": "basics/common_solver_opts.html#Memory-Optimizations-1",
    "page": "Common Solver Options",
    "title": "Memory Optimizations",
    "category": "section",
    "text": "calck: Turns on and off the internal ability for intermediate     interpolations (also known as intermediate density). Not the same as dense, which is post-solution interpolation. This defaults to dense || !isempty(saveat) ||  \"no custom callback is given\". This can be used to turn off interpolations (to save memory) if one isn\'t using interpolations when a custom callback is used. Another case where this may be used is to turn on interpolations for usage in the integrator interface even when interpolations are used nowhere else. Note that this is only required if the algorithm doesn\'t have a free or lazy interpolation (DP8()). If calck = false, saveat cannot be used. The rare keyword calck can be useful in event handling.\nalias_u0: allows the solver to alias the initial condition array that is contained in the problem struct. Defaults to false."
},

{
    "location": "basics/common_solver_opts.html#Miscellaneous-1",
    "page": "Common Solver Options",
    "title": "Miscellaneous",
    "category": "section",
    "text": "maxiters: Maximum number of iterations before stopping. Defaults to 1e5.\ncallback: Specifies a callback. Defaults to a callback function which performs the saving routine. For more information, see the Event Handling and Callback Functions manual page.\nisoutofdomain: Specifies a function isoutofdomain(u,p,t) where, when it returns true, it will reject the timestep. Disabled by default.\nunstable_check: Specifies a function unstable_check(dt,u,p,t) where, when it returns true, it will cause the solver to exit and throw a warning. Defaults to any(isnan,u), i.e. checking if any value is a NaN.\nverbose: Toggles whether warnings are thrown when the solver exits early. Defaults to true."
},

{
    "location": "basics/common_solver_opts.html#Progress-Monitoring-1",
    "page": "Common Solver Options",
    "title": "Progress Monitoring",
    "category": "section",
    "text": "These arguments control the usage of the progressbar in the Juno IDE.progress: Turns on/off the Juno progressbar. Default is false.\nprogress_steps: Numbers of steps between updates of the progress bar. Default is 1000.\nprogress_name: Controls the name of the progressbar. Default is the name of the problem type.\nprogress_message: Controls the message with the progressbar. Defaults to showing dt, t, the maximum of u."
},

{
    "location": "basics/common_solver_opts.html#Error-Calculations-1",
    "page": "Common Solver Options",
    "title": "Error Calculations",
    "category": "section",
    "text": "If you are using the test problems (ex: ODETestProblem), then the following options control the errors which are calculated:timeseries_errors: Turns on and off the calculation of errors at the steps which were taken, such as the l2 error. Default is true.\ndense_errors: Turns on and off the calculation of errors at the steps which require dense output and calculate the error at 100 evenly-spaced points throughout tspan. An example is the L2 error. Default is false."
},

{
    "location": "basics/common_solver_opts.html#Examples-1",
    "page": "Common Solver Options",
    "title": "Examples",
    "category": "section",
    "text": "The following lines are examples of how one could use the configuration of solve(). For these examples a 3-dimensional ODE problem is assumed, however the extention to other types is straightforward.solve(prob, AlgorithmName()) : The \"default\" setting, with a user-specifiedalgorithm (given by AlgorithmName()). All parameters get their default values.   This means that the solution is saved at the steps the Algorithm stops internally   and dense output is enabled if the chosen algorithm allows for it.All other integration parameters (e.g. stepsize) are chosen automatically.solve(prob, saveat = 0.01, abstol = 1e-9, reltol = 1e-9) : Standard settingfor accurate output at specified (and equidistant) time intervals, used for   e.g. Fourier Transform. The solution is given every 0.01 time units,   starting from tspan[1]. The solver used is Tsit5() since no keyword   alg_hits is given.solve(prob, maxiters = 1e7, progress = true, save_idxs = [1]) : Using longermaximum number of solver iterations can be useful when a given tspan is very   long. This example only saves the first of the variables of the system, either   to save size or because the user does not care about the others. Finally, with   progress = true you are enabling the progress bar, provided you are using   the Atom+Juno IDE set-up for your Julia."
},

{
    "location": "basics/solution.html#",
    "page": "Solution Handling",
    "title": "Solution Handling",
    "category": "page",
    "text": ""
},

{
    "location": "basics/solution.html#Solution-Handling-1",
    "page": "Solution Handling",
    "title": "Solution Handling",
    "category": "section",
    "text": ""
},

{
    "location": "basics/solution.html#Accessing-the-Values-1",
    "page": "Solution Handling",
    "title": "Accessing the Values",
    "category": "section",
    "text": "The solution type has a lot of built in functionality to help analysis. For example, it has an array interface for accessing the values. Internally, the solution type has two important fields:u which holds the Vector of values at each timestep\nt which holds the times of each timestep.Different solution types may add extra information as necessary, such as the derivative at each timestep du or the spatial discretization x, y, etc."
},

{
    "location": "basics/solution.html#Array-Interface-1",
    "page": "Solution Handling",
    "title": "Array Interface",
    "category": "section",
    "text": "Instead of working on the Vector{uType} directly, we can use the provided array interface.sol[i]to access the value at timestep i (if the timeseries was saved), andsol.t[i]to access the value of t at timestep i. For multi-dimensional systems, this will address first by component and lastly by time, and thussol[i,j]will be the ith component at timestep j. If the independent variables had shape (for example, was a matrix), then i is the linear index. We can also access solutions with shape:sol[i,j,k]gives the [i,j] component of the system at timestep k. The colon operator is supported, meaning thatsol[j,:]gives the timeseries for the jth component."
},

{
    "location": "basics/solution.html#Using-the-AbstractArray-Interface-1",
    "page": "Solution Handling",
    "title": "Using the AbstractArray Interface",
    "category": "section",
    "text": "The AbstractArray interface can be directly used. For example, for a vector system of variables sol[i,j] is a matrix with rows being the variables and columns being the timepoints. Operations like sol\' will transpose the solution type. Functionality written for AbstractArrays can directly use this. For example, the Base cov function computes correlations amongst columns, and thus:cov(sol)computes the correlation of the system state in time, whereascov(sol,2)computes the correlation between the variables. Similarly, mean(sol,2) is the mean of the variable in time, and var(sol,2) is the variance. Other statistical functions and packages which work on AbstractArray types will work on the solution type.At anytime, a true Array can be created using Array(sol)."
},

{
    "location": "basics/solution.html#Interpolations-1",
    "page": "Solution Handling",
    "title": "Interpolations",
    "category": "section",
    "text": "If the solver allows for dense output and dense=true was set for the solving (which is the default), then we can access the approximate value at a time t using the commandsol(t)Note that the interpolating function allows for t to be a vector and uses this to speed up the interpolation calculations. The full API for the interpolations issol(t,deriv=Val{0};idxs=nothing)The optional argument deriv lets you choose the number n derivative to solve the interpolation for, defaulting with n=0. Note that most of the derivatives have not yet been implemented (though it\'s not hard, it just has to be done by hand for each algorithm. Open an issue if there\'s a specific one you need). idxs allows you to choose the indices the interpolation should solve for. For example,sol(t,idxs=1:2:5)will return a Vector of length 3 which is the interpolated values at t for components 1, 3, and 5. idxs=nothing, the default, means it will return every component. In addition, we can dosol(t,idxs=1)and it will return a Number for the interpolation of the single value. Note that this interpolation only computes the values which are requested, and thus it\'s much faster on large systems to use this rather than computing the full interpolation and using only a few values.In addition, there is an inplace form:sol(out,t,deriv=Val{0};idxs=nothing)which will write the output to out. This allows one to use pre-allocated vectors for the output to improve the speed even more."
},

{
    "location": "basics/solution.html#Comprehensions-1",
    "page": "Solution Handling",
    "title": "Comprehensions",
    "category": "section",
    "text": "The solver interface also gives tools for using comprehensions over the solution. Using the tuples(sol) function, we can get a tuple for the output at each timestep. This allows one to do the following:[t+2u for (u,t) in tuples(sol)]One can use the extra components of the solution object as well as using zip. For example, say the solution type holds du, the derivative at each timestep. One can comprehend over the values using:[t+3u-du for (t,u,du) in zip(sol.t,sol.u,sol.du)]Note that the solution object acts as a vector in time, and so its length is the number of saved timepoints."
},

{
    "location": "basics/solution.html#Special-Fields-1",
    "page": "Solution Handling",
    "title": "Special Fields",
    "category": "section",
    "text": "The solution interface also includes some special fields. The problem object prob and the algorithm used to solve the problem alg are included in the solution. Additionally, the field dense is a boolean which states whether the interpolation functionality is available. Lastly, there is a mutable state tslocation which controls the plot recipe behavior. By default, tslocation=0. Its values have different meanings between partial and ordinary differential equations:tslocation=0  for non-spatial problems (ODEs) means that the plot recipe will plot the full solution. tslocation=i means that it will only plot the timepoint i.\ntslocation=0 for spatial problems (PDEs) means the plot recipe will plot the final timepoint. tslocation=i means that the plot recipe will plot the ith timepoint.What this means is that for ODEs, the plots will default to the full plot and PDEs will default to plotting the surface at the final timepoint. The iterator interface simply iterates the value of tslocation, and the animate function iterates the solution calling solve at each step."
},

{
    "location": "basics/solution.html#Return-Codes-(RetCodes)-1",
    "page": "Solution Handling",
    "title": "Return Codes (RetCodes)",
    "category": "section",
    "text": "The solution types have a retcode field which returns a symbol signifying the error state of the solution. The retcodes are as follows::Default: The solver did not set retcodes.\n:Success: The integration completed without erroring or the steady state solver from SteadyStateDiffEq found the steady state.\n:Terminated: The integration is terminated with terminate!(integrator). Note that this may occur by using TerminateSteadyState from the callback library DiffEqCallbacks.\n:MaxIters: The integration exited early because it reached its maximum number of iterations.\n:DtLessThanMin: The timestep method chose a stepsize which is smaller than the allowed minimum timestep, and exited early.\n:Unstable: The solver detected that the solution was unstable and exited early.\n:InitialFailure: The DAE solver could not find consistent initial conditions.\n:ConvergenceFailure: The internal implicit solvers failed to converge.\n:Failure: General uncategorized failures or errors."
},

{
    "location": "basics/solution.html#Problem-Specific-Features-1",
    "page": "Solution Handling",
    "title": "Problem-Specific Features",
    "category": "section",
    "text": "Extra fields for solutions of specific problems are specified in the appropriate problem definition page.  "
},

{
    "location": "basics/plot.html#",
    "page": "Plot Functions",
    "title": "Plot Functions",
    "category": "page",
    "text": ""
},

{
    "location": "basics/plot.html#Plot-Functions-1",
    "page": "Plot Functions",
    "title": "Plot Functions",
    "category": "section",
    "text": ""
},

{
    "location": "basics/plot.html#Standard-Plots-Using-the-Plot-Recipe-1",
    "page": "Plot Functions",
    "title": "Standard Plots Using the Plot Recipe",
    "category": "section",
    "text": "Plotting functionality is provided by recipes to Plots.jl. To plot solutions, simply call the plot(type) after importing Plots.jl and the plotter will generate appropriate plots.#]add Plots # You need to install Plots.jl before your first time using it!\nusing Plots\nplot(sol) # Plots the solutionMany of the types defined in the DiffEq universe, such as ODESolution, ConvergenceSimulation WorkPrecision, etc. have plot recipes to handle the default plotting behavior. Plots can be customized using all of the keyword arguments provided by Plots.jl. For example, we can change the plotting backend to the GR package and put a title on the plot by doing:gr()\nplot(sol,title=\"I Love DiffEqs!\")"
},

{
    "location": "basics/plot.html#Density-1",
    "page": "Plot Functions",
    "title": "Density",
    "category": "section",
    "text": "If the problem was solved with dense=true, then denseplot controls whether to use the dense function for generating the plot, and plotdensity is the number of evenly-spaced points (in time) to plot. For example:plot(sol,denseplot=false)means \"only plot the points which the solver stepped to\", while:plot(sol,plotdensity=1000)means to plot 1000 points using the dense function (since denseplot=true by default)."
},

{
    "location": "basics/plot.html#Choosing-Variables-1",
    "page": "Plot Functions",
    "title": "Choosing Variables",
    "category": "section",
    "text": "In the plot command, one can choose the variables to be plotted in each plot. The master form is:vars = [(f1,0,1), (f2,1,3), (f3,4,5)]which could be used to plot f1(var₀, var₁), f2(var₁, var₃), and f3(var₄, var₅), all on the same graph. (0 is considered to be time, or the independent variable). Functions f1, f2 and f3 should take in scalars and return a tuple. If no function is given, for example,vars = [(0,1), (1,3), (4,5)]this would mean \"plot var₁(t) vs t (time), var₃(var₁) vs var₁, and var₅(var₄) vs var₄ all on the same graph, putting the independent variables (t, var₁ and var₄) on the x-axis.\" While this can be used for everything, the following conveniences are provided:Everywhere in a tuple position where we only find an integer, this variable is plotted as a function of time.  For example, the list above is equivalent to:vars = [1, (1,3), (4,5)]andvars = [1, 3, 4]is the most concise way to plot the variables 1, 3, and 4 as a function of time.It is possible to omit the list if only one plot is wanted: (2,3) and 4 are respectively equivalent to [(2,3)] and [(0,4)].\nA tuple containing one or several lists will be expanded by associating corresponding elements of the lists with each other:vars = ([1,2,3], [4,5,6])is equivalent tovars = [(1,4), (2,5), (3,6)]andvars = (1, [2,3,4])is equivalent tovars = [(1,2), (1,3), (1,4)]Instead of using integers, one can use the symbols from a ParameterizedFunction. For example, vars=(:x,:y) will replace the symbols with the integer values for components :x and :y.\nn-dimensional groupings are allowed. For example, (1,2,3,4,5) would be a 5-dimensional plot between the associated variables."
},

{
    "location": "basics/plot.html#Complex-Numbers-and-High-Dimensional-Plots-1",
    "page": "Plot Functions",
    "title": "Complex Numbers and High Dimensional Plots",
    "category": "section",
    "text": "The recipe library DimensionalPlotRecipes.jl is provided for extra functionality on high dimensional numbers (complex numbers) and other high dimensional plots. See the README for more details on the extra controls that exist."
},

{
    "location": "basics/plot.html#Timespan-1",
    "page": "Plot Functions",
    "title": "Timespan",
    "category": "section",
    "text": "A plotting timespan can be chosen by the tspan argument in plot. For example:plot(sol,tspan=(0.0,40.0))only plots between t=0.0 and t=40.0. If denseplot=true these bounds will be respected exactly. Otherwise the first point inside and last point inside the interval will be plotted, i.e. no points outside the interval will be plotted."
},

{
    "location": "basics/plot.html#Example-1",
    "page": "Plot Functions",
    "title": "Example",
    "category": "section",
    "text": "using DifferentialEquations, Plots\nfunction lorenz(du,u,p,t)\n du[1] = p[1]*(u[2]-u[1])\n du[2] = u[1]*(p[2]-u[3]) - u[2]\n du[3] = u[1]*u[2] - p[3]*u[3]\nend\n\nu0 = [1., 5., 10.]\ntspan = (0., 100.)\np = (10.0,28.0,8/3)\nprob = ODEProblem(lorenz, u0, tspan,p)\nsol = solve(prob)\nxyzt = plot(sol, plotdensity=10000,lw=1.5)\nxy = plot(sol, plotdensity=10000, vars=(1,2))\nxz = plot(sol, plotdensity=10000, vars=(1,3))\nyz = plot(sol, plotdensity=10000, vars=(2,3))\nxyz = plot(sol, plotdensity=10000, vars=(1,2,3))\nplot(plot(xyzt,xyz),plot(xy, xz, yz, layout=(1,3),w=1), layout=(2,1))(Image: lorenz_plot)An example using the functions:f(x,y,z) = (sqrt(x^2+y^2+z^2),x)\nplot(sol,vars=(f,1,2,3))(Image: norm_plot)"
},

{
    "location": "basics/plot.html#Animations-1",
    "page": "Plot Functions",
    "title": "Animations",
    "category": "section",
    "text": "Using the iterator interface over the solutions, animations can also be generated via the animate(sol) command. One can choose the filename to save to via animate(sol,filename), while the frames per second fps and the density of steps to show every can be specified via keyword arguments. The rest of the arguments will be directly passed to the plot recipe to be handled as normal. For example, we can animate our solution with a larger line-width which saves every 4th frame via:#]add ImageMagick # You may need to install ImageMagick.jl before your first time using it!\n#using ImageMagick # Some installations require using ImageMagick for good animations\nanimate(sol,lw=3,every=4)Please see Plots.jl\'s documentation for more information on the available attributes."
},

{
    "location": "basics/plot.html#Plotting-Without-the-Plot-Recipe-1",
    "page": "Plot Functions",
    "title": "Plotting Without the Plot Recipe",
    "category": "section",
    "text": "What if you don\'t want to use Plots.jl? Odd choice, but that\'s okay! If the differential equation was described by a vector of values, then the solution object acts as an AbstractMatrix sol[i,j] for the ith variable at timepoint j. You can use this to plot solutions. For example, in PyPlot, Gadfly, GR, etc., you can do the following to plot the timeseries:plot(sol.t,sol\')since these plot along the columns, and sol\' has the timeseries along the column. Phase plots can be done similarly, for example:plot(sol[i,:],sol[j,:],sol[k,:])is a 3d phase plot between variables i, j, and k.Notice that this does not use the interpolation. When not using the plot recipe, the interpolation must be done manually. For example:n = 100 #number of timepoints\nts = range(0, stop=1, length=n)\nplot(sol(ts,idxs=i),sol(ts,idxs=j),sol(ts,idxs=k))is the phase space using values 0.001 apart in time."
},

{
    "location": "basics/integrator.html#",
    "page": "Integrator Interface",
    "title": "Integrator Interface",
    "category": "page",
    "text": ""
},

{
    "location": "basics/integrator.html#Integrator-Interface-1",
    "page": "Integrator Interface",
    "title": "Integrator Interface",
    "category": "section",
    "text": "The integrator interface gives one the ability to interactively step through the numerical solving of a differential equation. Through this interface, one can easily monitor results, modify the problem during a run, and dynamically continue solving as one sees fit."
},

{
    "location": "basics/integrator.html#DiffEqBase.step!",
    "page": "Integrator Interface",
    "title": "DiffEqBase.step!",
    "category": "function",
    "text": "step!(integ::DEIntegrator [, dt [, stop_at_tdt]])\n\nPerform one (successful) step on the integrator.\n\nAlternative, if a dt is given, then step! the integrator until there is a temporal difference ≥ dt in integ.t.  When true is passed to the optional third argument, the integrator advances exactly dt.\n\n\n\n\n\n"
},

{
    "location": "basics/integrator.html#DiffEqBase.check_error",
    "page": "Integrator Interface",
    "title": "DiffEqBase.check_error",
    "category": "function",
    "text": "check_error(integrator)\n\nCheck state of integrator and return one of the Return Codes\n\n\n\n\n\n"
},

{
    "location": "basics/integrator.html#DiffEqBase.check_error!",
    "page": "Integrator Interface",
    "title": "DiffEqBase.check_error!",
    "category": "function",
    "text": "check_error!(integrator)\n\nSame as check_error but also set solution\'s return code (integrator.sol.retcode) and run postamble!.\n\n\n\n\n\n"
},

{
    "location": "basics/integrator.html#Initialization-and-Stepping-1",
    "page": "Integrator Interface",
    "title": "Initialization and Stepping",
    "category": "section",
    "text": "To initialize an integrator, use the syntax:integrator = init(prob,alg;kwargs...)The keyword args which are accepted are the same Common Solver Options used by solve and the returned value is an integrator which satisfies typeof(integrator)<:DEIntegrator. One can manually choose to step via the step! command:step!(integrator)which will take one successful step. Additonally:step!(integrator,dt[,stop_at_tdt=false])passing a dt will make the integrator keep stepping until integrator.t+dt, and setting stop_at_tdt=true will add a tstop to force it to step to integrator.t+dtTo check whether or not the integration step was successful, you can call check_error(integrator) which returns one of the Return Codes (RetCodes).This type also implements an iterator interface, so one can step n times (or to the last tstop) using the take iterator:for i in take(integrator,n) endOne can loop to the end by using solve!(integrator) or using the iterator interface:for i in integrator endIn addition, some helper iterators are provided to help monitor the solution. For example, the tuples iterator lets you view the values:for (u,t) in tuples(integrator)\n  @show u,t\nendand the intervals iterator lets you view the full interval:for (tprev,uprev,u,t) in intervals(integrator)\n  @show tprev,t\nendAdditionally, you can make the iterator return specific time points via the TimeChoiceIterator:ts = range(0, stop=1, length=11)\nfor (u,t) in TimeChoiceIterator(integrator,ts)\n  @show u,t\nendLastly, one can dynamically control the \"endpoint\". The initialization simply makes prob.tspan[2] the last value of tstop, and many of the iterators are made to stop at the final tstop value. However, step! will always take a step, and one can dynamically add new values of tstops by modifiying the variable in the options field: add_tstop!(integrator,new_t).Finally, to solve to the last tstop, call solve!(integrator). Doing init and then solve! is equivalent to solve.DiffEqBase.step!\nDiffEqBase.check_error\nDiffEqBase.check_error!"
},

{
    "location": "basics/integrator.html#Handing-Integrators-1",
    "page": "Integrator Interface",
    "title": "Handing Integrators",
    "category": "section",
    "text": "The integrator<:DEIntegrator type holds all of the information for the intermediate solution of the differential equation. Useful fields are:t - time of the proposed step\nu - value at the proposed step\np - user-provided data\nopts - common solver options\nalg - the algorithm associated with the solution\nf - the function being solved\nsol - the current state of the solution\ntprev - the last timepoint\nuprev - the value at the last timepointThe p is the data which is provided by the user as a keyword arg in init. opts holds all of the common solver options, and can be mutated to change the solver characteristics. For example, to modify the absolute tolerance for the future timesteps, one can do:integrator.opts.abstol = 1e-9The sol field holds the current solution. This current solution includes the interpolation function if available, and thus integrator.sol(t) lets one interpolate efficiently over the whole current solution. Additionally, a a \"current interval interpolation function\" is provided on the integrator type via integrator(t). This uses only the solver information from the interval [tprev,t] to compute the interpolation, and is allowed to extrapolate beyond that interval."
},

{
    "location": "basics/integrator.html#DiffEqBase.set_t!",
    "page": "Integrator Interface",
    "title": "DiffEqBase.set_t!",
    "category": "function",
    "text": "set_t!(integrator::DEIntegrator, t::Real)\n\nSet current time point of the integrator to t.\n\n\n\n\n\n"
},

{
    "location": "basics/integrator.html#DiffEqBase.set_u!",
    "page": "Integrator Interface",
    "title": "DiffEqBase.set_u!",
    "category": "function",
    "text": "set_u!(integrator::DEIntegrator, u)\n\nSet current state of the integrator to u.\n\n\n\n\n\n"
},

{
    "location": "basics/integrator.html#DiffEqBase.set_ut!",
    "page": "Integrator Interface",
    "title": "DiffEqBase.set_ut!",
    "category": "function",
    "text": "set_ut!(integrator::DEIntegrator, u, t)\n\nSet current state of the integrator to u and t\n\n\n\n\n\n"
},

{
    "location": "basics/integrator.html#Note-about-mutating-1",
    "page": "Integrator Interface",
    "title": "Note about mutating",
    "category": "section",
    "text": "Be cautious: one should not directly mutate the t and u fields of the integrator. Doing so will destroy the accuracy of the interpolator and can harm certain algorithms. Instead if one wants to introduce discontinuous changes, one should use the Event Handling and Callback Functions. Modifications within a callback affect! surrounded by saves provides an error-free handling of the discontinuity.As low-level alternative to the callbacks, one can use set_t!, set_u! and set_ut! to mutate integrator states.  Note that certain integrators may not have efficient ways to modify u and t.  In such case, set_*! are as inefficient as reinit!.DiffEqBase.set_t!\nDiffEqBase.set_u!\nDiffEqBase.set_ut!"
},

{
    "location": "basics/integrator.html#Integrator-vs-Solution-1",
    "page": "Integrator Interface",
    "title": "Integrator vs Solution",
    "category": "section",
    "text": "The integrator and the solution have very different actions because they have very different meanings. The typeof(sol) <: DESolution type is a type with history: it stores all of the (requested) timepoints and interpolates/acts using the values closest in time. On the other hand, the typeof(integrator)<:DEIntegrator type is a local object. It only knows the times of the interval it currently spans, the current caches and values, and the current state of the solver (the current options, tolerances, etc.). These serve very different purposes:The integrator\'s interpolation can extrapolate, both forward and backward in in time. This is used to estimate events and is internally used for predictions.\nThe integrator is fully mutable upon iteration. This means that every time an iterator affect is used, it will take timesteps from the current time. This means that first(integrator)!=first(integrator) since the integrator will step once to evaluate the left and then step once more (not backtracking). This allows the iterator to keep dynamically stepping, though one should note that it may violate some immutablity assumptions commonly made about iterators.If one wants the solution object, then one can find it in integrator.sol."
},

{
    "location": "basics/integrator.html#Function-Interface-1",
    "page": "Integrator Interface",
    "title": "Function Interface",
    "category": "section",
    "text": "In addition to the type interface, a function interface is provided which allows for safe modifications of the integrator type, and allows for uniform usage throughout the ecosystem (for packages/algorithms which implement the functions). The following functions make up the interface:"
},

{
    "location": "basics/integrator.html#Saving-Controls-1",
    "page": "Integrator Interface",
    "title": "Saving Controls",
    "category": "section",
    "text": "savevalues!(integrator): Adds the current state to the sol."
},

{
    "location": "basics/integrator.html#Caches-1",
    "page": "Integrator Interface",
    "title": "Caches",
    "category": "section",
    "text": "get_tmp_cache(integrator): Returns a tuple of internal cache vectors which are safe to use as temporary arrays. This should be used for integrator interface and callbacks which need arrays to write into in order to be non-allocating. The length of the tuple is dependent on the method.\nfull_cache(integrator):  Returns an iterator over the cache arrays of the method. This can be used to change internal values as needed."
},

{
    "location": "basics/integrator.html#Stepping-Controls-1",
    "page": "Integrator Interface",
    "title": "Stepping Controls",
    "category": "section",
    "text": "u_modified!(integrator,bool): Bool which states whether a change to u occurred, allowing the solver to handle the discontinuity. By default, this is assumed to be true if a callback is used. This will result in the re-calculation of the derivative at t+dt, which is not necessary if the algorithm is FSAL and u does not experience a discontinuous change at the end of the interval. Thus if u is unmodified in a callback, a single call to the derivative calculation can be eliminated by u_modified!(integrator,false).\nget_proposed_dt(integrator):  Gets the proposed dt for the next timestep.\nset_proposed_dt!(integrator,dt):  Sets the proposed dt for the next timestep.\nset_proposed_dt!(integrator,integrator2):  Sets the timestepping of integrator to match that of integrator2. Note that due to PI control and step acceleration this is more than matching the factors in most cases.\nproposed_dt(integrator): Returns the dt of the proposed step.\nterminate!(integrator[, retcode = :Terminated]): Terminates the integrator by emptying tstops. This can be used in events and callbacks to immediately end the solution process.  Optionally, retcode may be specified (see: Return Codes (RetCodes)).\nchange_t_via_interpolation!(integrator,t,modify_save_endpoint=Val{false}): This option lets one modify the current t and changes all of the corresponding values using the local interpolation. If the current solution has already been saved, one can provide the optional value modify_save_endpoint to also modify the endpoint of sol in the same manner.\nadd_tstop!(integrator,t): Adds a tstop at time t.\nadd_saveat!(integrator,t): Adds a saveat time point at t."
},

{
    "location": "basics/integrator.html#Resizing-1",
    "page": "Integrator Interface",
    "title": "Resizing",
    "category": "section",
    "text": "resize!(integrator,k): Resizes the DE to a size k. This chops off the end of the array, or adds blank values at the end, depending on whether k>length(integrator.u).\ndeleteat!(integrator,idxs): Shrinks the ODE by deleting the idxs components.\naddat!(integrator,idxs): Grows the ODE by adding the idxs components. Must be contiguous indices.\nresize_non_user_cache!(integrator,k): Resizes the non-user facing caches to be compatible with a DE of size k. This includes resizing Jacobian caches. Note that in many cases, resize! simply resizes full_cache variables and then calls this function. This finer control is required for some AbstractArray operations.\ndeleteat_non_user_cache!(integrator,idxs): deleteat!s the non-user facing caches at indices idxs. This includes resizing Jacobian caches. Note that in many cases, deleteat! simply deleteat!s full_cache variables and then calls this function. This finer control is required for some AbstractArray operations.\naddat_non_user_cache!(integrator,idxs): addat!s the non-user facing caches at indices idxs. This includes resizing Jacobian caches. Note that in many cases, addat! simply addat!s full_cache variables and then calls this function. This finer control is required for some AbstractArray operations."
},

{
    "location": "basics/integrator.html#Reinit-1",
    "page": "Integrator Interface",
    "title": "Reinit",
    "category": "section",
    "text": "The reinit function lets you restart the integration at a new value. The full function is of the form:reinit!(integrator::ODEIntegrator,u0 = integrator.sol.prob.u0;\n  t0 = integrator.sol.prob.tspan[1], tf = integrator.sol.prob.tspan[2],\n  erase_sol = true,\n  tstops = integrator.opts.tstops_cache,\n  saveat = integrator.opts.saveat_cache,\n  d_discontinuities = integrator.opts.d_discontinuities_cache,\n  reset_dt = (integrator.dtcache == zero(integrator.dt)) && integrator.opts.adaptive,\n  reinit_callbacks = true, initialize_save = true,\n  reinit_cache = true)u0 is the value to start at. The starting time point and end point can be changed via t0 and tf. erase_sol allows one to start with no other values in the solution, or keep the previous solution. tstops, d_discontinuities, and saveat are reset as well, but can be ignored. reset_dt is a boolean for whether to reset the current value of dt using the automatic dt determination algorithm. reinit_callbacks is whether to run the callback initializations again (and initialize_save is for that). reinit_cache is whether to re-run the cache initialization function (i.e. resetting FSAL, not allocating vectors) which should usually be true for correctness.Additionally, once can access auto_dt_reset!(integrator::ODEIntegrator) which will run the auto dt initialization algorithm."
},

{
    "location": "basics/integrator.html#Misc-1",
    "page": "Integrator Interface",
    "title": "Misc",
    "category": "section",
    "text": "get_du(integrator): Returns the derivative at t.\nget_du!(out,integrator): Write the current derivative at t into out.\ncheck_error(integrator): Checks error conditions and updates the retcode."
},

{
    "location": "basics/integrator.html#Note-1",
    "page": "Integrator Interface",
    "title": "Note",
    "category": "section",
    "text": "Note that not all of these functions will be implemented for every algorithm. Some have hard limitations. For example, Sundials.jl cannot resize problems. When a function is not limited, an error will be thrown."
},

{
    "location": "basics/integrator.html#Additional-Options-1",
    "page": "Integrator Interface",
    "title": "Additional Options",
    "category": "section",
    "text": "The following options can additionally be specified in init (or be mutated in the opts) for further control of the integrator:advance_to_tstop: This makes step! continue to the next value in tstop.\nstop_at_next_tstop: This forces the iterators to stop at the next value of tstop.For example, if one wants to iterate but only stop at specific values, one can choose:integrator = init(prob,Tsit5();dt=1//2^(4),tstops=[0.5],advance_to_tstop=true)\nfor (u,t) in tuples(integrator)\n  @test t ∈ [0.5,1.0]\nendwhich will only enter the loop body at the values in tstops (here, prob.tspan[2]==1.0 and thus there are two values of tstops which are hit). Addtionally, one can solve! only to 0.5 via:integrator = init(prob,Tsit5();dt=1//2^(4),tstops=[0.5])\nintegrator.opts.stop_at_next_tstop = true\nsolve!(integrator)"
},

{
    "location": "basics/integrator.html#Plot-Recipe-1",
    "page": "Integrator Interface",
    "title": "Plot Recipe",
    "category": "section",
    "text": "Like the DESolution type, a plot recipe is provided for the DEIntegrator type. Since the DEIntegrator type is a local state type on the current interval, plot(integrator) returns the solution on the current interval. The same options for the plot recipe are provided as for sol, meaning one can choose variables via the vars keyword argument, or change the plotdensity / turn on/off denseplot.Additionally, since the integrator is an iterator, this can be used in the Plots.jl animate command to iteratively build an animation of the solution while solving the differential equation.For an example of manually chaining together the iterator interface and plotting, one should try the following:using DifferentialEquations, DiffEqProblemLibrary, Plots\n\n# Linear ODE which starts at 0.5 and solves from t=0.0 to t=1.0\nprob = ODEProblem((u,p,t)->1.01u,0.5,(0.0,1.0))\n\nusing Plots\nintegrator = init(prob,Tsit5();dt=1//2^(4),tstops=[0.5])\npyplot(show=true)\nplot(integrator)\nfor i in integrator\n  display(plot!(integrator,vars=(0,1),legend=false))\nend\nstep!(integrator); plot!(integrator,vars=(0,1),legend=false)\nsavefig(\"iteratorplot.png\")(Image: Iterator Plot)"
},

{
    "location": "basics/problem.html#",
    "page": "Problem interface",
    "title": "Problem interface",
    "category": "page",
    "text": ""
},

{
    "location": "basics/problem.html#Problem-interface-1",
    "page": "Problem interface",
    "title": "Problem interface",
    "category": "section",
    "text": ""
},

{
    "location": "basics/problem.html#Functional-and-Condensed-Problem-Inputs-1",
    "page": "Problem interface",
    "title": "Functional and Condensed Problem Inputs",
    "category": "section",
    "text": "Note that the initial condition can be written as a function of parameters and initial time:u0(p,t0)and be resolved before going to the solver. Additionally, the initial condition can be a distribution from Distributions.jl, in which case a sample initial condition will be taken each time init or solve is called.In addition, tspan supports the following forms. The single value form t is equivalent to (zero(t),t). The functional form is allowed:tspan(p)which outputs a tuple."
},

{
    "location": "basics/problem.html#Examples-1",
    "page": "Problem interface",
    "title": "Examples",
    "category": "section",
    "text": "prob = ODEProblem((u,p,t)->u,(p,t0)->p[1],(p)->(0.0,p[2]),(2.0,1.0))\nusing Distributions\nprob = ODEProblem((u,p,t)->u,(p,t)->Normal(p,1),(0.0,1.0),1.0)"
},

{
    "location": "basics/problem.html#Lower-Level-__init-and-__solve-1",
    "page": "Problem interface",
    "title": "Lower Level __init and __solve",
    "category": "section",
    "text": "At the high level, known problematic problems will emit warnings before entering the solver to better clarify the error to the user. The following cases are checked if the solver is adaptive:Integer times warn\nDual numbers must be in the initial conditions and timespans\nMeasurements.jl values must be in the initial conditions and timespansIf there is an exception to these rules, please file an issue. If one wants to go around the high level solve interface and its warnings, one can call __init or __solve instead."
},

{
    "location": "basics/problem.html#Modification-of-problem-types-1",
    "page": "Problem interface",
    "title": "Modification of problem types",
    "category": "section",
    "text": "Problem-related types in DifferentialEquations.jl are immutable.  This helps, e.g., parallel solvers to efficiently handle problem types.However, you may want to modify the problem after it is created.  For example, to simulate it for longer timespan.  It can be done by the remake function:prob1 = ODEProblem((u,p,t) -> u/2, 1.0, (0.0,1.0))\nprob2 = remake(prob1; tspan=(0.0,2.0))A general syntax of remake ismodified_problem = remake(original_problem;\n  field_1 = value_1,\n  field_2 = value_2,\n  ...\n)where field_N and value_N are renamed to appropriate field names and new desired values."
},

{
    "location": "basics/faq.html#",
    "page": "Frequently Asked Questions",
    "title": "Frequently Asked Questions",
    "category": "page",
    "text": ""
},

{
    "location": "basics/faq.html#Frequently-Asked-Questions-1",
    "page": "Frequently Asked Questions",
    "title": "Frequently Asked Questions",
    "category": "section",
    "text": "This page is a compilation of frequently asked questions and answers."
},

{
    "location": "basics/faq.html#Performance-1",
    "page": "Frequently Asked Questions",
    "title": "Performance",
    "category": "section",
    "text": ""
},

{
    "location": "basics/faq.html#Do-you-support-GPUs?-Multithreading?-Distributed-computation?-1",
    "page": "Frequently Asked Questions",
    "title": "Do you support GPUs? Multithreading? Distributed computation?",
    "category": "section",
    "text": "Yes. The *DiffEq.jl libraries (OrdinaryDiffEq.jl, StochasticDiffEq.jl, and DelayDiffEq.jl) are all written to be generic to the array and number types. This means they will adopt the implementation that is given by the array type. The in-place algorithms internally utilize Julia\'s broadcast (with some exceptions due to a Julia bug for now, see this issue) and Julia\'s mul! in-place matrix multiplication function. The out-of-place algorithms utilize standard arithmetical functions. Both additionally utilize the user\'s norm specified via the common interface options and, if a stiff solver, ForwardDiff/DiffEqDiffTools for the Jacobian calculation, and Base linear factorizations for the linear solve. For your type, you may likely need to give a better form of the norm, Jacobian, or linear solve calculations to fully utilize parallelism.GPUArrays.jl (CuArrays.jl), ArrayFire.jl, DistributedArrays.jl have been tested and work in various forms, where the last one is still not recommended for common use yet.The next question is whether it matters. Generally, your system has to be large for parallelism to matter. Using a multithreaded array for broadcast we find helpful around N>1000, though the Sundials manual says N>100,000. For high order Runge-Kutta methods it\'s likely lower than the Sundials estimate because of more operations packed into each internal step, but as always that will need more benchmarks to be precise and will depend on the problem being solved. GPUs generally require some intensive parallel operation in the user\'s f function to be viable, for example a matrix multiplication for a stencil computation in a PDE. If you\'re simply solving some ODE element-wise on a big array it likely won\'t do much or it will slow things down just due to how GPUs work. DistributedArrays require parallel linear solves to really matter, and thus are only recommended when you have a problem that cannot fit into memory or are using a stiff solver with a Krylov method for the linear solves."
},

{
    "location": "basics/faq.html#My-ODE-is-solving-really-slow...-what-do-I-do?-1",
    "page": "Frequently Asked Questions",
    "title": "My ODE is solving really slow... what do I do?",
    "category": "section",
    "text": "First, check for bugs. These solvers go through a ton of convergence tests and so if there\'s a solver issue, it\'s either just something to do with how numerical methods work or it\'s a user-error (generally the latter, though check the later part of the FAQ on normal numerical errors). User-errors in the f function causing a divergence of the solution is the most common reason for reported slow codes.If you have no bugs, great! The standard tricks for optimizing Julia code then apply. What you want to do first is make sure your function does not allocate. If your system is small (<=100 ODEs/SDEs/DDEs/DAEs?), then you should set your system up to use StaticArrays.jl. This is demonstrated in the ODE tutorial with static matrices. Static vectors/arrays are stack-allocated, and thus creating new arrays is free and the compiler doesn\'t have to heap-allocate any of the temporaries (that\'s the expensive part!). These have specialized super fast dispatches for arithmetic operations and extra things like LU-factorizations, and thus they are preferred when possible. However, they lose efficiency if they grow too large.For anything larger, you should use the in-place syntax f(du,u,p,t) and make sure that your function doesn\'t allocate. Assuming you know of a u0, you should be able to do:du = similar(u0)\n@time f(du,u0,p,t)and see close to zero allocations and close to zero memory allocated. If you see more, then you might have a type-instability or have temporary arrays. To find type-instabilities, you should do:@code_warntype f(du,u,p,t)and read the printout to see if there\'s any types that aren\'t inferred by the compiler, and fix them. If you have any global variables, you should make them const. As for allocations, some common things that allocate are:Array slicing, like u[1:5]. Instead, use @view u[1:5]\nMatrix multiplication with *. Instead of A*b, use A_mul_B!(c,A,b) for some pre-allocated cache vector c.\nNon-broadcasted expressions. Every expression on arrays should .= into another array, or it should be re-written to loop and do computations with scalar (or static array) values.For an example of optimizing a function resulting from a PDE discretization, see this blog post."
},

{
    "location": "basics/faq.html#The-stiff-solver-takes-forever-to-take-steps-for-my-PDE-discretization!-Why?-1",
    "page": "Frequently Asked Questions",
    "title": "The stiff solver takes forever to take steps for my PDE discretization! Why?",
    "category": "section",
    "text": "The solvers for stiff solvers require solving a nonlinear equation each step. In order to do so, they have to do a few Newton steps. By default, these methods assume that the Jacobian is dense, automatically calculate the Jacobian for you, and do a dense factorization. However, in many cases you may want to use alternatives that are more tuned for your problem.First of all, when available, it\'s recommended that you pass a function for computing your Jacobian. This is discussed in the performance overloads section. Jacobians are especially helpful for Rosenbrock methods.Secondly, if your Jacobian isn\'t dense, you shouldn\'t use a dense Jacobian! In the Sundials algorithm you can set linear_solver=:Band for banded Jacobians for example. More support is coming for this soon.But lastly, you shouldn\'t use a dense factorization for large sparse matrices. Instead, if you\'re using  a *DiffEq library you should specify a linear solver. For Sundials.jl, you should change the linear_solver option. See the ODE solve Sundials portion for details on that. Right now, Sundials.jl is the recommended method for stiff problems with large sparse Jacobians. linear_solver=:Band should be used if your Jacobian is banded and you can specify the band sizes. If you only know the Jacobian is sparse, linear_solver=:GMRES is a good option. Once again, a good reference for how to handle PDE discretizations can be found at this blog post."
},

{
    "location": "basics/faq.html#Complicated-Models-1",
    "page": "Frequently Asked Questions",
    "title": "Complicated Models",
    "category": "section",
    "text": ""
},

{
    "location": "basics/faq.html#Can-I-switch-my-ODE-function-in-the-middle-of-integration?-1",
    "page": "Frequently Asked Questions",
    "title": "Can I switch my ODE function in the middle of integration?",
    "category": "section",
    "text": "There are a few ways to do this. The simplest way is to just have a parameter to switch between the two. For example:function f(du,u,p,t)\n  if p == 0\n    du[1] = 2u[1]\n  else\n    du[1] = -2u[1]\n  end\n  du[2] = -u[2]\nendThen in a callback you can make the affect! function modify integrator.prob.p. For example, we can make it change when u[2]<0.5 via:condition(t,u,integrator) = u[2] - 0.5\naffect!(integrator) = integrator.prob.p = 1Then it will change betweeen the two ODE choices for du1 at that moment. Another way to do this is to make the ODE functions all be the same type via FunctionWrappers.jl, but that is unnecessary. With the way that modern processors work, there exists branch prediction and thus execution of a conditional is free if it\'s predictable which branch will be taken. In this case, almost every call to f takes the p==0 route until the callback, at which point it is almost always the else route. Therefore the processor will effectively get rid of the computational cost associated with this, so you\'re likely over-optimizing if you\'re going further (unless this change happens every step, but even then this is probably the cheapest part of the computation...)."
},

{
    "location": "basics/faq.html#Numerical-Error-1",
    "page": "Frequently Asked Questions",
    "title": "Numerical Error",
    "category": "section",
    "text": ""
},

{
    "location": "basics/faq.html#The-solver-doesn\'t-obey-physical-law-X-(e.g.-conservation-of-energy)-1",
    "page": "Frequently Asked Questions",
    "title": "The solver doesn\'t obey physical law X (e.g. conservation of energy)",
    "category": "section",
    "text": "Yes, this is because the numerical solution of the ODE is not the exact solution. There are a few ways that you can handle this problem. One way is to get a more exact solution. Thus instead ofsol = solve(prob,alg)usesol = solve(prob,alg,abstol=1e-10,reltol=1e-10)Of course, there\'s always a tradeoff between accuracy and efficiency, so play around to find out what\'s right for your problem.Another thing you can do is use a callback. There are some premade callbacks in the callback library which handle these sorts of things like projecting to manifolds and preserving positivity."
},

{
    "location": "basics/faq.html#The-symplectic-integrator-doesn\'t-conserve-energy?-1",
    "page": "Frequently Asked Questions",
    "title": "The symplectic integrator doesn\'t conserve energy?",
    "category": "section",
    "text": "Yes, symplectic integrators do not exactly conserve energy. It is a common misconception that they do. What symplectic integrators actually do is solve for a trajectory which rests on a symplectic manifold that is perturbed from the true solution\'s manifold by the truncation error. This means that symplectic integrators do not experience (very much) long time drift, but their orbit is not exactly the same as the true solution in phase space and thus you will see differences in energy that tend to look periodic. There is a small drift which grows linearly and is related to floating point error, but this drift is much less than standard methods. This is why symplectic methods are recommended for long time integration.For conserving energy, there are a few things you can do. First of all, the energy error is related to the integration error, so simply solving with higher accuracy will reduce the error. The results in the DiffEqBenchmarks show that using a DPRKN method with low tolerance can be a great choice. Another thing you can do is use the ManifoldProjection callback from the callback library."
},

{
    "location": "basics/faq.html#How-do-I-get-to-zero-error?-1",
    "page": "Frequently Asked Questions",
    "title": "How do I get to zero error?",
    "category": "section",
    "text": "You can\'t. For floating point numbers, you shouldn\'t use below abstol=1e-14 and reltol=1e-14. If you need lower than that, use arbitrary precision numbers like BigFloats or ArbFloats.jl."
},

{
    "location": "basics/faq.html#Autodifferentiation-and-Dual-Numbers-1",
    "page": "Frequently Asked Questions",
    "title": "Autodifferentiation and Dual Numbers",
    "category": "section",
    "text": ""
},

{
    "location": "basics/faq.html#Are-the-native-Julia-solvers-compatible-with-autodifferentiation?-1",
    "page": "Frequently Asked Questions",
    "title": "Are the native Julia solvers compatible with autodifferentiation?",
    "category": "section",
    "text": "Yes! Take a look at the sensitivity analysis page for more details.If the algorithm does not have differentiation of parameter-depedendent events,  then you simply need to make the initial condition have elements of Dual numbers.  If the algorithm uses Dual numbers, you need to make sure that time is also given by Dual numbers. To show this in action, let\'s say we want to find the Jacobian of solution of the Lotka-Volterra equation at t=10 with respect to the parameters.function func(du,u,p,t)\n  du[1] = p[1] * u[1] - p[2] * u[1]*u[2]\n  du[2] = -3 * u[2] + u[1]*u[2]\nend\nfunction f(p)\n  prob = ODEProblem(func,eltype(p).([1.0,1.0]),eltype(p).((0.0,10.0)),p)\n  # Lower tolerances to show the methods converge to the same value\n  solve(prob,Tsit5(),save_everystep=false,abstol=1e-12,reltol=1e-12)[end]\nendThis function takes in new parameters and spits out the solution at the end. We make the inital condition eltype(p).([1.0,1.0]) so that way it\'s typed to be Dual numbers whenever p is an array of Dual numbers, and we do the same for the timespan just to show what you\'d do if there was parameters-dependent events.  Then we can take the Jacobian via ForwardDiff.jl:using ForwardDiff\nForwardDiff.jacobian(f,[1.5,1.0])\n\n2×2 Array{Float64,2}:\n  2.20214   0.189782\n -6.2273   -0.700188and compare it to Calculus.jl:Calculus.jacobian(f,[1.5,1.0],:central)\n\n2×2 Array{Float64,2}:\n  2.20214   0.189782\n -6.2273   -0.700188"
},

{
    "location": "basics/faq.html#I-get-Dual-number-errors-when-I-solve-my-ODE-with-Rosenbrock-or-SDIRK-methods...?-1",
    "page": "Frequently Asked Questions",
    "title": "I get Dual number errors when I solve my ODE with Rosenbrock or SDIRK methods...?",
    "category": "section",
    "text": "This is because you\'re using a cache which is not compatible with autodifferentiaion via ForwardDiff.jl. For example, if we use the ODE function:using LinearAlgebra, OrdinaryDiffEq\nfunction foo(du, u, (A, tmp), t)\n    mul!(tmp, A, u)\n    @. du = u + tmp\n    nothing\nend\nprob = ODEProblem(foo, ones(5, 5), (0., 1.0), (ones(5,5), zeros(5,5)))\nsolve(prob, Rosenbrock23())Here we use a cached temporary array in order to avoid the allocations of matrix multiplication. When autodifferentiation occurs, the element type of u is Dual numbers, so A*u produces Dual numbers, so the error arises when it tries to write into tmp. There are two ways to avoid this. The first way, the easy way, is to just turn off autodifferentiation with the autodiff=false option in the solver. Every solver which uses autodifferentiation has this option. Thus we\'d solve this with:prob = ODEProblem(f,rand(4),(0.0,1.0))\nsol = solve(prob,Rosenbrock23(autodiff=false))and it will use a numerical differentiation fallback (DiffEqDiffTools.jl) to calculate Jacobians.We could use get_tmp and dualcache functions from DiffEqBase to solve this issue, e.g.,using LinearAlgebra, OrdinaryDiffEq\nusing DiffEqBase: get_tmp, dualcache\nfunction foo(du, u, (A, tmp), t)\n    tmp = DiffEqBase.get_tmp(tmp, u)\n    mul!(tmp, A, u)\n    @. du = u + tmp\n    nothing\nend\nprob = ODEProblem(foo, ones(5, 5), (0., 1.0), (ones(5,5), DiffEqBase.dualcache(zeros(5,5))))\nsolve(prob, TRBDF2())"
},

{
    "location": "basics/compatibility_chart.html#",
    "page": "Solver Compatibility Chart",
    "title": "Solver Compatibility Chart",
    "category": "page",
    "text": ""
},

{
    "location": "basics/compatibility_chart.html#Solver-Compatibility-Chart-1",
    "page": "Solver Compatibility Chart",
    "title": "Solver Compatibility Chart",
    "category": "section",
    "text": "This chart is for documenting the compatibility of the component solver packages to the common interface. An x means that the option is implemented or the add-on functionality will work with the given solver. A blank means that the option has not been implemented or that a given add-on has not been tested with a given package. If there are any errors in this chart, please file an issue or submit a pull-request.Option OrdinaryDiffEq.jl Sundials.jl ODE.jl ODEInterface.jl LSODA.jl StochasticDiffEq.jl DelayDiffEq.jl DASKR.jl DASSL.jl\nNonlinear Dense (continuous) output x x    x x x \nTolerance control x x x x x x x x x\nAdvanced stepsize control x 0  x 0 x x 0 \nMass Matrices^ x 0  x 0 x x 0 \nAnalytical Jacobians^† x x  x  x x x \nGeneral Performance Overloads^† x 0  0 0 x x 0 \ninternalnorm x 0 x 0 0 x x 0 \nInitial dt x x x x  x x x \nsave_everystep x x x x x x x x \nsaveat x x x x x x x x \ntstops x x  0  x x x \nd_discontinuities x   0  x x  \nisoutofdomain x  x   x x  \nAllows reverse time direction x x x x x x x  \nUnitful numbers x 0  0 0  x 0 \nArbitrary dimension arrays x x x x x x x x x\nComplex numbers p     x p  \nArbitrary precision x 0 x 0 0 x x 0 x\nApproxFun types x 0  0 0  x 0 \nProgress monitoring x     x x  \nIntegrator interface x x  0  x x  \nResizability x 0  0 0 x x 0 \nCache iterator x 0  0 0 x x 0 \nCan choose linear solvers x s    x x s x\nCan choose nonlinear solvers x 0  0 0 x x 0 x\nCan use out of place natively x 0 x 0 0 x x 0 x\nCan use inplace natively x x  x x x x x \nCompatible with DiffEqDevTools x x x x x x x x \nCompatible with ParameterizedFunctions x x x x x x x x \nContinuous Callbacks x x  x  x x  x\nDiscrete Callbacks x x  x  x x  \nMonte Carlo Simulations x x x x x x x x \nParameter Estimation x n n n n x x n x\nParameter Sensitivity Analysis x x x x x  x  \nPlotting and solution handling x x x x x x x x xx: Full compatibility\np: Partial compatibility, only in nonstiff methods unless the Jacobian is provided.\nn: General compatibility, but not compatible with routines which. require being able to autodifferentiate through the entire solver.\n0: Not possible. This is generally due to underlying inflexibility in a wrapped library.\ns: Special, Sundials has its own linear solver choices.\n^: Only stiff (implicit) methods.\n†: For packages with compatibility, no warning is given when a specific algorithm does not need to use this feature.All blank spaces are possible future additions."
},

{
    "location": "types/discrete_types.html#",
    "page": "Discrete Problems",
    "title": "Discrete Problems",
    "category": "page",
    "text": ""
},

{
    "location": "types/discrete_types.html#Discrete-Problems-1",
    "page": "Discrete Problems",
    "title": "Discrete Problems",
    "category": "section",
    "text": ""
},

{
    "location": "types/discrete_types.html#Mathematical-Specification-of-a-Discrete-Problem-1",
    "page": "Discrete Problems",
    "title": "Mathematical Specification of a Discrete Problem",
    "category": "section",
    "text": "To define an Discrete Problem, you simply need to give the function f and the initial condition u₀ which define a function map:u_n+1 = f(upt_n+1)f should be specified as f(u,p,t) (or in-place as f(du,u,p,t)), and u₀ should be an AbstractArray (or number) whose geometry matches the desired geometry of u. Note that we are not limited to numbers or vectors for u₀; one is allowed to provide u₀ as arbitrary matrices / higher dimension tensors as well. t_{n+1} is the current time at which the map is applied. For a FunctionMap with defaults, t_n = t0 + n*dt (with dt=1 being the default). For continuous-time Markov chains this is the time at which the change is occuring.Note that if the discrete solver is set to have scale_by_time=true, then the problem is interpreted as the map:u_n+1 = u_n + dt f(upt_n)"
},

{
    "location": "types/discrete_types.html#Problem-Type-1",
    "page": "Discrete Problems",
    "title": "Problem Type",
    "category": "section",
    "text": ""
},

{
    "location": "types/discrete_types.html#Constructors-1",
    "page": "Discrete Problems",
    "title": "Constructors",
    "category": "section",
    "text": "DiscreteProblem{isinplace}(f::ODEFunction,u0,tspan) : Defines the discrete problem with the specified functions.\nDiscreteProblem{isinplace}(f,u0,tspan) : Defines the discrete problem with the specified functions.\nDiscreteProblem{isinplace}(u0,tspan) : Defines the discrete problem with the identity map.For specifying Jacobians and mass matrices, see the DiffEqFunctions page."
},

{
    "location": "types/discrete_types.html#Fields-1",
    "page": "Discrete Problems",
    "title": "Fields",
    "category": "section",
    "text": "f: The function in the map.\nu0: The initial condition.\ntspan: The timespan for the problem.\ncallback: A callback to be applied to every solver which uses the problem. Defaults to a black CallbackSet, which will have no effect."
},

{
    "location": "types/discrete_types.html#Note-About-Timing-1",
    "page": "Discrete Problems",
    "title": "Note About Timing",
    "category": "section",
    "text": "Note that if no dt and not tstops is given, it\'s assumed that dt=1 and thus tspan=(0,n) will solve for n iterations. If in the solver dt is given, then the number of iterations will change. And if tstops is not empty, the solver will revert to the standard behavior of fixed timestep methods, which is \"step to each tstop\"."
},

{
    "location": "types/ode_types.html#",
    "page": "ODE Problems",
    "title": "ODE Problems",
    "category": "page",
    "text": ""
},

{
    "location": "types/ode_types.html#ODE-Problems-1",
    "page": "ODE Problems",
    "title": "ODE Problems",
    "category": "section",
    "text": ""
},

{
    "location": "types/ode_types.html#Mathematical-Specification-of-an-ODE-Problem-1",
    "page": "ODE Problems",
    "title": "Mathematical Specification of an ODE Problem",
    "category": "section",
    "text": "To define an ODE Problem, you simply need to give the function f and the initial condition u₀ which define an ODE:fracdudt = f(upt)f should be specified as f(u,p,t) (or in-place as f(du,u,p,t)), and u₀ should be an AbstractArray (or number) whose geometry matches the desired geometry of u. Note that we are not limited to numbers or vectors for u₀; one is allowed to provide u₀ as arbitrary matrices / higher dimension tensors as well."
},

{
    "location": "types/ode_types.html#Problem-Type-1",
    "page": "ODE Problems",
    "title": "Problem Type",
    "category": "section",
    "text": ""
},

{
    "location": "types/ode_types.html#Constructors-1",
    "page": "ODE Problems",
    "title": "Constructors",
    "category": "section",
    "text": "ODEProblem(f::ODEFunction,u0,tspan,callback=CallbackSet())\nODEProblem{isinplace}(f,u0,tspan,callback=CallbackSet()) : Defines the ODE with the specified functions. isinplace optionally sets whether the function is inplace or not. This is determined automatically, but not inferred.For specifying Jacobians and mass matrices, see the DiffEqFunctions page."
},

{
    "location": "types/ode_types.html#Fields-1",
    "page": "ODE Problems",
    "title": "Fields",
    "category": "section",
    "text": "f: The function in the ODE.\nu0: The initial condition.\ntspan: The timespan for the problem.\ncallback: A callback to be applied to every solver which uses the problem. Defaults to nothing."
},

{
    "location": "types/ode_types.html#Example-Problems-1",
    "page": "ODE Problems",
    "title": "Example Problems",
    "category": "section",
    "text": "Example problems can be found in DiffEqProblemLibrary.jl.To use a sample problem, such as prob_ode_linear, you can do something like:#] add DiffEqProblemLibrary\nusing DiffEqProblemLibrary\nprob = prob_ode_linear\nsol = solve(prob)DiffEqProblemLibrary.prob_ode_linear\nDiffEqProblemLibrary.prob_ode_2Dlinear\nDiffEqProblemLibrary.prob_ode_bigfloatlinear\nDiffEqProblemLibrary.prob_ode_bigfloat2Dlinear\nDiffEqProblemLibrary.prob_ode_large2Dlinear\nDiffEqProblemLibrary.prob_ode_2Dlinear_notinplace\nDiffEqProblemLibrary.prob_ode_threebody\nDiffEqProblemLibrary.prob_ode_pleides\nDiffEqProblemLibrary.prob_ode_vanderpol\nDiffEqProblemLibrary.prob_ode_vanderpol_stiff\nDiffEqProblemLibrary.prob_ode_rober\nDiffEqProblemLibrary.prob_ode_rigidbody"
},

{
    "location": "types/dynamical_types.html#",
    "page": "Dynamical, Hamiltonian and 2nd Order ODE Problems",
    "title": "Dynamical, Hamiltonian and 2nd Order ODE Problems",
    "category": "page",
    "text": ""
},

{
    "location": "types/dynamical_types.html#Dynamical,-Hamiltonian-and-2nd-Order-ODE-Problems-1",
    "page": "Dynamical, Hamiltonian and 2nd Order ODE Problems",
    "title": "Dynamical, Hamiltonian and 2nd Order ODE Problems",
    "category": "section",
    "text": "Dynamical ordinary differential equations, such as those arising from the definition of a Hamiltonian system or a second order ODE, have a special structure that can be utilized in the solution of the differential equation. On this page we describe how to define second order differential equations for their efficient numerical solution."
},

{
    "location": "types/dynamical_types.html#Mathematical-Specification-of-a-Dynamical-ODE-Problem-1",
    "page": "Dynamical, Hamiltonian and 2nd Order ODE Problems",
    "title": "Mathematical Specification of a Dynamical ODE Problem",
    "category": "section",
    "text": "These algorithms require a Partitioned ODE of the form:fracdvdt = f_1(ut) \nfracdudt = f_2(v) This is a Partitioned ODE partitioned into two groups, so the functions should be specified as f1(dv,v,u,p,t) and f2(du,v,u,p,t) (in the inplace form), where f1 is independent of v (unless specified by the solver), and f2 is independent of u and t. This includes discretizations arising from SecondOrderODEProblems where the velocity is not used in the acceleration function, and Hamiltonians where the potential is (or can be) time-dependent but the kinetic energy is only dependent on v.Note that some methods assume that the integral of f2 is a quadratic form. That means that f2=v\'*M*v, i.e. int f_2 = frac12 m v^2, giving du = v. This is equivalent to saying that the kinetic energy is related to v^2. The methods which require this assumption will lose accuracy if this assumption is violated. Methods listed make note of this requirement with \"Requires quadratic kinetic energy\"."
},

{
    "location": "types/dynamical_types.html#Constructor-1",
    "page": "Dynamical, Hamiltonian and 2nd Order ODE Problems",
    "title": "Constructor",
    "category": "section",
    "text": "DynamicalODEProblem(f::DynamicalODEFunction,v0,u0,tspan,callback=CallbackSet())\nDynamicalODEProblem{isinplace}(f1,f2,v0,u0,tspan,callback=CallbackSet())Defines the ODE with the specified functions. isinplace optionally sets whether the function is inplace or not. This is determined automatically, but not inferred."
},

{
    "location": "types/dynamical_types.html#Fields-1",
    "page": "Dynamical, Hamiltonian and 2nd Order ODE Problems",
    "title": "Fields",
    "category": "section",
    "text": "f1 and f2: The functions in the ODE.\nv0 and u0: The initial conditions.\ntspan: The timespan for the problem.\ncallback: A callback to be applied to every solver which uses the problem. Defaults to nothing."
},

{
    "location": "types/dynamical_types.html#Mathematical-Specification-of-a-2nd-Order-ODE-Problem-1",
    "page": "Dynamical, Hamiltonian and 2nd Order ODE Problems",
    "title": "Mathematical Specification of a 2nd Order ODE Problem",
    "category": "section",
    "text": "To define a 2nd Order ODE Problem, you simply need to give the function f and the initial condition u₀ which define an ODE:u = f(uupt)f should be specified as f(du,u,p,t) (or in-place as f(ddu,du,u,p,t)), and u₀ should be an AbstractArray (or number) whose geometry matches the desired geometry of u. Note that we are not limited to numbers or vectors for u₀; one is allowed to provide u₀ as arbitrary matrices / higher dimension tensors as well.From this form, a dynamical ODE:v = f(vupt) \nu = v is generated."
},

{
    "location": "types/dynamical_types.html#Constructors-1",
    "page": "Dynamical, Hamiltonian and 2nd Order ODE Problems",
    "title": "Constructors",
    "category": "section",
    "text": "SecondOrderODEProblem{isinplace}(f,du0,u0,tspan,callback=CallbackSet())Defines the ODE with the specified functions."
},

{
    "location": "types/dynamical_types.html#Fields-2",
    "page": "Dynamical, Hamiltonian and 2nd Order ODE Problems",
    "title": "Fields",
    "category": "section",
    "text": "f: The function for the second derivative.\ndu0: The initial derivative.\nu0: The initial condition.\ntspan: The timespan for the problem.\ncallback: A callback to be applied to every solver which uses the problem. Defaults to nothing."
},

{
    "location": "types/dynamical_types.html#Hamiltonian-Problems-1",
    "page": "Dynamical, Hamiltonian and 2nd Order ODE Problems",
    "title": "Hamiltonian Problems",
    "category": "section",
    "text": "HamiltonianProblems are provided by DiffEqPhysics.jl and provide an easy way to define equations of motion from the corresponding Hamiltonian. To define a HamiltonianProblem one only needs to specify the Hamiltonian:H(pq)and autodifferentiation (via ForwardDiff.jl) will create the appropriate equations."
},

{
    "location": "types/dynamical_types.html#Constructors-2",
    "page": "Dynamical, Hamiltonian and 2nd Order ODE Problems",
    "title": "Constructors",
    "category": "section",
    "text": "HamiltonianProblem{T}(H,p0,q0,tspan;kwargs...)"
},

{
    "location": "types/dynamical_types.html#Fields-3",
    "page": "Dynamical, Hamiltonian and 2nd Order ODE Problems",
    "title": "Fields",
    "category": "section",
    "text": "H: The Hamiltonian H(p,q,params) which returns a scalar.\np0: The initial momentums.\nq0: The initial positions.\ntspan: The timespan for the problem.\ncallback: A callback to be applied to every solver which uses the problem. Defaults to nothing."
},

{
    "location": "types/split_ode_types.html#",
    "page": "Split ODE Problems",
    "title": "Split ODE Problems",
    "category": "page",
    "text": ""
},

{
    "location": "types/split_ode_types.html#Split-ODE-Problems-1",
    "page": "Split ODE Problems",
    "title": "Split ODE Problems",
    "category": "section",
    "text": ""
},

{
    "location": "types/split_ode_types.html#Mathematical-Specification-of-a-Split-ODE-Problem-1",
    "page": "Split ODE Problems",
    "title": "Mathematical Specification of a Split ODE Problem",
    "category": "section",
    "text": "To define a SplitODEProblem, you simply need to give a two functions functions f_1 and f_2 along with an initial condition u₀ which define an ODE:fracdudt =  f_1(upt) + f_2(upt)f should be specified as f(u,p,t) (or in-place as f(du,u,p,t)), and u₀ should be an AbstractArray (or number) whose geometry matches the desired geometry of u. Note that we are not limited to numbers or vectors for u₀; one is allowed to provide u₀ as arbitrary matrices / higher dimension tensors as well.Many splits are at least partially linear. That is the equation:fracdudt =  Au + f_2(upt)For how to define a linear function A, see the documentation for the DiffEqOperators."
},

{
    "location": "types/split_ode_types.html#Constructors-1",
    "page": "Split ODE Problems",
    "title": "Constructors",
    "category": "section",
    "text": "SplitODEProblem(f::SplitFunction,u0,tspan,p=nothing;kwargs...)\nSplitODEProblem{isinplace}(f1,f2,u0,tspan,p=nothing;kwargs...)The isinplace parameter can be omitted and will be determined using the signature of f2. Note that both f1 and f2 should support the in-place style if isinplace is true or they should both support the out-of-place style if isinplace is false. You cannot mix up the two styles.Under the hood, a SplitODEProblem is just a regular ODEProblem whose f is a SplitFunction. Therefore you can solve a SplitODEProblem using the same solvers for ODEProblem. For solvers dedicated to split problems, see Split ODE Solvers.For specifying Jacobians and mass matrices, see the DiffEqFunctions page."
},

{
    "location": "types/split_ode_types.html#Fields-1",
    "page": "Split ODE Problems",
    "title": "Fields",
    "category": "section",
    "text": "f1, f2: The functions in the ODE.\nu0: The initial condition.\ntspan: The timespan for the problem.\np: The parameters for the problem.\ncallback: A callback to be applied to every solver which uses the problem. Defaults to nothing."
},

{
    "location": "types/steady_state_types.html#",
    "page": "Steady State Problems",
    "title": "Steady State Problems",
    "category": "page",
    "text": ""
},

{
    "location": "types/steady_state_types.html#Steady-State-Problems-1",
    "page": "Steady State Problems",
    "title": "Steady State Problems",
    "category": "section",
    "text": ""
},

{
    "location": "types/steady_state_types.html#Mathematical-Specification-of-a-Steady-State-Problem-1",
    "page": "Steady State Problems",
    "title": "Mathematical Specification of a Steady State Problem",
    "category": "section",
    "text": "To define an Steady State Problem, you simply need to give the function f which defines the ODE:fracdudt = f(upt)and an initial guess u₀ of where f(u,p,t)=0. f should be specified as f(u,p,t) (or in-place as f(du,u,p,t)), and u₀ should be an AbstractArray (or number) whose geometry matches the desired geometry of u. Note that we are not limited to numbers or vectors for u₀; one is allowed to provide u₀ as arbitrary matrices / higher dimension tensors as well.Note that for the steady-state to be defined, we must have that f is autonomous, that is f is independent of t. But the form which matches the standard ODE solver should still be used. The steady state solvers interpret the f by fixing t=0."
},

{
    "location": "types/steady_state_types.html#Problem-Type-1",
    "page": "Steady State Problems",
    "title": "Problem Type",
    "category": "section",
    "text": ""
},

{
    "location": "types/steady_state_types.html#Constructors-1",
    "page": "Steady State Problems",
    "title": "Constructors",
    "category": "section",
    "text": "SteadyStateProblem(f::ODEFunction,u0)\nSteadyStateProblem{isinplace}(f,u0)isinplace optionally sets whether the function is inplace or not. This is determined automatically, but not inferred. Additionally, the constructor from ODEProblems is provided:SteadyStateProblem(prob::ODEProblem)For specifying Jacobians and mass matrices, see the DiffEqFunctions page."
},

{
    "location": "types/steady_state_types.html#Fields-1",
    "page": "Steady State Problems",
    "title": "Fields",
    "category": "section",
    "text": "f: The function in the ODE.\nu0: The initial guess for the steady state."
},

{
    "location": "types/steady_state_types.html#Special-Solution-Fields-1",
    "page": "Steady State Problems",
    "title": "Special Solution Fields",
    "category": "section",
    "text": "The SteadyStateSolution type is different from the other DiffEq solutions because it does not have temporal information."
},

{
    "location": "types/bvp_types.html#",
    "page": "BVP Problems",
    "title": "BVP Problems",
    "category": "page",
    "text": ""
},

{
    "location": "types/bvp_types.html#BVP-Problems-1",
    "page": "BVP Problems",
    "title": "BVP Problems",
    "category": "section",
    "text": ""
},

{
    "location": "types/bvp_types.html#Mathematical-Specification-of-an-BVP-Problem-1",
    "page": "BVP Problems",
    "title": "Mathematical Specification of an BVP Problem",
    "category": "section",
    "text": "To define an BVP Problem, you simply need to give the function f and the initial condition u₀ which define an ODE:fracdudt = f(upt)along with an implicit function bc! which defines the residual equation, wherebc(upt) = 0is the manifold on which the solution must live. A common form for this is the two-point BVProblem where the manifold defines the solution at two points:u(t_0) = a\nu(t_f) = b"
},

{
    "location": "types/bvp_types.html#Problem-Type-1",
    "page": "BVP Problems",
    "title": "Problem Type",
    "category": "section",
    "text": ""
},

{
    "location": "types/bvp_types.html#Constructors-1",
    "page": "BVP Problems",
    "title": "Constructors",
    "category": "section",
    "text": "TwoPointBVProblem{isinplace}(f,bc!,u0,tspan)\nBVProblem{isinplace}(f,bc!,u0,tspan)For any BVP problem type, bc! is the inplace function:bc!(residual, u, p, t)where residual computed from the current u. u is an array of solution values where u[i] is at time t[i], while p are the parameters. For a TwoPointBVProblem, t = tspan. For the more general BVProblem, u can be all of the internal time points, and for shooting type methods u=sol the ODE solution. Note that all features of the ODESolution are present in this form. In both cases, the size of the residual matches the size of the initial condition."
},

{
    "location": "types/bvp_types.html#Fields-1",
    "page": "BVP Problems",
    "title": "Fields",
    "category": "section",
    "text": "f: The function for the ODE.\nbc: The boundary condition function.\nu0: The initial condition. Either the initial condition for the ODE as an initial value problem, or a Vector of values for u(t_i) for collocation methods\ntspan: The timespan for the problem."
},

{
    "location": "types/sde_types.html#",
    "page": "SDE Problems",
    "title": "SDE Problems",
    "category": "page",
    "text": ""
},

{
    "location": "types/sde_types.html#SDE-Problems-1",
    "page": "SDE Problems",
    "title": "SDE Problems",
    "category": "section",
    "text": ""
},

{
    "location": "types/sde_types.html#Mathematical-Specification-of-a-SDE-Problem-1",
    "page": "SDE Problems",
    "title": "Mathematical Specification of a SDE Problem",
    "category": "section",
    "text": "To define an SDE Problem, you simply need to give the forcing function f, the noise function g, and the initial condition u₀ which define an SDE:du = f(upt)dt + Σgᵢ(upt)dWⁱf and g should be specified as f(u,p,t) and  g(u,p,t) respectively, and u₀ should be an AbstractArray whose geometry matches the desired geometry of u. Note that we are not limited to numbers or vectors for u₀; one is allowed to provide u₀ as arbitrary matrices / higher dimension tensors as well. A vector of gs can also be defined to determine an SDE of higher Ito dimension."
},

{
    "location": "types/sde_types.html#Problem-Type-1",
    "page": "SDE Problems",
    "title": "Problem Type",
    "category": "section",
    "text": "Wraps the data which defines an SDE problemu = f(upt)dt + Σgᵢ(upt)dWⁱwith initial condition u0."
},

{
    "location": "types/sde_types.html#Constructors-1",
    "page": "SDE Problems",
    "title": "Constructors",
    "category": "section",
    "text": "SDEProblem(f::SDEFunction,g,u0,tspan,p=nothing;noise=WHITE_NOISE,noise_rate_prototype=nothing)\nSDEProblem{isinplace}(f,g,u0,tspan,p=nothing;noise=WHITE_NOISE,noise_rate_prototype=nothing) : Defines the SDE with the specified functions. The default noise is WHITE_NOISE. isinplace optionally sets whether the function is inplace or not. This is determined automatically, but not inferred.For specifying Jacobians and mass matrices, see the DiffEqFunctions page."
},

{
    "location": "types/sde_types.html#Fields-1",
    "page": "SDE Problems",
    "title": "Fields",
    "category": "section",
    "text": "f: The drift function in the SDE.\ng: The noise function in the SDE.\nu0: The initial condition.\ntspan: The timespan for the problem.\np: The optional parameters for the problem. Defaults to nothing.\nnoise: The noise process applied to the noise upon generation. Defaults to Gaussian white noise. For information on defining different noise processes, see the noise process documentation page\nnoise_rate_prototype: A prototype type instance for the noise rates, that is the output g. It can be any type which overloads A_mul_B! with itself being the middle argument. Commonly, this is a matrix or sparse matrix. If this is not given, it defaults to nothing, which means the problem should be interpreted as having diagonal noise.  \ncallback: A callback to be applied to every solver which uses the problem. Defaults to nothing."
},

{
    "location": "types/sde_types.html#Example-Problems-1",
    "page": "SDE Problems",
    "title": "Example Problems",
    "category": "section",
    "text": "Examples problems can be found in DiffEqProblemLibrary.jl.To use a sample problem, such as prob_sde_linear, you can do something like:#]add DiffEqProblemLibrary\nusing DiffEqProblemLibrary\nprob = prob_sde_linear\nsol = solve(prob)DiffEqProblemLibrary.prob_sde_linear\nDiffEqProblemLibrary.prob_sde_2Dlinear\nDiffEqProblemLibrary.prob_sde_wave\nDiffEqProblemLibrary.prob_sde_lorenz\nDiffEqProblemLibrary.prob_sde_cubic\nDiffEqProblemLibrary.prob_sde_additive\nDiffEqProblemLibrary.prob_sde_additivesystem"
},

{
    "location": "types/rode_types.html#",
    "page": "RODE Problems",
    "title": "RODE Problems",
    "category": "page",
    "text": ""
},

{
    "location": "types/rode_types.html#RODE-Problems-1",
    "page": "RODE Problems",
    "title": "RODE Problems",
    "category": "section",
    "text": ""
},

{
    "location": "types/rode_types.html#Mathematical-Specification-of-a-RODE-Problem-1",
    "page": "RODE Problems",
    "title": "Mathematical Specification of a RODE Problem",
    "category": "section",
    "text": "To define a RODE Problem, you simply need to give the function f and the initial condition u₀ which define an ODE:fracdudt = f(uptW(t))where W(t) is a random process. f should be specified as f(u,p,t,W) (or in-place as f(du,u,p,t,W)), and u₀ should be an AbstractArray (or number) whose geometry matches the desired geometry of u. Note that we are not limited to numbers or vectors for u₀; one is allowed to provide u₀ as arbitrary matrices / higher dimension tensors as well."
},

{
    "location": "types/rode_types.html#Constructors-1",
    "page": "RODE Problems",
    "title": "Constructors",
    "category": "section",
    "text": "RODEProblem(f::RODEFunction,u0,tspan,p=nothing;noise=WHITE_NOISE,rand_prototype=nothing,callback=nothing)\nRODEProblem{isinplace}(f,u0,tspan,p=nothing;noise=WHITE_NOISE,rand_prototype=nothing,callback=nothing,mass_matrix=I) : Defines the RODE with the specified functions. The default noise is WHITE_NOISE. isinplace optionally sets whether the function is inplace or not. This is determined automatically, but not inferred.For specifying Jacobians and mass matrices, see the DiffEqFunctions page."
},

{
    "location": "types/rode_types.html#Fields-1",
    "page": "RODE Problems",
    "title": "Fields",
    "category": "section",
    "text": "f: The drift function in the SDE.\nu0: The initial condition.\ntspan: The timespan for the problem.\np: The optional parameters for the problem. Defaults to nothing.\nnoise: The noise process applied to the noise upon generation. Defaults to Gaussian white noise. For information on defining different noise processes, see the noise process documentation page\nrand_prototype: A prototype type instance for the noise vector. It defaults to nothing, which means the problem should be interpreted as having a noise vector whose size matches u0.\ncallback: A callback to be applied to every solver which uses the problem. Defaults to nothing."
},

{
    "location": "types/dde_types.html#",
    "page": "DDE Problems",
    "title": "DDE Problems",
    "category": "page",
    "text": ""
},

{
    "location": "types/dde_types.html#DDE-Problems-1",
    "page": "DDE Problems",
    "title": "DDE Problems",
    "category": "section",
    "text": ""
},

{
    "location": "types/dde_types.html#Mathematical-Specification-of-a-DDE-Problem-1",
    "page": "DDE Problems",
    "title": "Mathematical Specification of a DDE Problem",
    "category": "section",
    "text": "To define a DDE Problem, you simply need to give the function f, the initial condition u_0 at time point t_0, and the history function h which together define a DDE:beginalign*\n    fracdudt = f(uhpt) qquad  (t geq t_0) \n    u(t_0) = u_0 \n    u(t) = h(t) qquad (t  t_0)\nendalign*f should be specified as f(u, h, p, t) (or in-place as f(du, u, h, p, t)), u_0 should be an AbstractArray (or number) whose geometry matches the desired geometry of u, and h should be specified as described below. The history function h is accessed for all delayed values. Note that we are not limited to numbers or vectors for u_0; one is allowed to provide u_0 as arbitrary matrices / higher dimension tensors as well."
},

{
    "location": "types/dde_types.html#Functional-Forms-of-the-History-Function-1",
    "page": "DDE Problems",
    "title": "Functional Forms of the History Function",
    "category": "section",
    "text": "The history function h can be called in the following ways:h(p, t): out-of-place calculation\nh(out, p, t): in-place calculation\nh(p, t, deriv::Type{Val{i}}): out-of-place calculation of the ith derivative\nh(out, p, t, deriv::Type{Val{i}}): in-place calculation of the ith derivative\nh(args...; idxs): calculation of h(args...) for indices idxsNote that a dispatch for the supplied history function of matching form is required for whichever function forms are used in the user derivative function f."
},

{
    "location": "types/dde_types.html#Declaring-Lags-1",
    "page": "DDE Problems",
    "title": "Declaring Lags",
    "category": "section",
    "text": "Lags are declared separately from their use. One can use any lag by simply using the interpolant of h at that point. However, one should use caution in order to achieve the best accuracy. When lags are declared, the solvers can more efficiently be more accurate and thus this is recommended."
},

{
    "location": "types/dde_types.html#Problem-Type-1",
    "page": "DDE Problems",
    "title": "Problem Type",
    "category": "section",
    "text": ""
},

{
    "location": "types/dde_types.html#Constructors-1",
    "page": "DDE Problems",
    "title": "Constructors",
    "category": "section",
    "text": "DDEProblem(f[, u0], h, tspan[, p]; <keyword arguments>)\nDDEProblem{isinplace}(f[, u0], h, tspan[, p]; <keyword arguments>)Parameter isinplace optionally sets whether the function is inplace or not. This is determined automatically, but not inferred.For specifying Jacobians and mass matrices, see the DiffEqFunctions page."
},

{
    "location": "types/dde_types.html#Arguments-1",
    "page": "DDE Problems",
    "title": "Arguments",
    "category": "section",
    "text": "f: The function in the DDE.\nu0: The initial condition. Defaults to the value h(p, first(tspan)) of the history function evaluated at the initial time point.\nh: The history function for the DDE before t0.\ntspan: The timespan for the problem.\np: The parameters with which function f is called. Defaults to nothing.\nconstant_lags: A collection of constant lags used by the history function h. Defaults to ().\ndependent_lags A tuple of functions (u, p, t) -> lag for the state-dependent lags used by the history function h. Defaults to ().\nneutral: If the DDE is neutral, i.e., if delays appear in derivative terms.\norder_discontinuity_t0: The order of the discontinuity at the initial time point. Defaults to 0 if an initial condition u0 is provided. Otherwise it is forced to be greater or equal than 1.\ncallback: A callback to be applied to every solver which uses the problem. Defaults to nothing."
},

{
    "location": "types/dae_types.html#",
    "page": "DAE Problems",
    "title": "DAE Problems",
    "category": "page",
    "text": ""
},

{
    "location": "types/dae_types.html#DAE-Problems-1",
    "page": "DAE Problems",
    "title": "DAE Problems",
    "category": "section",
    "text": ""
},

{
    "location": "types/dae_types.html#Mathematical-Specification-of-an-DAE-Problem-1",
    "page": "DAE Problems",
    "title": "Mathematical Specification of an DAE Problem",
    "category": "section",
    "text": "To define a DAE Problem, you simply need to give the function f and the initial condition u₀ which define an ODE:0 = f(duupt)f should be specified as f(du,u,p,t) (or in-place as f(resid,du,u,p,t)). Note that we are not limited to numbers or vectors for u₀; one is allowed to provide u₀ as arbitrary matrices / higher dimension tensors as well."
},

{
    "location": "types/dae_types.html#Problem-Type-1",
    "page": "DAE Problems",
    "title": "Problem Type",
    "category": "section",
    "text": ""
},

{
    "location": "types/dae_types.html#Constructors-1",
    "page": "DAE Problems",
    "title": "Constructors",
    "category": "section",
    "text": "DAEProblem(f::DAEFunction,du0,u0,tspan)\nDAEProblem{isinplace}(f,du0,u0,tspan) : Defines the DAE with the specified functions. isinplace optionally sets whether the function is inplace or not. This is determined automatically, but not inferred.For specifying Jacobians and mass matrices, see the DiffEqFunctions page."
},

{
    "location": "types/dae_types.html#Fields-1",
    "page": "DAE Problems",
    "title": "Fields",
    "category": "section",
    "text": "f: The function in the ODE.\ndu0: The initial condition for the derivative.\nu0: The initial condition.\ntspan: The timespan for the problem.\ncallback: A callback to be applied to every solver which uses the problem. Defaults to a black CallbackSet, which will have no effect.\ndifferential_vars: A logical array which declares which variables are the differential (non algebraic) vars (i.e. du\' is in the equations for this variable). Defaults to nothing. Some solvers may require this be set if an initial condition needs to be determined."
},

{
    "location": "types/dae_types.html#Example-Problems-1",
    "page": "DAE Problems",
    "title": "Example Problems",
    "category": "section",
    "text": "Examples problems can be found in DiffEqProblemLibrary.jl.To use a sample problem, such as prob_dae_resrob, you can do something like:#]add DiffEqProblemLibrary\nusing DiffEqProblemLibrary, Sundials\nprob = prob_dae_resrob\nsol = solve(prob,IDA())"
},

{
    "location": "types/jump_types.html#",
    "page": "Jump Problems",
    "title": "Jump Problems",
    "category": "page",
    "text": ""
},

{
    "location": "types/jump_types.html#Jump-Problems-1",
    "page": "Jump Problems",
    "title": "Jump Problems",
    "category": "section",
    "text": ""
},

{
    "location": "types/jump_types.html#Mathematical-Specification-of-an-problem-with-jumps-1",
    "page": "Jump Problems",
    "title": "Mathematical Specification of an problem with jumps",
    "category": "section",
    "text": "Jumps are defined as a Poisson process which changes states at some rate. When there are multiple possible jumps, the process is a compound Poisson process. On their own, a jump equation is continuous-time Markov Chain where the time to the next jump is exponentially distributed as calculated by the rate. This type of process, known in biology as \"Gillespie discrete stochastic simulations\" and modeled by the Chemical Master Equation (CME), is the same thing as adding jumps to a DiscreteProblem. However, any differential equation can be extended by jumps as well. For example, we have an ODE with jumps, denoted byfracdudt = f(upt) + Σ c_i(upt)dp_iwhere dp_i is a Poisson counter of rate lambda_i(upt). Extending a stochastic differential equation to have jumps is commonly known as a Jump Diffusion, and is denoted byfracdudt = f(upt) + Σgᵢ(ut)dWⁱ + Σ c_i(upt)dp_i"
},

{
    "location": "types/jump_types.html#Types-of-Jumps:-Regular,-Variable,-Constant-Rate-and-Mass-Action-1",
    "page": "Jump Problems",
    "title": "Types of Jumps: Regular, Variable, Constant Rate and Mass Action",
    "category": "section",
    "text": "A RegularJump is a set of jumps that do not make structural changes to the underlying equation. These kinds of jumps only change values of the dependent variable (u) and thus can be treated in an inexact manner. Other jumps, such as those which change the size of u, require exact handling which is also known as time-adaptive jumping. These can only be specified as a ConstantRateJump, MassActionJump, or a VariableRateJump.We denote a jump as variable rate if its rate function is dependent on values which may change between constant rate jumps. For example, if there are multiple jumps whose rates only change when one of them occur, than that set of jumps is a constant rate jump. If a jump\'s rate depends on the differential equation, time, or by some value which changes outside of any constant rate jump, then it is denoted as variable.A MassActionJump is a specialized representation for a collection of constant rate jumps that can each be interpreted as a standard mass action reaction. For systems comprised of many mass action reactions, using the MassActionJump type will offer improved performance. Note, only one MassActionJump should be defined per JumpProblem; it is then responsible for handling all mass action reaction type jumps. For systems with both mass action jumps and non-mass action jumps, one can create one MassActionJump to handle the mass action jumps, and create a number of ConstantRateJumps to handle the non-mass action jumps.RegularJumps are optimized for regular jumping algorithms like tau-leaping and hybrid algorithms. ConstantRateJumps and MassActionJumps are optimized for SSA algorithms. ConstantRateJumps, MassActionJumps and VariableRateJumps can be added to standard DiffEq algorithms since they are simply callbacks, while RegularJumps require special algorithms. "
},

{
    "location": "types/jump_types.html#Defining-a-Regular-Jump-1",
    "page": "Jump Problems",
    "title": "Defining a Regular Jump",
    "category": "section",
    "text": "The constructor for a RegularJump is:RegularJump(rate,c,c_prototype;mark_dist = nothing,constant_c = false)rate(out,u,p,t) is the function which computes the rate for every regular jump process\nc(dc,u,p,t,mark) is the current Stoichiometry matrix for each jump process\ndc is the cache array to be used for dc\nmark_dist is the distribution for the mark\nconstant_c denotes whether the Stoichiometry matrix c is constantdc is an n x m matrix, where n is the number of Poisson processes and m is the number of dependent variables (should match length(u)). rate is a vector equation which should compute the rates in to out which is a length n vector."
},

{
    "location": "types/jump_types.html#Defining-a-Constant-Rate-Jump-1",
    "page": "Jump Problems",
    "title": "Defining a Constant Rate Jump",
    "category": "section",
    "text": "The constructor for a ConstantRateJump is:ConstantRateJump(rate,affect!)rate(u,p,t) is a function which calculates the rate given the time and the state.\naffect!(integrator) is the effect on the equation, using the integrator interface."
},

{
    "location": "types/jump_types.html#Defining-a-Mass-Action-Jump-1",
    "page": "Jump Problems",
    "title": "Defining a Mass Action Jump",
    "category": "section",
    "text": "The constructor for a MassActionJump is:MassActionJump(rate_consts, reactant_stoich, net_stoich; scale_rates = true)rate_consts is a vector of the rate constants for each reaction.\nreactant_stoich is a vector whose kth entry is the reactant stoichiometry of the kth reaction. The reactant stoichiometry for an individual reaction is assumed to be represented as a vector of Pairs, mapping species id to stoichiometric coefficient.\nnet_stoich is assumed to have the same type as reactant_stoich; a vector whose kth entry is the net stoichiometry of the kth reaction. The net stoichiometry for an individual reaction is again represented as a vector of Pairs, mapping species id to the net change in the species when the reaction occurs.\nscale_rates is an optional parameter that specifies whether the rate constants correspond to stochastic rate constants in the sense used by Gillespie, and hence need to be rescaled. The default, scale_rates=true, corresponds to rescaling the passed in rate constants. See below.Notes for Mass Action JumpsWhen using MassActionJump the default behavior is to assume rate constants correspond to stochastic rate constants in the sense used by Gillespie (J. Comp. Phys., 1976, 22 (4)). This means that for a reaction such as 2A oversetkrightarrow B, the jump rate function constructed by MassActionJump would be k*A*(A-1)/2!. For a trimolecular reaction like 3A oversetkrightarrow B the rate function would be k*A*(A-1)*(A-2)/3!. To avoid having the reaction rates rescaled (by 1/2 and 1/6 for these two examples), one can pass the MassActionJump constructor the optional named parameter scale_rates=false, i.e. use\nMassActionJump(rates, reactant_stoich, net_stoich; scale_rates = false)\nZero order reactions can be passed as reactant_stoichs in one of two ways. Consider the varnothing oversetkrightarrow A reaction with rate k=1:\nk = [1.]\nreactant_stoich = [[0 => 1]]\nnet_stoich = [[1 => 1]]\njump = MassActionJump(k, reactant_stoich, net_stoich)\nAlternatively one can create an empty vector of pairs to represent the reaction:\nk = [1.]\nreactant_stoich = [Vector{Pair{Int,Int}}()]\nnet_stoich = [[1 => 1]]\njump = MassActionJump(k, reactant_stoich, net_stoich)\nFor performance reasons, it is recommended to order species indices in stoichiometry vectors from smallest to largest. That is \nreactant_stoich = [[1 => 2, 3 => 1, 4 => 2], [2 => 2, 3 => 2]]\nis preferred over\nreactant_stoich = [[3 => 1, 1 => 2, 4 = > 2], [3 => 2, 2 => 2]]"
},

{
    "location": "types/jump_types.html#Defining-a-Variable-Rate-Jump-1",
    "page": "Jump Problems",
    "title": "Defining a Variable Rate Jump",
    "category": "section",
    "text": "The constructor for a VariableRateJump is:VariableRateJump(rate,affect!;\n                   idxs = nothing,\n                   rootfind=true,\n                   save_positions=(true,true),\n                   interp_points=10,\n                   abstol=1e-12,reltol=0)Note that this is the same as defining a ContinuousCallback, except that instead of the condition function, you provide a rate(u,p,t) function for the rate at a given time and state."
},

{
    "location": "types/jump_types.html#Defining-a-Jump-Problem-1",
    "page": "Jump Problems",
    "title": "Defining a Jump Problem",
    "category": "section",
    "text": "To define a JumpProblem, you must first define the basic problem. This can be a DiscreteProblem if there is no differential equation, or an ODE/SDE/DDE/DAE if you would like to augment a differential equation with jumps. Denote this previously defined problem as prob. Then the constructor for the jump problem is:JumpProblem(prob,aggregator::Direct,jumps::JumpSet;\n            save_positions = typeof(prob) <: AbstractDiscreteProblem ? (false,true) : (true,true))The aggregator is the method for aggregating the constant jumps. These are defined below. jumps is a JumpSet which is just a gathering of jumps. Instead of passing a JumpSet, one may just pass a list of jumps themselves. For example:JumpProblem(prob,aggregator,jump1,jump2)and the internals will automatically build the JumpSet. save_positions is the save_positions argument built by the aggregation of the constant rate jumps.Note that a JumpProblem/JumpSet can only have 1 RegularJump (since a RegularJump itself describes multiple processes together). Similarly, it can only have one MassActionJump (since it also describes multiple processes together)."
},

{
    "location": "types/jump_types.html#Constant-Rate-Jump-Aggregators-1",
    "page": "Jump Problems",
    "title": "Constant Rate Jump Aggregators",
    "category": "section",
    "text": "Constant rate jump aggregators are the methods by which constant rate jumps, including MassActionJumps, are lumped together. This is required in all algorithms for both speed and accuracy. The current methods are:Direct: the Gillespie Direct method SSA.\nDirectCR: The Composition-Rejection Direct method of Slepoy et al. For large networks and linear chain-type networks it will often give better performance than Direct. (Requires dependency graph, see below.)\nDirectFW: the Gillespie Direct method SSA with FunctionWrappers. This aggregator uses a different internal storage format for collections of ConstantRateJumps. \nFRM: the Gillespie first reaction method SSA. Direct should generally offer better performance and be preferred to FRM.\nFRMFW: the Gillespie first reaction method SSA with FunctionWrappers.\nNRM: The Gibson-Bruck Next Reaction Method. For some reaction network structures  this may offer better performance than Direct (for example, large, linear chains of reactions). (Requires dependency graph, see below.) \nRSSA: The Rejection SSA (RSSA) method of Thanh et al. For very large reaction networks it often offers the best performance of all methods. (Requires dependency graph, see below.)\nSortingDirect: The Sorting Direct Method of McCollum et al. It will usually offer performance as good as Direct, and for some systems can offer substantially better performance. (Requires dependency graph, see below.)To pass the aggregator, pass the instantiation of the type. For example:JumpProblem(prob,Direct(),jump1,jump2)will build a problem where the constant rate jumps are solved using Gillespie\'s Direct SSA method."
},

{
    "location": "types/jump_types.html#Constant-Rate-Jump-Aggregators-Requiring-Dependency-Graphs-1",
    "page": "Jump Problems",
    "title": "Constant Rate Jump Aggregators Requiring Dependency Graphs",
    "category": "section",
    "text": "Italicized constant rate jump aggregators require the user to pass a dependency graph to JumpProblem. DirectCR, NRM and SortingDirect require a jump-jump dependency graph, passed through the named parameter dep_graph. i.e.JumpProblem(prob,DirectCR(),jump1,jump2; dep_graph=your_dependency_graph)For systems with only MassActionJumps, or those generated from a DiffEqBiological reaction_network, this graph will be auto-generated. Otherwise you must construct the dependency graph manually. Dependency graphs are represented as a Vector{Vector{Int}}, with the ith vector containing the indices of the jumps for which rates must be recalculated when the ith jump occurs.RSSA requires two different types of dependency graphs, passed through the following JumpProblem kwargs:vartojumps_map - A Vector{Vector{Int}} mapping each variable index, i, to a set of jump indices. The jump indices correspond to jumps with rate functions that depend on the value of u[i].\njumptovars_map - A Vector{Vector{Int}}  mapping each jump index to a set of variable indices. The corresponding variables are those that have their value, u[i], altered when the jump occurs.For systems generated from a DiffEqBiological reaction_network these will be auto-generated. Otherwise you must explicitly construct and pass in these mappings."
},

{
    "location": "types/jump_types.html#Recommendations-for-Constant-Rate-Jumps-1",
    "page": "Jump Problems",
    "title": "Recommendations for Constant Rate Jumps",
    "category": "section",
    "text": "For representing and aggregating constant rate jumps Use a MassActionJump to handle all jumps that can be represented as mass action reactions. This will generally offer the fastest performance. \nUse ConstantRateJumps for any remaining jumps.\nFor a small number of jumps, < ~10, Direct will often perform as well as the other aggregators.\nFor > ~10 jumps SortingDirect will often offer better performance than Direct.\nFor large number of jumps with sparse chain like structures and similar jump rates, for example continuous time random walks, DirectCR and then NRM often have the best performance.\nFor very large networks, with many updates per jump, RSSA will often substantially outperform the other methods.In general, for systems with sparse dependency graphs if Direct is slow, one of SortingDirect, DirectCR or RSSA will usually offer substantially better performance. See DiffEqBenchmarks.jl for benchmarks on several example networks."
},

{
    "location": "solvers/discrete_solve.html#",
    "page": "Discrete Solvers",
    "title": "Discrete Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "solvers/discrete_solve.html#Discrete-Solvers-1",
    "page": "Discrete Solvers",
    "title": "Discrete Solvers",
    "category": "section",
    "text": ""
},

{
    "location": "solvers/discrete_solve.html#DiscreteProblems-1",
    "page": "Discrete Solvers",
    "title": "DiscreteProblems",
    "category": "section",
    "text": "solve(prob::DiscreteProblem,alg;kwargs)Solves the discrete function map defined by prob using the algorithm alg. If no algorithm is given, a default algorithm will be chosen."
},

{
    "location": "solvers/discrete_solve.html#Recommended-Methods-1",
    "page": "Discrete Solvers",
    "title": "Recommended Methods",
    "category": "section",
    "text": "The implementation for solving discrete equations is the FunctionMap algorithm in OrdinaryDiffEq.jl. It allows the full common interface (including events/callbacks) to solve function maps, along with everything else like plot recipes, while completely ignoring the ODE functionality related to continuous equations (except for a tiny bit of initialization). However, the SimpleFunctionMap from SimpleDiffEq.jl can be more efficient if the mapping function is sufficiently cheap, but it doesn\'t have all of the extras like callbacks and saving support (but does have an integrator interface)."
},

{
    "location": "solvers/discrete_solve.html#Full-List-of-Methods-1",
    "page": "Discrete Solvers",
    "title": "Full List of Methods",
    "category": "section",
    "text": ""
},

{
    "location": "solvers/discrete_solve.html#OrdinaryDiffEq.jl-1",
    "page": "Discrete Solvers",
    "title": "OrdinaryDiffEq.jl",
    "category": "section",
    "text": "FunctionMap: A basic function map which implements the full common interface.OrdinaryDiffEq.jl also contains the FunctionMap algorithm which lets you  It has a piecewise constant interpolation and allows for all of the  callback/event handling capabilities (of course, with rootfind=false. If a  ContinuousCallback is given, it\'s always assumed rootfind=false).The constructor is:FunctionMap()\nFunctionMap{scale_by_time}()Every step is the updateu_n+1 = f(t_n+1u_n)If in addition scale_by_time is marked true (default is false),  then every step is the update:u_n+1 = u_n + dtf(t_n+1u_n)Notice that this is the same as updates from the Euler method, except in this case we assume that its a discrete change and thus the interpolation is piecewise constant."
},

{
    "location": "solvers/discrete_solve.html#SimpleDiffEq.jl-1",
    "page": "Discrete Solvers",
    "title": "SimpleDiffEq.jl",
    "category": "section",
    "text": "SimpleFunctionMap: A barebones implementation of a function map. Is optimally-efficient and has an integrator interface version, but does not support callbacks or saving controls."
},

{
    "location": "solvers/ode_solve.html#",
    "page": "ODE Solvers",
    "title": "ODE Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "solvers/ode_solve.html#ODE-Solvers-1",
    "page": "ODE Solvers",
    "title": "ODE Solvers",
    "category": "section",
    "text": "solve(prob::ODEProblem,alg;kwargs)Solves the ODE defined by prob using the algorithm alg. If no algorithm is given, a default algorithm will be chosen."
},

{
    "location": "solvers/ode_solve.html#Recommended-Methods-1",
    "page": "ODE Solvers",
    "title": "Recommended Methods",
    "category": "section",
    "text": "It is suggested that you try choosing an algorithm using the alg_hints keyword argument. However, in some cases you may want something specific, or you may just be curious. This guide is to help you choose the right algorithm."
},

{
    "location": "solvers/ode_solve.html#Unknown-Stiffness-Problems-1",
    "page": "ODE Solvers",
    "title": "Unknown Stiffness Problems",
    "category": "section",
    "text": "When the stiffness of the problem is unknown, it is recommended you use a stiffness detection and auto-switching algorithm. These methods are multi-paradigm and allow for efficient solution of both stiff and non-stiff problems. The cost for auto-switching is very minimal but the choices are restrained and so they are a good go-to method when applicable.For default tolerances, AutoTsit5(Rosenbrock23()) is a good choice. For lower tolerances, using AutoVern7 or AutoVern9 with Rodas4, KenCarp4, or Rodas5 can all be good choices depending on the problem. For very large systems (>1000 ODEs?), consider using lsoda."
},

{
    "location": "solvers/ode_solve.html#Non-Stiff-Problems-1",
    "page": "ODE Solvers",
    "title": "Non-Stiff Problems",
    "category": "section",
    "text": "For non-stiff problems, the native OrdinaryDiffEq.jl algorithms are vastly more efficient than the other choices. For most non-stiff problems, we recommend Tsit5. When more robust error control is required, BS5 is a good choice. If at moderate tolerances an the interpolation error is very important, consider the OwrenZen5 method. For fast solving at higher tolerances, we recommend BS3 (or OwrenZen3 if the interpolation error is important). For high accuracy but with the range of Float64 (~1e-8-1e-12), we recommend Vern6, Vern7, or Vern8 as efficient choices.For high accuracy non-stiff solving (BigFloat and tolerances like <1e-12), we recommend the Vern9 method. If a high-order method is needed with a high order interpolant, then you should choose Vern9 which is Order 9 with an Order 9 interpolant. If you need extremely high accuracy (<1e-30?) and do not need an interpolant, try the Feagin12 or Feagin14 methods. Note that the Feagin methods are the only high-order optimized methods which do not include a high-order interpolant (they do include a 3rd order Hermite interpolation if needed). Note that these high order RK methods are more robust than the high order Adams-Bashforth methods to discontinuities and achieve very high precision, and are much more efficient than the extrapolation methods. However, the VCABM method can be a good choice for high accuracy when the system of equations is very large (>1,000 ODEs?), the function calculation is very expensive, or the solution is very smooth.If strict error bounds are needed, then adaptive methods with defect controls are required. Defect controls use an error measurement on the interpolating polynomial to make the error estimate better capture the error over the full interval. For medium accuracy calculations, RK4 is a good choice."
},

{
    "location": "solvers/ode_solve.html#Stiff-Problems-1",
    "page": "ODE Solvers",
    "title": "Stiff Problems",
    "category": "section",
    "text": "For stiff problems at high tolerances (>1e-2?) it is recommended that you use Rosenbrock23 or TRBDF2. These are robust to oscillations and massive stiffness is needed, though are only efficient when low accuracy is needed. Rosenbrock23 is more efficient for small systems where re-evaluating and re-factorizing the Jacobian is not too costly, and for sufficiently large systems TRBDF2 will be more efficient. ABDF2 can be the most efficient the largest systems or most expensive f.At medium tolerances (>1e-8?) it is recommended you use Rodas5, Rodas4P (the former is more efficient but the later is more reliable), Kvaerno5, or KenCarp4. As native DifferentialEquations.jl solvers, many Julia numeric types (such as BigFloats, ArbFloats, or DecFP) will work. When the equation is defined via the @ode_def macro, these will be the most efficient.For faster solving at low tolerances (<1e-9) but when Vector{Float64} is used, use radau.For asymptotically large systems of ODEs (N>1000?) where f is very costly and the complex eigenvalues are minimal (low oscillations), in that case CVODE_BDF will be the most efficient but requires Vector{Float64}. CVODE_BDF will also do surprisingly well if the solution is smooth. However, this method can be less stiff than other methods and stuff may fail at low accuracy situations. Another good choice for this regime is lsoda."
},

{
    "location": "solvers/ode_solve.html#Special-Properties-of-Stiff-Integrators-1",
    "page": "ODE Solvers",
    "title": "Special Properties of Stiff Integrators",
    "category": "section",
    "text": "ImplicitMidpoint is a symmetric and symplectic integrator. Trapezoid is a symmetric (almost symplectic) integrator with adaptive timestepping. ImplicitEuler is an extension to the common algorithm with adaptive timestepping and efficient quasi-Newton Jacobian re-usage which is fully strong-stability presurving (SSP) for hyperbolic PDEs.Notice that Rodas4 loses accuracy on discretizations of nonlinear parabolic PDEs, and thus it\'s suggested you replace it with Rodas4P in those situations which is 3rd order. ROS3P is only third order and achieves 3rd order on such problems and can thus be more efficient in this case."
},

{
    "location": "solvers/ode_solve.html#Translations-from-MATLAB/Python/R-1",
    "page": "ODE Solvers",
    "title": "Translations from MATLAB/Python/R",
    "category": "section",
    "text": "For users familiar with MATLAB/Python/R, good translations of the standard library methods are as follows:ode23 –> BS3()\node45/dopri5 –> DP5(), though in most cases Tsit5() is more efficient\node23s –> Rosenbrock23(), though in most cases Rodas4() is more efficient\node113 –> VCABM(), though in many cases Vern7() is more efficient\ndop853 –> DP8(), though in most cases Vern7() is more efficient\node15s/vode –> QNDF(), though in many cases CVODE_BDF(), Rodas4() or radau() are more efficient\node23t –> Trapezoid() for efficiency and GenericTrapezoid() for robustness\node23tb –> TRBDF2\nlsoda –> lsoda() (requires ]add LSODA; using LSODA)\node15i –> IDA(), though in many cases Rodas4() can handle the DAE and is significantly more efficient"
},

{
    "location": "solvers/ode_solve.html#Full-List-of-Methods-1",
    "page": "ODE Solvers",
    "title": "Full List of Methods",
    "category": "section",
    "text": ""
},

{
    "location": "solvers/ode_solve.html#OrdinaryDiffEq.jl-for-Non-Stiff-Equations-1",
    "page": "ODE Solvers",
    "title": "OrdinaryDiffEq.jl for Non-Stiff Equations",
    "category": "section",
    "text": "Unless otherwise specified, the OrdinaryDiffEq algorithms all come with a 3rd order Hermite polynomial interpolation. The algorithms denoted as having a \"free\" interpolation means that no extra steps are required for the interpolation. For the non-free higher order interpolating functions, the extra steps are computed lazily (i.e. not during the solve).The OrdinaryDiffEq.jl algorithms achieve the highest performance for non-stiff equations while being the most generic: accepting the most Julia-based types, allow for sophisticated event handling, etc. On stiff ODEs these algorithms again consistently among the top. OrdinaryDiffEq.jl is recommended for most ODE problems."
},

{
    "location": "solvers/ode_solve.html#Explicit-Runge-Kutta-Methods-1",
    "page": "ODE Solvers",
    "title": "Explicit Runge-Kutta Methods",
    "category": "section",
    "text": "Euler- The canonical forward Euler method. Fixed timestep only.\nMidpoint - The second order midpoint method. Uses embedded Euler method for adaptivity.\nHeun - The second order Heun\'s method. Uses embedded Euler method for adaptivity.\nRalston - The optimized second order midpoint method. Uses embedded Euler. method for adaptivity.\nRK4 - The canonical Runge-Kutta Order 4 method. Uses a defect control for adaptive stepping using maximum error over the whole interval.\nBS3 - Bogacki-Shampine 3/2 method.\nOwrenZen3 - Owren-Zennaro optimized interpolantion 3/2 method (free 3th order interpolant).\nOwrenZen4 - Owren-Zennaro optimized interpolantion 4/3 method (free 4th order interpolant).\nOwrenZen5 - Owren-Zennaro optimized interpolantion 5/4 method (free 5th order interpolant).\nDP5 - Dormand-Prince\'s 5/4 Runge-Kutta method. (free 4th order interpolant).\nTsit5 - Tsitouras 5/4 Runge-Kutta method. (free 4th order interpolant).\nAnas5(w) - 4th order Runge-Kutta method designed for periodic problems. Requires a periodicity estimate w which when accurate the method becomes 5th order (and is otherwise 4th order with less error for better estimates).\nTanYam7 - Tanaka-Yamashita 7 Runge-Kutta method.\nDP8 - Hairer\'s 8/5/3 adaption of the Dormand-Prince Runge-Kutta method. (7th order interpolant).\nTsitPap8 - Tsitouras-Papakostas 8/7 Runge-Kutta method.\nFeagin10 - Feagin\'s 10th-order Runge-Kutta method.\nFeagin12 - Feagin\'s 12th-order Runge-Kutta method.\nFeagin14 - Feagin\'s 14th-order Runge-Kutta method.Example usage:alg = Tsit5()\nsolve(prob,alg)  Additionally, the following algorithms have a lazy interpolant:BS5 - Bogacki-Shampine 5/4 Runge-Kutta method. (lazy 5th order interpolant).\nVern6 - Verner\'s \"Most Efficient\" 6/5 Runge-Kutta method. (lazy 6th order interpolant).\nVern7 - Verner\'s \"Most Efficient\" 7/6 Runge-Kutta method. (lazy 7th order interpolant).\nVern8 - Verner\'s \"Most Efficient\" 8/7 Runge-Kutta method. (lazy 8th order interpolant)\nVern9 - Verner\'s \"Most Efficient\" 9/8 Runge-Kutta method. (lazy 9th order interpolant)These methods require a few extra steps in order to compute the high order interpolation, but these steps are only taken when the interpolation is used. These methods when lazy assume that the parameter vector p will be unchanged between the moment of the interval solving and the interpolation. If p is changed in a ContinuousCallback, or in a DiscreteCallback and the continuous solution is used after the full solution, then set lazy=false.Example:solve(prob,Vern7()) # lazy by default\nsolve(prob,Vern7(lazy=false))"
},

{
    "location": "solvers/ode_solve.html#Parallel-Explicit-Runge-Kutta-Methods-1",
    "page": "ODE Solvers",
    "title": "Parallel Explicit Runge-Kutta Methods",
    "category": "section",
    "text": "KuttaPRK2p5 - A 5 parallel, 2 processor explicit Runge-Kutta method of 5th order.These methods utilize multithreading on the f calls to parallelize the problem. This requires that simultaneous calls to f are thread-safe."
},

{
    "location": "solvers/ode_solve.html#Explicit-Strong-Stability-Preserving-Runge-Kutta-Methods-for-Hyperbolic-PDEs-(Conservation-Laws)-1",
    "page": "ODE Solvers",
    "title": "Explicit Strong-Stability Preserving Runge-Kutta Methods for Hyperbolic PDEs (Conservation Laws)",
    "category": "section",
    "text": "SSPRK22 - The two-stage, second order strong stability preserving (SSP) method of Shu and Osher (SSP coefficient 1, free 2nd order SSP interpolant). Fixed timestep only.\nSSPRK33 - The three-stage, third order strong stability preserving (SSP) method of Shu and Osher (SSP coefficient 1, free 2nd order SSP interpolant). Fixed timestep only.\nSSPRK53 - The five-stage, third order strong stability preserving (SSP) method of Ruuth (SSP coefficient 2.65, free 3rd order Hermite interpolant). Fixed timestep only.\nSSPRK63 - The six-stage, third order strong stability preserving (SSP) method of Ruuth (SSP coefficient 3.518, free 3rd order Hermite interpolant). Fixed timestep only.\nSSPRK73 - The seven-stage, third order strong stability preserving (SSP) method of Ruuth (SSP coefficient 4.2879, free 3rd order Hermite interpolant). Fixed timestep only.\nSSPRK83 - The eight-stage, third order strong stability preserving (SSP) method of Ruuth (SSP coefficient 5.107, free 3rd order Hermite interpolant). Fixed timestep only.\nSSPRK432 - A  3/2 adaptive strong stability preserving (SSP) method with five stages (SSP coefficient 2, free 2nd order SSP interpolant).\nSSPRK932 - A  3/2 adaptive strong stability preserving (SSP) method with nine stages (SSP coefficient 6, free 3rd order Hermite interpolant).\nSSPRK54 - The five-stage, fourth order strong stability preserving (SSP) method of Spiteri and Ruuth (SSP coefficient 1.508, 3rd order Hermite interpolant). Fixed timestep only.\nSSPRK104 - The ten-stage, fourth order strong stability preserving method of Ketcheson (SSP coefficient 6, free 3rd order Hermite interpolant). Fixed timestep only.\nSSPRKMSVS32 - 3-stage, 2nd order SSP-optimal linear multistep method. (SSP coefficent 0.5, 3rd order Hermite interpolant). Fixed timestep only.\nSSPRKMSVS43 - 4-stage, 3rd order SSP-optimal linear multistep method. (SSP coefficent 0.33, 3rd order Hermite interpolant). Fixed timestep only.The SSP coefficients of the methods can be queried as ssp_coefficient(alg). All explicit SSP methods take two optional arguments SSPXY(stage_limiter!, step_limiter!), where stage_limiter! and step_limiter are functions taking arguments of the form limiter!(u, f, t). Here, u is the new solution value (updated inplace) after an explicit Euler stage / the whole time step , f the time derivative function (semidiscretisation for PDEs), and t the current time. These limiters can be used to enforce physical constraints, e.g. the positivity preserving limiters of Zhang and Shu (Zhang, Xiangxiong, and Chi-Wang Shu. \"Maximum-principle-satisfying and positivity-preserving high-order schemes for conservation laws: survey and new developments.\" Proceedings of the Royal Society of London A: Mathematical, Physical and Engineering Sciences. The Royal Society, 2011.)."
},

{
    "location": "solvers/ode_solve.html#Low-Storage-Methods-1",
    "page": "ODE Solvers",
    "title": "Low-Storage Methods",
    "category": "section",
    "text": "ORK256 - 5-stage, second order low-storage method for wave propogation equations. Fixed timestep only.\nSSPRK53_2N1 and SSPRK53_2N2 - 5-stage, third order low-storage methods with large SSP coefficients. (SSP coefficient 2.18 and 2.15, free 3rd order Hermite interpolant). Fixed timestep only.\nCarpenterKennedy2N54 - The five-stage, fourth order low-storage method of Carpenter and Kennedy (free 3rd order Hermite interpolant). Fixed timestep only. Designed for hyperbolic PDEs (stability properties).\nNDBLSRK124 - 12-stage, fourth order low-storage method with optimized stability regions for advection-dominated problems. Fixed timestep only.\nNDBLSRK134 - 13-stage, fourth order low-storage method with optimized stability regions for advection-dominated problems. Fixed timestep only.\nNDBLSRK144 - 14-stage, fourth order low-storage method with optimized stability regions for advection-dominated problems. Fixed timestep only.\nCFRLDDRK64 - 6-stage, fourth order low-storage, low-dissipation, low-dispersion scheme. Fixed timestep only.\nTSLDDRK74 - 7-stage, fourth order low-storage low-dissipation, low-dispersion scheme with maximal accuracy and stability limit along the imaginary axes. Fixed timestep only.\nDGLDDRK73_C - 7-stage, third order low-storage low-dissipation, low-dispersion scheme for discontinuous Galerkin space discretizations applied to wave propagation problems, optimized for PDE discretizations when maximum spatial step is small due to geometric features of computational domain. Fixed timestep only.\nDGLDDRK84_C - 8-stage, fourth order low-storage low-dissipation, low-dispersion scheme for discontinuous Galerkin space discretizations applied to wave propagation problems, optimized for PDE discretizations when maximum spatial step is small due to geometric features of computational domain. Fixed timestep only.\nDGLDDRK84_F - 8-stage, fourth order low-storage low-dissipation, low-dispersion scheme for discontinuous Galerkin space discretizations applied to wave propagation problems, optimized for PDE discretizations when the maximum spatial step size is not constrained. Fixed timestep only.\nHSLDDRK64 - 6-stage, fourth order low-stage, low-dissipation, low-dispersion scheme. Fixed timestep only.\nRK46NL - 6-stage, fourth order low-stage, low-dissipation, low-dispersion scheme. Fixed timestep only.\nParsaniKetchesonDeconinck3S32 - 3-stage, second order (3S) low-storage scheme, optimised for for the spectral difference method applied to wave propagation problems.\nParsaniKetchesonDeconinck3S82 - 8-stage, second order (3S) low-storage scheme, optimised for for the spectral difference method applied to wave propagation problems.\nParsaniKetchesonDeconinck3S53 - 5-stage, third order (3S) low-storage scheme, optimised for for the spectral difference method applied to wave propagation problems.\nParsaniKetchesonDeconinck3S173 - 17-stage, third order (3S) low-storage scheme, optimised for for the spectral difference method applied to wave propagation problems.\nParsaniKetchesonDeconinck3S94 - 9-stage, fourth order (3S) low-storage scheme, optimised for for the spectral difference method applied to wave propagation problems.\nParsaniKetchesonDeconinck3S184 - 18-stage, fourth order (3S) low-storage scheme, optimised for for the spectral difference method applied to wave propagation problems.\nParsaniKetchesonDeconinck3S105 - 10-stage, fifth order (3S) low-storage scheme, optimised for for the spectral difference method applied to wave propagation problems.\nParsaniKetchesonDeconinck3S205 - 20-stage, fifth order (3S) low-storage scheme, optimised for for the spectral difference method applied to wave propagation problems.\nCKLLSRK43_2 - 4-stage, third order low-storage scheme, optimised for compressible Navier–Stokes equations..\nCKLLSRK54_3C - 5-stage, fourth order low-storage scheme, optimised for compressible Navier–Stokes equations.\nCKLLSRK95_4S - 9-stage, fifth order low-storage scheme, optimised for compressible Navier–Stokes equations.\nCKLLSRK95_4C - 9-stage, fifth order low-storage scheme, optimised for compressible Navier–Stokes equations.\nCKLLSRK95_4M - 9-stage, fifth order low-storage scheme, optimised for compressible Navier–Stokes equations.\nCKLLSRK54_3C_3R - 5-stage, fourth order low-storage scheme, optimised for compressible Navier–Stokes equations.\nCKLLSRK54_3M_3R - 5-stage, fourth order low-storage scheme, optimised for compressible Navier–Stokes equations.\nCKLLSRK54_3N_3R - 5-stage, fourth order low-storage scheme, optimised for compressible Navier–Stokes equations.\nCKLLSRK85_4C_3R - 8-stage, fifth order low-storage scheme, optimised for compressible Navier–Stokes equations.\nCKLLSRK85_4M_3R - 8-stage, fifth order low-storage scheme, optimised for compressible Navier–Stokes equations.\nCKLLSRK85_4P_3R - 8-stage, fifth order low-storage scheme, optimised for compressible Navier–Stokes equations.\nCKLLSRK54_3N_4R - 5-stage, fourth order low-storage scheme, optimised for compressible Navier–Stokes equations.\nCKLLSRK54_3M_4R - 5-stage, fourth order low-storage scheme, optimised for compressible Navier–Stokes equations.\nCKLLSRK65_4M_4R - 6-stage, fifth order low-storage scheme, optimised for compressible Navier–Stokes equations.\nCKLLSRK85_4FM_4R - 8-stage, fifth order low-storage scheme, optimised for compressible Navier–Stokes equations.\nCKLLSRK75_4M_5R - 7-stage, fifth order low-storage scheme, optimised for compressible Navier–Stokes equations.NOTE: All the 2N Methods (ORK256, CarpenterKennedy2N54, NDBLSRK124, NDBLSRK134, NDBLSRK144, DGLDDRK73_C, DGLDDRK84_C, DGLDDRK84_F and HSLDDRK64) work on the basic principle of being able to perform step S1 = S1 + F(S2) in just 2 registers. Certain optimizations have been done to achieve this theoritical limit (when alias_u0 is set) but have a limitation that du should always be on the left hand side (assignments only) in the implementation.Example - This is an invalid implementation for 2N methods:function f(du,u,p,t)\n  du[1] = u[1] * u[2]\n  du[2] = du[1] * u[2] # du appears on the RHS\nendIf you don\'t wish to have the optimization and have to use du on the RHS, please set the keyword argument williamson_condition to false in the algorithm (by default it is set to true). In this case 3 registers worth memory would be needed instead.Example :alg = CarpenterKennedy2N54(;williamson_condition=false)So the above implementation of f becomes valid."
},

{
    "location": "solvers/ode_solve.html#Parallelized-Explicit-Extrapolation-Methods-1",
    "page": "ODE Solvers",
    "title": "Parallelized Explicit Extrapolation Methods",
    "category": "section",
    "text": "The following are adaptive order, adaptive step size extrapolation methods:AitkenNevillie - Euler extrapolation using Aitken-Neville with the Romberg Sequence.\nExtrapolationMidpointDeuflhard - Midpoint extrapolation using Barycentric coordinates\nExtrapolationMidpointHairerWanner - Midpoint extrapolation using Barycentric coordinates, following Hairer\'s ODEX in the adaptivity behavior.These methods have arguments for max_order, min_order, and init_order on the adaptive order algorithm. threading denotes whether to automatically multithread the f evaluations, allowing for a high degree of within-method parallelism. The defaults are:max_order=10\nmin_order=1 except for ExtrapolationMidpointHairerWanner it\'s 2.\ninit_order=5\nthreading=trueAdditionally, the ExtrapolationMidpointDeuflhard and ExtrapolationMidpointHairerWanner methods have the additional argument:sequence: the step-number sequences, also called the subdividingsequence. Possible values are :harmonic, :romberg or :bulirsch. Default  is :harmonic.To override, utilize the keyword arguments. For example:alg = ExtrapolationMidpointDeuflhard(max_order=7,min_order=4,init_order=4,sequence=:bulirsch,threading=false)\nsolve(prob,alg)Note that the order that is referred to is the extrapolation order. For AitkenNevillie this is the order of the method, for the others an extrapolation order of n gives an order 2(n+1) method."
},

{
    "location": "solvers/ode_solve.html#Explicit-Multistep-Methods-1",
    "page": "ODE Solvers",
    "title": "Explicit Multistep Methods",
    "category": "section",
    "text": "Methods using the approximation at more than one previous mesh point to determine the approximation at the next point are called multistep methods. These methods tend to be more efficient as the size of the system or the cost of f increases."
},

{
    "location": "solvers/ode_solve.html#Adams-Bashforth-Explicit-Methods-1",
    "page": "ODE Solvers",
    "title": "Adams-Bashforth Explicit Methods",
    "category": "section",
    "text": "These methods require a choice of dt.AB3 - The 3-step third order multistep method. Ralston\'s Second Order Method is used to calculate starting values.\nAB4 - The 4-step fourth order multistep method. Runge-Kutta method of order 4 is used to calculate starting values.  \nAB5 - The 5-step fifth order multistep method. Runge-Kutta method of order 4 is used to calculate starting values.  \nABM32 - It is third order method. In ABM32, AB3 works as predictor and Adams Moulton 2-steps method works as Corrector. Ralston\'s Second Order Method is used to calculate starting values.  \nABM43 - It is fourth order method. In ABM43, AB4 works as predictor and Adams Moulton 3-steps method works as Corrector. Runge-Kutta method of order 4 is used to calculate starting values.  \nABM54 - It is fifth order method. In ABM54, AB5 works as predictor and Adams Moulton 4-steps method works as Corrector. Runge-Kutta method of order 4 is used to calculate starting values."
},

{
    "location": "solvers/ode_solve.html#Adaptive-step-size-Adams-explicit-Methods-1",
    "page": "ODE Solvers",
    "title": "Adaptive step size Adams explicit Methods",
    "category": "section",
    "text": "VCAB3 - The 3rd order Adams method. Bogacki-Shampine 3/2 method is used to calculate starting values.  \nVCAB4 - The 4th order Adams method. Runge-Kutta 4 is used to calculate starting values.  \nVCAB5 - The 5th order Adams method. Runge-Kutta 4 is used to calculate starting values.\nVCABM3 - The 3rd order Adams-Moulton method. Bogacki-Shampine 3/2 method is used to calculate starting values.  \nVCABM4 - The 4th order Adams-Moulton method. Runge-Kutta 4 is used to calculate starting values.  \nVCABM5 - The 5th order Adams-Moulton method. Runge-Kutta 4 is used to calculate starting values.\nVCABM - An adaptive order adaptive time Adams Moulton method. It uses an order adaptivity algorithm is derived from Shampine\'s DDEABM.\nAN5 - An adaptive 5th order fixed-leading coefficient Adams method in Nordsieck form.\nJVODE_Adams - An adaptive time adaptive order fixed-leading coefficient Adams method in Nordsieck form. The order adaptivity algorithm is derived from Sundials\' CVODE_Adams. In development."
},

{
    "location": "solvers/ode_solve.html#OrdinaryDiffEq.jl-for-Stiff-Equations-1",
    "page": "ODE Solvers",
    "title": "OrdinaryDiffEq.jl for Stiff Equations",
    "category": "section",
    "text": ""
},

{
    "location": "solvers/ode_solve.html#SDIRK-Methods-1",
    "page": "ODE Solvers",
    "title": "SDIRK Methods",
    "category": "section",
    "text": "ImplicitEuler - A 1st order implicit solver. A-B-L-stable. Adaptive timestepping through a divided differences estimate via memory. Strong-stability preserving (SSP).\nImplicitMidpoint - A second order A-stable symplectic and symmetric implicit solver. Good for highly stiff equations which need symplectic integration.\nTrapezoid - A second order A-stable symmetric ESDIRK method. \"Almost symplectic\" without numerical dampening. Also known as Crank-Nicolson when applied to PDEs. Adaptive timestepping via divided differences on the memory. Good for highly stiff equations which are non-oscillatory.\nTRBDF2 - A second order A-B-L-S-stable one-step ESDIRK method. Includes stiffness-robust error estimates for accurate adaptive timestepping, smoothed derivatives for highly stiff and oscillatory problems.\nGenericImplicitEuler - A 1st order A-B-L-stable implicit solver with adaptive timestepping through a divided differences estimate via memory. Strong-stability preserving (SSP). Uses an external nonlinear solver. Defaults to trust region dogleg with full Newton, making it more robust to numerical instability at the cost of being less efficient.\nGenericTrapezoid - A second order A-stable symplectic implicit solver. Also known as Crank-Nicolson when applied to PDEs. Adaptive timestepping via divided differences on the memory. Good for highly stiff equations which are non-oscillatory. Uses an external nonlinear solver. Defaults to trust region dogleg with full Newton, making it more robust to numerical instability at the cost of being less efficient.\nSDIRK2 - An A-B-L stable 2nd order SDIRK method\nKvaerno3 - An A-L stable stiffly-accurate 3rd order ESDIRK method\nKenCarp3 - An A-L stable stiffly-accurate 3rd order ESDIRK method with splitting\nCash4 - An A-L stable 4th order SDIRK method\nHairer4 - An A-L stable 4rd order SDIRK method\nHairer42 - An A-L stable 4rd order SDIRK method\nKvaerno4 - An A-L stable stiffly-accurate 4rd order ESDIRK method\nKenCarp4 - An A-L stable stiffly-accurate 4rd order ESDIRK method with splitting\nKvaerno5 - An A-L stable stiffly-accurate 5rd order ESDIRK method\nKenCarp5 - An A-L stable stiffly-accurate 5rd order ESDIRK method with splitting"
},

{
    "location": "solvers/ode_solve.html#Fully-Implicit-Runge-Kutta-Methods-(FIRK)-1",
    "page": "ODE Solvers",
    "title": "Fully-Implicit Runge-Kutta Methods (FIRK)",
    "category": "section",
    "text": "RadauIIA5 - An A-B-L stable fully implicit Runge-Kutta method with internal tableau complex basis transform for efficiency."
},

{
    "location": "solvers/ode_solve.html#Parallel-Diagonally-Implicit-Runge-Kutta-Methods-1",
    "page": "ODE Solvers",
    "title": "Parallel Diagonally Implicit Runge-Kutta Methods",
    "category": "section",
    "text": "PDIRK44 - A 2 processor 4th order diagonally non-adaptive implicit method.These methods also have option nlsolve same as SDIRK methods. These methods also need f to be thread safe. It parallelises the nlsolve calls inside the method.  "
},

{
    "location": "solvers/ode_solve.html#Rosenbrock-Methods-1",
    "page": "ODE Solvers",
    "title": "Rosenbrock Methods",
    "category": "section",
    "text": "ROS3P - 3rd order A-stable and stiffly stable Rosenbrock method. Keeps high accuracy on discretizations of nonlinear parabolic PDEs.\nRodas3 - 3rd order A-stable and stiffly stable Rosenbrock method.\nRosShamp4- An A-stable 4th order Rosenbrock method.\nVeldd4 - A 4th order D-stable Rosenbrock method.\nVelds4 - A 4th order A-stable Rosenbrock method.\nGRK4T - An efficient 4th order Rosenbrock method.\nGRK4A - An A-stable 4th order Rosenbrock method. Essentially \"anti-L-stable\" but efficient.\nRos4LStab - A 4th order L-stable Rosenbrock method.\nRodas4 - A 4th order A-stable stiffly stable Rosenbrock method with a stiff-aware 3rd order interpolant\nRodas42 - A 4th order A-stable stiffly stable Rosenbrock method with a stiff-aware 3rd order interpolant\nRodas4P - A 4th order A-stable stiffly stable Rosenbrock method with a stiff-aware 3rd order interpolant. 4th order on linear parabolic problems and 3rd order accurate on nonlinear parabolic problems (as opposed to lower if not corrected).\nRodas5 - A 5th order A-stable stiffly stable Rosenbrock method. Currently has a Hermite interpolant because its stiff-aware 3rd order interpolant is not yet implemented."
},

{
    "location": "solvers/ode_solve.html#Rosenbrock-W-Methods-1",
    "page": "ODE Solvers",
    "title": "Rosenbrock-W Methods",
    "category": "section",
    "text": "Rosenbrock23 - An Order 2/3 L-Stable Rosenbrock-W method which is good for very stiff equations with oscillations at low tolerances. 2nd order stiff-aware interpolation.\nRosenbrock32 - An Order 3/2 A-Stable Rosenbrock-W method which is good for mildy stiff equations without oscillations at low tolerances. Note that this method is prone to instability in the presence of oscillations, so use with caution. 2nd order stiff-aware interpolation.\nRosenbrockW6S4OS - A 4th order L-stable Rosenbrock-W method (fixed step only).\nROS34PW1a - A 4th order L-stable Rosenbrock-W method.\nROS34PW1b - A 4th order L-stable Rosenbrock-W method.\nROS34PW2 - A 4th order stiffy accurate Rosenbrock-W method for PDAEs.\nROS34PW3 - A 4th order strongly A-stable (Rinf~0.63) Rosenbrock-W method."
},

{
    "location": "solvers/ode_solve.html#Stabilized-Explicit-Methods-1",
    "page": "ODE Solvers",
    "title": "Stabilized Explicit Methods",
    "category": "section",
    "text": "ROCK2 - Second order stabilized Runge-Kutta method. Exhibits high stability for real eigenvalues and is smoothened to allow for moderate sized complex eigenvalues.\nROCK4 - Fourth order stabilized Runge-Kutta method. Exhibits high stability for real eigenvalues and is smoothened to allow for moderate sized complex eigenvalues.\nRKC - Second order stabilized Runge-Kutta method. Exhibits high stability for real eigenvalues and is smoothened to allow for moderate sized complex eigenvalues.\nSERK2v2 - Second order stabilized extrapolated Runge-Kutta method. Exhibits high stability for real eigenvalues and is smoothened to allow for moderate sized complex eigenvalues.\nESERK5 - Fifth order stabilized extrapolated Runge-Kutta method. Exhibits high stability for real eigenvalues and is smoothened to allow for moderate sized complex eigenvalues.ROCK methods offer a min_stages and max_stages functionality. SERK methods derive higher orders by Aitken-Neville algorithm. SERK2v2 is defaulted to Predictive control but has option of PI control."
},

{
    "location": "solvers/ode_solve.html#Parallelized-Implicit-Extrapolation-Methods-1",
    "page": "ODE Solvers",
    "title": "Parallelized Implicit Extrapolation Methods",
    "category": "section",
    "text": "The following are adaptive order, adaptive step size extrapolation methods:ImplicitEulerExtrapolation - Extrapolation of implicit Euler method with Romberg sequence. Similar to Hairer\'s SEULEX.\nImplicitDeuflhardExtrapolation - Midpoint extrapolation using Barycentric coordinates\nImplicitHairerWannerExtrapolation - Midpoint extrapolation using Barycentric coordinates, following Hairer\'s SODEX in the adaptivity behavior.These methods have arguments for max_order, min_order, and init_order on the adaptive order algorithm. threading denotes whether to automatically multithread the f evaluations and J/W instantiations+factorizations, allowing for a high degree of within-method parallelism. The defaults are:max_order=10\nmin_order=1 except for ImplicitHairerWannerExtrapolation it\'s 2.\ninit_order=5\nthreading=trueAdditionally, the ImplicitDeuflhardExtrapolation and ImplicitHairerWannerExtrapolation methods have the additional argument:sequence: the step-number sequences, also called the subdividingsequence. Possible values are :harmonic, :romberg or :bulirsch. Default  is :harmonic.To override, utilize the keyword arguments. For example:alg = ImplicitEulerExtrapolation(max_order=7,min_order=4,init_order=4,sequence=:bulirsch)\nsolve(prob,alg)Note that the order that is referred to is the extrapolation order. For ImplicitEulerExtrapolation this is the order of the method, for the others an extrapolation order of n gives an order 2(n+1) method."
},

{
    "location": "solvers/ode_solve.html#Parallelized-DIRK-Methods-1",
    "page": "ODE Solvers",
    "title": "Parallelized DIRK Methods",
    "category": "section",
    "text": "These methods parallelize the J/W instantiation and factorization, making them efficient on small highly stiff ODEs. Has an option threading=true to turn on/off multithreading.PDIRK44: a 4th order 2-processor DIRK method."
},

{
    "location": "solvers/ode_solve.html#Exponential-Methods-for-Linear-and-Affine-Problems-1",
    "page": "ODE Solvers",
    "title": "Exponential Methods for Linear and Affine Problems",
    "category": "section",
    "text": "LinearExponential - Exact solution formula for linear, time-independent problems. Expects the right hand side function to be a AbstractDiffEqOperator.Options:krylov - symbol. One of\n:off (default) - cache the operator beforehand. Requires Matrix(A) method defined for the operator A.\n:simple - uses simple Krylov approximations with fixed subspace size m.\n:adaptive - uses adaptive Krylov approximations with internal timestepping.\nm - integer, default: 30. Controls the size of Krylov subsapce if krylov=:simple, and the initial subspace size if krylov=:adaptive.\niop - integer, default: 0. If not zero, determines the length of the incomplete orthogonalization procedure (IOP) [1]. Note that if the linear operator/jacobian is hermitian, then the Lanczos algorithm will always be used and the IOP setting is ignored."
},

{
    "location": "solvers/ode_solve.html#Exponential-Runge-Kutta-Methods-1",
    "page": "ODE Solvers",
    "title": "Exponential Runge-Kutta Methods",
    "category": "section",
    "text": "These methods are all fixed timestepping only.LawsonEuler - First order exponential Euler scheme.\nNorsettEuler - First order exponential-RK scheme. Alias: ETD1.\nETD2 - Second order Exponential Time Differencing method (in development).\nETDRK2 - 2nd order exponential-RK scheme.\nETDRK3 - 3rd order exponential-RK scheme.\nETDRK4 - 4th order exponential-RK scheme.\nHochOst4 - 4th order exponential-RK scheme with stiff order 4.The methods are intended for semilinear problems constructed by SplitODEProblem or SplitODEFunction. They can also be used for a general nonlinear problem, in which case the jacobian of the right hand side is used as the linear operator in each time step.Except for ETD2, all methods come with these options, which can be set in the methods\' constructor:krylov - boolean, default: false. Determines whether Krylov approximation or operator caching is used, the latter only available for semilinear problems.\nm - integer, default: 30. Controls the size of Krylov subsapce.\niop - integer, default: 0. If not zero, determines the length of the incomplete orthogonalization procedure (IOP) [1]. Note that if the linear operator/jacobian is hermitian, then the Lanczos algorithm will always be used and the IOP setting is ignored.\nautodiff and chunksize: autodiff control if problem is not semilinear and explicit jacobian is not given. See Extra Options for more details."
},

{
    "location": "solvers/ode_solve.html#Adaptive-Exponential-Rosenbrock-Methods-1",
    "page": "ODE Solvers",
    "title": "Adaptive Exponential Rosenbrock Methods",
    "category": "section",
    "text": "Exprb32 - 3rd order adaptive Exponential-Rosenbrock scheme.\nExprb43 - 4th order adaptive Exponential-Rosenbrock scheme.The exponential rosenbrock methods cannot be applied to semilinear problems. Options for the solvers are the same as Exponential Runge-Kutta Methods except that Krylov approximation is always used."
},

{
    "location": "solvers/ode_solve.html#Exponential-Propagation-Iterative-Runge-Kutta-Methods-(EPIRK)-1",
    "page": "ODE Solvers",
    "title": "Exponential Propagation Iterative Runge-Kutta Methods (EPIRK)",
    "category": "section",
    "text": "These methods are all fixed timestepping only.Exp4 - 4th order EPIRK scheme.\nEPIRK4s3A - 4th order EPIRK scheme with stiff order 4.\nEPIRK4s3B - 4th order EPIRK scheme with stiff order 4.\nEPIRK5P1 - 5th order EPIRK scheme.\nEPIRK5P2 - 5th order EPIRK scheme.\nEPIRK5s3 - 5th order \"horizontal\" EPIRK scheme with stiff order 5. Broken.\nEXPRB53s3- 5th order EPIRK scheme with stiff order 5.Options:adaptive_krylov - boolean, default: true. Determines if the adaptive Krylov algorithm with timestepping of Neisen & Wright is used.\nm - integer, default: 30. Controls the size of Krylov subsapce, or the size for the first step if adaptive_krylov=true.\niop - integer, default: 0. If not zero, determines the length of the incomplete orthogonalization procedure (IOP) [1]. Note that if the linear operator/jacobian is hermitian, then the Lanczos algorithm will always be used and the IOP setting is ignored.\nautodiff and chunksize: autodiff control if problem is not semilinear and explicit jacobian is not given. See Extra Options for more details.It should be noted that many of the methods are still at an experimental stage of development, and thus should be used with caution."
},

{
    "location": "solvers/ode_solve.html#Multistep-Methods-1",
    "page": "ODE Solvers",
    "title": "Multistep Methods",
    "category": "section",
    "text": "Quasi-constant stepping is the time stepping strategy which matches the classic GEAR, LSODE,  and ode15s integrators. The variable-coefficient methods match the ideas of the classic EPISODE integrator and early VODE designs. The Fixed Leading Coefficient (FLC) methods match the behavior of the classic VODE and Sundials CVODE integrator.QNDF1 - An adaptive order 1 quasi-constant timestep L-stable numerical differentiation function (NDF) method. Optional parameter kappa defaults to Shampine\'s accuracy-optimal -0.1850.\nQBDF1 - An adaptive order 1 L-stable BDF method. This is equivalent to implicit Euler but using the BDF error estimator.\nABDF2 - An adaptive order 2 L-stable fixed leading coefficient multistep BDF method.\nQNDF - An adaptive order quasi-constant timestep NDF method. Utilizes Shampine\'s accuracy-optimal kappa values as defaults (has a keyword argument for a tuple of kappa coefficients).\nQBDF - An adaptive order quasi-constant timestep BDF method.\nJVODE_BDF - An adaptive time adaptive order fixed-leading coefficient BDF method in Nordsieck form. In development.\nMEBDF2 - The second order Modified Extended BDF method, which has improved stability properties over the standard BDF. Fixed timestep only."
},

{
    "location": "solvers/ode_solve.html#Implicit-Strong-Stability-Preserving-Runge-Kutta-Methods-for-Hyperbolic-PDEs-(Conservation-Laws)-1",
    "page": "ODE Solvers",
    "title": "Implicit Strong-Stability Preserving Runge-Kutta Methods for Hyperbolic PDEs (Conservation Laws)",
    "category": "section",
    "text": "SSPSDIRK2 - A second order A-L stable symplectic SDIRK method with the strong stability preserving (SSP) property (SSP coefficient 2). Fixed timestep only."
},

{
    "location": "solvers/ode_solve.html#Extra-Options-1",
    "page": "ODE Solvers",
    "title": "Extra Options",
    "category": "section",
    "text": "All of the Rosenbrock and SDIRK methods allow for specification of linsolve: the linear solver which is used. For more information on specifying the linear solver, see the manual page on solver specification.The following methods allow for specification of nlsolve: the nonlinear solver which is used:GenericImplicitEuler\nGenericTrapezoidNote that performance overload information (Jacobians etc.) are not used in this mode. This can control autodifferentiation of the Jacobian as well. For more information on specifying the nonlinear solver, see the manual page on solver specification.Additionally, the Rosenbrock and SDIRK methods have differentiation controls. In each of these, autodiff can be set to turn on/off autodifferentiation, and chunk_size can be used to set the chunksize of the Dual  numbers (see the documentation for ForwardDiff.jl for details). In addition, the Rosenbrock and SDIRK methods can set diff_type, which is the type of numerical differentiation that is used (when autodifferentiation is disabled). The choices are Val{:central}, Val{:forward} or Val{:complex}.Examples:sol = solve(prob,Rosenbrock23()) # Standard, uses autodiff\nsol = solve(prob,Rosenbrock23(chunk_size=10)) # Autodiff with chunksize of 10\nsol = solve(prob,Rosenbrock23(autodiff=false)) # Numerical differentiation with central differencing\nsol = solve(prob,Rosenbrock23(autodiff=false,diff_type=Val{:forward})) # Numerical differentiation with forward differencing"
},

{
    "location": "solvers/ode_solve.html#Tableau-Method-1",
    "page": "ODE Solvers",
    "title": "Tableau Method",
    "category": "section",
    "text": "Additionally, there is the tableau method:ExplicitRK - A general Runge-Kutta solver which takes in a tableau. Can be adaptive. Tableaus are specified via the keyword argument tab=tableau. The default tableau is for Dormand-Prince 4/5. Other supplied tableaus can be found in the Supplied Tableaus section.Example usage:alg = ExplicitRK(tableau=constructDormandPrince())\nsolve(prob,alg)"
},

{
    "location": "solvers/ode_solve.html#CompositeAlgorithm-1",
    "page": "ODE Solvers",
    "title": "CompositeAlgorithm",
    "category": "section",
    "text": "One unique feature of OrdinaryDiffEq.jl is the CompositeAlgorithm, which allows you to, with very minimal overhead, design a multimethod which switches between chosen algorithms as needed. The syntax is CompositeAlgorithm(algtup,choice_function) where algtup is a tuple of OrdinaryDiffEq.jl algorithms, and choice_function is a function which declares which method to use in the following step. For example, we can design a multimethod which uses Tsit5() but switches to Vern7() whenever dt is too small:choice_function(integrator) = (Int(integrator.dt<0.001) + 1)\nalg_switch = CompositeAlgorithm((Tsit5(),Vern7()),choice_function)The choice_function takes in an integrator and thus all of the features available in the Integrator Interface can be used in the choice function.A helper algorithm was created for building 2-method automatic switching for stiffness detection algorithms. This is the AutoSwitch algorithm with the following options:AutoSwitch(nonstiffalg::nAlg, stiffalg::sAlg;\n           maxstiffstep=10, maxnonstiffstep=3,\n           nonstifftol::T=9//10, stifftol::T=9//10,\n           dtfac=2.0, stiffalgfirst=false)The nonstiffalg must have an appropriate stiffness estimate built into the method. The stiffalg can receive its estimate from the Jacobian calculation. maxstiffstep is the number of stiffness detects before switching to the stiff algorithm and maxnonstiffstep is vice versa. nonstifftol and stifftol are the tolerances associated with the stiffness comparison against the stability region. Decreasing stifftol makes switching to the non-stiff algorithm less likely. Decreasing nonstifftol makes switching to the stiff algorithm more likely. dtfac is the factor that dt is changed when switching: multiplied when going from non-stiff to stiff and divided when going stiff to non-stiff. stiffalgfirst denotes whether the first step should use the stiff algorithm."
},

{
    "location": "solvers/ode_solve.html#Pre-Built-Stiffness-Detecting-and-Auto-Switching-Algorithms-1",
    "page": "ODE Solvers",
    "title": "Pre-Built Stiffness Detecting and Auto-Switching Algorithms",
    "category": "section",
    "text": "These methods require a Autoalg(stiffalg) to be chosen as the method to switch to when the ODE is stiff. It can be any of the OrdinaryDiffEq.jl one-step stiff methods and has all of the arguments of the AutoSwitch algorithm.AutoTsit5 - Tsit5 with automated switching.\nAutoDP5 - DP5 with automated switching.\nAutoVern6 - Vern6 with automated switching.\nAutoVern7 - Vern7 with automated switching.\nAutoVern8 - Vern8 with automated switching.\nAutoVern9 - Vern9 with automated switching.Example:tsidas_alg = AutoTsit5(Rodas5())\nsol = solve(prob,tsidas_alg)\n\ntsidas_alg = AutoTsit5(Rodas5(),nonstifftol = 11/10)Is the Tsit5 method with automatic switching to Rodas5."
},

{
    "location": "solvers/ode_solve.html#Sundials.jl-1",
    "page": "ODE Solvers",
    "title": "Sundials.jl",
    "category": "section",
    "text": "Note that this setup is not automatically included with DifferentialEquations.jl. To use the following algorithms, you must install and use Sundials.jl:]add Sundials\nusing SundialsThe Sundials suite is built around multistep methods. These methods are more efficient than other methods when the cost of the function calculations is really high, but for less costly functions the cost of nurturing the timestep overweighs the benefits. However, the BDF method is a classic method for stiff equations and \"generally works\".CVODE_BDF - CVode Backward Differentiation Formula (BDF) solver.\nCVODE_Adams - CVode Adams-Moulton solver.\nARKODE - Explicit and ESDIRK Runge-Kutta methods of orders 2-8 depending on choice of options.The Sundials algorithms all come with a 3rd order Hermite polynomial interpolation. Note that the constructors for the Sundials algorithms take two main arguments:method - This is the method for solving the implicit equation. For BDF this defaults to :Newton while for Adams this defaults to :Functional. These choices match the recommended pairing in the Sundials.jl manual. However, note that using the :Newton method may take less iterations but requires more memory than the :Function iteration approach.\nlinearsolver - This is the linear solver which is used in the :Newton method.The choices for the linear solver are::Dense - A dense linear solver.\n:Band - A solver specialized for banded Jacobians. If used, you must set the position of the upper and lower non-zero diagonals via jac_upper and jac_lower.\n:Diagonal - This method is specialized for diagonal Jacobians.\n:GMRES - A GMRES method. Recommended first choice Krylov method\n:BCG - A Biconjugate gradient method.\n:PCG - A preconditioned conjugate gradient method. Only for symmetric linear systems.\n:TFQMR - A TFQMR method.\n:KLU - A sparse factorization method. Requires that the user specifies a Jacobian. The Jacobian must be set as a sparse matrix in the ODEProblem type.Example:CVODE_BDF() # BDF method using Newton + Dense solver\nCVODE_BDF(method=:Functional) # BDF method using Functional iterations\nCVODE_BDF(linear_solver=:Band,jac_upper=3,jac_lower=3) # Banded solver with nonzero diagonals 3 up and 3 down\nCVODE_BDF(linear_solver=:BCG) # Biconjugate gradient method                                   The main options for ARKODE are the choice between explicit and implicit and the method order, given via:ARKODE(Sundials.Explicit()) # Solve with explicit tableau of default order 4\nARKODE(Sundials.Implicit(),order = 3) # Solve with explicit tableau of order 3The order choices for explicit are 2 through 8 and for implicit 3 through 5. Specific methods can also be set through the etable and itable options for explicit and implicit tableaus respectively. The available tableaus are:etable:HEUN_EULER_2_1_2: 2nd order Heun\'s method\nBOGACKI_SHAMPINE_4_2_3:\nARK324L2SA_ERK_4_2_3: explicit portion of Kennedy and Carpenter\'s 3rd order method\nZONNEVELD_5_3_4: 4th order explicit method\nARK436L2SA_ERK_6_3_4: explicit portion of Kennedy and Carpenter\'s 4th order method\nSAYFY_ABURUB_6_3_4: 4th order explicit method\nCASH_KARP_6_4_5: 5th order explicit method\nFEHLBERG_6_4_5: Fehlberg\'s classic 5th order method\nDORMAND_PRINCE_7_4_5: the classic 5th order Dormand-Prince method\nARK548L2SA_ERK_8_4_5: explicit portion of Kennedy and Carpenter\'s 5th order method\nVERNER_8_5_6: Verner\'s classic 5th order method\nFEHLBERG_13_7_8: Fehlberg\'s 8th order methoditable:SDIRK_2_1_2: An A-B-stable 2nd order SDIRK method\nBILLINGTON_3_3_2: A second order method with a 3rd order error predictor of less stability\nTRBDF2_3_3_2: The classic TR-BDF2 method\nKVAERNO_4_2_3: an L-stable 3rd order ESDIRK method\nARK324L2SA_DIRK_4_2_3: implicit portion of Kennedy and Carpenter\'s 3th order method\nCASH_5_2_4: Cash\'s 4th order L-stable SDIRK method\nCASH_5_3_4: Cash\'s 2nd 4th order L-stable SDIRK method\nSDIRK_5_3_4: Hairer\'s 4th order SDIRK method\nKVAERNO_5_3_4: Kvaerno\'s 4th order ESDIRK method\nARK436L2SA_DIRK_6_3_4: implicit portion of Kennedy and Carpenter\'s 4th order method\nKVAERNO_7_4_5: Kvaerno\'s 5th order ESDIRK method\nARK548L2SA_DIRK_8_4_5: implicit portion of Kennedy and Carpenter\'s 5th order methodThese can be set for example via:ARKODE(Sundials.Explicit(),etable = Sundials.DORMAND_PRINCE_7_4_5)\nARKODE(Sundials.Implicit(),itable = Sundials.KVAERNO_4_2_3)All of the additional options are available. The full constructor is:CVODE_BDF(;method=:Newton,linear_solver=:Dense,\n          jac_upper=0,jac_lower=0,\n          stored_upper = jac_upper + jac_lower,\n          non_zero=0,krylov_dim=0,\n          stability_limit_detect=false,\n          max_hnil_warns = 10,\n          max_order = 5,\n          max_error_test_failures = 7,\n          max_nonlinear_iters = 3,\n          max_convergence_failures = 10)\n\nCVODE_Adams(;method=:Functional,linear_solver=:None,\n            jac_upper=0,jac_lower=0,\n            stored_upper = jac_upper + jac_lower,\n            krylov_dim=0,\n            stability_limit_detect=false,\n            max_hnil_warns = 10,\n            max_order = 12,\n            max_error_test_failures = 7,\n            max_nonlinear_iters = 3,\n            max_convergence_failures = 10)\n\nARKODE(stiffness=Sundials.Implicit();\n      method=:Newton,linear_solver=:Dense,\n      jac_upper=0,jac_lower=0,stored_upper = jac_upper+jac_lower,\n      non_zero=0,krylov_dim=0,\n      max_hnil_warns = 10,\n      max_error_test_failures = 7,\n      max_nonlinear_iters = 3,\n      max_convergence_failures = 10,\n      predictor_method = 0,\n      nonlinear_convergence_coefficient = 0.1,\n      dense_order = 3,\n      order = 4,\n      set_optimal_params = false,\n      crdown = 0.3,\n      dgmax = 0.2,\n      rdiv = 2.3,\n      msbp = 20,\n      adaptivity_method = 0\n      )See the CVODE manual and the ARKODE manual for details on the additional options."
},

{
    "location": "solvers/ode_solve.html#ODEInterface.jl-1",
    "page": "ODE Solvers",
    "title": "ODEInterface.jl",
    "category": "section",
    "text": "The ODEInterface algorithms are the classic Fortran algorithms. While the non-stiff algorithms are superseded by the more featured and higher performance Julia implementations from OrdinaryDiffEq.jl, the stiff solvers such as radau are some of the most efficient methods available (but are restricted for use on arrays of Float64).Note that this setup is not automatically included with DifferentialEquations.jl. To use the following algorithms, you must install and use ODEInterfaceDiffEq.jl:]add ODEInterfaceDiffEq\nusing ODEInterfaceDiffEqdopri5 - Hairer\'s classic implementation of the Dormand-Prince 4/5 method.\ndop853 - Explicit Runge-Kutta 8(5,3) by Dormand-Prince.\nodex - GBS extrapolation-algorithm based on the midpoint rule.\nseulex - Extrapolation-algorithm based on the linear implicit Euler method.\nradau - Implicit Runge-Kutta (Radau IIA) of variable order between 5 and 13.\nradau5 - Implicit Runge-Kutta method (Radau IIA) of order 5.\nrodas - Rosenbrock 4(3) method.\nddeabm - Adams-Bashforth-Moulton Predictor-Corrector method (order between 1 and 12)\nddebdf - Backward Differentiation Formula (orders between 1 and 5)Note that while the output only has a linear interpolation, a higher order interpolation is used for intermediate dense output for saveat and for event handling."
},

{
    "location": "solvers/ode_solve.html#LSODA.jl-1",
    "page": "ODE Solvers",
    "title": "LSODA.jl",
    "category": "section",
    "text": "This setup provides a wrapper to the algorithm LSODA, a well-known method which uses switching to solve both stiff and non-stiff equations.lsoda - The LSODA wrapper algorithm.Note that this setup is not automatically included with DifferentialEquaitons.jl. To use the following algorithms, you must install and use LSODA.jl:]add LSODA\nusing LSODA"
},

{
    "location": "solvers/ode_solve.html#SimpleDiffEq.jl-1",
    "page": "ODE Solvers",
    "title": "SimpleDiffEq.jl",
    "category": "section",
    "text": "This setup provides access to simplified versions of a few ODE solvers. They mostly exist for experimentation, but offer shorter compile times. They have limitations compared to OrdinaryDiffEq.jl and are not generally faster.SimpleTsit5 - A fixed timestep integrator form of Tsit5. Not compatible with events.\nSimpleATsit5 - An adaptive Tsit5 with an interpolation in its simplest form. Not compatible with events.\nGPUSimpleATsit5 - A version of SimpleATsit5 without the integrator interface. Only allows solve.Note that this setup is not automatically included with DifferentialEquaitons.jl. To use the following algorithms, you must install and use SimpleDiffEq.jl:]add SimpleDiffEq\nusing SimpleDiffEq"
},

{
    "location": "solvers/ode_solve.html#ODE.jl-1",
    "page": "ODE Solvers",
    "title": "ODE.jl",
    "category": "section",
    "text": "Note that this setup is not automatically included with DifferentialEquaitons.jl. To use the following algorithms, you must install and use ODE.jl:]add ODE\nusing ODEode23 - Bogacki-Shampine\'s order 2/3 Runge-Kutta  method\node45 - A Dormand-Prince order 4/5 Runge-Kutta method\node23s - A modified Rosenbrock order 2/3 method due to Shampine\node78 - A Fehlburg order 7/8 Runge-Kutta method\node4 - The classic Runge-Kutta order 4 method\node4ms - A fixed-step, fixed order Adams-Bashforth-Moulton method†\node4s - A 4th order Rosenbrock method due to Shampine†: Does not step to the interval endpoint. This can cause issues with discontinuity detection, and discrete variables need to be updated appropriately."
},

{
    "location": "solvers/ode_solve.html#MATLABDiffEq.jl-1",
    "page": "ODE Solvers",
    "title": "MATLABDiffEq.jl",
    "category": "section",
    "text": "These algorithms require that the problem was defined using a ParameterizedFunction via the @ode_def macro. Note that this setup is not automatically included with DifferentialEquaitons.jl. To use the following algorithms, you must install and use MATLABDiffEq.jl:]add https://github.com/JuliaDiffEq/MATLABDiffEq.jl\nusing MATLABDiffEqThis requires a licensed MATLAB installation. The available methods are:ode23\node45\node113\node23s\node23t\node23tb\node15s\node15iFor more information on these algorithms, see the MATLAB documentation."
},

{
    "location": "solvers/ode_solve.html#GeometricIntegrators.jl-1",
    "page": "ODE Solvers",
    "title": "GeometricIntegrators.jl",
    "category": "section",
    "text": ""
},

{
    "location": "solvers/ode_solve.html#Note:-This-package-currently-segfaults-on-non-Linux-Julia-v1.0!-1",
    "page": "ODE Solvers",
    "title": "Note: This package currently segfaults on non-Linux Julia v1.0!",
    "category": "section",
    "text": "GeometricIntegrators.jl is a set of fixed timestep algorithms written in Julia. Note that this setup is not automatically included with DifferentialEquaitons.jl. To use the following algorithms, you must install and use GeometricIntegratorsDiffEq.jl:]add https://github.com/JuliaDiffEq/GeometricIntegratorsDiffEq.jl\nusing GeometricIntegratorsDiffEqGIEuler - 1st order Euler method\nGIMidpoint - 2nd order explicit midpoint method\nGIHeun - 2nd order Heun\'s method\nGIKutta - 3rd order Kutta\'s method\nGIERK4 - standard 4th order Runge-Kutta\nGIERK438 - 4th order Runge-Kutta, 3/8\'s rule\nGIImplicitEuler - 1st order implicit Euler method\nGIImplicitMidpoint - 2nd order implicit midpoint method\nGIRadIIA2 - 2-stage order 3 Radau-IIA\nGIRadIIA3 - 3-stage order 5 Radau-IIA\nGISRK3 - 3-stage order 4 symmetric Runge-Kutta method\nGIGLRK(s) - Gauss-Legendre Runge-Kutta method of order 2sNote that all of these methods require the user supplies dt."
},

{
    "location": "solvers/ode_solve.html#BridgeDiffEq.jl-1",
    "page": "ODE Solvers",
    "title": "BridgeDiffEq.jl",
    "category": "section",
    "text": "Bridge.jl is a set of fixed timestep algorithms written in Julia. These methods are made and optimized for out-of-place functions on immutable (static vector) types. Note that this setup is not automatically included with DifferentialEquaitons.jl. To use the following algorithms, you must install and use BridgeDiffEq.jl:]add https://github.com/JuliaDiffEq/BridgeDiffEq.jl\nusing BridgeDiffEqBridgeR3 - 3rd order Ralston method\nBridgeBS3 - 3rd order Bogacki-Shampine method"
},

{
    "location": "solvers/ode_solve.html#TaylorIntegration.jl-1",
    "page": "ODE Solvers",
    "title": "TaylorIntegration.jl",
    "category": "section",
    "text": "TaylorIntegration.jl is a pure-Julia implementation of an adaptive order Taylor series method for high accuracy integration of ODEs. These methods are optimized when the absolute tolerance is required to be very low. Note that this setup is not automatically included with DifferentialEquaitons.jl. To use the following algorithms, you must install and use TaylorIntegration.jl:]add TaylorIntegration\nusing TaylorIntegrationTaylorMethod(order) - Taylor integration method with maximal order (required)Note: this method is much faster if you put @taylorize on your derivative function!"
},

{
    "location": "solvers/ode_solve.html#QuDiffEq.jl-1",
    "page": "ODE Solvers",
    "title": "QuDiffEq.jl",
    "category": "section",
    "text": "QuDiffEq.jl is a pacakge for solving differential equations using quantum algorithm. It makes use of the Yao framework for simulating quantum circuits.Note that this setup is not automatically included with DifferentialEquaitons.jl. To use the following algorithms, you must install and use QuDiffEq.jl:]add https://github.com/QuantumBFS/QuDiffEq.jl\nusing QuDiffEqQuLDE(k) - Algorithm based on truncated Taylor series. The method linearizes a system of non-linear differential equations and solves the resultant by means of a quantum circuit. k selects the order in the Taylor series aprroximation (for the quantum circuit).\nQuNLDE(k,ϵ)- Algorithm uses forward Euler to solve quadratc differential equations. k selects the order in the Taylor series aprroximation (for the quantum circuit). ϵ sets the precision for Hamiltonian evolution."
},

{
    "location": "solvers/ode_solve.html#NeuralNetDiffEq.jl-1",
    "page": "ODE Solvers",
    "title": "NeuralNetDiffEq.jl",
    "category": "section",
    "text": "This method trains a neural network using Flux.jl to approximate the solution of the ODE. Currently this method isn\'t competitive but it is a fun curiosity that will be improved with future integration with Zygote.Note that this setup is not automatically included with DifferentialEquaitons.jl. To use the following algorithms, you must install and use NeuralNetDiffEq.jl:]add NeuralNetDiffEq\nusing NeuralNetDiffEqnnode(chain,opt=ADAM(0.1)) - Defines a neural network solver which utilizes a Flux.jl chain under the hood which must be supplied by the user. Defaults to using the ADAM optimization method, but the user can pass any Flux.jl optimizer."
},

{
    "location": "solvers/ode_solve.html#List-of-Supplied-Tableaus-1",
    "page": "ODE Solvers",
    "title": "List of Supplied Tableaus",
    "category": "section",
    "text": "A large variety of tableaus have been supplied by default via DiffEqDevTools.jl. The list of tableaus can be found in the developer docs. To use them, note you must install the library:]add DiffEqDevTools\nusing DiffEqDevToolsFor the most useful and common algorithms, a hand-optimized version is supplied in OrdinaryDiffEq.jl which is recommended for general uses (i.e. use DP5 instead of ExplicitRK with tableau=constructDormandPrince()). However, these serve as a good method for comparing between tableaus and understanding the pros/cons of the methods. Implemented are every published tableau (that I know exists). Note that user-defined tableaus also are accepted. To see how to define a tableau, checkout the premade tableau source code. Tableau docstrings should have appropriate citations (if not, file an issue).Plot recipes are provided which will plot the stability region for a given tableau.[1]: Koskela, A. (2015). Approximating the matrix exponential of an advection-diffusion operator using the incomplete orthogonalization method. In Numerical Mathematics and Advanced Applications-ENUMATH 2013 (pp. 345-353). Springer, Cham."
},

{
    "location": "solvers/dynamical_solve.html#",
    "page": "Dynamical, Hamiltonian, and 2nd Order ODE Solvers",
    "title": "Dynamical, Hamiltonian, and 2nd Order ODE Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "solvers/dynamical_solve.html#Dynamical,-Hamiltonian,-and-2nd-Order-ODE-Solvers-1",
    "page": "Dynamical, Hamiltonian, and 2nd Order ODE Solvers",
    "title": "Dynamical, Hamiltonian, and 2nd Order ODE Solvers",
    "category": "section",
    "text": "Dynamical ODEs, such as those arising from Hamiltonians or second order ordinary differential equations, give rise to a special structure that can be specialized on in the solver for more efficiency. These algorithms require an ODE defined in the following ways:DynamicalODEProblem{isinplace}(f1,f2,u0,v0,tspan;kwargs...)\nSecondOrderODEProblem{isinplace}(f,du0,u0,tspan;kwargs...)\nHamiltonianProblem{T}(H,q0,p0,tspan;kwargs...)These correspond to partitioned equations of motion:fracdvdt = f_1(tu) \nfracdudt = f_2(v) The functions should be specified as f1(dv,v,u,p,t) and f2(du,v,u,p,t) (in the inplace form), where f1 is independent of v (unless specified by the solver), and f2 is independent of t and u. This includes discretizations arising from SecondOrderODEProblems where the velocity is not used in the acceleration function, and Hamiltonians where the potential is (or can be) time-dependent but the kinetic energy is only dependent on v.Note that some methods assume that the integral of f2 is a quadratic form. That means that f2=v\'*M*v, i.e. int f_2 = frac12 m v^2, giving du = v. This is equivalent to saying that the kinetic energy is related to v^2. The methods which require this assumption will lose accuracy if this assumption is violated. Methods listed below make note of this requirement with \"Requires quadratic kinetic energy\"."
},

{
    "location": "solvers/dynamical_solve.html#Recommendations-1",
    "page": "Dynamical, Hamiltonian, and 2nd Order ODE Solvers",
    "title": "Recommendations",
    "category": "section",
    "text": "When energy conservation is required, use a symplectic method. Otherwise the Runge-Kutta-Nyström methods will be more efficient. Energy is mostly conserved by Runge-Kutta-Nyström methods, but is not conserved for long time integrations. Thus it is suggested that for shorter integrations you use Runge-Kutta-Nyström methods as well.As a go-to method for efficiency, DPRKN6 is a good choice. DPRKN12 is a good choice when high accuracy, like tol<1e-10 is necessary. However, DPRKN6 is the only Runge-Kutta-Nyström method with a higher order interpolant (all default to order 3 Hermite, whereas DPRKN6 is order 6th interpolant) and thus in cases where interpolation matters (ex: event handling) one should use DPRKN6. For very smooth problems with expensive acceleration function evaluations, IRKN4 can be a good choice as it minimizes the number of evaluations.For symplectic methods, higher order algorithms are the most efficient when higher accuracy is needed, and when less accuracy is needed lower order methods do better. Optimized efficiency methods take more steps and thus have more force calculations for the same order, but have smaller error. Thus the \"optimized efficiency\" algorithms are recommended if your force calculation is not too sufficiency large, while the other methods are recommend when force calculations are really large (for example, like in MD simulations VelocityVerlet is very popular since it only requires one force calculation per timestep). A good go-to method would be McAte5, and a good high order choice is KahanLi8."
},

{
    "location": "solvers/dynamical_solve.html#Standard-ODE-Integrators-1",
    "page": "Dynamical, Hamiltonian, and 2nd Order ODE Solvers",
    "title": "Standard ODE Integrators",
    "category": "section",
    "text": "The standard ODE integrators will work on Dynamical ODE problems via a transformation to a first-order ODE. See the ODE solvers page for more details."
},

{
    "location": "solvers/dynamical_solve.html#Specialized-OrdinaryDiffEq.jl-Integrators-1",
    "page": "Dynamical, Hamiltonian, and 2nd Order ODE Solvers",
    "title": "Specialized OrdinaryDiffEq.jl Integrators",
    "category": "section",
    "text": "Unless otherwise specified, the OrdinaryDiffEq algorithms all come with a 3rd order Hermite polynomial interpolation. The algorithms denoted as having a \"free\" interpolation means that no extra steps are required for the interpolation. For the non-free higher order interpolating functions, the extra steps are computed lazily (i.e. not during the solve)."
},

{
    "location": "solvers/dynamical_solve.html#Runge-Kutta-Nyström-Integrators-1",
    "page": "Dynamical, Hamiltonian, and 2nd Order ODE Solvers",
    "title": "Runge-Kutta-Nyström Integrators",
    "category": "section",
    "text": "Nystrom4: 4th order explicit Runge-Kutta-Nyström method. Allows acceleration to depend on velocity. Fixed timestep only.\nIRKN3: 4th order explicit two-step Runge-Kutta-Nyström method. Fixed timestep only.\nIRKN4: 4th order explicit two-step Runge-Kutta-Nyström method. Can be more efficient for smooth problems. Fixed timestep only.\nERKN4: 4th order Runge-Kutta-Nyström method which is integrates the periodic properties of the harmonic oscillator exactly. Gets extra efficiency on periodic problems.\nERKN5: 5th order Runge-Kutta-Nyström method which is integrates the periodic properties of the harmonic oscillator exactly. Gets extra efficiency on periodic problems.\nNystrom4VelocityIndependent: 4th order explicit Runge-Kutta-Nyström method. Fixed timestep only.\nNystrom5VelocityIndependent: 5th order explicit Runge-Kutta-Nyström method. Fixed timestep only.\nDPRKN6: 6th order explicit adaptive Runge-Kutta-Nyström method. Free 6th order interpolant.\nDPRKN8: 8th order explicit adaptive Runge-Kutta-Nyström method.\nDPRKN12: 12th order explicit adaptive Runge-Kutta-Nyström method."
},

{
    "location": "solvers/dynamical_solve.html#Symplectic-Integrators-1",
    "page": "Dynamical, Hamiltonian, and 2nd Order ODE Solvers",
    "title": "Symplectic Integrators",
    "category": "section",
    "text": "Note that all symplectic integrators are fixed timestep only.SymplecticEuler: First order explicit symplectic integrator\nVelocityVerlet: 2nd order explicit symplectic integrator.\nVerletLeapfrog: 2nd order explicit symplectic integrator.\nPseudoVerletLeapfrog: 2nd order explicit symplectic integrator.\nMcAte2: Optimized efficiency 2nd order explicit symplectic integrator.\nRuth3: 3rd order explicit symplectic integrator.\nMcAte3: Optimized efficiency 3rd order explicit symplectic integrator.\nCandyRoz4: 4th order explicit symplectic integrator.\nMcAte4: 4th order explicit symplectic integrator. Requires quadratic kinetic energy.\nCalvoSanz4: Optimized efficiency 4th order explicit symplectic integrator.\nMcAte42: 4th order explicit symplectic integrator. (Broken)\nMcAte5: Optimized efficiency 5th order explicit symplectic integrator. Requires quadratic kinetic energy\nYoshida6: 6th order explicit symplectic integrator.\nKahanLi6: Optimized efficiency 6th order explicit symplectic integrator.\nMcAte8: 8th order explicit symplectic integrator.\nKahanLi8: Optimized efficiency 8th order explicit symplectic integrator.\nSofSpa10: 10th order explicit symplectic integrator."
},

{
    "location": "solvers/dynamical_solve.html#GeometricIntegrators.jl-1",
    "page": "Dynamical, Hamiltonian, and 2nd Order ODE Solvers",
    "title": "GeometricIntegrators.jl",
    "category": "section",
    "text": "GeometricIntegrators.jl is a set of fixed timestep algorithms written in Julia. Note that this setup is not automatically included with DifferentialEquaitons.jl. To use the following algorithms, you must install and use GeometricIntegratorsDiffEq.jl:Pkg.clone(\"https://github.com/JuliaDiffEq/GeometricIntegratorsDiffEq.jl\")\nusing GeometricIntegratorsDiffEqGISymplecticEulerA - First order explicit symplectic Euler A\nGISymplecticEulerB - First order explicit symplectic Euler B\nGILobattoIIIAIIIB2 - Second order Gauss-Labatto-IIIA-IIIB\nGILobattoIIIBIIIA2 - Second order Gauss-Labatto-IIIB-IIIA"
},

{
    "location": "solvers/split_ode_solve.html#",
    "page": "Split ODE Solvers",
    "title": "Split ODE Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "solvers/split_ode_solve.html#Split-ODE-Solvers-1",
    "page": "Split ODE Solvers",
    "title": "Split ODE Solvers",
    "category": "section",
    "text": "The solvers which are available for a SplitODEProblem depend on the input linearity and number of components. Each solver has functional form (or many) that it allows."
},

{
    "location": "solvers/split_ode_solve.html#Implicit-Explicit-(IMEX)-ODE-1",
    "page": "Split ODE Solvers",
    "title": "Implicit-Explicit (IMEX) ODE",
    "category": "section",
    "text": "The Implicit-Explicit (IMEX) ODE is a SplitODEProblem with two functions:fracdudt =  f_1(tu) + f_2(tu)where the first function is the stiff part and the second function is the non-stiff part (implicit integration on f1, explicit integration on f2)."
},

{
    "location": "solvers/split_ode_solve.html#Recommended-Methods-1",
    "page": "Split ODE Solvers",
    "title": "Recommended Methods",
    "category": "section",
    "text": "The recommended method in most cases is KenCarp4. In cases of extreme stiffness or for high tolerances, KenCarp3 can be a good choice. The ARKODE methods are generally inefficient and diverge unless the options are tweaked to match the problem, though for large enough PDEs the ARKODE method with linear_solver=:GMRES is a good choice."
},

{
    "location": "solvers/split_ode_solve.html#OrdinaryDiffEq.jl-1",
    "page": "Split ODE Solvers",
    "title": "OrdinaryDiffEq.jl",
    "category": "section",
    "text": "SplitEuler: 1st order fully explicit method. Used for testing accuracy of splits.\nIMEXEuler : 1st order explicit Euler mixed with implicit Euler. Fixed time step only.\nCNAB2: Crank-Nicholson Adams Bashforth Order 2. Fixed time step only.\nCNLF: Crank-Nicholson Leapfrog of Order 2. Fixed time step only.\nSBDF2 : 2nd order IMEX BDF method. Fixed time step only.\nSBDF3 : 3rd order IMEX BDF method. Fixed time step only. In development.\nSBDF4 : 4th order IMEX BDF method. Fixed time step only. In development.\nKenCarp3: An A-L stable stiffly-accurate 3rd order ESDIRK method.\nKenCarp4: An A-L stable stiffly-accurate 4rd order ESDIRK method.\nKenCarp5: An A-L stable stiffly-accurate 5rd order ESDIRK method."
},

{
    "location": "solvers/split_ode_solve.html#Sundials.jl-1",
    "page": "Split ODE Solvers",
    "title": "Sundials.jl",
    "category": "section",
    "text": "ARKODE: An additive Runge-Kutta method. Order between 3rd and 5th. For a list of available options, please see its ODE solver page"
},

{
    "location": "solvers/split_ode_solve.html#Semilinear-ODE-1",
    "page": "Split ODE Solvers",
    "title": "Semilinear ODE",
    "category": "section",
    "text": "The Semilinear ODE is a SplitODEProblem with one linear operator and one nonlinear function:fracdudt =  Au + f(tu)See the documentation page for DiffEqOperator for details about how to define linear operators from a matrix or finite difference discretization of derivative operators.The appropriate algorithms for this form are:"
},

{
    "location": "solvers/split_ode_solve.html#OrdinaryDiffEq.jl-2",
    "page": "Split ODE Solvers",
    "title": "OrdinaryDiffEq.jl",
    "category": "section",
    "text": "GenericIIF1 - First order Implicit Integrating Factor method. Fixed timestepping only. Doesn\'t support Krylov approximation.\nGenericIIF2 - Second order Implicit Integrating Factor method. Fixed timestepping only. Doesn\'t support Krylov approximation.\nLawsonEuler - First order exponential Euler scheme. Fixed timestepping only.\nNorsettEuler - First order exponential-RK scheme. Fixed timestepping only. Alias: ETD1.\nETD2 - Second order Exponential Time Differencing method (in development). Fixed timestepping only. Doesn\'t support Krylov approximation.\nETDRK2 - 2nd order exponential-RK scheme. Fixed timestepping only.\nETDRK3 - 3rd order exponential-RK scheme. Fixed timestepping only.\nETDRK4 - 4th order exponential-RK scheme. Fixed timestepping only.\nHochOst4 - 4th order exponential-RK scheme with stiff order 4. Fixed timestepping only.\nExprb32 - 3rd order adaptive Exponential Rosenbrock scheme (in development).\nExprb43 - 4th order adaptive Exponential Rosenbrock scheme (in development).Note that the generic algorithms GenericIIF1 and GenericIIF2 allow for a choice of nlsolve.By default, the exponential methods cache matrix functions such as exp(dt*A) to accelerate the time stepping for small systems. For large systems, using Krylov-based versions of the methods can allow for lazy calculation of exp(dt*A)*v and similar entities, and thus improve performance.To tell a solver to use Krylov methods, pass krylov=true to its constructor. You can also manually set the size of the Krylov subspace by setting the m parameter, which defaults to 30. For exampleLawsonEuler(krylov=true, m=50)constructs a Lawson-Euler method which uses a size-50 Krylov subspace. Note that m only sets an upper bound to the Krylov subspace size. If a convergence criterion is met (determined by the reltol of the integrator), \"happy breakdown\" will occur and the Krylov subspace will only be constructed partially.For more advanced control over the Krylov algorithms, you can change the length of the incomplete orthogonalization procedure (IOP) [1] by setting the iop parameter in the constructor. By default, IOP is turned off and full Arnoldi iteration is used. Note that if the linear operator is hermitian, then the Lanczos algorithm will always be used and IOP setting is ignored.[1]: Koskela, A. (2015). Approximating the matrix exponential of an advection-diffusion operator using the incomplete orthogonalization method. In Numerical Mathematics and Advanced Applications-ENUMATH 2013 (pp. 345-353). Springer, Cham."
},

{
    "location": "solvers/steady_state_solve.html#",
    "page": "Steady State Solvers",
    "title": "Steady State Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "solvers/steady_state_solve.html#Steady-State-Solvers-1",
    "page": "Steady State Solvers",
    "title": "Steady State Solvers",
    "category": "section",
    "text": "solve(prob::SteadyStateProblem,alg;kwargs)Solves for the steady states in the problem defined by prob using the algorithm alg. If no algorithm is given, a default algorithm will be chosen."
},

{
    "location": "solvers/steady_state_solve.html#Recommended-Methods-1",
    "page": "Steady State Solvers",
    "title": "Recommended Methods",
    "category": "section",
    "text": "DynamicSS is a good choice if you think you may have multiple steady states or a bad initial guess. SSRootfind can be faster if you have a good initial guess. For DynamicSS, in many cases an adaptive stiff solver, like a Rosenbrock method (Rodas5 or CVODE_BDF), is a good way to allow for very large time steps as the steady state approaches. Note that if you use CVODE_BDF you may need to give a starting dt via dt=....."
},

{
    "location": "solvers/steady_state_solve.html#Full-List-of-Methods-1",
    "page": "Steady State Solvers",
    "title": "Full List of Methods",
    "category": "section",
    "text": ""
},

{
    "location": "solvers/steady_state_solve.html#SteadyStateDiffEq.jl-1",
    "page": "Steady State Solvers",
    "title": "SteadyStateDiffEq.jl",
    "category": "section",
    "text": "SSRootfind : Uses a rootfinding algorithm to find a steady state. Defaults to using NLsolve.jl. A different algorithm can be specified via the nlsolve keyword argument.\nDynamicSS : Uses an ODE solver to find the steady state. Automatically terminates when close to the steady state. DynamicSS(alg;abstol=1e-8,reltol=1e-6,tspan=Inf) requires that an ODE algorithm is given as the first argument.  The absolute and relative tolerances specify the termination conditions on the derivative\'s closeness to zero.  This internally uses the TerminateSteadyState callback from the Callback Library.  The simulated time for which given ODE is solved can be limited by tspan.  If tspan is a number, it is equivalent to passing (zero(tspan), tspan).Example usage:sol = solve(prob,SSRootfind())\nsol = solve(prob,DynamicSS(Tsit5()))\nusing Sundials\nsol = solve(prob,DynamicSS(CVODE_BDF()),dt=1.0)"
},

{
    "location": "solvers/bvp_solve.html#",
    "page": "BVP Solvers",
    "title": "BVP Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "solvers/bvp_solve.html#BVP-Solvers-1",
    "page": "BVP Solvers",
    "title": "BVP Solvers",
    "category": "section",
    "text": ""
},

{
    "location": "solvers/bvp_solve.html#Recomended-Methods-1",
    "page": "BVP Solvers",
    "title": "Recomended Methods",
    "category": "section",
    "text": "GeneralMIRK4 is a good well-rounded method when the problem is not too large. It uses highly stable trust region methods to solve a 4th order fully implicit Runge-Kutta scheme. As an alternative on general BVProblems, the Shooting method paired with an appropriate integrator for the IVP, such as Shooting(Tsit5()), is a flexible and efficient option. This will allow one to combine callbacks/event handling with the BVP solver, and the high order interpolations can be used to define complex boundary conditions. However, Shooting methods can in some cases be prone to sensitivity of the boundary condition.When the problem is a large two-point boundary value problem that is sensitive to the boundary conditions, MIRK4 utilizes a sparse Jacobian to greatly improve the efficiency."
},

{
    "location": "solvers/bvp_solve.html#Full-List-of-Methods-1",
    "page": "BVP Solvers",
    "title": "Full List of Methods",
    "category": "section",
    "text": ""
},

{
    "location": "solvers/bvp_solve.html#BoundaryValueDiffEq.jl-1",
    "page": "BVP Solvers",
    "title": "BoundaryValueDiffEq.jl",
    "category": "section",
    "text": "Shooting - A wrapper over initial value problem solvers.\nGeneralMIRK4 - A 4th order collocation method using an implicit Runge-Kutta tableau solved using a trust region dogleg method from NLsolve.jl.\nMIRK4 - A 4th order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian. Compatible only with two-point boundary value problems."
},

{
    "location": "solvers/jump_solve.html#",
    "page": "Jump Problem Solvers",
    "title": "Jump Problem Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "solvers/jump_solve.html#Jump-Problem-Solvers-1",
    "page": "Jump Problem Solvers",
    "title": "Jump Problem Solvers",
    "category": "section",
    "text": "solve(prob::JumpProblem,alg;kwargs)"
},

{
    "location": "solvers/jump_solve.html#Recommended-Methods-1",
    "page": "Jump Problem Solvers",
    "title": "Recommended Methods",
    "category": "section",
    "text": "A JumpProblem(prob,aggregator,jumps...) come in two forms. The first major form is if it does not have a RegularJump. In this case, it can be solved with any integrator on  prob. However, in the case of a pure JumpProblem (a JumpProblem over a  DiscreteProblem), there are special algorithms available.  The SSAStepper() is an efficient streamlined algorithm for running the  aggregator version of the SSA for pure ConstantRateJump and/or MassActionJump problems. However, it is not compatible with event handling. If events are necessary, then FunctionMap does well.If there is a RegularJump, then specific methods must be used. The current recommended method is SimpleTauLeaping."
},

{
    "location": "solvers/jump_solve.html#Special-Methods-for-Pure-Jump-Problems-1",
    "page": "Jump Problem Solvers",
    "title": "Special Methods for Pure Jump Problems",
    "category": "section",
    "text": "If you are using jumps with a differential equations, use the same methods as in the case of the differential equation solving. However, the following algorithms are optimized for pure jump problems."
},

{
    "location": "solvers/jump_solve.html#DiffEqJump.jl-1",
    "page": "Jump Problem Solvers",
    "title": "DiffEqJump.jl",
    "category": "section",
    "text": "SSAStepper: a stepping algorithm for pure ConstantRateJump and/or MassActionJump JumpProblems. Does not support event handling, but does support saving controls like saveat."
},

{
    "location": "solvers/jump_solve.html#RegularJump-Compatible-Methods-1",
    "page": "Jump Problem Solvers",
    "title": "RegularJump Compatible Methods",
    "category": "section",
    "text": ""
},

{
    "location": "solvers/jump_solve.html#DiffEqJump.jl-2",
    "page": "Jump Problem Solvers",
    "title": "DiffEqJump.jl",
    "category": "section",
    "text": "SimpleTauLeaping: a tau-leaping algorithm for pure RegularJump JumpProblems. Requires a choice of dt.\nRegularSSA: a version of SSA for pure RegularJump JumpProblems."
},

{
    "location": "solvers/sde_solve.html#",
    "page": "SDE Solvers",
    "title": "SDE Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "solvers/sde_solve.html#SDE-Solvers-1",
    "page": "SDE Solvers",
    "title": "SDE Solvers",
    "category": "section",
    "text": ""
},

{
    "location": "solvers/sde_solve.html#Recommended-Methods-1",
    "page": "SDE Solvers",
    "title": "Recommended Methods",
    "category": "section",
    "text": "For most Ito diagonal and scalar noise problems where a good amount of accuracy is required and mild stiffness may be an issue, the SOSRI algorithm should do well. If the problem has additive noise, then SOSRA will be the optimal algorithm. At low tolerances (<1e-4?) SRA3 will be more efficient, though SOSRA is more robust to stiffness. For commutative noise, RKMilCommute is a strong order 1.0 method which utilizes the commutivity property to greatly speed up the Wiktorsson approximation and can choose between Ito and Stratonovich. For non-commutative noise, difficult problems usually require adaptive time stepping in order to be efficient. In this case, LambaEM and LambaEulerHeun are adaptive and handle general non-diagonal problems (for Ito and Stratonovich interpretations respectively). If adaptivity isn\'t necessary, the EM and EulerHeun are good choices (for Ito and Stratonovich interpretations respectively).For stiff problems with additive noise, the high order adaptive method SKenCarp is highly preferred and will solve problems with similar efficiency as ODEs. If possible, stiff problems should be converted to make use of this additive noise solver. If the noise term is large/stiff, then the split-step methods are required in order for the implicit methods to be stable. For Ito in this case, use ISSEM and for Stratonovich use ISSEulerHeun. These two methods can handle any noise form.If the noise term is not too large, for stiff problems with diagonal noise, ImplicitRKMil is the most efficient method and can choose between Ito and Stratonovich. For each of the theta methods, the parameter theta can be chosen. The default is theta=1/2 which will not dampen numerical oscillations and thus is symmetric (and almost symplectic) and will lead to less error when noise is sufficiently small. However, theta=1/2 is not L-stable in the drift term, and thus one can receive more stability (L-stability in the drift term) with theta=1, but with a tradeoff of error efficiency in the low noise case. In addition, the option symplectic=true will turns these methods into an implicit Midpoint extension which is symplectic in distribution but has an accuracy tradeoff."
},

{
    "location": "solvers/sde_solve.html#Mass-Matrices-and-Stochastic-DAEs-1",
    "page": "SDE Solvers",
    "title": "Mass Matrices and Stochastic DAEs",
    "category": "section",
    "text": "The stiff methods can solve stochastic equations with mass matrices (including stochastic DAEs written in mass matrix form) when either symplectic=true or theta=1. These methods interpret the mass matrix equation as:Mu = f(tu)dt + Mg(tu)dW_ti.e. with no mass matrix inversion applied to the g term. Thus these methods apply noise per dependent variable instead of on the combinations of the dependent variables and this is designed for phenomenological noise on the dependent variables (like multiplicative or additive noise)"
},

{
    "location": "solvers/sde_solve.html#Special-Noise-Forms-1",
    "page": "SDE Solvers",
    "title": "Special Noise Forms",
    "category": "section",
    "text": "Some solvers are for specialized forms of noise. Diagonal noise is the default setup. Non-diagonal noise is specified via setting noise_rate_prototype to a matrix in the SDEProblem type. A special form of non-diagonal noise, commutative noise, occurs when the noise satisfies the following condition:sum_i=1^d g_ij_1(tx) fracpartial g_kj_2(tx)partial x_i = sum_i=1^d g_ij_2(tx) fracpartial g_kj_1(tx)partial x_ifor every j_1j_2 and k. Additive noise is when g(tu)=g(t), i.e. is independent of u. Multiplicative noise is g_i(tu)=a_i u."
},

{
    "location": "solvers/sde_solve.html#Special-Keyword-Arguments-1",
    "page": "SDE Solvers",
    "title": "Special Keyword Arguments",
    "category": "section",
    "text": "save_noise: Determines whether the values of W are saved whenever the timeseries is saved. Defaults to true.\ndelta: The delta adaptivity parameter for the natural error estimator. Determines the balance between drift and diffusion error. For more details, see the publication.\nseed: Sets the seed for the random number generator. This overrides any seed set in the SDEProblem."
},

{
    "location": "solvers/sde_solve.html#Full-List-of-Methods-1",
    "page": "SDE Solvers",
    "title": "Full List of Methods",
    "category": "section",
    "text": ""
},

{
    "location": "solvers/sde_solve.html#StochasticDiffEq.jl-1",
    "page": "SDE Solvers",
    "title": "StochasticDiffEq.jl",
    "category": "section",
    "text": "Each of the StochasticDiffEq.jl solvers come with a linear interpolation. Orders are given in terms of strong order."
},

{
    "location": "solvers/sde_solve.html#Nonstiff-Methods-1",
    "page": "SDE Solvers",
    "title": "Nonstiff Methods",
    "category": "section",
    "text": "EM- The Euler-Maruyama method. Strong Order 0.5 in the Ito sense. Has an optional argument split=true for controlling step splitting. When splitting is enabled, the stability with large diffusion eigenvalues is improved. Can handle all forms of noise, including non-diagonal, scalar, and colored noise. Fixed time step only.†\nLambaEM- A modified Euler-Maruyama method with adaptive time stepping with an error estimator based on Lamba and Rackauckas. Has an optional argument split=true for controlling step splitting. When splitting is enabled, the stability with   large diffusion eigenvalues is improved. Strong Order 0.5 in the Ito sense. Can handle all forms of noise, including non-diagonal, scalar, and colored noise.†\nEulerHeun - The Euler-Heun method. Strong Order 0.5 in the Stratonovich sense. Can handle all forms of noise, including non-diagonal, scalar, and colored noise. Fixed time step only.†\nLambaEulerHeun - A modified Euler-Heun method with adaptive time stepping with an error estimator based on Lamba due to Rackauckas. Strong order 0.5 in the Stratonovich sense. Can handle all forms of noise, including non-diagonal, scalar, and colored noise.†\nRKMil - An explicit Runge-Kutta discretization of the strong order 1.0 Milstein method. Defaults to solving the Ito problem, but RKMil(interpretation=:Stratonovich) makes it solve the Stratonovich problem. Only handles scalar and diagonal noise.†\nRKMilCommute - An explicit Runge-Kutta discretization of the strong order 1.0 Milstein method for commutative noise problems. Defaults to solving the Ito problem, but RKMilCommute(interpretation=:Stratonovich) makes it solve the Stratonovich problem. Uses a 1.5/2.0 error estimate for adaptive time stepping.†\nWangLi3SMil_A - fixed step-size explicit 3-stage Milstein methods for Ito problem with strong and weak order 1.0\nWangLi3SMil_B - fixed step-size explicit 3-stage Milstein methods for Ito problem with strong and weak order 1.0\nWangLi3SMil_C - fixed step-size explicit 3-stage Milstein methods for Ito problem with strong and weak order 1.0\nWangLi3SMil_D - fixed step-size explicit 3-stage Milstein methods for Ito problem with strong and weak order 1.0\nWangLi3SMil_E - fixed step-size explicit 3-stage Milstein methods for Ito problem with strong and weak order 1.0\nWangLi3SMil_F - fixed step-size explicit 3-stage Milstein methods for Ito problem with strong and weak order 1.0\nSRA - Adaptive strong order 1.5 methods for additive Ito and Stratonovich SDEs. Default tableau is for SRA1. Can handle diagonal, non-diagonal and scalar additive noise.\nSRI - Adaptive strong order 1.5 methods for diagonal/scalar Ito SDEs. Default tableau is for SRIW1.\nSRIW1 - Adaptive strong order 1.5 and weak order 2.0 for diagonal/scalar Ito SDEs.†\nSRIW2 - Adaptive strong order 1.5 and weak order 3.0 for diagonal/scalar Ito SDEs.†\nSOSRI - Stability-optimized adaptive strong order 1.5 and weak order 2.0 for diagonal/scalar Ito SDEs. Stable at high tolerances and robust to stiffness.†\nSOSRI2 - Stability-optimized adaptive strong order 1.5 and weak order 2.0 for diagonal/scalar Ito SDEs. Stable at high tolerances and robust to stiffness.†\nSRA1 - Adaptive strong order 1.5 for additive Ito and Stratonovich SDEs with weak order 2. Can handle diagonal, non-diagonal, and scalar additive noise.†\nSRA2 - Adaptive strong order 1.5 for additive Ito and Stratonovich SDEs with weak order 2. Can handle diagonal, non-diagonal, and scalar additive noise.†\nSRA3 - Adaptive strong order 1.5 for additive Ito and Stratonovich SDEs with weak order 3. Can handle non-diagonal and scalar additive noise.†\nSOSRA - A stability-optimized adaptive SRA. Strong order 1.5 for additive Ito and Stratonovich SDEs with weak order 2. Can handle diagonal, non-diagonal, and scalar additive noise. Stable at high tolerances and robust to stiffness.†\nSOSRA2 - A stability-optimized adaptive SRA. Strong order 1.5 for additive Ito and Stratonovich SDEs with weak order 2. Can handle diagonal, non-diagonal, and scalar additive noise. Stable at high tolerances and robust to stiffness.†Example usage:sol = solve(prob,SRIW1())3-stage Milstein Methods WangLi3SMil_A, WangLi3SMil_B, WangLi3SMil_D, WangLi3SMil_E and WangLi3SMil_F are currently implemented for 1-dimensional and diagonal noise only."
},

{
    "location": "solvers/sde_solve.html#Tableau-Controls-1",
    "page": "SDE Solvers",
    "title": "Tableau Controls",
    "category": "section",
    "text": "For SRA and SRI, the following option is allowed:tableau: The tableau for an :SRA or :SRI algorithm. Defaults to SRIW1 or SRA1."
},

{
    "location": "solvers/sde_solve.html#S-ROCK-Methods-1",
    "page": "SDE Solvers",
    "title": "S-ROCK Methods",
    "category": "section",
    "text": "SROCK1 - is a fixed step size stabilized explicit method for stiff problems. Defaults to solving th Ito problem but SROCK1(interpretation=:Stratonovich) can make it solve the Stratonovich problem. Strong order of convergence is 0.5 and weak order 1, but is optimised to get order 1 in case os scalar/diagonal noise.\nSROCKEM - is fixed step Euler-Mayurama with first order ROCK stabilization thus can handle stiff problems. Only for Ito problems. Defaults to strong and weak order 1.0, but can solve with weak order 0.5 as SROCKEM(strong_order_1=false). This method can handle 1-dimensional, diagonal and multi-dimensional noise.\nSROCK2 - is a weak second order and strong first order fixed step stabilized method for stiff Ito problems.This method can handle 1-dimensional, diagonal and multi-dimensional noise.\nSKSROCK - is fixed step stabilized explicit method for stiff Ito problems. Strong order 0.5 and weak order 1. This method has a better stability domain then SROCK1. Also it allows special post-processing techniques in case of ergodic dynamical systems, in the context of ergodic Brownian dynamics, to achieve order 2 accuracy. SKSROCK(;post_processing=true) will make use of post processing. By default it doesn\'t use post processing. Post processing is optional and under development. The rest of the method is completely functional and can handle 1-dimensional, diagonal and multi-dimensional noise.  \nTangXiaoSROCK2 - is a fixed step size stabilized expicit method for stiff problems. Only for Ito problems. Weak order of 2 and strog order of 1. Has 5 versions with different stability domains which can be used as TangXiaoSROCK2(version_num=i) where i is 1-5. Under Development."
},

{
    "location": "solvers/sde_solve.html#Stiff-Methods-1",
    "page": "SDE Solvers",
    "title": "Stiff Methods",
    "category": "section",
    "text": "ImplicitEM - An order 0.5 Ito drift-implicit method. This is a theta method which defaults to theta=1/2 or the Trapezoid method on the drift term. This method defaults to symplectic=false, but when true and theta=1/2 this is the implicit Midpoint method on the drift term and is symplectic in distribution. Can handle all forms of noise, including non-diagonal, scalar, and colored noise. Uses a 1.0/1.5 heuristic for adaptive time stepping.\nImplicitEulerHeun - An order 0.5 Stratonovich drift-implicit method. This is a theta method which defaults to theta=1/2 or the Trapezoid method on the drift term. This method defaults to symplectic=false, but when true and theta=1/2 this is the implicit Midpoint method on the drift term and is symplectic in distribution. Can handle all forms of noise, including non-diagonal, scalar, and colored noise. Uses a 1.0/1.5 heuristic for adaptive time stepping.\nImplicitRKMil - An order 1.0 drift-implicit method. This is a theta method which defaults to theta=1/2 or the Trapezoid method on the drift term. Defaults to solving the Ito problem, but ImplicitRKMil(interpretation=:Stratonovich) makes it solve the Stratonovich problem. This method defaults to symplectic=false, but when true and theta=1/2 this is the implicit Midpoint method on the drift term and is symplectic in distribution. Handles diagonal and scalar noise. Uses a 1.5/2.0 heuristic for adaptive time stepping.\nISSEM - An order 0.5 split-step Ito implicit method. It is fully implicit, meaning it can handle stiffness in the noise term. This is a theta method which defaults to theta=1 or the Trapezoid method on the drift term. This method defaults to symplectic=false, but when true and theta=1/2 this is the implicit Midpoint method on the drift term and is symplectic in distribution. Can handle all forms of noise, including non-diagonal, scalar, and colored noise. Uses a 1.0/1.5 heuristic for adaptive time stepping.\nISSEulerHeun - An order 0.5 split-step Stratonovich implicit method. It is fully implicit, meaning it can handle stiffness in the noise term. This is a theta method which defaults to theta=1 or the Trapezoid method on the drift term. This method defaults to symplectic=false, but when true and theta=1/2 this is the implicit Midpoint method on the drift term and is symplectic in distribution. Can handle all forms of noise, including non-diagonal,Q scalar, and colored noise. Uses a 1.0/1.5 heuristic for adaptive time stepping.\nSKenCarp - Adaptive L-stable drfit-implicit strong order 1.5 for additive Ito and Stratonovich SDEs with weak order 2. Can handle diagonal, non-diagonal and scalar additive noise.*†"
},

{
    "location": "solvers/sde_solve.html#Derivative-Based-Methods-1",
    "page": "SDE Solvers",
    "title": "Derivative-Based Methods",
    "category": "section",
    "text": "The following methods require analytic derivatives of the diffusion term.PCEuler - The predictor corrector euler method. Strong Order 0.5 in the Ito sense. Requires the ggprime function, which is defined as\n  textggprime^k(tx) = sum_j=1^m sum_i=1^d g_ij(tx) fracpartial g_kj(tx)partial x_i\nThis can also be understood more intuitively in vector/matrix form as,\ntextggprime(tx) = sum_j=1^m barmathcalJvec g^(j)(tx) vec g^(j)(tx)\nwhere vec g^(j) is the noise vector for the j\'th noise channel and barmathcalJ is the Jacobian of the j\'th   noise vector.\nThe default settings for the drift implicitness is theta=0.5 and the diffusion implicitness is eta=0.5.  "
},

{
    "location": "solvers/sde_solve.html#StochasticCompositeAlgorithm-1",
    "page": "SDE Solvers",
    "title": "StochasticCompositeAlgorithm",
    "category": "section",
    "text": "One unique feature of StochasticDiffEq.jl is the StochasticCompositeAlgorithm, which allows you to, with very minimal overhead, design a multimethod which switches between chosen algorithms as needed. The syntax is StochasticCompositeAlgorithm(algtup,choice_function) where algtup is a tuple of StochasticDiffEq.jl algorithms, and choice_function is a function which declares which method to use in the following step. For example, we can design a multimethod which uses EM() but switches to RKMil() whenever dt is too small:choice_function(integrator) = (Int(integrator.dt<0.001) + 1)\nalg_switch = StochasticCompositeAlgorithm((EM(),RKMil()),choice_function)The choice_function takes in an integrator and thus all of the features available in the Integrator Interface can be used in the choice function."
},

{
    "location": "solvers/sde_solve.html#SimpleDiffEq.jl-1",
    "page": "SDE Solvers",
    "title": "SimpleDiffEq.jl",
    "category": "section",
    "text": "This setup provides access to simplified versions of a few SDE solvers. They mostly exist for experimentation, but offer shorter compile times. They have limitations compared to StochasticDiffEq.jl.SimpleEM - A fixed timestep solve method for Euler-Maruyama. Only works with non-colored Gaussian noise.Note that this setup is not automatically included with DifferentialEquaitons.jl. To use the following algorithms, you must install and use SimpleDiffEq.jl:]add SimpleDiffEq\nusing SimpleDiffEq"
},

{
    "location": "solvers/sde_solve.html#BridgeDiffEq.jl-1",
    "page": "SDE Solvers",
    "title": "BridgeDiffEq.jl",
    "category": "section",
    "text": "Bridge.jl is a set of fixed timestep algorithms written in Julia. These methods are made and optimized for out-of-place functions on immutable (static vector) types. Note that this setup is not automatically included with DifferentialEquaitons.jl. To use the following algorithms, you must install and use BridgeDiffEq.jl:Pkg.clone(\"https://github.com/JuliaDiffEq/BridgeDiffEq.jl\")\nusing BridgeDiffEqBridgeEuler - Strong order 0.5 Euler-Maruyama method for Ito equations.†\nBridgeHeun - Strong order 0.5 Euler-Heun method for Stratonovich equations.†\nBridgeSRK - Strong order 1.0 derivative-free stochastic Runge-Kutta method for scalar (<:Number) Ito equations.†"
},

{
    "location": "solvers/sde_solve.html#Notes-1",
    "page": "SDE Solvers",
    "title": "Notes",
    "category": "section",
    "text": "†: Does not step to the interval endpoint. This can cause issues with discontinuity detection, and discrete variables need to be updated appropriately.*:  Note that although SKenCarp uses the same table as KenCarp3, solving a ODE problem using SKenCarp by setting g(du,u,p,t) = du .= 0 will take much more steps than KenCarp3 because error estimator of SKenCarp is different (because of noise terms) and default value of qmax (maximum permissible ratio of relaxing/tightening dt for adaptive steps) is smaller for StochasticDiffEq algorithms."
},

{
    "location": "solvers/rode_solve.html#",
    "page": "RODE Solvers",
    "title": "RODE Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "solvers/rode_solve.html#RODE-Solvers-1",
    "page": "RODE Solvers",
    "title": "RODE Solvers",
    "category": "section",
    "text": ""
},

{
    "location": "solvers/rode_solve.html#Recommended-Methods-1",
    "page": "RODE Solvers",
    "title": "Recommended Methods",
    "category": "section",
    "text": "Currently, the only implemented method is the RandomEM method in StochasticDiffEq.jl. It is strong order alpha for a alpha-Holder continuous noise process."
},

{
    "location": "solvers/rode_solve.html#Full-List-of-Methods-1",
    "page": "RODE Solvers",
    "title": "Full List of Methods",
    "category": "section",
    "text": ""
},

{
    "location": "solvers/rode_solve.html#StochasticDiffEq.jl-1",
    "page": "RODE Solvers",
    "title": "StochasticDiffEq.jl",
    "category": "section",
    "text": "Each of the StochasticDiffEq.jl solvers come with a linear interpolation.RandomEM- The Euler-Maruyama method for RODEs. Strong order matching Holder continuity.Example usage:sol = solve(prob,RandomEM())"
},

{
    "location": "solvers/dde_solve.html#",
    "page": "DDE Solvers",
    "title": "DDE Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "solvers/dde_solve.html#DDE-Solvers-1",
    "page": "DDE Solvers",
    "title": "DDE Solvers",
    "category": "section",
    "text": "solve(prob::AbstractDDEProblem, alg; kwargs)Solves the DDE defined by prob using the algorithm alg. If no algorithm is given, a default algorithm will be chosen."
},

{
    "location": "solvers/dde_solve.html#Recommended-Methods-1",
    "page": "DDE Solvers",
    "title": "Recommended Methods",
    "category": "section",
    "text": "The recommended method for DDE problems are the MethodOfSteps algorithms. These are constructed from an OrdinaryDiffEq.jl algorithm as follows:MethodOfSteps(alg; constrained = false, fpsolve = FPFunctional(; max_iter = 10))where alg is an OrdinaryDiffEq.jl algorithm. Most algorithms should work."
},

{
    "location": "solvers/dde_solve.html#Nonstiff-DDEs-1",
    "page": "DDE Solvers",
    "title": "Nonstiff DDEs",
    "category": "section",
    "text": "The standard algorithm choice is MethodOfSteps(Tsit5()). This is a highly efficient FSAL 5th order algorithm with free interpolants which should handle most problems. For fast solving where non-strict error control is needed, choosing MethodOfSteps(BS3()) can do well. Using BS3 is similar to the MATLAB dde23. For algorithms where strict error control is needed, it is recommended that one uses MethodOfSteps(Vern6()). Benchmarks show that going to higher order methods like MethodOfSteps(DP8()) may not be beneficial."
},

{
    "location": "solvers/dde_solve.html#Stiff-DDEs-and-Differential-Algebraic-Delay-Equations-(DADEs)-1",
    "page": "DDE Solvers",
    "title": "Stiff DDEs and Differential-Algebraic Delay Equations (DADEs)",
    "category": "section",
    "text": "For stiff DDEs, the SDIRK and Rosenbrock methods are very efficient as they will reuse the Jacobian in the unconstrained stepping iterations. One should choose from the methods which have stiff-aware interpolants for better stability. MethodOfSteps(Rosenbrock23()) is a good low order method choice. Additionally, the Rodas methods like MethodOfSteps(Rodas4()) are good choices because of their higher order stiff-aware interpolant.Additionally, DADEs can be solved by specifying the problem in mass matrix form. The Rosenbrock methods are good choices in these situations."
},

{
    "location": "solvers/dde_solve.html#Lag-Handling-1",
    "page": "DDE Solvers",
    "title": "Lag Handling",
    "category": "section",
    "text": "Lags are declared separately from their use. One can use any lag by simply using the interpolant of h at that point. However, one should use caution in order to achieve the best accuracy. When lags are declared, the solvers can more efficiently be more accurate. Constant delays are propagated until the order is higher than the order of the integrator. If state-dependent delays are declared, the algorithm will detect discontinuities arising from these delays and adjust the step size such that these discontinuities are included in the mesh, if steps are rejected. This way, all discontinuities are treated exactly.If there are undeclared lags, the discontinuities due to delays are not tracked. In this case, one should only use residual control methods like MethodOfSteps(RK4()), which is the current best choice, as these will step more accurately. Still, residual control is an error-prone method. We recommend setting the tolerances lower in order to get accurate results, though this may be costly since it will use a rejection-based approach to adapt to the delay discontinuities."
},

{
    "location": "solvers/dde_solve.html#Special-Keyword-Arguments-1",
    "page": "DDE Solvers",
    "title": "Special Keyword Arguments",
    "category": "section",
    "text": "discontinuity_interp_points - Number of interpolation points used to track discontinuities arising from dependent delays. Defaults to 10. Only relevant if dependent delays are declared.\ndiscontinuity_abstol and discontinuity_reltol - These are absolute and relative tolerances used by the check whether the time point at the beginning of the current step is a discontinuity arising from dependent delays. Defaults to 1/10^12 and 0. Only relevant if dependent delays are declared."
},

{
    "location": "solvers/dde_solve.html#Note-1",
    "page": "DDE Solvers",
    "title": "Note",
    "category": "section",
    "text": "If the method is having trouble, one may want to adjust the fixed-point iteration. Decreasing the absolute tolerance and the relative tolerance by specifying the keyword arguments abstol and reltol when solving the DDE problem, and increasing the maximal number of iterations by specifying the keyword argument max_iter in the MethodOfSteps algorithm, can help ensure that the steps are correct. If the problem still is not correctly converging, one should lower dtmax. For problems with only constant delays, in the worst case scenario, one may need to set constrained = true which will constrain timesteps to at most the size of the minimal lag and hence forces more stability at the cost of smaller timesteps."
},

{
    "location": "solvers/dae_solve.html#",
    "page": "DAE Solvers",
    "title": "DAE Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "solvers/dae_solve.html#DAE-Solvers-1",
    "page": "DAE Solvers",
    "title": "DAE Solvers",
    "category": "section",
    "text": ""
},

{
    "location": "solvers/dae_solve.html#Recomended-Methods-1",
    "page": "DAE Solvers",
    "title": "Recomended Methods",
    "category": "section",
    "text": "For medium to low accuracy DAEs in mass matrix form, the Rodas4 and Rodas42 methods are good choices which will get good efficiency. The OrdinaryDiffEq.jl methods are also the only methods which allow for Julia-defined number types. For high accuracy (error <1e-7) on problems of Vector{Float64} defined in mass matrix form, radau is an efficient method.If the problem cannot be defined in mass matrix form, the recommended method for performance is IDA from the Sundials.jl package if you are solving problems with Float64. It\'s a very well-optimized method, and allows you to have a little bit of control over the linear solver to better tailor it to your problem. A similar algorithm is daskr. Which one is more efficient is problem-dependent."
},

{
    "location": "solvers/dae_solve.html#Full-List-of-Methods-1",
    "page": "DAE Solvers",
    "title": "Full List of Methods",
    "category": "section",
    "text": ""
},

{
    "location": "solvers/dae_solve.html#OrdinaryDiffEq.jl-1",
    "page": "DAE Solvers",
    "title": "OrdinaryDiffEq.jl",
    "category": "section",
    "text": "These methods require the DAE to be an ODEProblem in mass matrix form. For extra options for the solvers, see the ODE solver page."
},

{
    "location": "solvers/dae_solve.html#Rosenbrock-Methods-1",
    "page": "DAE Solvers",
    "title": "Rosenbrock Methods",
    "category": "section",
    "text": "ROS3P - 3rd order A-stable and stiffly stable Rosenbrock method. Keeps high accuracy on discretizations of nonlinear parabolic PDEs.\nRodas3 - 3rd order A-stable and stiffly stable Rosenbrock method.\nRosShamp4- An A-stable 4th order Rosenbrock method.\nVeldd4 - A 4th order D-stable Rosenbrock method.\nVelds4 - A 4th order A-stable Rosenbrock method.\nGRK4T - An efficient 4th order Rosenbrock method.\nGRK4A - An A-stable 4th order Rosenbrock method. Essentially \"anti-L-stable\" but efficient.\nRos4LStab - A 4th order L-stable Rosenbrock method.\nRodas4 - A 4th order A-stable stiffly stable Rosenbrock method with a stiff-aware 3rd order interpolant\nRodas42 - A 4th order A-stable stiffly stable Rosenbrock method with a stiff-aware 3rd order interpolant\nRodas4P - A 4th order A-stable stiffly stable Rosenbrock method with a stiff-aware 3rd order interpolant. 4th order on linear parabolic problems and 3rd order accurate on nonlinear parabolic problems (as opposed to lower if not corrected).\nRodas5 - A 5th order A-stable stiffly stable Rosenbrock method. Currently has a Hermite interpolant because its stiff-aware 3rd order interpolant is not yet implemented."
},

{
    "location": "solvers/dae_solve.html#Rosenbrock-W-Methods-1",
    "page": "DAE Solvers",
    "title": "Rosenbrock-W Methods",
    "category": "section",
    "text": "Rosenbrock23 - An Order 2/3 L-Stable Rosenbrock-W method which is good for very stiff equations with oscillations at low tolerances. 2nd order stiff-aware interpolation.\nRosenbrock32 - An Order 3/2 A-Stable Rosenbrock-W method which is good for mildy stiff equations without oscillations at low tolerances. Note that this method is prone to instability in the presence of oscillations, so use with caution. 2nd order stiff-aware interpolation.\nRosenbrockW6S4OS - A 4th order L-stable Rosenbrock-W method (fixed step only).\nROS34PW1a - A 4th order L-stable Rosenbrock-W method.\nROS34PW1b - A 4th order L-stable Rosenbrock-W method.\nROS34PW2 - A 4th order stiffy accurate Rosenbrock-W method for PDAEs.\nROS34PW3 - A 4th order strongly A-stable (Rinf~0.63) Rosenbrock-W method."
},

{
    "location": "solvers/dae_solve.html#SDIRK-Methods-1",
    "page": "DAE Solvers",
    "title": "SDIRK Methods",
    "category": "section",
    "text": "SDIRK MethodsImplicitEuler - Stage order 1. A-B-L-stable. Adaptive timestepping through a divided differences estimate via memory. Strong-stability presurving (SSP).\nImplicitMidpoint - Stage order 1. Symplectic. Good for when symplectic integration is required."
},

{
    "location": "solvers/dae_solve.html#Sundials.jl-1",
    "page": "DAE Solvers",
    "title": "Sundials.jl",
    "category": "section",
    "text": "Note that this setup is not automatically included with DifferentialEquations.jl. To use the following algorithms, you must install and use Sundials.jl:]add Sundials\nusing SundialsIDA - This is the IDA method from the Sundials.jl package.Note that the constructors for the Sundials algorithms take a main argument:linearsolver - This is the linear solver which is used in the Newton iterations. The choices are:\n:Dense - A dense linear solver.\n:Band - A solver specialized for banded Jacobians. If used, you must set the position of the upper and lower non-zero diagonals via jac_upper and jac_lower.\n:GMRES - A GMRES method. Recommended first choice Krylov method\n:BCG - A Biconjugate gradient method.\n:PCG - A preconditioned conjugate gradient method. Only for symmetric linear systems.\n:TFQMR - A TFQMR method.\n:KLU - A sparse factorization method. Requires that the user specifies a Jacobian. The Jacobian must be set as a sparse matrix in the ODEProblem type.Example:IDA() # Newton + Dense solver\nIDA(linear_solver=:Band,jac_upper=3,jac_lower=3) # Banded solver with nonzero diagonals 3 up and 3 down\nIDA(linear_solver=:BCG) # Biconjugate gradient method                                   All of the additional options are available. The constructor is:IDA(;linear_solver=:Dense,jac_upper=0,jac_lower=0,krylov_dim=0,\n    max_order = 5,\n    max_error_test_failures = 7,\n    max_nonlinear_iters = 3,\n    nonlinear_convergence_coefficient = 0.33,\n    nonlinear_convergence_coefficient_ic = 0.0033,\n    max_num_steps_ic = 5,\n    max_num_jacs_ic = 4,\n    max_num_iters_ic = 10,\n    max_num_backs_ic = 100,\n    use_linesearch_ic = true,\n    max_convergence_failures = 10,\n    init_all = false)See the Sundials manual for details on the additional options. The option init_all controls the initial condition consistancy routine. If the initial conditions are inconsistant (i.e. they do not satisfy the implicit equation), init_all=false means that the algebraic variables and derivatives will be modified in order to satisfy the DAE. If init_all=true, all initial conditions will be modified to satify the DAE."
},

{
    "location": "solvers/dae_solve.html#DASKR.jl-1",
    "page": "DAE Solvers",
    "title": "DASKR.jl",
    "category": "section",
    "text": "DASKR.jl is not automatically included by DifferentialEquations.jl. To use this algorithm, you will need to install and use the package:]add DASKR\nusing DASKRdaskr - This is a wrapper for the well-known DASKR algorithm.All additional options are available. The constructor is:function daskr(;linear_solver=:Dense,\n                  jac_upper=0,jac_lower=0,max_order = 5,\n                  non_negativity_enforcement = 0,\n                  non_negativity_enforcement_array = nothing,\n                  max_krylov_iters = nothing,\n                  num_krylov_vectors = nothing,\n                  max_number_krylov_restarts = 5,\n                  krylov_convergence_test_constant = 0.05,\n                  exclude_algebraic_errors = false)Choices for the linear solver are::Dense\n:Banded\n:SPIGMR, a Krylov method"
},

{
    "location": "solvers/dae_solve.html#DASSL.jl-1",
    "page": "DAE Solvers",
    "title": "DASSL.jl",
    "category": "section",
    "text": "dassl - A native Julia implementation of the DASSL algorithm."
},

{
    "location": "solvers/dae_solve.html#ODEInterfaceDiffEq.jl-1",
    "page": "DAE Solvers",
    "title": "ODEInterfaceDiffEq.jl",
    "category": "section",
    "text": "These methods require the DAE to be an ODEProblem in mass matrix form. For extra options for the solvers, see the ODE solver page.seulex - Extrapolation-algorithm based on the linear implicit Euler method.\nradau - Implicit Runge-Kutta (Radau IIA) of variable order between 5 and 13.\nradau5 - Implicit Runge-Kutta method (Radau IIA) of order 5.\nrodas - Rosenbrock 4(3) method."
},

{
    "location": "solvers/benchmarks.html#",
    "page": "Solver Benchmarks",
    "title": "Solver Benchmarks",
    "category": "page",
    "text": ""
},

{
    "location": "solvers/benchmarks.html#Solver-Benchmarks-1",
    "page": "Solver Benchmarks",
    "title": "Solver Benchmarks",
    "category": "section",
    "text": "Benchmarks for the solvers can be found at DiffEqBenchmarks.jl. Many different problems are tested. However, if you would like additional problems to be benchmarked, please open an issue or PR at the DiffEqBenchmarks.jl repository with the code that defines the DEProblem."
},

{
    "location": "features/performance_overloads.html#",
    "page": "DiffEqFunctions (Jacobians, Gradients, etc.) and Jacobian Types",
    "title": "DiffEqFunctions (Jacobians, Gradients, etc.) and Jacobian Types",
    "category": "page",
    "text": ""
},

{
    "location": "features/performance_overloads.html#DiffEqFunctions-(Jacobians,-Gradients,-etc.)-and-Jacobian-Types-1",
    "page": "DiffEqFunctions (Jacobians, Gradients, etc.) and Jacobian Types",
    "title": "DiffEqFunctions (Jacobians, Gradients, etc.) and Jacobian Types",
    "category": "section",
    "text": "The DiffEq ecosystem provides an extensive interface for declaring extra functions associated with the differential equation\'s data. In traditional libraries there is usually only one option: the Jacobian. However, we allow for a large array of pre-computed functions to speed up the calculations. This is offered via the DiffEqFunction types which can be passed to the problems."
},

{
    "location": "features/performance_overloads.html#Function-Type-Definitions-1",
    "page": "DiffEqFunctions (Jacobians, Gradients, etc.) and Jacobian Types",
    "title": "Function Type Definitions",
    "category": "section",
    "text": ""
},

{
    "location": "features/performance_overloads.html#Function-Choice-Definitions-1",
    "page": "DiffEqFunctions (Jacobians, Gradients, etc.) and Jacobian Types",
    "title": "Function Choice Definitions",
    "category": "section",
    "text": "The full interface available to the solvers is as follows:jac: The Jacobian of the differential equation with respect to the state variable u at a time t with parameters p.\ntgrad: The gradient of the differential equation with respect to t at state u with parameters p.\nparamjac: The Jacobian of the differential equation with respect to p at state u at time t.\nanalytic: Defines an analytical solution using u0 at time t with p which will cause the solvers to return errors. Used for testing.\ninvW: The inverse of M - gamma*J where J is the jac.\ninvW_t: The inverse of M/gamma - J where J is the jac.\nggprime: See the definition in the SDEProblem page."
},

{
    "location": "features/performance_overloads.html#ODEFunction-1",
    "page": "DiffEqFunctions (Jacobians, Gradients, etc.) and Jacobian Types",
    "title": "ODEFunction",
    "category": "section",
    "text": "function ODEFunction{iip,recompile}(f;\n                 mass_matrix=I,\n                 analytic=nothing, # (u0,p,t)\n                 tgrad=nothing, # (dT,u,p,t) or (u,p,t)\n                 jac=nothing, # (J,u,p,t) or (u,p,t)\n                 jac_prototype=nothing, # Type for the Jacobian\n                 invW=nothing, # (iW,u,p,t) or (u,p,t)\n                 invW_t=nothing, # (iW,u,p,t) or (u,p,t)\n                 paramjac = nothing, # (pJ,u,p,t) or (u,p,t)\n                 syms = nothing) # collection of names for variables"
},

{
    "location": "features/performance_overloads.html#DynamicalODEFunction-1",
    "page": "DiffEqFunctions (Jacobians, Gradients, etc.) and Jacobian Types",
    "title": "DynamicalODEFunction",
    "category": "section",
    "text": "DynamicalODEFunction{iip,recompile}(f1, # (du,u,v,p,t) or (u,v,p,t)\n                                    f2; # (du,u,v,p,t) or (u,v,p,t)\n                                    mass_matrix=(I,I), # Mass matrix for each partition\n                                    analytic=nothing)"
},

{
    "location": "features/performance_overloads.html#SplitFunction-1",
    "page": "DiffEqFunctions (Jacobians, Gradients, etc.) and Jacobian Types",
    "title": "SplitFunction",
    "category": "section",
    "text": "SplitFunction{iip,recompile}(f1, # ODEFunction\n                        f2; # ODEFunction\n                        mass_matrix=I,\n                        _func_cache=nothing, # This is a cache used in f = f1+f2\n                        analytic=nothing)"
},

{
    "location": "features/performance_overloads.html#SDEFunction-1",
    "page": "DiffEqFunctions (Jacobians, Gradients, etc.) and Jacobian Types",
    "title": "SDEFunction",
    "category": "section",
    "text": "function SDEFunction{iip,recompile}(f,g;\n                 mass_matrix=I,\n                 analytic=nothing,\n                 tgrad=nothing,\n                 jac=nothing,\n                 jac_prototype=nothing,\n                 invW=nothing,\n                 invW_t=nothing,\n                 paramjac = nothing,\n                 ggprime = nothing,\n                 syms = nothing)"
},

{
    "location": "features/performance_overloads.html#SplitSDEFunction-1",
    "page": "DiffEqFunctions (Jacobians, Gradients, etc.) and Jacobian Types",
    "title": "SplitSDEFunction",
    "category": "section",
    "text": "SplitSDEFunction{iip,recompile}(f1, # ODEFunction\n                           f2, # ODEFunction\n                           g;\n                           mass_matrix=I,\n                           _func_cache=nothing,\n                           analytic=nothing)"
},

{
    "location": "features/performance_overloads.html#RODEFunction-1",
    "page": "DiffEqFunctions (Jacobians, Gradients, etc.) and Jacobian Types",
    "title": "RODEFunction",
    "category": "section",
    "text": "function RODEFunction{iip,recompile}(f;\n                 mass_matrix=I,\n                 analytic=nothing,\n                 tgrad=nothing,\n                 jac=nothing,\n                 jac_prototype=nothing,\n                 invW=nothing,\n                 invW_t=nothing,\n                 paramjac = nothing,\n                 syms = nothing)"
},

{
    "location": "features/performance_overloads.html#DAEFunction-1",
    "page": "DiffEqFunctions (Jacobians, Gradients, etc.) and Jacobian Types",
    "title": "DAEFunction",
    "category": "section",
    "text": "function DAEFunction{iip,recompile}(f;\n                 mass_matrix=I,\n                 analytic=nothing,\n                 tgrad=nothing,\n                 jac=nothing, # (J,du,u,p,gamma,t) or (du,u,p,gamma,t)\n                 jac_prototype=nothing,\n                 invW=nothing,\n                 invW_t=nothing,\n                 paramjac = nothing,\n                 syms = nothing)Note that the Jacobian of a DAE is defined as gamma*dG/d(du) + dG/du where gamma is given by the solver."
},

{
    "location": "features/performance_overloads.html#DDEFunction-1",
    "page": "DiffEqFunctions (Jacobians, Gradients, etc.) and Jacobian Types",
    "title": "DDEFunction",
    "category": "section",
    "text": "function DDEFunction{iip,recompile}(f;\n                 mass_matrix=I,\n                 analytic=nothing,\n                 tgrad=nothing,\n                 jac=nothing,\n                 jac_prototype=nothing,\n                 invW=nothing,\n                 invW_t=nothing,\n                 paramjac = nothing,\n                 syms = nothing)"
},

{
    "location": "features/performance_overloads.html#Inplace-Specification-and-No-Recompile-Mode-1",
    "page": "DiffEqFunctions (Jacobians, Gradients, etc.) and Jacobian Types",
    "title": "Inplace Specification and No-Recompile Mode",
    "category": "section",
    "text": "Each DiffEqFunction type can be called with an \"is inplace\" (iip) choice.ODEFunction(f)\nODEFunction{iip}(f)which is a boolean for whether the function is in the inplace form (mutating to change the first value). This is automatically determined using the methods table but note that for full type-inferrability of the DEProblem this iip-ness should be specified.Additionally, the functions are fully specialized to reduce the runtimes. If one would instead like to not specialize on the functions to reduce compile time, then one can set recompile to false.ODEFunction{iip,false}(f)This makes the ODE solver compilation independent of the function and so changing the function will not cause recompilation. One can change the default value by changing the const RECOMPILE_BY_DEFAULT = true to false in the DiffEqBase.jl source code."
},

{
    "location": "features/performance_overloads.html#Specifying-Jacobian-Types-1",
    "page": "DiffEqFunctions (Jacobians, Gradients, etc.) and Jacobian Types",
    "title": "Specifying Jacobian Types",
    "category": "section",
    "text": "The jac field of an inplace style DiffEqFunction has the signature jac(J,u,p,t), which updates the jacobian J in-place. The intended type for J can sometimes be inferred (e.g. when it is just a dense Matrix), but not in general. To supply the type information, you can provide a jac_prototype in the function\'s constructor.The following example creates an inplace ODEFunction whose jacobian is a Diagonal:using LinearAlgebra\nf = (du,u,p,t) -> du .= t .* u\njac = (J,u,p,t) -> (J[1,1] = t; J[2,2] = t; J)\njp = Diagonal(zeros(2))\nfun = ODEFunction(f; jac=jac, jac_prototype=jp)Note that the integrators will always make a deep copy of fun.jac_prototype, so there\'s no worry of aliasing.In general the jacobian prototype can be anything that has mul! defined, in particular sparse matrices or custom lazy types that support mul!. A special case is when the jac_prototype is a AbstractDiffEqLinearOperator, in which case you do not need to supply jac as it is automatically set to update_coefficients!. Refer to the DiffEqOperators section for more information on setting up time/parameter dependent operators."
},

{
    "location": "features/performance_overloads.html#Examples-1",
    "page": "DiffEqFunctions (Jacobians, Gradients, etc.) and Jacobian Types",
    "title": "Examples",
    "category": "section",
    "text": ""
},

{
    "location": "features/performance_overloads.html#Declaring-Explicit-Jacobians-for-ODEs-1",
    "page": "DiffEqFunctions (Jacobians, Gradients, etc.) and Jacobian Types",
    "title": "Declaring Explicit Jacobians for ODEs",
    "category": "section",
    "text": "The most standard case, declaring a function for a Jacobian is done by overloading the function f(du,u,p,t) with an in-place updating function for the Jacobian: f_jac(J,u,p,t) where the value type is used for dispatch. For example, take the LotkaVolterra model:function f(du,u,p,t)\n  du[1] = 2.0 * u[1] - 1.2 * u[1]*u[2]\n  du[2] = -3 * u[2] + u[1]*u[2]\nendTo declare the Jacobian we simply add the dispatch:function f_jac(J,u,p,t)\n  J[1,1] = 2.0 - 1.2 * u[2]\n  J[1,2] = -1.2 * u[1]\n  J[2,1] = 1 * u[2]\n  J[2,2] = -3 + u[1]\n  nothing\nendThen we can supply the Jacobian with our ODE as:ff = ODEFunction(f;jac=f_jac)and use this in an ODEProblem:prob = ODEProblem(ff,ones(2),(0.0,10.0))"
},

{
    "location": "features/performance_overloads.html#Declaring-Explicit-Jacobians-for-DAEs-1",
    "page": "DiffEqFunctions (Jacobians, Gradients, etc.) and Jacobian Types",
    "title": "Declaring Explicit Jacobians for DAEs",
    "category": "section",
    "text": "For fully implicit ODEs (DAEProblems), a slightly different Jacobian function is necessary. For the DAEG(duupt) = resThe Jacobian should be given in the form gamma*dG/d(du) + dG/du where gamma is given by the solver. This means that the signature is:f(::Type{Val{:jac}},J,du,u,p,gamma,t)For example, for the equationfunction testjac(res,du,u,p,t)\n  res[1] = du[1] - 2.0 * u[1] + 1.2 * u[1]*u[2]\n  res[2] = du[2] -3 * u[2] - u[1]*u[2]\nendwe would define the Jacobian as:function testjac(::Type{Val{:jac}},J,du,u,p,gamma,t)\n  J[1,1] = gamma - 2.0 + 1.2 * u[2]\n  J[1,2] = 1.2 * u[1]\n  J[2,1] = - 1 * u[2]\n  J[2,2] = gamma - 3 - u[1]\n  nothing\nend"
},

{
    "location": "features/performance_overloads.html#Symbolically-Calculating-the-Functions-1",
    "page": "DiffEqFunctions (Jacobians, Gradients, etc.) and Jacobian Types",
    "title": "Symbolically Calculating the Functions",
    "category": "section",
    "text": "ParameterizedFunctions.jl automatically calculates as many of these functions as possible and generates the ODEFunction using SymEngine. Thus, for good performance with the least work, it is one can try ParameterizedFunctions.jl.Additionally, an up-and-coming effort in the JuliaDiffEq ecosystem is ModelingToolkit.jl for performing these calculations more generically."
},

{
    "location": "features/diffeq_arrays.html#",
    "page": "DiffEq-Specific Array Types",
    "title": "DiffEq-Specific Array Types",
    "category": "page",
    "text": ""
},

{
    "location": "features/diffeq_arrays.html#DiffEq-Specific-Array-Types-1",
    "page": "DiffEq-Specific Array Types",
    "title": "DiffEq-Specific Array Types",
    "category": "section",
    "text": "In many cases, a standard array may not be enough to fully hold the data for a model. Many of the solvers in DifferentialEquations.jl (only the native Julia methods) allow you to solve problems on AbstractArray types which allow you to extend the meaning of an array. This page describes some of the AbstractArray types which can be helpful for modeling differential equations problems."
},

{
    "location": "features/diffeq_arrays.html#ArrayPartitions-1",
    "page": "DiffEq-Specific Array Types",
    "title": "ArrayPartitions",
    "category": "section",
    "text": "ArrayPartitions in DiffEq are used for heterogeneous arrays. For example, PartitionedODEProblem solvers use them internally to turn the separate parts into a single array. You can construct an ArrayPartition using RecursiveArrayTools.jl:using RecursiveArrayTools\nA = ArrayPartition(x::AbstractArray...)where is x a list of arrays. The resulting A will act like a single array, and its broadcast will be type stable, allowing for it to be used inside of the native Julia DiffEq solvers in an efficient way. This is a good way to generate an array which has different units for different parts, or different amounts of precision."
},

{
    "location": "features/diffeq_arrays.html#Usage-1",
    "page": "DiffEq-Specific Array Types",
    "title": "Usage",
    "category": "section",
    "text": "An ArrayPartition acts like a single array. A[i] indexes through the first array, then the second, etc. all linearly. But A.x is where the arrays are stored. Thus forusing RecursiveArrayTools\nA = ArrayPartition(y,z)We would have A.x[1]==y and A.x[2]==z. Broadcasting like f.(A) is efficient."
},

{
    "location": "features/diffeq_arrays.html#Example:-Dynamics-Equations-1",
    "page": "DiffEq-Specific Array Types",
    "title": "Example: Dynamics Equations",
    "category": "section",
    "text": "In this example we will show using heterogeneous units in dynamics equations. Our arrays will be:using Unitful, RecursiveArrayTools, DiffEqBase, OrdinaryDiffEq\nusing LinearAlgebra\n\nr0 = [1131.340, -2282.343, 6672.423]u\"km\"\nv0 = [-5.64305, 4.30333, 2.42879]u\"km/s\"\nΔt = 86400.0*365u\"s\"\nμ = 398600.4418u\"km^3/s^2\"\nrv0 = ArrayPartition(r0,v0)Here, r0 is the initial positions, and v0 are the initial velocities. rv0 is the ArrayPartition initial condition. We now write our update function in terms of the ArrayPartition:function f(dy, y, μ, t)\n    r = norm(y.x[1])\n    dy.x[1] .= y.x[2]\n    dy.x[2] .= -μ .* y.x[1] / r^3\nendNotice that y.x[1] is the r part of y, and y.x[2] is the v part of y. Using this kind of indexing is type stable, even though the array itself is heterogeneous. Note that one can also use things like 2y or y.+x and the broadcasting will be efficient.Now to solve our equations, we do the same thing as always in DiffEq:prob = ODEProblem(f, rv0, (0.0u\"s\", Δt), μ)\nsol = solve(prob, Vern8())"
},

{
    "location": "features/diffeq_arrays.html#MultiScaleArrays-1",
    "page": "DiffEq-Specific Array Types",
    "title": "MultiScaleArrays",
    "category": "section",
    "text": "The multi-scale modeling functionality is provided by MultiScaleArrays.jl. It allows for designing a multi-scale model as an extension of an array, which in turn can be directly used in the native Julia solvers of DifferentialEquations.jl.For more information, please see the MultiScaleArrays.jl README."
},

{
    "location": "features/diffeq_arrays.html#DEDataArrays-1",
    "page": "DiffEq-Specific Array Types",
    "title": "DEDataArrays",
    "category": "section",
    "text": "The DEDataArray{T} type allows one to add other \"non-continuous\" variables to an array, which can be useful in many modeling situations involving lots of events. To define an DEDataArray, make a type which subtypes DEDataArray{T} with a field x for the \"array of continuous variables\" for which you would like the differential equation to treat directly. The other fields are treated as \"discrete variables\". For example:mutable struct MyDataArray{T,1} <: DEDataArray{T,1}\n    x::Array{T,1}\n    a::T\n    b::Symbol\nendIn this example, our resultant array is a SimType, and its data which is presented to the differential equation solver will be the array x. Any array which the differential equation solver can use is allowed to be made as the field x, including other DEDataArrays. Other than that, you can add whatever fields you please, and let them be whatever type you please.These extra fields are carried along in the differential equation solver that the user can use in their f equation and modify via callbacks. For example, inside of a an update function, it is safe to do:function f(du,u,p,t)\n  u.a = t\nendto update the discrete variables (unless the algorithm notes that it does not step to the endpoint, in which case a callback must be used to update appropriately.)Note that the aliases DEDataVector and DEDataMatrix cover the one and two dimensional cases."
},

{
    "location": "features/diffeq_arrays.html#Example:-A-Control-Problem-1",
    "page": "DiffEq-Specific Array Types",
    "title": "Example: A Control Problem",
    "category": "section",
    "text": "In this example we will use a DEDataArray to solve a problem where control parameters change at various timepoints. First we will buildmutable struct SimType{T} <: DEDataVector{T}\n    x::Array{T,1}\n    f1::T\nendas our DEDataVector. It has an extra field f1 which we will use as our control variable. Our ODE function will use this field as follows:function f(du,u,p,t)\n    du[1] = -0.5*u[1] + u.f1\n    du[2] = -0.5*u[2]\nendNow we will setup our control mechanism. It will be a simple setup which uses set timepoints at which we will change f1. At t=5.0 we will want to increase the value of f1, and at t=8.0 we will want to decrease the value of f1. Using the DiscreteCallback interface, we code these conditions as follows:const tstop1 = [5.]\nconst tstop2 = [8.]\n\n\nfunction condition(u,t,integrator)\n  t in tstop1\nend\n\nfunction condition2(u,t,integrator)\n  t in tstop2\nendNow we have to apply an effect when these conditions are reached. When condition is hit (at t=5.0), we will increase f1 to 1.5. When condition2 is reached, we will decrease f1 to -1.5. This is done via the functions:function affect!(integrator)\n  for c in full_cache(integrator)\n    c.f1 = 1.5\n  end\nend\n\nfunction affect2!(integrator)\n  for c in full_cache(integrator)\n    c.f1 = -1.5\n  end\nendNotice that we have to loop through the full_cache array (provided by the integrator interface) to ensure that all internal caches are also updated. With these functions we can build our callbacks:save_positions = (true,true)\n\ncb = DiscreteCallback(condition, affect!, save_positions=save_positions)\n\nsave_positions = (false,true)\n\ncb2 = DiscreteCallback(condition2, affect2!, save_positions=save_positions)\n\ncbs = CallbackSet(cb,cb2)Now we define our initial condition. We will start at [10.0;10.0] with f1=0.0.u0 = SimType([10.0;10.0], 0.0)\nprob = ODEProblem(f,u0,(0.0,10.0))Lastly we solve the problem. Note that we must pass tstop values of 5.0 and 8.0 to ensure the solver hits those timepoints exactly:const tstop = [5.;8.]\nsol = solve(prob,Tsit5(),callback = cbs, tstops=tstop)(Image: data_array_plot)It\'s clear from the plot how the controls affected the outcome."
},

{
    "location": "features/diffeq_arrays.html#Data-Arrays-vs-ParameterizedFunctions-1",
    "page": "DiffEq-Specific Array Types",
    "title": "Data Arrays vs ParameterizedFunctions",
    "category": "section",
    "text": "The reason for using a DEDataArray is because the solution will then save the control parameters. For example, we can see what the control parameter was at every timepoint by checking:[sol[i].f1 for i in 1:length(sol)]A similar solution can be achieved using a ParameterizedFunction. We could have instead created our function as:function f(du,u,p,t)\n    du[1] = -0.5*u[1] + p\n    du[2] = -0.5*u[2]\nend\nu0 = SimType([10.0;10.0], 0.0)\np = 0.0\nprob = ODEProblem(f,u0,(0.0,10.0),p)\nconst tstop = [5.;8.]\nsol = solve(prob,Tsit5(),callback = cbs, tstops=tstop)where we now change the callbacks to changing the parameter:function affect!(integrator)\n  integrator.p = 1.5\nend\n\nfunction affect2!(integrator)\n  integrator.p = -1.5\nendThis will also solve the equation and get a similar result. It will also be slightly faster in some cases. However, if the equation is solved in this manner, there will be no record of what the parameter was at each timepoint. That is the tradeoff between DEDataArrays and ParameterizedFunctions."
},

{
    "location": "features/diffeq_operator.html#",
    "page": "DiffEqOperators",
    "title": "DiffEqOperators",
    "category": "page",
    "text": ""
},

{
    "location": "features/diffeq_operator.html#DiffEqOperators-1",
    "page": "DiffEqOperators",
    "title": "DiffEqOperators",
    "category": "section",
    "text": "The AbstractDiffEqOperator interface is an interface for declaring parts of a differential equation as linear or affine. This then allows the solvers to exploit linearity to achieve maximal performance."
},

{
    "location": "features/diffeq_operator.html#Note:-this-functionality-is-under-heavy-development.-Interfaces-may-change.-Treat-this-as-experimental-until-the-end-of-Summer-2019.-1",
    "page": "DiffEqOperators",
    "title": "Note: this functionality is under heavy development. Interfaces may change. Treat this as experimental until the end of Summer 2019.",
    "category": "section",
    "text": ""
},

{
    "location": "features/diffeq_operator.html#Using-DiffEqOperators-1",
    "page": "DiffEqOperators",
    "title": "Using DiffEqOperators",
    "category": "section",
    "text": "AbstractDiffEqOperators act like functions. When defined, A has function calls A(u,p,t) and A(du,u,p,t) that act like A*u. These operators update via a function update_coefficients!(A,u,p,t)."
},

{
    "location": "features/diffeq_operator.html#Constructors-1",
    "page": "DiffEqOperators",
    "title": "Constructors",
    "category": "section",
    "text": ""
},

{
    "location": "features/diffeq_operator.html#Wrapping-an-Array:-DiffEqArrayOperator-1",
    "page": "DiffEqOperators",
    "title": "Wrapping an Array: DiffEqArrayOperator",
    "category": "section",
    "text": "DiffEqArrayOperator is for defining an operator directly from an array. The operator is of the formalpha(t)A(upt)for some scalar α and time plus possibly state dependent A. The constructor is:DiffEqArrayOperator(A::AbstractMatrix{T},α=1.0,\n                             update_func = DEFAULT_UPDATE_FUNC)A is the operator array. α is the scalar coefficient. If α is a function α(t), then it will update the coefficient as necessary. update_func(A,u,p,t)  is the function called by update_coefficients!(A,u,p,t) (along with α if it\'s  a function). If left as its default, then update_func is trivial which signifies A is a constant."
},

{
    "location": "features/diffeq_operator.html#AffineDiffEqOperator-1",
    "page": "DiffEqOperators",
    "title": "AffineDiffEqOperator",
    "category": "section",
    "text": "For As = (A1,A2,...,An) and Bs = (B1,B2,...,Bn) where each of the Ai and Bi are DiffEqLinearOperators, the following constructor:function AffineDiffEqOperator{T}(As,Bs,u_cache=nothing)builds an operator L = (A1 + A2 + ... An)*u + B1 + B2 + ... + Bn. u_cache is for designating a type of internal cache for non-allocating evaluation of L(du,u,p,t). If not given, the function L(du,u,p,t) is not available. Note that in solves which exploit this structure, this function call is not necessary. It\'s only used as the fallback in ODE solvers which were not developed for this structure."
},

{
    "location": "features/diffeq_operator.html#Formal-Properties-of-DiffEqOperators-1",
    "page": "DiffEqOperators",
    "title": "Formal Properties of DiffEqOperators",
    "category": "section",
    "text": "These are the formal properties that an AbstractDiffEqOperator should obey for it to work in the solvers."
},

{
    "location": "features/diffeq_operator.html#AbstractDiffEqOperator-Interface-Description-1",
    "page": "DiffEqOperators",
    "title": "AbstractDiffEqOperator Interface Description",
    "category": "section",
    "text": "Function call and multiplication: L(du,u,p,t) for inplace and du = L(u,p,t) for out-of-place, meaning L*u and A_mul_B!.\nIf the operator is not a constant, update it with (u,p,t). A mutating form, i.e. update_coefficients!(A,u,p,t) that changes the internal coefficients, and a out-of-place form B = update_coefficients(A,u,p,t).\nis_constant(A) trait for whether the operator is constant or not."
},

{
    "location": "features/diffeq_operator.html#AbstractDiffEqLinearOpeartor-Interface-Description-1",
    "page": "DiffEqOperators",
    "title": "AbstractDiffEqLinearOpeartor Interface Description",
    "category": "section",
    "text": "AbstractDiffEqLinearOperator <: AbstractDiffEqOperator\nCan absorb under multiplication by a scalar. In all algorithms things like dt*L show up all the time, so the linear operator must be able to absorb such constants.\nis_constant(A) trait for whether the operator is constant or not.\nOptional: diagonal, symmetric, etc traits from LinearMaps.jl.\nOptional: expm(A). Required for simple exponential integration.\nOptional: expmv(A,u,t) = expm(t*A)*u and expmv!(v,A::DiffEqOperator,u,t) Required for sparse-saving exponential integration.\nOptional: factorizations. A_ldiv_B, factorize et. al. This is only required for algorithms which use the factorization of the operator (Crank-Nicholson), and only for when the default linear solve is used."
},

{
    "location": "features/noise_process.html#",
    "page": "Noise Processes",
    "title": "Noise Processes",
    "category": "page",
    "text": ""
},

{
    "location": "features/noise_process.html#Noise-Processes-1",
    "page": "Noise Processes",
    "title": "Noise Processes",
    "category": "section",
    "text": "Noise processes are essential in continuous stochastic modeling. The NoiseProcess types are distributionally-exact, meaning they are not solutions of stochastic differential equations and instead are directly generated according to their analytical distributions. These processes are used as the noise term in the SDE and RODE solvers. Additionally, the noise processes themselves can be simulated and solved using the DiffEq common interface (including the Monte Carlo interface).This page first describes how to use noise processes in SDEs, and analyze/simulate them directly noise processes. Then it describes the standard noise processes which are available. Processes like WienerProcess, CorrelatedWienerProcess, GeometricBrownianMotionProcess, BrownianBridgeProcess and OrnsteinUhlenbeckProcess are pre-defined. Then it is shown how one can define the distributions for a new NoiseProcess.In addition to the NoiseProcess type, more general AbstractNoiseProcesses are defined. The NoiseGrid allows you to define a noise process from a set of pre-calculated points (the \"normal\" way). The NoiseApproximation allows you to define a new noise process as the solution to some stochastic differential equation. While these methods are only approximate, they are more general and allow the user to easily define their own colored noise to use in simulations.The NoiseWrapper allows one to wrap a NoiseProcess from a previous simulation to re-use it in a new simulation in a way that follows the same stochastic trajectory (even if different points are hit, for example solving with a smaller dt) in a distributionally-exact manner. It is demonstrated how the NoiseWrapper can be used to wrap the NoiseProcess of one SDE/RODE solution in order to re-use the same noise process in another simulation.Lastly, the NoiseFunction allows you to use any function of time as the noise process. Together, this functionality allows you to define any colored noise process and use this efficiently and accurately in your simulations."
},

{
    "location": "features/noise_process.html#Using-Noise-Processes-1",
    "page": "Noise Processes",
    "title": "Using Noise Processes",
    "category": "section",
    "text": ""
},

{
    "location": "features/noise_process.html#Passing-a-Noise-Process-to-a-Problem-Type-1",
    "page": "Noise Processes",
    "title": "Passing a Noise Process to a Problem Type",
    "category": "section",
    "text": "AbstractNoiseProcesses can be passed directly to the problem types to replace the standard Wiener process (Brownian motion) with your choice of noise. To do this, simply construct the noise and pass it to the noise keyword argument:μ = 1.0\nσ = 2.0\nW = GeometricBrownianMotionProcess(μ,σ,0.0,1.0,1.0)\n# ...\n# Define f,g,u0,tspan for a SDEProblem\n# ...\nprob = SDEProblem(f,g,u0,tspan,noise=W)"
},

{
    "location": "features/noise_process.html#Basic-Interface-1",
    "page": "Noise Processes",
    "title": "Basic Interface",
    "category": "section",
    "text": "The NoiseProcess acts like a DiffEq solution. For some noise process W, you can get its ith timepoint like W[i] and the associated time W.t[i]. If the NoiseProcess has a bridging distribution defined, it can be interpolated to arbitrary time points using W(t). Note that every interpolated value is saved to the NoiseProcess so that way it can stay distributionally correct. A plot recipe is provided which plots the timeseries."
},

{
    "location": "features/noise_process.html#Direct-Simulation-of-the-Noise-Process-1",
    "page": "Noise Processes",
    "title": "Direct Simulation of the Noise Process",
    "category": "section",
    "text": "Since the NoiseProcess types are distribution-exact and do not require the stochastic differential equation solvers, many times one would like to directly simulate trajectories from these proecesses. The NoiseProcess has a NoiseProcessProblem type:NoiseProblem(noise,tspan)for which solve works. For example, we can simulate a distributionally-exact Geometric Brownian Motion solution by:μ = 1.0\nσ = 2.0\nW = GeometricBrownianMotionProcess(μ,σ,0.0,1.0,1.0)\nprob = NoiseProblem(W,(0.0,1.0))\nsol = solve(prob;dt=0.1)solve requires the dt is given, the solution it returns is a NoiseProcess which has stepped through the timespan. Because this follows the common interface, all of the normal functionality works. For example, we can use the Monte Carlo functionality as follows:monte_prob = MonteCarloProblem(prob)\nsol = solve(monte_prob;dt=0.1,num_monte=100)simulates 100 Geometric Brownian Motions."
},

{
    "location": "features/noise_process.html#Direct-Interface-1",
    "page": "Noise Processes",
    "title": "Direct Interface",
    "category": "section",
    "text": "Most of the time, a NoiseProcess is received from the solution of a stochastic or random differential equation, in which case sol.W gives the NoiseProcess and it is already defined along some timeseries. In other cases, NoiseProcess types are directly simulated (see below). However, NoiseProcess types can also be directly acted on. The basic functionality is given by calculate_step! to calculate a future time point, and accept_step! to accept the step. If steps are rejected, the Rejection Sampling with Memory algorithm is applied to keep the solution distributionally exact. This kind of stepping is done via:W = WienerProcess(0.0,1.0,1.0)\ndt = 0.1\nW.dt = dt\nsetup_next_step!(W)\nfor i in 1:10\n  accept_step!(W,dt)\nend"
},

{
    "location": "features/noise_process.html#Noise-Process-Types-1",
    "page": "Noise Processes",
    "title": "Noise Process Types",
    "category": "section",
    "text": "This section describes the available NoiseProcess types. Note that all keyword arguments are splatted into the NoiseProcess constructor, and thus options like reset are available on the pre-built processes."
},

{
    "location": "features/noise_process.html#Wiener-Process-(White-Noise)-1",
    "page": "Noise Processes",
    "title": "Wiener Process (White Noise)",
    "category": "section",
    "text": "The WienerProcess, also known as Gaussian white noise, Brownian motion, or the noise in the Langevin equation, is the stationary process with distribution N(0,t). The constructor is:WienerProcess(t0,W0,Z0=nothing;kwargs...)\nWienerProcess!(t0,W0,Z0=nothing;kwargs...)"
},

{
    "location": "features/noise_process.html#Correlated-Noise-1",
    "page": "Noise Processes",
    "title": "Correlated Noise",
    "category": "section",
    "text": "One can define a CorrelatedWienerProcess which is a Wiener process with correlations between the Wiener processes. The constructor is:CorrelatedWienerProcess(Γ,t0,W0,Z0=nothing;kwargs...)\nCorrelatedWienerProcess!(Γ,t0,W0,Z0=nothing;kwargs...)where Γ is the constant covariance matrix."
},

{
    "location": "features/noise_process.html#Geometric-Brownian-Motion-1",
    "page": "Noise Processes",
    "title": "Geometric Brownian Motion",
    "category": "section",
    "text": "A GeometricBrownianMotion process is a Wiener process with constant drift μ and constant diffusion σ. I.e. this is the solution of the stochastic differential equationdX_t = mu X_t dt + sigma X_t dW_tThe GeometricBrownianMotionProcess is distribution exact (meaning, not a numerical solution of the stochastic differential equation, and instead follows the exact distribution properties). It can be back interpolated exactly as well. The constructor is:GeometricBrownianMotionProcess(μ,σ,t0,W0,Z0=nothing;kwargs...)\nGeometricBrownianMotionProcess!(μ,σ,t0,W0,Z0=nothing;kwargs...)"
},

{
    "location": "features/noise_process.html#Brownian-Bridge-1",
    "page": "Noise Processes",
    "title": "Brownian Bridge",
    "category": "section",
    "text": "A BrownianBridge process is a Wiener process with a pre-defined start and end value. This process is distribution exact and back be back interpolated exactly as well. The constructor is:BrownianBridge(t0,tend,W0,Wend,Z0=nothing,Zend=nothing;kwargs...)\nBrownianBridge!(t0,tend,W0,Wend,Z0=nothing,Zend=nothing;kwargs...)where W(t0)=W₀, W(tend)=Wend, and likewise for the Z process if defined."
},

{
    "location": "features/noise_process.html#Ornstein-Uhlenbeck-1",
    "page": "Noise Processes",
    "title": "Ornstein-Uhlenbeck",
    "category": "section",
    "text": "One can define a Ornstein-Uhlenbeck process which is a Wiener process defined by the stochastic differential equationdX_t = theta (mu - X_t) dt + sigma dW_tThe OrnsteinUhlenbeckProcess is distribution exact (meaning, not a numerical solution of the stochastic differential equation, and instead follows the exact distribution properties). The constructor is:OrnsteinUhlenbeckProcess(Θ,μ,σ,t0,W0,Z0=nothing;kwargs...)\nOrnsteinUhlenbeckProcess!(Θ,μ,σ,t0,W0,Z0=nothing;kwargs...)"
},

{
    "location": "features/noise_process.html#Direct-Construction-of-a-NoiseProcess-1",
    "page": "Noise Processes",
    "title": "Direct Construction of a NoiseProcess",
    "category": "section",
    "text": "A NoiseProcess is a type defined asNoiseProcess(t0,W0,Z0,dist,bridge;\n             iip=DiffEqBase.isinplace(dist,3),\n             rswm = RSWM(),save_everystep=true,\n             rng = Xorshifts.Xoroshiro128Plus(rand(UInt64)),\n             reset = true, reseed = true)t0 is the first timepoint\nW0 is the first value of the process.\nZ0 is the first value of the psudo-process. This is necessary for higher order algorithms. If it\'s not needed, set to nothing.\ndist the distribution for the steps over time.\nbridge the bridging distribution. Optional, but required for adaptivity and interpolating at new values.\nsave_everystep whether to save every step of the Brownian timeseries.\nrng the local RNG used for generating the random numbers.\nreset whether to reset the process with each solve.\nreseed whether to reseed the process with each solve.The signature for the dist isdist!(rand_vec,W,dt,rng)for inplace functions, andrand_vec = dist(W,dt,rng)otherwise. The signature for bridge isbridge!(rand_vec,W,W0,Wh,q,h,rng)and the out of place syntax isrand_vec = bridge!(W,W0,Wh,q,h,rng)Here, W is the noise process, W0 is the left side of the current interval, Wh is the right side of the current interval, h is the interval length, and q is the proportion from the left where the interpolation is occuring."
},

{
    "location": "features/noise_process.html#Direct-Construction-Example-1",
    "page": "Noise Processes",
    "title": "Direct Construction Example",
    "category": "section",
    "text": "The easiest way to show how to directly construct a NoiseProcess is by example. Here we will show how to directly construct a NoiseProcess which generates Gaussian white noise.This is the noise process which uses randn!. A special dispatch is added for complex numbers for (randn()+im*randn())/sqrt(2). This function is DiffEqBase.wiener_randn (or with ! respectively).The first function that must be defined is the noise distribution. This is how to generate W(t+dt) given that we know W(x) for xt₀t. For Gaussian white noise, we know thatW(dt)  N(0dt)for W(0)=0 which defines the stepping distribution. Thus its noise distribution function is:@inline function WHITE_NOISE_DIST(W,dt,rng)\n  if typeof(W.dW) <: AbstractArray && !(typeof(W.dW) <: SArray)\n    return @fastmath sqrt(abs(dt))*wiener_randn(rng,W.dW)\n  else\n    return @fastmath sqrt(abs(dt))*wiener_randn(rng,typeof(W.dW))\n  end\nendfor the out of place versions, and for the inplace versionsfunction INPLACE_WHITE_NOISE_DIST(rand_vec,W,dt,rng)\n  wiener_randn!(rng,rand_vec)\n  sqrtabsdt = @fastmath sqrt(abs(dt))\n  @. rand_vec *= sqrtabsdt\nendOptionally, we can provide a bridging distribution. This is the distribution of W(qh) for q01 given that we know W(0)=0 and W(h)=Wₕ. For Brownian motion, this is known as the Brownian Bridge, and is well known to have the distribution:W(qh)  N(qWₕ(1-q)qh)Thus we have the out-of-place and in-place versions as:function WHITE_NOISE_BRIDGE(W,W0,Wh,q,h,rng)\n  if typeof(W.dW) <: AbstractArray\n    return @fastmath sqrt((1-q)*q*abs(h))*wiener_randn(rng,W.dW)+q*Wh\n  else\n    return @fastmath sqrt((1-q)*q*abs(h))*wiener_randn(rng,typeof(W.dW))+q*Wh\n  end\nend\nfunction INPLACE_WHITE_NOISE_BRIDGE(rand_vec,W,W0,Wh,q,h,rng)\n  wiener_randn!(rng,rand_vec)\n  #rand_vec .= sqrt((1.-q).*q.*abs(h)).*rand_vec.+q.*Wh\n  sqrtcoeff = @fastmath sqrt((1-q)*q*abs(h))\n  @. rand_vec = sqrtcoeff*rand_vec+q*Wh\nendThese functions are then placed in a noise process:NoiseProcess(t0,W0,Z0,WHITE_NOISE_DIST,WHITE_NOISE_BRIDGE;kwargs)\nNoiseProcess(t0,W0,Z0,INPLACE_WHITE_NOISE_DIST,INPLACE_WHITE_NOISE_BRIDGE;kwargs)Notice that we can optionally provide an alternative adaptive algorithm for the timestepping rejections. RSWM() defaults to the Rejection Sampling with Memory 3 algorithm (RSwM3).Note that the standard constructors are simply:WienerProcess(t0,W0,Z0=nothing) = NoiseProcess(t0,W0,Z0,WHITE_NOISE_DIST,WHITE_NOISE_BRIDGE;kwargs)\nWienerProcess!(t0,W0,Z0=nothing) = NoiseProcess(t0,W0,Z0,INPLACE_WHITE_NOISE_DIST,INPLACE_WHITE_NOISE_BRIDGE;kwargs)These will generate a Wiener process, which can be stepped with step!(W,dt), and interpolated as W(t)."
},

{
    "location": "features/noise_process.html#Non-Standard-Noise-Processes-1",
    "page": "Noise Processes",
    "title": "Non-Standard Noise Processes",
    "category": "section",
    "text": "In addition to the mathematically-defined noise processes above, there exists more generic functionality for building noise processes from other noise processes, from arbitrary functions, from arrays, and from approximations of stochastic differential equations."
},

{
    "location": "features/noise_process.html#NoiseWrapper-1",
    "page": "Noise Processes",
    "title": "NoiseWrapper",
    "category": "section",
    "text": "This produces a new noise process from an old one, which will use its interpolation to generate the noise. This allows you to re-use a previous noise process not just with the same timesteps, but also with new (adaptive) timesteps as well. Thus this is very good for doing Multi-level Monte Carlo schemes and strong convergence testing.To wrap a noise process, simply use:NoiseWrapper(W::NoiseProcess;reset=true)"
},

{
    "location": "features/noise_process.html#NoiseFunction-1",
    "page": "Noise Processes",
    "title": "NoiseFunction",
    "category": "section",
    "text": "This allows you to use any arbitrary function W(t) as a NoiseProcess. This will use the function lazily, only caching values required to minimize function calls, but not store the entire noise array. This requires an initial time point t0 in the domain of W. A second function is needed if the desired SDE algorithm requires multiple processes.NoiseFunction(t0,W,Z=nothing;noise_prototype=W(t0),reset=true)Additionally, one can use an in-place function W(out1,out2,t) for more efficient generation of the arrays for multi-dimensional processes. When the in-place version is used without a dispatch for the out-of-place version, the noise_prototype needs to be set."
},

{
    "location": "features/noise_process.html#NoiseGrid-1",
    "page": "Noise Processes",
    "title": "NoiseGrid",
    "category": "section",
    "text": "A noise grid builds a noise process from arrays of points. For example, you can generate your desired noise process as an array W with timepoints t, and use the constructor:NoiseGrid(t,W,Z=nothing;reset=true)to build the associated noise process. This process comes with a linear interpolation of the given points, and thus the grid does not have to match the grid of integration. Thus this can be used for adaptive solutions as well. However, one must make note that the fidelity of the noise process is linked to how fine the noise grid is determined: if the noise grid is sparse on points compared to the integration, then its distributional properties may be slightly perturbed by the linear interpolation. Thus its suggested that the grid size at least approximately match the number of time steps in the integration to ensure accuracy.For a one-dimensional process, W should be an AbstractVector of Numbers. For multi-dimensional processes, W should be an AbstractVector of the noise_prototype."
},

{
    "location": "features/noise_process.html#NoiseApproximation-1",
    "page": "Noise Processes",
    "title": "NoiseApproximation",
    "category": "section",
    "text": "In many cases, one would like to define a noise process directly by a stochastic differential equation which does not have an analytical solution. Of course, this will not be distributionally-exact and how well the properties match depends on how well the differential equation is integrated, but in many cases this can be used as a good approximation when other methods are much more difficult.A NoiseApproximation is defined by a DEIntegrator. The constructor for a NoiseApproximation is:NoiseApproximation(source1::DEIntegrator,source2::Union{DEIntegrator,Nothing}=nothing;reset=true)The DEIntegrator should have a final time point of integration far enough such that it will not halt during the integration. For ease of use, you can use a final time point as Inf. Note that the time points do not have to match the time points of the future integration since the interpolant of the SDE solution will be used. Thus the limiting factor is error tolerance and not hitting specific points."
},

{
    "location": "features/noise_process.html#Examples-Using-Non-Standard-Noise-Processes-1",
    "page": "Noise Processes",
    "title": "Examples Using Non-Standard Noise Processes",
    "category": "section",
    "text": ""
},

{
    "location": "features/noise_process.html#NoiseGrid-2",
    "page": "Noise Processes",
    "title": "NoiseGrid",
    "category": "section",
    "text": "In this example, we will show you how to define your own version of Brownian motion using an array of pre-calculated points. In normal usage you should use WienerProcess instead since this will have distributionally-exact interpolations while the noise grid uses linear interpolations, but this is a nice example of the workflow.To define a NoiseGrid you need to have a set of time points and a set of values for the process. Let\'s define a Brownian motion on (0.0,1.0) with a dt=0.001. To do this,dt = 0.001\nt = 0:dt:1\nbrownian_values = cumsum([0;[sqrt(dt)*randn() for i in 1:length(t)-1]])Now we build the NoiseGrid using these values:W = NoiseGrid(t,brownian_values)We can then pass W as the noise argument of an SDEProblem to use it in an SDE."
},

{
    "location": "features/noise_process.html#NoiseWrapper-Example-1",
    "page": "Noise Processes",
    "title": "NoiseWrapper Example",
    "category": "section",
    "text": "In this example, we will solve an SDE three times:First to generate a noise process\nSecond with the same timesteps to show the values are the same\nThird with half-sized timstepsFirst we will generate a noise process by solving an SDE:using StochasticDiffEq,  DiffEqBase, DiffEqNoiseProcess\nf1(t,u) = 1.01u\ng1(t,u) = 1.01u\ndt = 1//2^(4)\nprob1 = SDEProblem(f1,g1,1.0,(0.0,1.0))\nsol1 = solve(prob1,EM(),dt=dt,save_noise = true)Now we wrap the noise into a NoiseWrapper and solve the same problem:W2 = NoiseWrapper(sol1.W)\nprob1 = SDEProblem(f1,g1,1.0,(0.0,1.0),noise=W2)\nsol2 = solve(prob1,EM(),dt=dt)We can test@test sol1.u ≈ sol2.uto see that the values are essentially equal. Now we can use the same process to solve the same trajectory with a smaller dt:W3 = NoiseWrapper(sol1.W)\nprob2 = SDEProblem(f1,g1,1.0,(0.0,1.0),noise=W3)\n\ndt = 1//2^(5)\nsol3 = solve(prob2,EM(),dt=dt)We can plot the results to see what this looks like:using Plots\nplot(sol1)\nplot!(sol2)\nplot!(sol3)(Image: noise_process)In this plot, sol2 covers up sol1 because they hit essentially the same values. You can see that sol3 its similar to the others, because it\'s using the same underlying noise process just sampled much finer.To double check, we see that:plot(sol1.W)\nplot!(sol2.W)\nplot!(sol3.W)(Image: coupled_wiener)the coupled Wiener processes coincide at every other time point, and the intermediate timepoints were calculated according to a Brownian bridge."
},

{
    "location": "features/noise_process.html#Adaptive-NoiseWrapper-Example-1",
    "page": "Noise Processes",
    "title": "Adaptive NoiseWrapper Example",
    "category": "section",
    "text": "Here we will show that the same noise can be used with the adaptive methods using the NoiseWrapper. SRI and SRIW1 use slightly different error estimators, and thus give slightly different stepping behavior. We can see how they solve the same 2D SDE differently by using the noise wrapper:prob = SDEProblem(f1,g1,ones(2),(0.0,1.0))\nsol4 = solve(prob,SRI(),abstol=1e-8, save_noise = true)\n\nW2 = NoiseWrapper(sol4.W)\nprob2 = SDEProblem(f1,g1,ones(2),(0.0,1.0),noise=W2)\nsol5 = solve(prob2,SRIW1(),abstol=1e-8)\n\nusing Plots\nplot(sol4)\nplot!(sol5)(Image: SRI_SRIW1_diff)"
},

{
    "location": "features/noise_process.html#NoiseApproximation-Example-1",
    "page": "Noise Processes",
    "title": "NoiseApproximation Example",
    "category": "section",
    "text": "In this example we will show how to use the NoiseApproximation in order to build our own Geometric Brownian Motion from its stochastic differential equation definition. In normal usage, you should use the GeometricBrownianMotionProcess instead since that is more efficient and distributionally-exact.First, let\'s define the SDEProblem. Here will use a timespan (0.0,Inf) so that way the noise can be used over an indefinite integral.const μ = 1.5\nconst σ = 1.2\nf(t,u) = μ*u\ng(t,u) = σ*u\nprob = SDEProblem(f,g,1.0,(0.0,Inf))Now we build the noise process by building the integrator and sending that integrator to the NoiseApproximation constructor:integrator = init(prob,SRIW1())\nW = NoiseApproximation(integrator)We can use this noise process like any other noise process. For example, we can now build a geometric Brownian motion whose noise process is colored noise that itself is a geometric Brownian motion:prob = SDEProblem(f,g,1.0,(0.0,Inf),noise=W)The possibilities are endless."
},

{
    "location": "features/noise_process.html#NoiseFunction-Example-1",
    "page": "Noise Processes",
    "title": "NoiseFunction Example",
    "category": "section",
    "text": "The NoiseFunction is pretty simple: pass a function. As a silly example, we can use exp as a noise process by doing:f(t) = exp(t)\nW = NoiseFunction(0.0,f)If it\'s multi-dimensional and an in-place function is used, the noise_prototype must be given. For example:f(out,t) = (out.=exp(t))\nW = NoiseFunction(0.0,f,noise_prototype=rand(4))This allows you to put arbitrarily weird noise into SDEs and RODEs. Have fun."
},

{
    "location": "features/linear_nonlinear.html#",
    "page": "Specifying (Non)Linear Solvers",
    "title": "Specifying (Non)Linear Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "features/linear_nonlinear.html#Specifying-(Non)Linear-Solvers-1",
    "page": "Specifying (Non)Linear Solvers",
    "title": "Specifying (Non)Linear Solvers",
    "category": "section",
    "text": "One of the key features of DifferentialEquations.jl is its flexibility. Keeping with this trend, many of the native Julia solvers provided by DifferentialEquations.jl allow you to choose the method for linear and nonlinear solving. This section details how to make that choice."
},

{
    "location": "features/linear_nonlinear.html#Linear-Solvers:-linsolve-Specification-1",
    "page": "Specifying (Non)Linear Solvers",
    "title": "Linear Solvers: linsolve Specification",
    "category": "section",
    "text": "For differential equation integrators which use linear solvers, an argument to the method linsolve determines the linear solver which is used. The signature is:linsolve! = linsolve(Val{:init},f,x;kwargs...)\nlinsolve!(x,A,b,matrix_updated=false;kwargs...)This is an in-place function which updates x by solving Ax=b. The user should specify the function linsolve(Val{:init},f,x) which returns a linsolve! function. The setting matrix_updated determines whether the matrix A has changed from the last call. This can be used to smartly cache factorizations.Note that linsolve! needs to accept splatted keyword arguments. The possible arguments passed to the linear solver are as follows:Pl, a pre-specified left preconditioner which utilizes the internal adaptive norm estimates\nPr, a pre-specified right preconditioner which utilizes the internal adaptive norm estimates\ntol, a linear solver tolerance specified from the ODE solver\'s implicit handling"
},

{
    "location": "features/linear_nonlinear.html#Pre-Built-Linear-Solver-Choices-1",
    "page": "Specifying (Non)Linear Solvers",
    "title": "Pre-Built Linear Solver Choices",
    "category": "section",
    "text": "The following choices of pre-built linear solvers exist:DefaultLinSolve\nLinSolveFactorize\nLinSolveGPUFactorize\nLinSolveGMRES\nLinSolveCG\nLinSolveBiCGStabl\nLinSolveChebyshev\nLinSolveMINRES\nLinSolveIterativeSolvers"
},

{
    "location": "features/linear_nonlinear.html#DefaultLinSolve-1",
    "page": "Specifying (Non)Linear Solvers",
    "title": "DefaultLinSolve",
    "category": "section",
    "text": "The default linear solver is DefaultLinSolve. This method is adaptive, and automatically chooses an LU factorization choose for dense and sparse arrays, and is compatible with GPU-based arrays. When the Jacobian is an AbstractDiffEqOperator, i.e. is matrix-free, DefaultLinSolve defaults to using a gmres iterative solver."
},

{
    "location": "features/linear_nonlinear.html#Basic-linsolve-method-choice:-Factorization-by-LinSolveFactorize-1",
    "page": "Specifying (Non)Linear Solvers",
    "title": "Basic linsolve method choice: Factorization by LinSolveFactorize",
    "category": "section",
    "text": "The easiest way to specify a linsolve is by a factorization function which generates a type on which \\ (or A_ldiv_B!) is called.  This is done through the helper function LinSolveFactorize which makes the appropriate function. For example, the  Rosenbrock23 takes in a linsolve function, which we can choose to be a  QR-factorization by:Rosenbrock23(linsolve=LinSolveFactorize(qrfact!))LinSolveFactorize takes in a function which returns an object that can \\. Direct methods like qrfact! will automatically cache the factorization, making it efficient for small dense problems.However, for large sparse problems, you can let \\ be an iterative method. For example, using PETSc.jl, we can define our factorization function to be:linsolve = LinSolveFactorize((A) -> KSP(A, ksp_type=\"gmres\", ksp_rtol=1e-6))This function creates a KSP type which makes \\ perform the GMRES iterative method provided by PETSc.jl. Thus if we pass this function into the algorithm as the factorization method, all internal linear solves will happen by PETSc.jl."
},

{
    "location": "features/linear_nonlinear.html#GPU-offloading-of-factorization-with-LinSolveGPUFactorize-1",
    "page": "Specifying (Non)Linear Solvers",
    "title": "GPU offloading of factorization with LinSolveGPUFactorize",
    "category": "section",
    "text": "If one has a problem with a sufficiently large Jacobian (~100x100) and a sufficiently powerful GPU, it can make sense to offload the factorization and backpropogation steps to the GPU. For this, the LinSolveGPUFactorize linear solver is provided. It works similarly to LinSolveFactorize, but the matrix is automatically sent to the GPU as a CuArray and the ldiv! is performed against a CUDA QR factorization of the matrix.Note that this method requires that you have done using CuArrays in your script. A working installation of CuArrays.jl is required, which requires an installation of CUDA Toolkit."
},

{
    "location": "features/linear_nonlinear.html#IterativeSolvers.jl-Based-Methods-1",
    "page": "Specifying (Non)Linear Solvers",
    "title": "IterativeSolvers.jl-Based Methods",
    "category": "section",
    "text": "The signature for LinSolveIterativeSolvers is:LinSolveIterativeSolvers(generate_iterator,args...;\n                         Pl=IterativeSolvers.Identity(),\n                         Pr=IterativeSolvers.Identity(),\n                         kwargs...)where Pl is the left preconditioner, Pr is the right preconditioner, and the other args... and kwargs... are passed into the iterative solver chosen in generate_iterator which designates the construction of an iterator from IterativeSolvers.jl. For example, using gmres_iterable! would make a version that uses IterativeSolvers.gmres. The following are aliases to common choices:LinSolveGMRES – GMRES\nLinSolveCG – CG (Conjugate Gradient)\nLinSolveBiCGStabl – BiCGStabl Stabilized Bi-Conjugate Gradient\nLinSolveChebyshev – Chebyshev\nLinSolveMINRES – MINRESwhich all have the same arguments as LinSolveIterativeSolvers except with generate_iterator pre-specified."
},

{
    "location": "features/linear_nonlinear.html#Implementing-Your-Own-LinSolve:-How-LinSolveFactorize-Was-Created-1",
    "page": "Specifying (Non)Linear Solvers",
    "title": "Implementing Your Own LinSolve: How LinSolveFactorize Was Created",
    "category": "section",
    "text": "In order to make your own linsolve functions, let\'s look at how the LinSolveFactorize function is created. For example, for an LU-Factorization, we would like to use lufact! to do our linear solving. We can directly write this as:function linsolve!(::Type{Val{:init}},f,u0,kwargs...)\n  function _linsolve!(x,A,b,update_matrix=false,kwargs...)\n    _A = lufact!(A)\n    ldiv!(x,_A,b)\n  end\nendThis initialization function returns a linear solving function that always computes the LU-factorization and then does the solving. This method works fine and you can pass it to the methods likeRosenbrock23(linsolve=linsolve!)and it will work, but this method does not cache _A, the factorization. This means that, even if A has not changed, it will re-factorize the matrix.To change this, we can instead create a call-overloaded type. The generalized form of this is:mutable struct LinSolveFactorize{F}\n  factorization::F\n  A\nend\nLinSolveFactorize(factorization) = LinSolveFactorize(factorization,nothing)\nfunction (p::LinSolveFactorize)(x,A,b,matrix_updated=false)\n  if matrix_updated\n    p.A = p.factorization(A)\n  end\n  A_ldiv_B!(x,p.A,b)\nend\nfunction (p::LinSolveFactorize)(::Type{Val{:init}},f,u0_prototype)\n  LinSolveFactorize(p.factorization,nothing)\nend\nlinsolve = LinSolveFactorize(lufact!)LinSolveFactorize is a type which holds the factorization method and the pre-factorized matrix. When linsolve is passed to the ODE/SDE/etc. solver, it will use the function linsolve(Val{:init},f,u0_prototype) to create a LinSolveFactorize object which holds the factorization method and a cache for holding a factorized matrix. Thenfunction (p::LinSolveFactorize)(x,A,b,matrix_updated=false)\n  if matrix_updated\n    p.A = p.factorization(A)\n  end\n  A_ldiv_B!(x,p.A,b)\nendis what\'s used in the solver\'s internal loop. If matrix_updated is true, it will re-compute the factorization. Otherwise it just solves the linear system with the cached factorization. This general idea of using a call-overloaded type can be employed to do many other things."
},

{
    "location": "features/linear_nonlinear.html#Nonlinear-Solvers:-nlsolve-Specification-1",
    "page": "Specifying (Non)Linear Solvers",
    "title": "Nonlinear Solvers: nlsolve Specification",
    "category": "section",
    "text": "Nonlinear solvers can be chosen via the nlsolve option. Most algorithms use nonlinear solvers that are specialized for implicit ODE solvers. There are three pre-built nlsolves:NLNewton(): It is a modified Newton iteration solver, and it is the default nlsolve for most of the implicit ODE solvers. It converges the fastest, but requires more memory usage and linear system solve.\nNLAnderson(n::Int): It is an Anderson acceleration solver. It converges faster than NLFunctional but slower than NLNewton. It does not require to solve a linear system. In development.\nNLFunctional(): It is a functional (Picard) iteration solver. It converges the slowest, but requires the least amount of memory.One can specify a nonlinear solver byImplicitEuler(nlsolve = NLFunctional())"
},

{
    "location": "features/linear_nonlinear.html#Nonlinear-Solvers-for-Generic-Implicit-ODE-Solvers-1",
    "page": "Specifying (Non)Linear Solvers",
    "title": "Nonlinear Solvers for Generic Implicit ODE Solvers",
    "category": "section",
    "text": "For ODE solvers with names that begin with Generic, they take more generic nlsolve. An nlsolve function should have two dispatches:nlsolve(Val{init},f,u0_prototype) : Does an initialization phase. Returns a type init_f for later use in the solver. u0_prototype is the expected type for the initial condition u0.\nnlsolve(init_f,u0) : Solves for the root units the initialized f and the initial condition u0. Returns the zeros of the equation."
},

{
    "location": "features/linear_nonlinear.html#Basic-nlsolve-method:-NLSOLVEJL_SETUP-1",
    "page": "Specifying (Non)Linear Solvers",
    "title": "Basic nlsolve method: NLSOLVEJL_SETUP",
    "category": "section",
    "text": "By default, a basic nonlinear solver setup is given as NLSOLVEJL_SETUP. For example, the default nlsolve in GenericTrapezoid isGenericTrapezoid(nlsolve=NLSOLVEJL_SETUP())This will use NLsolve.jl with autodifferentiation to solve the nonlinear systems. NLSOLVEJL_SETUP has two options:chunk_size : The autodifferentiation chunk size. Integer. Defaults to ForwardDiff.jl\'s auto-detection.\nautodiff : Whether to use autodifferentiation. Defaults to true.For example, to turn off autodifferentiation, useTrapezoid(nlsolve=NLSOLVEJL_SETUP(autodiff=false))"
},

{
    "location": "features/linear_nonlinear.html#How-NLSOLVEJL_SETUP-Was-Created-1",
    "page": "Specifying (Non)Linear Solvers",
    "title": "How NLSOLVEJL_SETUP Was Created",
    "category": "section",
    "text": "To create a nonlinear solver, you need to define the two functions. Here we use a call-overloaded type so that way we can hold the chunk size and autodifferentiation information.struct NLSOLVEJL_SETUP{CS,AD} end\nNLSOLVEJL_SETUP(;chunk_size=0,autodiff=true) = NLSOLVEJL_SETUP{chunk_size,autodiff}()The solver function just calls NLsolve and returns the zeros(p::NLSOLVEJL_SETUP)(f,u0) = (res=NLsolve.nlsolve(f,u0); res.zero)while the initialization function has a different initialization for autodifferentiation or not:function (p::NLSOLVEJL_SETUP{CS,AD})(::Type{Val{:init}},f,u0_prototype) where {CS,AD}\n  if AD\n    return non_autodiff_setup(f,u0_prototype)\n  else\n    return autodiff_setup(f,u0_prototype,Val{determine_chunksize(initial_x,CS)})\n  end\nendWe need to declare the get_chunksize trait for the solver:get_chunksize(x::NLSOLVEJL_SETUP{CS,AD}) where {CS,AD} = CSThe initialization functions are directly for NLsolve. See the NLsolve.jl docs for the types of inputs it expects to see. This does exactly that:function autodiff_setup(f!, initial_x::Vector,chunk_size::Type{Val{CS}}) where CS\n\n    permf! = (fx, x) -> f!(x, fx)\n\n    fx2 = copy(initial_x)\n    jac_cfg = ForwardDiff.JacobianConfig{CS}(initial_x, initial_x)\n    g! = (x, gx) -> ForwardDiff.jacobian!(gx, permf!, fx2, x, jac_cfg)\n\n    fg! = (x, fx, gx) -> begin\n        jac_res = DiffBase.DiffResult(fx, gx)\n        ForwardDiff.jacobian!(jac_res, permf!, fx2, x, jac_cfg)\n        DiffBase.value(jac_res)\n    end\n\n    return DifferentiableMultivariateFunction(f!, g!, fg!)\nend\n\nfunction non_autodiff_setup(f!, initial_x::Vector)\n  DifferentiableMultivariateFunction(f!)\nend"
},

{
    "location": "features/callback_functions.html#",
    "page": "Event Handling and Callback Functions",
    "title": "Event Handling and Callback Functions",
    "category": "page",
    "text": ""
},

{
    "location": "features/callback_functions.html#Event-Handling-and-Callback-Functions-1",
    "page": "Event Handling and Callback Functions",
    "title": "Event Handling and Callback Functions",
    "category": "section",
    "text": ""
},

{
    "location": "features/callback_functions.html#Introduction-to-Callback-Functions-1",
    "page": "Event Handling and Callback Functions",
    "title": "Introduction to Callback Functions",
    "category": "section",
    "text": "DifferentialEquations.jl allows for using callback functions to inject user code into the solver algorithms. It allows for safely and accurately applying events and discontinuities. Multiple callbacks can be chained together, and these callback types can be used to build libraries of extension behavior."
},

{
    "location": "features/callback_functions.html#The-Callback-Types-1",
    "page": "Event Handling and Callback Functions",
    "title": "The Callback Types",
    "category": "section",
    "text": "The callback types are defined as follows. There are two callback types: the ContinuousCallback and the DiscreteCallback. The ContinuousCallback is applied when a continuous condition function hits zero. This type of callback implements what is known in other problem solving environments as an Event. A DiscreteCallback is applied when its condition function is true."
},

{
    "location": "features/callback_functions.html#ContinuousCallbacks-1",
    "page": "Event Handling and Callback Functions",
    "title": "ContinuousCallbacks",
    "category": "section",
    "text": "ContinuousCallback(condition,affect!,affect_neg!=affect!;\n                   rootfind = true,\n                   initialize = (c,u,t,integrator) -> nothing,\n                   save_positions = (true,true),\n                   interp_points=10,\n                   abstol=1e-9,reltol=0\n                   idxs=nothing)The arguments are defined as follows:condition: This is a function condition(u,t,integrator) for declaring when the callback should be used. A callback is initiated if the condition hits 0 within the time interval. See the Integrator Interface documentation for information about integrator.\naffect!: This is the function affect!(integrator) where one is allowed to modify the current state of the integrator. If you do not pass an affect_neg! function, it is called when condition is found to be 0 (at a root) and the cross is either an upcrossing (from negative to positive) or a downcrossing (from positive to negative). You need to explicitly pass nothing as the affect_neg! argument if it should only be called at upcrossings, e.g. ContinuousCallback(condition, affect!, nothing). For more information on what can be done, see the Integrator Interface manual page. Modifications to u are safe in this function.\naffect_neg!: This is the function affect_neg!(integrator) where one is allowed to modify the current state of the integrator. This is called when condition is found to be 0 (at a root) and the cross is an downcrossing (from positive to negative). For more information on what can be done, see the Integrator Interface manual page. Modifications to u are safe in this function.\nrootfind: This is a boolean for whether to rootfind the event location. If this is set to true, the solution will be backtracked to the point where condition==0. Otherwise the systems and the affect! will occur at t+dt.\ninterp_points: The number of interpolated points to check the condition. The condition is found by checking whether any interpolation point / endpoint has a different sign. If interp_points=0, then conditions will only be noticed if the sign of condition is different at t than at t+dt. This behavior is not robust when the solution is oscillatory, and thus it\'s recommended that one use some interpolation points (they\'re cheap to compute!). 0 within the time interval.\nsave_positions: Boolean tuple for whether to save before and after the affect!. This saving will occur just before and after the event, only at event times, and does not depend on options like saveat, save_everystep, etc. (i.e. if saveat=[1.0,2.0,3.0], this can still add a save point at 2.1 if true). For discontinuous changes like a modification to u to be handled correctly (without error), one should set save_positions=(true,true).\nidxs: The components which will be interpolated into the condition. Defaults to nothing which means u will be all components.\ninitialize: This is a function (c,u,t,integrator) which can be used to initialize the state of the callback c. It should modify the argument c and the return is ignored.Additionally, keyword arguments for abstol and reltol can be used to specify a tolerance from zero for the rootfinder: if the starting condition is less than the tolerance from zero, then no root will be detected. This is to stop repeat events happening just after a previously rootfound event. The default has abstol=1e-14 and reltol=0."
},

{
    "location": "features/callback_functions.html#DiscreteCallback-1",
    "page": "Event Handling and Callback Functions",
    "title": "DiscreteCallback",
    "category": "section",
    "text": "DiscreteCallback(condition,affect!;\n                 save_positions=(true,true),\n                 initialize = (c,u,t,integrator) -> nothing)condition: This is a function condition(u,t,integrator) for declaring when the callback should be used. A callback is initiated if the condition evaluates to true. See the Integrator Interface documentation for information about integrator.\naffect!: This is the function affect!(integrator) where one is allowed to modify the current state of the integrator. For more information on what can be done, see the Integrator Interface manual page.\nsave_positions: Boolean tuple for whether to save before and after the affect!. This saving will occur just before and after the event, only at event times, and does not depend on options like saveat, save_everystep, etc. (i.e. if saveat=[1.0,2.0,3.0], this can still add a save point at 2.1 if true). For discontinuous changes like a modification to u to be handled correctly (without error), one should set save_positions=(true,true).\ninitialize: This is a function (c,u,t,integrator) which can be used to initialize the state of the callback c. It should modify the argument c and the return is ignored."
},

{
    "location": "features/callback_functions.html#CallbackSet-1",
    "page": "Event Handling and Callback Functions",
    "title": "CallbackSet",
    "category": "section",
    "text": "Multiple callbacks can be chained together to form a CallbackSet. A CallbackSet is constructed by passing the constructor ContinuousCallback, DiscreteCallback,  VectorContinuousCallback or other CallbackSet instances:CallbackSet(cb1,cb2,cb3)You can pass as many callbacks as you like. When the solvers encounter multiple callbacks, the following rules apply:ContinuousCallbacks and VectorContinuousCallbacks are applied before DiscreteCallbacks. (This is because they often implement event-finding that will backtrack the timestep to smaller than dt).\nFor ContinuousCallbacks and VectorContinuousCallbacks, the event times are found by rootfinding and only the first ContinuousCallback or VectorContinuousCallback affect is applied.\nThe DiscreteCallbacks are then applied in order. Note that the ordering only matters for the conditions: if a previous callback modifies u in such a way that the next callback no longer evaluates condition to true, its affect will not be applied."
},

{
    "location": "features/callback_functions.html#VectorContinuousCallback-1",
    "page": "Event Handling and Callback Functions",
    "title": "VectorContinuousCallback",
    "category": "section",
    "text": "VectorContinuousCallback(condition,affect!,len;\n                   initialize = INITIALIZE_DEFAULT,\n                   idxs = nothing,\n                   rootfind=true,\n                   save_positions=(true,true),\n                   affect_neg! = affect!,\n                   interp_points=10,\n                   abstol=10eps(),reltol=0)VectorContinuousCallback is also a subtype of AbstractContinuousCallback. CallbackSet is not feasible when you have a large number of callbacks, as it doesn\'t scale well. For this reason, we have VectorContinuousCallback - it allows you to have a single callback for multiple events. condition - This is a function condition(out, u, t, integrator) which should save the condition value in the array out at the right index. Maximum index of out should be specified in the len property of callback. So this way you can have a chain of len events, which would cause the ith event to trigger when out[i] = 0.\naffect! - This is a function affect!(integrator, event_index) which lets you modify integrator and it tells you about which event occured using event_idx i.e. gives you index i for which out[i] came out to be zero.\nlen - Number of callbacks chained. This is compulsory to be specified.Rest of the fields have the same meaning as ContinuousCallback."
},

{
    "location": "features/callback_functions.html#Using-Callbacks-1",
    "page": "Event Handling and Callback Functions",
    "title": "Using Callbacks",
    "category": "section",
    "text": "The callback type is then sent to the solver (or the integrator) via the callback keyword argument:sol = solve(prob,alg,callback=cb)You can supply nothing, a single DiscreteCallback or ContinuousCallback, or a CallbackSet."
},

{
    "location": "features/callback_functions.html#Note-About-Saving-1",
    "page": "Event Handling and Callback Functions",
    "title": "Note About Saving",
    "category": "section",
    "text": "When a callback is supplied, the default saving behavior is turned off. This is because otherwise events would \"double save\" one of the values. To re-enable the standard saving behavior, one must have the first save_positions value be true for at least one callback."
},

{
    "location": "features/callback_functions.html#Modifying-the-Stepping-Within-A-Callback-1",
    "page": "Event Handling and Callback Functions",
    "title": "Modifying the Stepping Within A Callback",
    "category": "section",
    "text": "A common issue with callbacks is that they cause a large discontinuous change, and so it may be wise to pull down dt after such a change. To control the timestepping from a callback, please see the timestepping controls in the integrator interface. Specifically, set_proposed_dt! is used to set the next stepsize, and terminate! can be used to cause the simulation to stop."
},

{
    "location": "features/callback_functions.html#DiscreteCallback-Examples-1",
    "page": "Event Handling and Callback Functions",
    "title": "DiscreteCallback Examples",
    "category": "section",
    "text": ""
},

{
    "location": "features/callback_functions.html#Example-1:-AutoAbstol-1",
    "page": "Event Handling and Callback Functions",
    "title": "Example 1: AutoAbstol",
    "category": "section",
    "text": "MATLAB\'s Simulink has the option for an automatic absolute tolerance. In this example we will implement a callback which will add this behavior to any JuliaDiffEq solver which implments the integrator and callback interface.The algorithm is as follows. The default value is set to start at 1e-6, though we will give the user an option for this choice. Then as the simulation progresses, at each step the absolute tolerance is set to the maximum value that has been reached so far times the relative tolerance. This is the behavior that we will implement in affect!.Since the effect is supposed to occur every timestep, we use the trivial condition:condition = function (u,t,integrator)\n    true\nendwhich always returns true. For our effect we will overload the call on a type. This type will have a value for the current maximum. By doing it this way, we can store the current state for the running maximum. The code is as follows:mutable struct AutoAbstolAffect{T}\n  curmax::T\nend\n# Now make `affect!` for this:\nfunction (p::AutoAbstolAffect)(integrator)\n  p.curmax = max(p.curmax,integrator.u)\n  integrator.opts.abstol = p.curmax * integrator.opts.reltol\n  u_modified!(integrator,false)\nendThis makes affect!(integrator) use an internal mutating value curmax to update the absolute tolerance of the integrator as the algorithm states.Lastly, we can wrap it in a nice little constructor:function AutoAbstol(save=true;init_curmax=1e-6)\n  affect! = AutoAbstolAffect(init_curmax)\n  condtion = (u,t,integrator) -> true\n  save_positions = (save,false)\n  DiscreteCallback(condtion,affect!,save_positions=save_positions)\nendThis creates the DiscreteCallback from the affect! and condition functions that we implemented. Nowcb = AutoAbstol(save=true;init_curmax=1e-6)returns the callback that we created. We can then solve an equation using this by simply passing it with the callback keyword argument. Using the integrator interface rather than the solve interface, we can step through one by one to watch the absolute tolerance increase:integrator = init(prob,BS3(),callback=cb)\nat1 = integrator.opts.abstol\nstep!(integrator)\nat2 = integrator.opts.abstol\n@test at1 < at2\nstep!(integrator)\nat3 = integrator.opts.abstol\n@test at2 < at3Note that this example is contained in DiffEqCallbacks.jl, a library of useful callbacks for JuliaDiffEq solvers."
},

{
    "location": "features/callback_functions.html#Example-2:-A-Control-Problem-1",
    "page": "Event Handling and Callback Functions",
    "title": "Example 2: A Control Problem",
    "category": "section",
    "text": "Another example of a DiscreteCallback is the control problem demonstrated on the DiffEq-specific arrays page."
},

{
    "location": "features/callback_functions.html#ContinuousCallback-Examples-1",
    "page": "Event Handling and Callback Functions",
    "title": "ContinuousCallback Examples",
    "category": "section",
    "text": ""
},

{
    "location": "features/callback_functions.html#Example-1:-Bouncing-Ball-1",
    "page": "Event Handling and Callback Functions",
    "title": "Example 1: Bouncing Ball",
    "category": "section",
    "text": "Let\'s look at the bouncing ball. Let the first variable y is the height which changes by v the velocity, where the velocity is always changing at -g which is the gravitational constant. This is the equation:function f(du,u,p,t)\n  du[1] = u[2]\n  du[2] = -p\nendAll we have to do in order to specify the event is to have a function which should always be positive with an event occurring at 0. For now at least that\'s how it\'s specified. If a generalization is needed we can talk about this (but it needs to be \"root-findable\"). For here it\'s clear that we just want to check if the ball\'s height ever hits zero:function condition(u,t,integrator) # Event when event_f(u,t) == 0\n  u[1]\nendNotice that here we used the values u instead of the value from the integrator. This is because the values u,t will be appropriately modified at the interpolation points, allowing for the rootfinding behavior to occur.Now we have to say what to do when the event occurs. In this case we just flip the velocity (the second variable)function affect!(integrator)\n  integrator.u[2] = -integrator.u[2]\nendThe callback is thus specified by:cb = ContinuousCallback(condition,affect!)Then you can solve and plot:u0 = [50.0,0.0]\ntspan = (0.0,15.0)\np = 9.8\nprob = ODEProblem(f,u0,tspan,p)\nsol = solve(prob,Tsit5(),callback=cb)\nplot(sol)(Image: BallBounce)As you can see from the resulting image, DifferentialEquations.jl is smart enough to use the interpolation to hone in on the time of the event and apply the event back at the correct time. Thus one does not have to worry about the adaptive timestepping \"overshooting\" the event as this is handled for you. Notice that the event macro will save the value(s) at the discontinuity.The callback is robust to having multiple discontinuities occur. For example, we can integrate for long time periods and get the desired behavior:u0 = [50.0,0.0]\ntspan = (0.0,100.0)\nprob = ODEProblem(f,u0,tspan,p)\nsol = solve(prob,Tsit5(),callback=cb)\nplot(sol,plotdensity=10000)(Image: bounce_long)"
},

{
    "location": "features/callback_functions.html#Example-2:-Terminating-an-Integration-1",
    "page": "Event Handling and Callback Functions",
    "title": "Example 2: Terminating an Integration",
    "category": "section",
    "text": "In many cases you might want to terminate an integration when some condition is satisfied. To terminate an integration, use terminate!(integrator) as the affect! in a callback.In this example we will solve the differential equation:u0 = [1.,0.]\nfunction fun2(du,u,p,t)\n   du[2] = -u[1]\n   du[1] = u[2]\nend\ntspan = (0.0,10.0)\nprob = ODEProblem(fun2,u0,tspan)which has cosine and -sine as the solutions respectively. We wish to solve until the sine part, u[2] becomes positive. There are two things we may be looking for.A DiscreteCallback will cause this to halt at the first step such that the condition is satisfied. For example, we could use:condition(u,t,integrator) = u[2]>0\naffect!(integrator) = terminate!(integrator)\ncb = DiscreteCallback(condition,affect!)\nsol = solve(prob,Tsit5(),callback=cb)(Image: discrete_terminate)However, in many cases we wish to halt exactly at the point of time that the condition is satisfied. To do that, we use a continuous callback. The condition must thus be a function which is zero at the point we want to halt. Thus we use the following:condition(u,t,integrator) = u[2]\naffect!(integrator) = terminate!(integrator)\ncb = ContinuousCallback(condition,affect!)\nsol = solve(prob,Tsit5(),callback=cb)(Image: simple_terminate)Note that this uses rootfinding to approximate the \"exact\" moment of the crossing. Analytically we know the value is pi, and here the integration terminates atsol.t[end] # 3.1415902502224307Using a more accurate integration increases the accuracy of this prediction:sol = solve(prob,Vern8(),callback=cb,reltol=1e-12,abstol=1e-12)\nsol.t[end] # 3.1415926535896035\n#π = 3.141592653589703...Now say we wish to find the when the first period is over, i.e. we want to ignore the upcrossing and only stop on the downcrossing. We do this by ignoring the affect! and only passing an affect! for the second:condition(u,t,integrator) = u[2]\naffect!(integrator) = terminate!(integrator)\ncb = ContinuousCallback(condition,nothing,affect!)\nsol = solve(prob,Tsit5(),callback=cb)(Image: downcrossing_terminate)Notice that passing only one affect! is the same as ContinuousCallback(condition,affect!,affect!), i.e. both upcrossings and downcrossings will activate the event. Using ContinuousCallback(condition,affect!,nothing)will thus be the same as above because the first event is an upcrossing."
},

{
    "location": "features/callback_functions.html#Example-3:-Growing-Cell-Population-1",
    "page": "Event Handling and Callback Functions",
    "title": "Example 3: Growing Cell Population",
    "category": "section",
    "text": "Another interesting issue is with models of changing sizes. The ability to handle such events is a unique feature of DifferentialEquations.jl! The problem we would like to tackle here is a cell population. We start with 1 cell with a protein X which increases linearly with time with rate parameter α. Since we are going to be changing the size of the population, we write the model in the general form:const α = 0.3\nfunction f(du,u,p,t)\n  for i in 1:length(u)\n    du[i] = α*u[i]\n  end\nendOur model is that, whenever the protein X gets to a concentration of 1, it triggers a cell division. So we check to see if any concentrations hit 1:function condition(u,t,integrator) # Event when event_f(u,t) == 0\n  1-maximum(u)\nendAgain, recall that this function finds events as when condition==0, so 1-maximum(u) is positive until a cell has a concentration of X which is 1, which then triggers the event. At the event, we have that the cell splits into two cells, giving a random amount of protein to each one. We can do this by resizing the cache (adding 1 to the length of all of the caches) and setting the values of these two cells at the time of the event:function affect!(integrator)\n  u = integrator.u\n  resize!(integrator,length(u)+1)\n  maxidx = findmax(u)[2]\n  Θ = rand()\n  u[maxidx] = Θ\n  u[end] = 1-Θ\n  nothing\nendAs noted in the Integrator Interface, resize!(integrator,length(integrator.u)+1) is used to change the length of all of the internal caches (which includes u) to be their current length + 1, growing the ODE system. Then the following code sets the new protein concentrations. Now we can solve:callback = ContinuousCallback(condition,affect!)\nu0 = [0.2]\ntspan = (0.0,10.0)\nprob = ODEProblem(f,u0,tspan)\nsol = solve(prob,callback=callback)The plot recipes do not have a way of handling the changing size, but we can plot from the solution object directly. For example, let\'s make a plot of how many cells there are at each time. Since these are discrete values, we calculate and plot them directly:plot(sol.t,map((x)->length(x),sol[:]),lw=3,\n     ylabel=\"Number of Cells\",xlabel=\"Time\")(Image: NumberOfCells)Now let\'s check-in on a cell. We can still use the interpolation to get a nice plot of the concentration of cell 1 over time. This is done with the command:ts = range(0, stop=10, length=100)\nplot(ts,map((x)->x[1],sol.(ts)),lw=3,\n     ylabel=\"Amount of X in Cell 1\",xlabel=\"Time\")(Image: Cell1)Notice that every time it hits 1 the cell divides, giving cell 1 a random amount of X which then grows until the next division.Note that one macro which was not shown in this example is deleteat! on the caches. For example, to delete the second cell, we could use:deleteat!(integrator,2)This allows you to build sophisticated models of populations with births and deaths."
},

{
    "location": "features/callback_functions.html#VectorContinuousCallback-Example-1",
    "page": "Event Handling and Callback Functions",
    "title": "VectorContinuousCallback Example",
    "category": "section",
    "text": ""
},

{
    "location": "features/callback_functions.html#Example-1:-Bouncing-Ball-with-multiple-walls-1",
    "page": "Event Handling and Callback Functions",
    "title": "Example 1: Bouncing Ball with multiple walls",
    "category": "section",
    "text": "This is similar to the above Bouncing Ball example, but now we have two more vertical walls, at x = 0 and x = 10.0. We have our ODEFunction as - function f(du,u,p,t)\n  du[1] = u[2]\n  du[2] = -p\n  du[3] = u[4]\n  du[4] = 0.0\nendwhere u[1] denotes y-coordinate, u[2] denotes velocity in y-direction, u[3] denotes x-coordinate and u[4] denotes velocity in x-direction. We will make a VectorContinuousCallback of length 2 - one for x axis collision, one for walls parallel to y axis. function condition(out,u,t,integrator) # Event when event_f(u,t) == 0\n  out[1] = u[1]\n  out[2] = (u[3] - 10.0)u[3]\nend\n\nfunction affect!(integrator, idx)\n  if idx == 1\n    integrator.u[2] = -0.9integrator.u[2]\n  elseif idx == 2\n    integrator.u[4] = -0.9integrator.u[4]\n  end\nend\n\ncb = VectorContinuousCallback(condition,affect!,2)It is evident that out[2] will be zero when u[3] (x-coordinate) is either 0.0 or 10.0. And when that happens, we flip the velocity with some coefficient of restitution (0.9).Completeting rest of the code-u0 = [50.0,0.0,0.0,2.0]\ntspan = (0.0,15.0)\np = 9.8\nprob = ODEProblem(f,u0,tspan,p)\nsol = solve(prob,Tsit5(),callback=cb,dt=1e-3,adaptive=false)\nplot(sol,vars=(1,3))And you get the following output: (Image: Cell1)"
},

{
    "location": "features/callback_library.html#",
    "page": "Callback Library",
    "title": "Callback Library",
    "category": "page",
    "text": ""
},

{
    "location": "features/callback_library.html#Callback-Library-1",
    "page": "Callback Library",
    "title": "Callback Library",
    "category": "section",
    "text": "DiffEqCallbacks.jl provides a library of various helpful callbacks which can be used with any component solver which implements the callback interface. It adds the following callbacks which are available to users of DifferentialEquations.jl."
},

{
    "location": "features/callback_library.html#Manifold-Conservation-and-Projection-1",
    "page": "Callback Library",
    "title": "Manifold Conservation and Projection",
    "category": "section",
    "text": "In many cases, you may want to declare a manifold on which a solution lives. Mathematically, a manifold M is defined by a function g as the set of points where g(u)=0. An embedded manifold can be a lower dimensional object which constrains the solution. For example, g(u)=E(u)-C where E is the energy of the system in state u, meaning that the energy must be constant (energy preservation). Thus by defining the manifold the solution should live on, you can retain desired properties of the solution.It is a consequence of convergence proofs both in the deterministic and stochastic cases that post-step projection to manifolds keep the same convergence rate (stochastic requires a truncation in the proof, details details), thus any algorithm can be easily extended to conserve properties. If the solution is supposed to live on a specific manifold or conserve such property, this guarantees the conservation law without modifying the convergence properties."
},

{
    "location": "features/callback_library.html#Constructor-1",
    "page": "Callback Library",
    "title": "Constructor",
    "category": "section",
    "text": "ManifoldProjection(g; nlsolve=NLSOLVEJL_SETUP(), save=true, autonomous=numargs(g)==2, nlopts=Dict{Symbol,Any}())g: The residual function for the manifold. This is an inplace function of form g(u, resid) or g(t, u, resid) which writes to the residual the difference from the manifold components.\nnlsolve: A nonlinear solver as defined in the nlsolve format.\nsave: Whether to do the save after the callback is applied. Standard saving is unchanged.\nautonomous: Whether g is an autonomous function of the form g(u, resid).\nnlopts: Optional arguments to nonlinear solver which can be any of the NLsolve keywords."
},

{
    "location": "features/callback_library.html#Example-1",
    "page": "Callback Library",
    "title": "Example",
    "category": "section",
    "text": "Here we solve the harmonic oscillator:u0 = ones(2)\nfunction f(du,u,p,t)\n  du[1] = u[2]\n  du[2] = -u[1]\nend\nprob = ODEProblem(f,u0,(0.0,100.0))However, this problem is supposed to conserve energy, and thus we define our manifold to conserve the sum of squares:function g(resid,u,p,t)\n  resid[1] = u[2]^2 + u[1]^2 - 2\n  resid[2] = 0\nendTo build the callback, we just callcb = ManifoldProjection(g)Using this callback, the Runge-Kutta method Vern7 conserves energy. Note that the standard saving occurs after the step and before the callback, and thus we set  save_everystep=false to turn off all standard saving and let the callback save after the projection is applied.sol = solve(prob,Vern7(),save_everystep=false,callback=cb)\n@test sol[end][1]^2 + sol[end][2]^2 ≈ 2(Image: manifold_projection)"
},

{
    "location": "features/callback_library.html#Saveat-Warning-1",
    "page": "Callback Library",
    "title": "Saveat Warning",
    "category": "section",
    "text": "Note that the ManifoldProjection callback modifies the endpoints of the integration intervals and thus breaks assumptions of internal interpolations. Because of this, the values for given by saveat will not be order-matching. However, the interpolation error can be proportional to the change by the projection, so if the projection is making small changes then one is still safe. However, if there are large changes from each projection, you should consider only saving at stopping/projection times. To do this, set tstops to the same values as saveat. There is a performance hit by doing so because now the integrator is forced to stop at every saving point, but this is guerenteed to match the order of the integrator even with the ManifoldProjection."
},

{
    "location": "features/callback_library.html#AutoAbstol-1",
    "page": "Callback Library",
    "title": "AutoAbstol",
    "category": "section",
    "text": "Many problem solving environments such as MATLAB provide a way to automatically adapt the absolute tolerance to the problem. This helps the solvers automatically \"learn\" what appropriate limits are. Via the callback interface, DiffEqCallbacks.jl implements a callback AutoAbstol which has the same behavior as the MATLAB implementation, that is the absolute tolerance starts and at each iteration it is set to the maximum value that the state has thus far reached times the relative tolerance. If init_curmax is zero, then the initial value is determined by the abstol of the solver. Otherwise this is the initial value for the current maximum abstol.To generate the callback, use the constructor:AutoAbstol(save=true;init_curmax=0.0)"
},

{
    "location": "features/callback_library.html#PositiveDomain-1",
    "page": "Callback Library",
    "title": "PositiveDomain",
    "category": "section",
    "text": "Especially in biology and other natural sciences, a desired property of dynamical systems is the positive invariance of the positive cone, i.e. non-negativity of variables at time t_0 ensures their non-negativity at times t geq t_0 for which the solution is defined. However, even if a system satisfies this property mathematically it can be difficult for ODE solvers to ensure it numerically, as these MATLAB examples show.In order to deal with this problem one can specify isoutofdomain=(u,p,t) -> any(x -> x < 0, u) as additional solver option, which will reject any step that leads to non-negative values and reduce the next time step. However, since this approach only rejects steps and hence calculations might be repeated multiple times until a step is accepted, it can be computationally expensive.Another approach is taken by a PositiveDomain callback in DiffEqCallbacks.jl, which is inspired by Shampine\'s et al. paper about non-negative ODE solutions. It reduces the next step by a certain scale factor until the extrapolated value at the next time point is non-negative with a certain tolerance. Extrapolations are cheap to compute but might be inaccurate, so if a time step is changed it is additionally reduced by a safety factor of 0.9. Since extrapolated values are only non-negative up to a certain tolerance and in addition actual calculations might lead to negative values, also any negative values at the current time point are set to 0. Hence by this callback non-negative values at any time point are ensured in a computationally cheap way, but the quality of the solution depends on how accurately extrapolations approximate next time steps.Please note that the system should be defined also outside the positive domain, since even with these approaches negative variables might occur during the calculations. Moreover, one should follow Shampine\'s et. al. advice and set the derivative x_i of a negative component x_i to max 0 f_i(x t), where t denotes the current time point with state vector x and f_i is the i-th component of function f in an ODE system x = f(x t)."
},

{
    "location": "features/callback_library.html#Constructor-2",
    "page": "Callback Library",
    "title": "Constructor",
    "category": "section",
    "text": "PositiveDomain(u=nothing; save=true, abstol=nothing, scalefactor=nothing)u: A prototype of the state vector of the integrator. A copy of it is saved and extrapolated values are written to it. If it is not specified every application of the callback allocates a new copy of the state vector.\nsave: Whether to do the standard saving (applied after the callback).\nabstol: Tolerance up to which negative extrapolated values are accepted. Element-wise tolerances are allowed. If it is not specified every application of the callback uses the current absolute tolerances of the integrator.\nscalefactor: Factor by which an unaccepted time step is reduced. If it is not specified time steps are halved."
},

{
    "location": "features/callback_library.html#GeneralDomain-1",
    "page": "Callback Library",
    "title": "GeneralDomain",
    "category": "section",
    "text": "A GeneralDomain callback in DiffEqCallbacks.jl generalizes the concept of a PositiveDomain callback to arbitrary domains. Domains are specified by in-place functions g(u, resid) or g(t, u, resid) that calculate residuals of a state vector u at time t relative to that domain. As for PositiveDomain, steps are accepted if residuals of the extrapolated values at the next time step are below a certain tolerance. Moreover, this callback is automatically coupled with a ManifoldProjection that keeps all calculated state vectors close to the desired domain, but in contrast to a PositiveDomain callback the nonlinear solver in a ManifoldProjection can not guarantee that all state vectors of the solution are actually inside the domain. Thus a PositiveDomain callback should in general be preferred."
},

{
    "location": "features/callback_library.html#Constructor-3",
    "page": "Callback Library",
    "title": "Constructor",
    "category": "section",
    "text": "function GeneralDomain(g, u=nothing; nlsolve=NLSOLVEJL_SETUP(), save=true,\n                       abstol=nothing, scalefactor=nothing, autonomous=numargs(g)==2,\n                       nlopts=Dict(:ftol => 10*eps()))g: The residual function for the domain. This is an inplace function of form g(resid, u, p, t) which writes to the residual the difference from the domain.\nu: A prototype of the state vector of the integrator and the residuals. Two copies of it are saved, and extrapolated values and residuals are written to them. If it is not specified every application of the callback allocates two new copies of the state vector.\nnlsolve: A nonlinear solver as defined in the nlsolve format which is passed to a ManifoldProjection.\nsave: Whether to do the standard saving (applied after the callback).\nabstol: Tolerance up to which residuals are accepted. Element-wise tolerances are allowed. If it is not specified every application of the callback uses the current absolute tolerances of the integrator.\nscalefactor: Factor by which an unaccepted time step is reduced. If it is not specified time steps are halved.\nautonomous: Whether g is an autonomous function of the form g(u, resid).\nnlopts: Optional arguments to nonlinear solver of a ManifoldProjection which can be any of the NLsolve keywords. The default value of ftol = 10*eps() ensures that convergence is only declared if the infinite norm of residuals is very small and hence the state vector is very close to the domain."
},

{
    "location": "features/callback_library.html#Stepsize-Limiters-1",
    "page": "Callback Library",
    "title": "Stepsize Limiters",
    "category": "section",
    "text": "In many cases there is a known maximal stepsize for which the computation is stable and produces correct results. For example, in hyperbolic PDEs one normally needs to ensure that the stepsize stays below some Delta t_FE determined by the CFL condition. For nonlinear hyperbolic PDEs this limit can be a function dtFE(u,p,t) which changes throughout the computation. The stepsize limiter lets you pass a function which will adaptively limit the stepsizes to match these constraints."
},

{
    "location": "features/callback_library.html#Constructor-4",
    "page": "Callback Library",
    "title": "Constructor",
    "category": "section",
    "text": "StepsizeLimiter(dtFE;safety_factor=9//10,max_step=false,cached_dtcache=0.0)dtFE: The function for the maximal timestep, called as dtFE(u,p,t) using the previous values of u, p, and t.\nsafety_factor: The factor below the true maximum that will be stepped to which defaults to 9//10.\nmax_step: Makes every step equal to safety_factor*dtFE(u,p,t) when the solver is set to adaptive=false.\ncached_dtcache: Should be set to match the type for time when not using Float64 values."
},

{
    "location": "features/callback_library.html#FunctionCallingCallback-1",
    "page": "Callback Library",
    "title": "FunctionCallingCallback",
    "category": "section",
    "text": "The function calling callback lets you define a function func(u,t,integrator) which gets calls at the time points of interest. The constructor is:  FunctionCallingCallback(func;\n                 funcat=Vector{Float64}(),\n                 func_everystep=isempty(funcat),\n                 func_start = true,\n                 tdir=1)func(u, t, integrator) is the function to be called.\nfuncat values that the function is sure to be evaluated at.\nfunc_everystep whether to call the function after each integrator step.\nfunc_start whether the function is called the initial condition.\ntdir should be sign(tspan[end]-tspan[1]). It defaults to 1 and should   be adapted if tspan[1] > tspan[end]."
},

{
    "location": "features/callback_library.html#SavingCallback-1",
    "page": "Callback Library",
    "title": "SavingCallback",
    "category": "section",
    "text": "The saving callback lets you define a function save_func(u, t, integrator) which returns quantities of interest that shall be saved."
},

{
    "location": "features/callback_library.html#Constructor-5",
    "page": "Callback Library",
    "title": "Constructor",
    "category": "section",
    "text": "SavingCallback(save_func, saved_values::SavedValues;\n               saveat=Vector{eltype(saved_values.t)}(),\n               save_everystep=isempty(saveat),\n               tdir=1)save_func(u, t, integrator) returns the quantities which shall be saved. Note that this should allocate the output (not as a view to u).\nsaved_values::SavedValues is the types that save_func will return, i.e. save_func(u, t, integrator)::savevalType. It\'s specified via SavedValues(typeof(t),savevalType), i.e. give the type for time and the type that save_func will output (or higher compatible type).\nsaveat mimicks saveat in solve from solve.\nsave_everystep mimicks save_everystep from solve.\nsave_start mimicks save_start from solve.\ntdir should be sign(tspan[end]-tspan[1]). It defaults to 1 and should be adapted if tspan[1] > tspan[end].The outputted values are saved into saved_values. Time points are found via saved_values.t and the values are saved_values.saveval."
},

{
    "location": "features/callback_library.html#Example-2",
    "page": "Callback Library",
    "title": "Example",
    "category": "section",
    "text": "In this example we will solve a matrix equation and at each step save a tuple of values which contains the current trace and the norm of the matrix. We build the SavedValues cache to use Float64 for time and Tuple{Float64,Float64} for the saved values, and then call the solver with the callback.using DiffEqCallbacks, OrdinaryDiffEq, LinearAlgebra\nprob = ODEProblem((du,u,p,t) -> du .= u, rand(4,4), (0.0,1.0))\nsaved_values = SavedValues(Float64, Tuple{Float64,Float64})\ncb = SavingCallback((u,t,integrator)->(tr(u),norm(u)), saved_values)\nsol = solve(prob, Tsit5(), callback=cb)\n\nprint(saved_values.saveval)\n#=\nTuple{Float64,Float64}[(2.23186, 2.49102), (2.46675, 2.75318), (3.16138, 3.52847), (4.42011, 4.93337), (6.06683, 6.77129)]\n=#Note that the values are retrieved from the cache as .saveval, and the time points are found as .t. If we want to control the saved times, we use saveat in the callback. The save controls like saveat act analogously to how they act in the solve function.saved_values = SavedValues(Float64, Tuple{Float64,Float64})\ncb = SavingCallback((u,t,integrator)->(tr(u),norm(u)), saved_values, saveat=0.0:0.1:1.0)\nsol = solve(prob, Tsit5(), callback=cb)\nprint(saved_values.saveval)\nprint(saved_values.t)\n\n#=\nTuple{Float64,Float64}[(2.23186, 2.49102), (2.46659, 2.753), (2.726, 3.04254), (3.0127, 3.36253), \n(3.32955, 3.71617), (3.67972, 4.107), (4.06672, 4.53893), (4.49442, 5.0163), (4.9671, 5.54387), \n(5.48949, 6.12692), (6.06683, 6.77129)]\n[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]\n=#"
},

{
    "location": "features/callback_library.html#IterativeCallback-1",
    "page": "Callback Library",
    "title": "IterativeCallback",
    "category": "section",
    "text": "IterativeCallback is a callback to be used to iteratively apply some effect. For example, if given the first effect at t₁, you can define t₂ to apply the next effect.A IterativeCallback is constructed as follows:function IterativeCallback(time_choice, user_affect!,tType = Float64;\n                           initialize = DiffEqBase.INITIALIZE_DEFAULT,\n                           initial_affect = false, kwargs...)where time_choice(integrator) determines the time of the next callback and user_affect! is the effect applied to the integrator at the stopping points. If nothing is returned for the time choice then the iterator ends."
},

{
    "location": "features/callback_library.html#PeriodicCallback-1",
    "page": "Callback Library",
    "title": "PeriodicCallback",
    "category": "section",
    "text": "PeriodicCallback can be used when a function should be called periodically in terms of integration time (as opposed to wall time), i.e. at t = tspan[1], t = tspan[1] + Δt, t = tspan[1] + 2Δt, and so on. This callback can, for example, be used to model a digital controller for an analog system, running at a fixed rate."
},

{
    "location": "features/callback_library.html#Constructor-6",
    "page": "Callback Library",
    "title": "Constructor",
    "category": "section",
    "text": "PeriodicCallback(f, Δt::Number; kwargs...)where f is the function to be called periodically, Δt is the period, and kwargs are keyword arguments accepted by the DiscreteCallback constructor (see the DiscreteCallback section)."
},

{
    "location": "features/callback_library.html#TerminateSteadyState-1",
    "page": "Callback Library",
    "title": "TerminateSteadyState",
    "category": "section",
    "text": "TerminateSteadyState can be used to solve the problem for the steady-state by running the solver until the derivatives of the problem converge to 0 or tspan[2] is reached. This is an alternative approach to root finding (see the Steady State Solvers section). The constructor of this callback is:TerminateSteadyState(abstol = 1e-8, reltol = 1e-6, test = allDerivPass)where abstol and reltol are the absolute and relative tolerance, respectively. These tolerances may be specified as scalars or as arrays of the same length as the states of the problem. test represents the function that evaluates the condition for termination. The default condition is that all derivatives should become smaller than abstol and the states times reltol. The user can pass any other function to implement a different termination condition. Such function should take four arguments: integrator (see Integrator Interface for details), abstol and reltol."
},

{
    "location": "features/ensemble.html#",
    "page": "Parallel Ensemble Simulations",
    "title": "Parallel Ensemble Simulations",
    "category": "page",
    "text": ""
},

{
    "location": "features/ensemble.html#Parallel-Ensemble-Simulations-1",
    "page": "Parallel Ensemble Simulations",
    "title": "Parallel Ensemble Simulations",
    "category": "section",
    "text": "Performing Monte Carlo simulations, solving with a predetermined set of initial conditions, and GPU-parallelizing a parameter search all fall under the ensemble simulation interface. This interface allows one to declare a template DEProblem to parallelize, how to tweak the template in trajectories many trajectories, solve each in parallel batches, reduce the solutions down to specific answers, and compute summary statistics on the results."
},

{
    "location": "features/ensemble.html#Performing-an-Ensemble-Simulation-1",
    "page": "Parallel Ensemble Simulations",
    "title": "Performing an Ensemble Simulation",
    "category": "section",
    "text": ""
},

{
    "location": "features/ensemble.html#Building-a-Problem-1",
    "page": "Parallel Ensemble Simulations",
    "title": "Building a Problem",
    "category": "section",
    "text": "To perform a simulation on an ensemble of trajectories, define a EnsembleProblem. The constructor is:EnsembleProblem(prob::DEProblem;\n                output_func = (sol,i) -> (sol,false),\n                prob_func= (prob,i,repeat)->(prob),\n                reduction = (u,data,I)->(append!(u,data),false),\n                u_init = [])prob_func: The function by which the problem is to be modified. prob is the problem, i is the unique id 1:trajectories for the problem, and repeat is for if the iteration of the repeat. At first it\'s 0, but if rerun was true this will be 1, 2, etc. counting the number of times problem i has been repeated.\noutput_func: The function determines what is saved from the solution to the output array. Defaults to saving the solution itself. The output is (out,rerun) where out is the output and rerun is a boolean which designates whether to rerun\nreduction: This function determines how to reduce the data in each batch. Defaults to appending the data from the batches. The second part of the output determines whether the simulation has converged. If true, the simulation will exit early. By default, this is always false.One can specify a function prob_func which changes the problem. For example:function prob_func(prob,i,repeat)\n  @. prob.u0 = randn()*prob.u0\n  prob\nendmodifies the initial condition for all of the problems by a standard normal random number (a different random number per simulation). Notice that since problem types are immutable, it uses .=. Otherwise, one can just create a new problem type:function prob_func(prob,i,repeat)\n  @. prob.u0 = u0_arr[i]\n  prob\nendIf your function is a ParameterizedFunction, you can do similar modifications to prob.f to perform a parameter search. The output_func is a reduction function. It\'s arguments are the generated solution and the unique index for the run. For example, if we wish to only save the 2nd coordinate at the end of each solution, we can do:output_func(sol,i) = (sol[end,2],false)Thus the ensemble simulation would return as its data an array which is the end value of the 2nd dependent variable for each of the runs."
},

{
    "location": "features/ensemble.html#Solving-the-Problem-1",
    "page": "Parallel Ensemble Simulations",
    "title": "Solving the Problem",
    "category": "section",
    "text": "sim = solve(prob,alg,ensemblealg,kwargs...)The keyword arguments take in the arguments for the common solver interface and will pass them to the differential equation solver. The ensemblealg is optional, and will default to EnsembleThreads(). The special keyword arguments to note are:trajectories: The number of simulations to run. This argument is required.\nbatch_size : The size of the batches on which the reductions are applies. Defaults to trajectories.\npmap_batch_size: The size of the pmap batches. Default is  batch_size÷100 > 0 ? batch_size÷100 : 1"
},

{
    "location": "features/ensemble.html#EnsembleAlgorithms-1",
    "page": "Parallel Ensemble Simulations",
    "title": "EnsembleAlgorithms",
    "category": "section",
    "text": "The choice of ensemble algorithm allows for control over how the multiple trajectories are handled. Currently, the ensemble algorithm types are:EnsembleSerial() - No parallelism\nEnsembleThreads() - This uses multithreading. It\'s local (single computer, shared memory) parallelism only. Fastest when the trajectories are quick.\nEnsembleDistributed() - The default. Uses pmap internally. It will use as many processors as you have Julia processes. To add more processes, use addprocs(n). See Julia\'s documentation for more details. Recommended for the case when each trajectory calculation isn\'t \"too quick\" (at least about a millisecond each?).\nEnsembleSplitThreads() - This uses threading on each process, splitting the problem into nprocs() even parts. This is for solving many quick trajectories on a multi-node machine. It\'s recommended you have one process on each node.For example, EnsembleThreads() is invoked by:solve(ensembleprob,alg,EnsembleThreads();trajectories=1000)Additionally, an experimental GPU-based ensembling method is provided by DiffEqGPU.jl. Add and import that package to get EnsembleGPUArray()."
},

{
    "location": "features/ensemble.html#Solution-Type-1",
    "page": "Parallel Ensemble Simulations",
    "title": "Solution Type",
    "category": "section",
    "text": "The resulting type is a EnsembleSimulation, which includes the array of solutions."
},

{
    "location": "features/ensemble.html#Plot-Recipe-1",
    "page": "Parallel Ensemble Simulations",
    "title": "Plot Recipe",
    "category": "section",
    "text": "There is a plot recipe for a AbstractEnsembleSimulation which composes all of the plot recipes for the component solutions. The keyword arguments are passed along. A useful argument to use is linealpha which will change the transparency of the plots. An additional argument is idxs which allows you to choose which components of the solution to plot. For example, if the differential equation is a vector of 9 values, idxs=1:2:9 will plot only the solutions of the odd components. An other additional argument is zcolors which allows you to pass a zcolor for each series. For details about zcolor see the  documentation for Plots.jl."
},

{
    "location": "features/ensemble.html#Analyzing-an-Ensemble-Experiment-1",
    "page": "Parallel Ensemble Simulations",
    "title": "Analyzing an Ensemble Experiment",
    "category": "section",
    "text": "Analysis tools are included for generating summary statistics and summary plots for a EnsembleSimulation."
},

{
    "location": "features/ensemble.html#Time-steps-vs-time-points-1",
    "page": "Parallel Ensemble Simulations",
    "title": "Time steps vs time points",
    "category": "section",
    "text": "For the summary statistics, there are two types. You can either summarize by time steps or by time points. Summarizing by time steps assumes that the time steps are all the same time point, i.e. the integrator used a fixed dt or the values were saved using saveat. Summarizing by time points requires interpolating the solution."
},

{
    "location": "features/ensemble.html#Analysis-at-a-time-step-or-time-point-1",
    "page": "Parallel Ensemble Simulations",
    "title": "Analysis at a time step or time point",
    "category": "section",
    "text": "get_timestep(sim,i) # Returns an iterator of each simulation at time step i\nget_timepoint(sim,t) # Returns an iterator of each simulation at time point t\ncomponentwise_vectors_timestep(sim,i) # Returns a vector of each simulation at time step i\ncomponentwise_vectors_timepoint(sim,t) # Returns a vector of each simulation at time point t"
},

{
    "location": "features/ensemble.html#Summary-Statistics-1",
    "page": "Parallel Ensemble Simulations",
    "title": "Summary Statistics",
    "category": "section",
    "text": ""
},

{
    "location": "features/ensemble.html#Single-Time-Statistics-1",
    "page": "Parallel Ensemble Simulations",
    "title": "Single Time Statistics",
    "category": "section",
    "text": "The available functions for time steps are:timestep_mean(sim,i) # Computes the mean of each component at time step i\ntimestep_median(sim,i) # Computes the median of each component at time step i\ntimestep_quantile(sim,q,i) # Computes the quantile q of each component at time step i\ntimestep_meanvar(sim,i)  # Computes the mean and variance of each component at time step i\ntimestep_meancov(sim,i,j) # Computes the mean at i and j, and the covariance, for each component\ntimestep_meancor(sim,i,j) # Computes the mean at i and j, and the correlation, for each component\ntimestep_weighted_meancov(sim,W,i,j) # Computes the mean at i and j, and the weighted covariance W, for each componentThe available functions for time points are:timepoint_mean(sim,t) # Computes the mean of each component at time t\ntimepoint_median(sim,t) # Computes the median of each component at time t\ntimepoint_quantile(sim,q,t) # Computes the quantile q of each component at time t\ntimepoint_meanvar(sim,t) # Computes the mean and variance of each component at time t\ntimepoint_meancov(sim,t1,t2) # Computes the mean at t1 and t2, the covariance, for each component\ntimepoint_meancor(sim,t1,t2) # Computes the mean at t1 and t2, the correlation, for each component\ntimepoint_weighted_meancov(sim,W,t1,t2) # Computes the mean at t1 and t2, the weighted covariance W, for each component"
},

{
    "location": "features/ensemble.html#Full-Timeseries-Statistics-1",
    "page": "Parallel Ensemble Simulations",
    "title": "Full Timeseries Statistics",
    "category": "section",
    "text": "Additionally, the following functions are provided for analyzing the full timeseries. The mean and meanvar versions return a DiffEqArray which can be directly plotted. The meancov and meancor return a matrix of tuples, where the tuples are the (mean_t1,mean_t2,cov or cor).The available functions for the time steps are:timeseries_steps_mean(sim) # Computes the mean at each time step\ntimeseries_steps_median(sim) # Computes the median at each time step\ntimeseries_steps_quantile(sim,q) # Computes the quantile q at each time step\ntimeseries_steps_meanvar(sim) # Computes the mean and variance at each time step\ntimeseries_steps_meancov(sim) # Computes the covariance matrix and means at each time step\ntimeseries_steps_meancor(sim) # Computes the correlation matrix and means at each time step\ntimeseries_steps_weighted_meancov(sim) # Computes the weighted covariance matrix and means at each time stepThe available functions for the time points are:timeseries_point_mean(sim,ts) # Computes the mean at each time point in ts\ntimeseries_point_median(sim,ts) # Computes the median at each time point in ts\ntimeseries_point_quantile(sim,q,ts) # Computes the quantile q at each time point in ts\ntimeseries_point_meanvar(sim,ts) # Computes the mean and variance at each time point in ts\ntimeseries_point_meancov(sim,ts) # Computes the covariance matrix and means at each time point in ts\ntimeseries_point_meancor(sim,ts) # Computes the correlation matrix and means at each time point in ts\ntimeseries_point_weighted_meancov(sim,ts) # Computes the weighted covariance matrix and means at each time point in ts"
},

{
    "location": "features/ensemble.html#EnsembleSummary-1",
    "page": "Parallel Ensemble Simulations",
    "title": "EnsembleSummary",
    "category": "section",
    "text": "The EnsembleSummary type is included to help with analyzing the general summary statistics. Two constructors are provided:EnsembleSummary(sim;quantile=[0.05,0.95])\nEnsembleSummary(sim,ts;quantile=[0.05,0.95])The first produces a (mean,var) summary at each time step. As with the summary statistics, this assumes that the time steps are all the same. The second produces a (mean,var) summary at each time point t in ts. This requires the ability to interpolate the solution. Quantile is used to determine the qlow and qhigh quantiles at each timepoint. It defaults to the 5% and 95% quantiles."
},

{
    "location": "features/ensemble.html#Plot-Recipe-2",
    "page": "Parallel Ensemble Simulations",
    "title": "Plot Recipe",
    "category": "section",
    "text": "The EnsembleSummary comes with a plot recipe for visualizing the summary statistics. The extra keyword arguments are:idxs: the solution components to plot. Defaults to plotting all components.\nerror_style: The style for plotting the error. Defaults to ribbon. Other choices are :bars for error bars and :none for no error bars.\nci_type : Defaults to :quantile which has (qlow,qhigh) quantiles whose limits were determined when constructing the EnsembleSummary. Gaussian CI 1.96*(standard error of the mean) can be set using ci_type=:SEM.One useful argument is fillalpha which controls the transparency of the ribbon around the mean."
},

{
    "location": "features/ensemble.html#Example-1:-Solving-an-ODE-With-Different-Initial-Conditions-1",
    "page": "Parallel Ensemble Simulations",
    "title": "Example 1: Solving an ODE With Different Initial Conditions",
    "category": "section",
    "text": ""
},

{
    "location": "features/ensemble.html#Random-Initial-Conditions-1",
    "page": "Parallel Ensemble Simulations",
    "title": "Random Initial Conditions",
    "category": "section",
    "text": "Let\'s test the sensitivity of the linear ODE to its initial condition. To do this, we would like to solve the linear ODE 100 times and plot what the trajectories look like. Let\'s start by opening up some extra processes so that way the computation will be parallelized. This will use pmap as default, which means that the required functions must be made available to all processes. This can be achieved with @everywhere macro:addprocs()\n@everywhere using DifferentialEquationsNow let\'s define the linear ODE which is our base problem:# Linear ODE which starts at 0.5 and solves from t=0.0 to t=1.0\nprob = ODEProblem((u,p,t)->1.01u,0.5,(0.0,1.0))For our ensemble simulation, we would like to change the initial condition around. This is done through the prob_func. This function takes in the base problem and modifies it to create the new problem that the trajectory actually solves. Here we will take the base problem, multiply the initial condition by a rand(), and use that for calculating the trajectory:@everywhere function prob_func(prob,i,repeat)\n  ODEProblem(prob.f,rand()*prob.u0,prob.tspan)\nendNow we build and solve the EnsembleProblem with this base problem and prob_func:ensemble_prob = EnsembleProblem(prob,prob_func=prob_func)\nsim = solve(ensemble_prob,Tsit5(),trajectories=100)We can use the plot recipe to plot what the 100 ODEs look like:using Plots\nplotly()\nplot(sim,linealpha=0.4)(Image: monte_carlo_plot)We note that if we wanted to find out what the initial condition was for a given trajectory, we can retrieve it from the solution. sim[i] returns the ith solution object. sim[i].prob is the problem that specific trajectory solved, and sim[i].prob.u0 would then be the initial condition used in the ith trajectory.Note: If the problem has callbacks, the functions for the condition and affect! must be named functions (not anonymous functions)."
},

{
    "location": "features/ensemble.html#Using-multithreading-1",
    "page": "Parallel Ensemble Simulations",
    "title": "Using multithreading",
    "category": "section",
    "text": "The previous ensemble simulation can also be parallelized using a multithreading approach, which will make use of the different cores within a single computer. Because the memory is shared across the different threads, it is not necessary to use the @everywhere macro. Instead, the same problem can be implemented simply as:using DifferentialEquations\nprob = ODEProblem((u,p,t)->1.01u,0.5,(0.0,1.0))\nfunction prob_func(prob,i,repeat)\n  ODEProblem(prob.f,rand()*prob.u0,prob.tspan)\nend\nensemble_prob = EnsembleProblem(prob,prob_func=prob_func)\nsim = solve(ensemble_prob,Tsit5(),EnsembleThreads(),trajectories=100)The number of threads to be used has to be defined outside of Julia, in the environmental variable JULIA_NUM_THREADS (see Julia\'s documentation for details)."
},

{
    "location": "features/ensemble.html#Pre-Determined-Initial-Conditions-1",
    "page": "Parallel Ensemble Simulations",
    "title": "Pre-Determined Initial Conditions",
    "category": "section",
    "text": "In many cases, you may already know what initial conditions you want to use. This can be specified by the i argument of the prob_func. This i is the unique index of each trajectory. So, if we have trajectories=100, then we have i as some index in 1:100, and it\'s different for each trajectory.So, if we wanted to use a grid of evenly spaced initial conditions from 0 to 1, we could simply index the linspace type:initial_conditions = range(0, stop=1, length=100)\nfunction prob_func(prob,i,repeat)\n  prob.u0 = initial_conditions[i]\n  prob\nend"
},

{
    "location": "features/ensemble.html#Example-2:-Solving-an-SDE-with-Different-Parameters-1",
    "page": "Parallel Ensemble Simulations",
    "title": "Example 2: Solving an SDE with Different Parameters",
    "category": "section",
    "text": "Let\'s solve the same SDE but with varying parameters. Let\'s create a Lotka-Volterra  system with multiplicative noise. Our Lotka-Volterra system will have as its  drift component:function f(du,u,p,t)\n  du[1] = p[1] * u[1] - p[2] * u[1]*u[2]\n  du[2] = -3 * u[2] + u[1]*u[2]\nendFor our noise function we will use multiplicative noise:function g(du,u,p,t)\n  du[1] = p[3]*u[1]\n  du[2] = p[4]*u[2]\nendNow we build the SDE with these functions:p = [1.5,1.0,0.1,0.1]\nprob = SDEProblem(f,g,[1.0,1.0],(0.0,10.0),p)This is the base problem for our study. What would like to do with this experiment is keep the same parameters in the deterministic component each time, but very the parameters for the amount of noise using 0.3rand(2) as our parameters. Once again, we do this with a prob_func, and here we modify the parameters in prob.p:function prob_func(prob,i,repeat)\n  prob.p[3:4] = 0.3rand(2)\n  prob\nendNow we solve the problem 10 times and plot all of the trajectories in phase space:ensemble_prob = EnsembleProblem(prob,prob_func=prob_func)\nsim = solve(ensemble_prob,SRIW1(),trajectories=10)\nusing Plots; plotly()\nusing Plots; plot(sim,linealpha=0.6,color=:blue,vars=(0,1),title=\"Phase Space Plot\")\nplot!(sim,linealpha=0.6,color=:red,vars=(0,2),title=\"Phase Space Plot\")(Image: monte_lotka_blue)We can then summarize this information with the mean/variance bounds using a EnsembleSummary plot. We will take the mean/quantile at every 0.1 time units and directly plot the summary:summ = EnsembleSummary(sim,0:0.1:10)\npyplot() # Note that plotly does not support ribbon plots\nplot(summ,fillalpha=0.5)(Image: monte_carlo_quantile)Note that here we used the quantile bounds, which default to [0.05,0.95] in the EnsembleSummary constructor. We can change to standard error of the mean bounds using ci_type=:SEM in the plot recipe."
},

{
    "location": "features/ensemble.html#Example-3:-Using-the-Reduction-to-Halt-When-Estimator-is-Within-Tolerance-1",
    "page": "Parallel Ensemble Simulations",
    "title": "Example 3: Using the Reduction to Halt When Estimator is Within Tolerance",
    "category": "section",
    "text": "In this problem we will solve the equation just as many times as needed to get the standard error of the mean for the final time point below our tolerance 0.5. Since we only care about the endpoint, we can tell the output_func to discard the rest of the data.function output_func(sol,i)\n  last(sol)\nendOur prob_func will simply randomize the initial condition:# Linear ODE which starts at 0.5 and solves from t=0.0 to t=1.0\nprob = ODEProblem((u,p,t)->1.01u,0.5,(0.0,1.0))\n\nfunction prob_func(prob,i,repeat)\n  ODEProblem(prob.f,rand()*prob.u0,prob.tspan)\nendOur reduction function will append the data from the current batch to the previous batch, and declare convergence if the standard error of the mean is calculated as sufficiently small:function reduction(u,batch,I)\n  u = append!(u,batch)\n  u,((var(u)/sqrt(last(I)))/mean(u)<0.5) ? true : false\nendThen we can define and solve the problem:prob2 = EnsembleProblem(prob,prob_func=prob_func,output_func=output_func,reduction=reduction,u_init=Vector{Float64}())\nsim = solve(prob2,Tsit5(),trajectories=10000,batch_size=20)Since batch_size=20, this means that every 20 simulations, it will take this batch, append the results to the previous batch, calculate (var(u)/sqrt(last(I)))/mean(u), and if that\'s small enough, exit the simulation. In this case, the simulation exits only after 20 simulations (i.e. after calculating the first batch). This can save a lot of time!In addition to saving time by checking convergence, we can save memory by reducing between batches. For example, say we only care about the mean at the end once again. Instead of saving the solution at the end for each trajectory, we can instead save the running summation of the endpoints:function reduction(u,batch,I)\n  u+sum(batch),false\nend\nprob2 = EnsembleProblem(prob,prob_func=prob_func,output_func=output_func,reduction=reduction,u_init=0.0)\nsim2 = solve(prob2,Tsit5(),trajectories=100,batch_size=20)this will sum up the endpoints after every 20 solutions, and save the running sum. The final result will have sim2.u as simply a number, and thus sim2.u/100 would be the mean."
},

{
    "location": "features/ensemble.html#Example-4:-Using-the-Analysis-Tools-1",
    "page": "Parallel Ensemble Simulations",
    "title": "Example 4: Using the Analysis Tools",
    "category": "section",
    "text": "In this example we will show how to analyze a EnsembleSolution. First, let\'s generate a 10 solution Monte Carlo experiment. For our problem we will use a 4x2 system of linear stochastic differential equations:function f(du,u,p,t)\n  for i = 1:length(u)\n    du[i] = 1.01*u[i]\n  end\nend\nfunction σ(du,u,p,t)\n  for i in 1:length(u)\n    du[i] = .87*u[i]\n  end\nend\nprob = SDEProblem(f,σ,ones(4,2)/2,(0.0,1.0)) #prob_sde_2DlinearTo solve this 10 times, we use the EnsembleProblem constructor and solve with trajectories=10. Since we wish to compare values at the timesteps, we need to make sure the steps all hit the same times. Thus we set adaptive=false and explicitly give a dt.prob2 = EnsembleProblem(prob)\nsim = solve(prob2,SRIW1(),dt=1//2^(3),trajectories=10,adaptive=false)Note that if you don\'t do the timeseries_steps calculations, this code is compatible with adaptive timestepping. Using adaptivity is usually more efficient!We can compute the mean and the variance at the 3rd timestep using:m,v = timestep_meanvar(sim,3)or we can compute the mean and the variance at the t=0.5 using:m,v = timepoint_meanvar(sim,0.5)We can get a series for the mean and the variance at each time step using:m_series,v_series = timeseries_steps_meanvar(sim)or at chosen values of t:ts = 0:0.1:1\nm_series = timeseries_point_mean(sim,ts)Note that these mean and variance series can be directly plotted. We can compute covariance matrices similarly:timeseries_steps_meancov(sim) # Use the time steps, assume fixed dt\ntimeseries_point_meancov(sim,0:1//2^(3):1,0:1//2^(3):1) # Use time points, interpolateFor general analysis, we can build a EnsembleSummary type.summ = EnsembleSummary(sim)will summarize at each time step, whilesumm = EnsembleSummary(sim,0.0:0.1:1.0)will summarize at the 0.1 time points using the interpolations. To visualize the results we can plot it. Since there are 8 components to the differential equation, this can get messy, so let\'s only plot the 3rd component:plot(summ;idxs=3)(Image: monte_ribbon)We can change to errorbars instead of ribbons and plot two different indices:plot(summ;idxs=(3,5),error_style=:bars)(Image: monte_bars)Or we can simply plot the mean of every component over time:plot(summ;error_style=:none)(Image: monte_means)"
},

{
    "location": "features/io.html#",
    "page": "I/O: Saving and Loading Solution Data",
    "title": "I/O: Saving and Loading Solution Data",
    "category": "page",
    "text": ""
},

{
    "location": "features/io.html#I/O:-Saving-and-Loading-Solution-Data-1",
    "page": "I/O: Saving and Loading Solution Data",
    "title": "I/O: Saving and Loading Solution Data",
    "category": "section",
    "text": "The ability to save and load solutions is important for handling large datasets and analyzing the results over multiple Julia sessions. This page explains the existing functionality for doing so."
},

{
    "location": "features/io.html#Tabular-Data:-IterableTables-1",
    "page": "I/O: Saving and Loading Solution Data",
    "title": "Tabular Data: IterableTables",
    "category": "section",
    "text": "An interface to IterableTables.jl is provided. This IterableTables link allows you to use a solution type as the data source to convert to other tabular data formats. For example, let\'s solve a 4x2 system of ODEs:f_2dlinear = (du,u,p,t) -> du.=1.01u\nprob = ODEProblem(f_2dlinear,rand(2,2),(0.0,1.0))\nsol1 =solve(prob,Euler();dt=1//2^(4))then we can convert this to a dataframe using DataFrame:using IterableTables, DataFrames\ndf = DataFrame(sol1)\n\n# Result\n17×5 DataFrames.DataFrame\n│ Row │ timestamp │ value 1  │ value 2  │ value 3  │ value 4  │\n├─────┼───────────┼──────────┼──────────┼──────────┼──────────┤\n│ 1   │ 0.0       │ 0.110435 │ 0.569561 │ 0.918336 │ 0.508044 │\n│ 2   │ 0.0625    │ 0.117406 │ 0.605515 │ 0.976306 │ 0.540114 │\n│ 3   │ 0.125     │ 0.124817 │ 0.643738 │ 1.03794  │ 0.574208 │\n│ 4   │ 0.1875    │ 0.132696 │ 0.684374 │ 1.10345  │ 0.610455 │\n│ 5   │ 0.25      │ 0.141073 │ 0.727575 │ 1.17311  │ 0.64899  │\n│ 6   │ 0.3125    │ 0.149978 │ 0.773503 │ 1.24716  │ 0.689958 │\n│ 7   │ 0.375     │ 0.159445 │ 0.822331 │ 1.32589  │ 0.733511 │\n│ 8   │ 0.4375    │ 0.16951  │ 0.87424  │ 1.40959  │ 0.779814 │\n│ 9   │ 0.5       │ 0.18021  │ 0.929427 │ 1.49857  │ 0.82904  │\n│ 10  │ 0.5625    │ 0.191586 │ 0.988097 │ 1.59316  │ 0.881373 │\n│ 11  │ 0.625     │ 0.20368  │ 1.05047  │ 1.69373  │ 0.93701  │\n│ 12  │ 0.6875    │ 0.216537 │ 1.11678  │ 1.80065  │ 0.996159 │\n│ 13  │ 0.75      │ 0.230206 │ 1.18728  │ 1.91432  │ 1.05904  │\n│ 14  │ 0.8125    │ 0.244738 │ 1.26222  │ 2.03516  │ 1.12589  │\n│ 15  │ 0.875     │ 0.260187 │ 1.3419   │ 2.16363  │ 1.19697  │\n│ 16  │ 0.9375    │ 0.276611 │ 1.42661  │ 2.30021  │ 1.27252  │\n│ 17  │ 1.0       │ 0.294072 │ 1.51667  │ 2.44541  │ 1.35285  │If a ParameterizedFunction is used, the output will use the variable names:using ParameterizedFunctions\n\nf = @ode_def begin\n  dx = a*x - b*x*y\n  dy = -3y + x*y\nend a b\n\nprob = ODEProblem(f,[1.0,1.0],(0.0,1.0),[1.5,1.0])\nsol2 =solve(prob,Tsit5())\n\ndf = DataFrame(sol2)\n\n7×3 DataFrames.DataFrame\n│ Row │ timestamp │ x       │ y        │\n├─────┼───────────┼─────────┼──────────┤\n│ 1   │ 0.0       │ 1.0     │ 1.0      │\n│ 2   │ 0.0776085 │ 1.04549 │ 0.857668 │\n│ 3   │ 0.232645  │ 1.17587 │ 0.63946  │\n│ 4   │ 0.429118  │ 1.41968 │ 0.456996 │\n│ 5   │ 0.679082  │ 1.87672 │ 0.324733 │\n│ 6   │ 0.944406  │ 2.58825 │ 0.263362 │\n│ 7   │ 1.0       │ 2.77285 │ 0.25871  │Additionally, this data can be saved to a CSV:using CSV\nCSV.write(\"out.csv\",df)For more information on using the IterableTables interface and other output formats, see IterableTables.jl."
},

{
    "location": "features/io.html#JLD2-and-BSON.jl-1",
    "page": "I/O: Saving and Loading Solution Data",
    "title": "JLD2 and BSON.jl",
    "category": "section",
    "text": "JLD2.jl and BSON.jl will work with the full solution type if you bring the required functions back into scope before loading. For eaxmple, if we save the solution:using OrdinaryDiffEq, JLD2\nf(u,p,t) = 1.01*u\nu0=1/2\ntspan = (0.0,1.0)\nprob = ODEProblem(f,u0,tspan)\nsol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)\n@save \"out.jld2\" solthen we can get the full solution type back, interpolations and all, if we load the dependent functions first:using JLD2\nusing OrdinaryDiffEq\nf(u,p,t) = 1.01*u\nJLD2.@load \"out.jld2\" solThe example with BSON.jl is:using OrdinaryDiffEq\nf_2dlinear = (du,u,p,t) -> du.=1.01u\nprob = ODEProblem(f_2dlinear,rand(2,2),(0.0,1.0))\nsol1 =solve(prob,Euler();dt=1//2^(4))\n\nusing BSON\nbson(\"test.bson\",Dict(:sol1=>sol1))\n\n# New session\nusing OrdinaryDiffEq\nusing BSON\nBSON.load(\"test.bson\")If you load it without the DE function then for some algorithms the interpolation may not work, and for all algorithms you\'ll need at least a solver package or DiffEqBase.jl in scope in order for the solution interface (plot recipes, array indexing, etc.) to work. If none of these are put into scope, the solution type will still load and hold all of the values (so sol.u and sol.t will work), but none of the interface will be available."
},

{
    "location": "features/io.html#JLD-1",
    "page": "I/O: Saving and Loading Solution Data",
    "title": "JLD",
    "category": "section",
    "text": "Don\'t use JLD. It\'s dead. Julia types can be saved via JLD.jl. However, they cannot save types which have functions, which means that the solution type is currently not compatible with JLD.using JLD\nJLD.save(\"out.jld\",\"sol\",sol)"
},

{
    "location": "features/low_dep.html#",
    "page": "Low Dependency Usage",
    "title": "Low Dependency Usage",
    "category": "page",
    "text": ""
},

{
    "location": "features/low_dep.html#Low-Dependency-Usage-1",
    "page": "Low Dependency Usage",
    "title": "Low Dependency Usage",
    "category": "section",
    "text": "DifferentialEquations.jl is a large library containing the functionality of many different solver and addon packages. However in many cases you may want to cut down on the size of the dependency and only use the parts of the the library which are essential to your application. This is possible due to JuliaDiffEq\'s modular package structure."
},

{
    "location": "features/low_dep.html#Common-Example:-Using-only-OrdinaryDiffEq.jl-1",
    "page": "Low Dependency Usage",
    "title": "Common Example: Using only OrdinaryDiffEq.jl",
    "category": "section",
    "text": "One common example is using only the ODE solvers OrdinaryDiffEq.jl. The solvers all reexport DiffEqBase.jl (which holds the problem and solution types) and so OrdinaryDiffEq.jl is all that\'s needed. Thus replacingusing DifferentialEquationswithusing OrdinaryDiffEqwill work if these are the only features you are using."
},

{
    "location": "features/low_dep.html#Generalizing-the-Idea-1",
    "page": "Low Dependency Usage",
    "title": "Generalizing the Idea",
    "category": "section",
    "text": "In general, you will always need DiffEqBase.jl, since it defines all of the fundamental types, but the solvers will automatically reexport it. For solvers, you typically only need that solver package. So DiffEqBase+Sundials, DiffEqBase+LSODA, etc. will get you the common interface with that specific solver setup. DiffEqBase.jl is a very lightweight dependency, so there is no issue here! For PDEs, you normally need DiffEqBase+DiffEqPDEBase in addition to the solver package.For the addon packages, you will normally need DiffEqBase, the solver package you choose, and the addon package. So for example, for parameter estimation you would likely want DiffEqBase+OrdinaryDiffEq+DiffEqParamEstim. If you arne\'t sure which package a specific command is from, they using @which. For example, from the parameter estimation docs we have:using DifferentialEquations\nfunction f(du,u,p,t)\n  dx = p[1]*u[1] - u[1]*u[2]\n  dy = -3*u[2] + u[1]*u[2]\nend\n\nu0 = [1.0;1.0]\ntspan = (0.0,10.0)\np = [1.5]\nprob = ODEProblem(f,u0,tspan,p)\nsol = solve(prob,Tsit5())\nt = collect(range(0, stop=10, length=200))\nrandomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])\nusing RecursiveArrayTools\ndata = convert(Array,randomized)\ncost_function = build_loss_objective(prob,t,data,Tsit5(),maxiters=10000)If we wanted to know where build_loss_objective came from, we can do:@which build_loss_objective(prob,t,data,Tsit5(),maxiters=10000)\n\n(::DiffEqParamEstim.#kw##build_loss_objective)(::Array{Any,1}, ::DiffEqParamEstim.#build_loss_objective, prob::DiffEqBase.DEProblem, t, data, alg)This says it\'s in the DiffEqParamEstim.jl package. Thus in this case, we could have doneusing OrdinaryDiffEq, DiffEqParamEstiminstead of the full using DifferentialEquations. Note that due to the way Julia dependencies work, any internal function in the package will work. The only dependencies you need to explicitly using are the functions you are specifically calling. Thus this method can be used to determine all of the DiffEq packages you are using."
},

{
    "location": "features/progress_bar.html#",
    "page": "Juno Progress Bar Integration",
    "title": "Juno Progress Bar Integration",
    "category": "page",
    "text": ""
},

{
    "location": "features/progress_bar.html#Juno-Progress-Bar-Integration-1",
    "page": "Juno Progress Bar Integration",
    "title": "Juno Progress Bar Integration",
    "category": "section",
    "text": "DifferentialEquations.jl integrates with the Juno progress bar in order to make long calculations more manageable. By default this feature is off for ODE and SDE solvers, but can be turned on via the keyword argument progressbar=true. The progress bar updates every progress_steps timesteps, which has a default value of 1000. Note that making this value really low could cause a performance hit, though from some basic testing it seems that with updates of at least 1000 steps on number (the fastest problems) there\'s no discernable performance degradation, giving a high upper bound.Note that the progressbar also includes a time estimate. This time-estimate is provided by linear extrapolation for how long it has taken to get to what percentage. For adaptive timestepping methods this should only be used as a rough estimate since the timesteps may (and will) change. By scrolling over the progressbar one will also see the current timestep. This can be used to track the solution\'s progress and find tough locations for the solvers."
},

{
    "location": "analysis/parameterized_functions.html#",
    "page": "ParameterizedFunctions",
    "title": "ParameterizedFunctions",
    "category": "page",
    "text": ""
},

{
    "location": "analysis/parameterized_functions.html#ParameterizedFunctions-1",
    "page": "ParameterizedFunctions",
    "title": "ParameterizedFunctions",
    "category": "section",
    "text": ""
},

{
    "location": "analysis/parameterized_functions.html#Installation-1",
    "page": "ParameterizedFunctions",
    "title": "Installation",
    "category": "section",
    "text": "This functionality does not come standard with DifferentialEquations.jl. To use this functionality, you must install ParameterizedFunctions.jl:]add ParameterizedFunctions\nusing ParameterizedFunctions"
},

{
    "location": "analysis/parameterized_functions.html#Function-Definition-Macros-1",
    "page": "ParameterizedFunctions",
    "title": "Function Definition Macros",
    "category": "section",
    "text": "DifferentialEquations.jl provides a set of macros for more easily and legibly defining your differential equations. It exploits the standard notation for mathematically writing differential equations and the notation for \"punching differential equations into the computer\"; effectively doing the translation step for you. This is best shown by an example. Say we want to solve the ROBER model. Using the @ode_def macro from ParameterizedFunctions.jl, we can do this by writing:using ParameterizedFunctions\nf = @ode_def begin\n  dy₁ = -k₁*y₁+k₃*y₂*y₃\n  dy₂ =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃\n  dy₃ =  k₂*y₂^2\nend k₁ k₂ k₃This looks just like pseudocode! The macro will expand this to the \"standard form\", i.e. the ugly computer form:function f(du,u,p,t)\n  du[1] = -p[1]*u[1] + p[3]*u[2]*u[3]\n  du[2] = p[1]*u[1] - p[2]*u[2]^2 - p[3]*u[2]*u[3]\n  du[3] = p[2]*u[2]^2\nendNote that one doesn\'t need to use numbered variables: DifferentialEquations.jl will number the variables for you. For example, the following defines the function for the Lotka-Volterra model, with full Unicode support to boot:f = @ode_def begin\n  d🐁  = α*🐁  - β*🐁*🐈\n  d🐈 = -γ*🐈 + δ*🐁*🐈\nend α β γ δ"
},

{
    "location": "analysis/parameterized_functions.html#Limitations-1",
    "page": "ParameterizedFunctions",
    "title": "Limitations",
    "category": "section",
    "text": "The macro is a Domain-Specific Language (DSL) and thus has different internal semantics than standard Julia functions. In particular:Control sequences and conditionals (while, for, if) will not work in the macro.\nIntermediate calculations (likes that don\'t start with d_) are incompatible with the Jacobian etc. calculations.\nThe macro has to use t for the independent variable."
},

{
    "location": "analysis/parameterized_functions.html#Extra-Optimizations-1",
    "page": "ParameterizedFunctions",
    "title": "Extra Optimizations",
    "category": "section",
    "text": "Because the ParameterizedFunction defined by the macro holds the definition at a symbolic level, optimizations are provided by SymEngine. Using the symbolic calculator, in-place functions for many things such as Jacobians, Hessians, etc. are symbolically pre-computed. In addition, functions for the inverse Jacobian, Hessian, etc. are also pre-computed. In addition, parameter gradients and Jacobians are also used.Normally these will be computed fast enough that the user doesn\'t have to worry. However, in some cases you may want to restrict the number of functions (or get rid of a warning). For more information, please see the ParameterizedFunctions.jl documentation."
},

{
    "location": "analysis/parameter_estimation.html#",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "Parameter Estimation and Bayesian Analysis",
    "category": "page",
    "text": ""
},

{
    "location": "analysis/parameter_estimation.html#Parameter-Estimation-and-Bayesian-Analysis-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "Parameter Estimation and Bayesian Analysis",
    "category": "section",
    "text": "Parameter estimation for ODE models, also known as dynamic data analysis, is provided by the DiffEq suite."
},

{
    "location": "analysis/parameter_estimation.html#Installation-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "Installation",
    "category": "section",
    "text": "This functionality does not come standard with DifferentialEquations.jl. To use this functionality, you must install DiffEqParamEstim.jl:]add DiffEqParamEstim\nusing DiffEqParamEstimFor the Bayesian, methods, you must install DiffEqBayes.jl:]add DiffEqBayes\nusing DiffEqBayes"
},

{
    "location": "analysis/parameter_estimation.html#Recommended-Methods-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "Recommended Methods",
    "category": "section",
    "text": "The recommended method is to use build_loss_objective with the optimizer of your choice. This method can thus be paired with global optimizers from packages like BlackBoxOptim.jl or NLopt.jl which can be much less prone to finding local minima than local optimization methods. Also, it allows the user to define the cost function in the way they choose as a function loss(sol), and thus can fit using any cost function on the solution, making it applicable to fitting non-temporal data and other types of problems. Also, build_loss_objective works for all of the DEProblem types, allowing it to optimize parameters on ODEs, SDEs, DDEs, DAEs, etc.However, this method requires repeated solution of the differential equation. If the data is temporal data, the most efficient method is the two_stage_method which does not require repeated solutions but is not as accurate. Usage of the two_stage_method should have a post-processing step which refines using a method like build_loss_objective."
},

{
    "location": "analysis/parameter_estimation.html#Optimization-Based-Methods-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "Optimization-Based Methods",
    "category": "section",
    "text": ""
},

{
    "location": "analysis/parameter_estimation.html#Nonlinear-Regression-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "Nonlinear Regression",
    "category": "section",
    "text": "build_loss_objective builds an objective function to be used with Optim.jl and MathProgBase-associated solvers like NLopt.function build_loss_objective(prob::DEProblem,alg,loss_func\n                              regularization=nothing;\n                              mpg_autodiff = false,\n                              verbose_opt = false,\n                              verbose_steps = 100,\n                              prob_generator = (prob,p)->remake(prob,p=p),\n                              kwargs...)The first argument is the DEProblem to solve, and next is the alg to use. The alg must match the problem type, which can be any DEProblem (ODEs, SDEs, DAEs, DDEs, etc.). regularization defaults to nothing which has no regularization function. One can also choose verbose_opt and verbose_steps, which, in the optimization routines, will print the steps and the values at the steps every verbose_steps steps. mpg_autodiff uses autodifferentiation to define the derivative for the MathProgBase solver. The extra keyword arguments are passed to the differential equation solver."
},

{
    "location": "analysis/parameter_estimation.html#Multiple-Shooting-objective-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "Multiple Shooting objective",
    "category": "section",
    "text": "Multiple Shooting is generally used in Boundary Value Problems (BVP) and is more robust than the regular objective function used in these problems. It proceeds as follows:Divide up the time span into short time periods and solve the equationwith the current parameters which here consist of both, the parameters of the   differential equations and also the initial values for the short time periods.This objective additionally involves a discontinuity error term that imposeshigher cost if the end of the solution of one time period doesn\'t match the   beginning of the next one.Merge the solutions from the shorter intervals and then calculate the loss.For consistency multiple_shooting_objective takes exactly the same arguments as build_loss_objective. It also has the option for discontinuity_error as a keyword argument which assigns weight to the error occurring due to the discontinuity that arises from the breaking up of the time span."
},

{
    "location": "analysis/parameter_estimation.html#Two-Stage-method-(Non-Parametric-Collocation)-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "Two Stage method (Non-Parametric Collocation)",
    "category": "section",
    "text": "The two-stage method is a collocation method for estimating parameters without requiring repeated solving of the differential equation. It does so by determining a smoothed estimated trajectory of the data (local quadratic polynomial fit by least squares) and optimizing the derivative function and the data\'s timepoints to match the derivatives of the smoothed trajectory. This method has less accuracy than other methods but is much faster, and is a good method to try first to get in the general \"good parameter\" region, to then finish using one of the other methods.function two_stage_method(prob::DEProblem,tpoints,data;kernel= :Epanechnikov,\n                          loss_func = L2DistLoss,mpg_autodiff = false,\n                          verbose = false,verbose_steps = 100)"
},

{
    "location": "analysis/parameter_estimation.html#The-Loss-Function-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "The Loss Function",
    "category": "section",
    "text": "loss_func(sol)is a function which reduces the problem\'s solution to a scalar which the optimizer will try to minimize. While this is very flexible, two convenience routines are included for fitting to data with standard cost functions:L2Loss(t,data;differ_weight=nothing,data_weight=nothing,\n              colloc_grad=nothing,dudt=nothing)where t is the set of timepoints which the data is found at, and data are the values that are known where each column corresponds to measures of the values of the system. L2Loss is an optimized version of the L2-distance. The data_weight is a scalar or vector of weights for the loss function which must match the size of the data. Note that minimization of a weighted L2Loss is equivalent to maximum likelihood estimation of a heteroskedastic Normally distributed likelihood. differ_weight allows one to add a weight on the first differencing terms sol[i+1]-sol[i] against the data first differences. This smooths out the loss term and can make it easier to fit strong solutions of stochastic models, but is zero (nothing) by default. Additionally, colloc_grad allows one to give a matrix of the collocation gradients for the data. This is used to add an interpolation derivative term, like the two-stage method. A convenience function colloc_grad(t,data) returns a collocation gradient from a 3rd order spline calculated by Dierckx.jl, which can be used as the colloc_grad. Note that, with a collocation gradient and regularization, this loss is equivalent to a 4DVAR.Additionally, we include a more flexible log-likelihood approach:LogLikeLoss(t,distributions,diff_distributions=nothing)In this case, there are two forms. The simple case is where distributions[i,j] is the likelihood distributions from a UnivariateDistribution from Distributions.jl, where it corresponds to the likelihood at t[i] for component j. The second case is where distributions[i] is a MultivariateDistribution which corresponds to the likelihood at t[i] over the vector of components. This likelihood function then calculates the negative of the total loglikelihood over time as its objective value (negative since optimizers generally find minimums, and thus this corresponds to maximum likelihood estimation). The third term, diff_distributions, acts similarly but allows putting a distribution on the first difference terms sol[i+1]-sol[i].Note that these distributions can be generated via fit_mle on some dataset against some chosen distribution type."
},

{
    "location": "analysis/parameter_estimation.html#Note-About-Loss-Functions-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "Note About Loss Functions",
    "category": "section",
    "text": "For parameter estimation problems, it\'s not uncommon for the optimizers to hit unstable regions of parameter space. This causes warnings that the solver exited early, and the built-in loss functions like L2Loss automatically handle this. However, if using a user-supplied loss function, you should make sure it\'s robust to these issues. One common pattern is to apply infinite loss when the integration is not successful. Using the retcodes, this can be done via:function my_loss_function(sol)\n   tot_loss = 0.0\n   if any((s.retcode != :Success for s in sol))\n     tot_loss = Inf\n   else\n     # calculation for the loss here\n   end\n   tot_loss\nend"
},

{
    "location": "analysis/parameter_estimation.html#Note-on-First-Differencing-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "Note on First Differencing",
    "category": "section",
    "text": "L2Loss(t,data,differ_weight=0.3,data_weight=0.7)First differencing incorporates the differences of data points at consecutive time points which adds more information about the trajectory in the loss function. Adding first differencing is helpful in cases where the L2Loss alone leads to non-identifiable parameters but adding a first differencing term makes it more identifiable. This can be noted on stochastic differential equation models, where this aims to capture the autocorrelation and therefore helps us avoid getting the same stationary distribution despite different trajectories and thus wrong parameter estimates."
},

{
    "location": "analysis/parameter_estimation.html#The-Regularization-Function-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "The Regularization Function",
    "category": "section",
    "text": "The regularization can be any function of p, the parameter vector:regularization(p)The Regularization helper function builds a regularization using a penalty function penalty from PenaltyFunctions.jl:Regularization(λ,penalty=L2Penalty())The regularization defaults to L2 if no penalty function is specified. λ is the weight parameter for the addition of the regularization term."
},

{
    "location": "analysis/parameter_estimation.html#The-Problem-Generator-Function-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "The Problem Generator Function",
    "category": "section",
    "text": "The argument prob_generator allows one to specify a function for generating new problems from a given parameter set. By default, this just builds a new problem which fixes the element types in a way that\'s autodifferentiation compatible and adds the new parameter vector p. For example, the code for this is:prob_generator = (prob,p) -> remake(prob,u0=convert.(eltype(p),prob.u0),p=p)Then the new problem with these new values is returned.One can use this to change the meaning of the parameters using this function. For example, if one instead wanted to optimize the initial conditions for a function without parameters, you could change this to:prob_generator = (prob,p) -> remake(prob.f,u0=p)which simply uses p as the initial condition in the initial value problem."
},

{
    "location": "analysis/parameter_estimation.html#LeastSquaresOptim.jl-objective-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "LeastSquaresOptim.jl objective",
    "category": "section",
    "text": "build_lsoptim_objective builds an objective function to be used with LeastSquaresOptim.jl.build_lsoptim_objective(prob,tspan,t,data;\n                        prob_generator = (prob,p) -> remake(prob,u0=convert.(eltype(p),prob.u0),p=p),\n                        kwargs...)The arguments are the same as build_loss_objective."
},

{
    "location": "analysis/parameter_estimation.html#lm_fit-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "lm_fit",
    "category": "section",
    "text": "lm_fit is a function for fitting the parameters of an ODE using the Levenberg-Marquardt algorithm. This algorithm is really bad and thus not recommended since, for example, the Optim.jl algorithms on an L2 loss are more performant and robust. However, this is provided for completeness as most other differential equation libraries use an LM-based algorithm, so this allows one to test the increased effectiveness of not using LM.lm_fit(prob::DEProblem,tspan,t,data,p0;\n       prob_generator = (prob,p) -> remake(prob,u0=convert.(eltype(p),prob.u0),p=p),\n       kwargs...)The arguments are similar to before, but with p0 being the initial conditions for the parameters and the kwargs as the args passed to the LsqFit curve_fit function (which is used for the LM solver). This returns the fitted parameters."
},

{
    "location": "analysis/parameter_estimation.html#MAP-estimate-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "MAP estimate",
    "category": "section",
    "text": "You can also add a prior option to build_loss_objective and multiple_shooting_objective that essentially turns it into MAP by multiplying the loglikelihood (the cost) by the prior. The option was added as a keyword argument priors, it can take in either an array of univariate distributions for each of the parameters or a multivariate distribution.ms_obj = multiple_shooting_objective(ms_prob,Tsit5(),L2Loss(t,data);priors=priors,discontinuity_weight=1.0,abstol=1e-12,reltol=1e-12)"
},

{
    "location": "analysis/parameter_estimation.html#Bayesian-Methods-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "Bayesian Methods",
    "category": "section",
    "text": "The following methods require the DiffEqBayes.jl]add DiffEqBayes\nusing DiffEqBayes"
},

{
    "location": "analysis/parameter_estimation.html#stan_inference-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "stan_inference",
    "category": "section",
    "text": "stan_inference(prob::ODEProblem,t,data,priors = nothing;alg=:rk45,\n               num_samples=1000, num_warmup=1000, reltol=1e-3,\n               abstol=1e-6, maxiter=Int(1e5),likelihood=Normal,\n               vars=(StanODEData(),InverseGamma(2,3)))stan_inference uses Stan.jl to perform the Bayesian inference. The Stan installation process is required to use this function. The input requires that the function is defined by a ParameterizedFunction with the @ode_def macro. t is the array of time and data is the array where the first dimension (columns) corresponds to the array of system values. priors is an array of prior distributions for each parameter, specified via a Distributions.jl type. alg is a choice between :rk45 and :bdf, the two internal integrators of Stan. num_samples is the number of samples to take per chain, and num_warmup is the number of MCMC warmup steps. abstol and reltol are the keyword arguments for the internal integrator. likelihood is the likelihood distribution to use with the arguments from vars, and vars is a tuple of priors for the distributions of the likelihood hyperparameters. The special value StanODEData() in this tuple denotes the position that the ODE solution takes in the likelihood\'s parameter list."
},

{
    "location": "analysis/parameter_estimation.html#turing_inference-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "turing_inference",
    "category": "section",
    "text": "function turing_inference(prob::DiffEqBase.DEProblem,alg,t,data,priors; \n                              likelihood_dist_priors, likelihood, num_samples=1000, \n                              sampler = Turing.NUTS(num_samples, 0.65), syms, kwargs...)turing_inference uses Turing.jl to perform its parameter inference. prob can be any DEProblem with a corresponding alg choice. t is the array of time points and data is the set of observations for the differential equation system at time point t[i] (or higher dimensional). priors is an array of prior distributions for each parameter, specified via a Distributions.jl type. num_samples is the number of samples per MCMC chain. The extra kwargs are given to the internal differential equation solver."
},

{
    "location": "analysis/parameter_estimation.html#dynamichmc_inference-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "dynamichmc_inference",
    "category": "section",
    "text": "dynamichmc_inference(prob::DEProblem,alg,t,data,priors,transformations;\n                      σ = 0.01,ϵ=0.001,initial=Float64[])dynamichmc_inference uses DynamicHMC.jl to  perform the bayesian parameter estimation. prob can be any DEProblem, data is the set  of observations for our model which is to be used in the Bayesian Inference process. priors represent the  choice of prior distributions for the parameters to be determined, passed as an array of Distributions.jl distributions. t is the array of time points. transformations  is an array of Tranformations imposed for constraining the  parameter values to specific domains. initial values for the parameters can be passed, if not passed the means of the  priors are used. ϵ can be used as a kwarg to pass the initial step size for the NUTS algorithm.      "
},

{
    "location": "analysis/parameter_estimation.html#abc_inference-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "abc_inference",
    "category": "section",
    "text": "abc_inference(prob::DEProblem, alg, t, data, priors; ϵ=0.001,\n     distancefunction = euclidean, ABCalgorithm = ABCSMC, progress = false,\n     num_samples = 500, maxiterations = 10^5, kwargs...)abc_inference uses ApproxBayes.jl which uses Approximate Bayesian Computation (ABC) to perform its parameter inference. prob can be any DEProblem with a corresponding alg choice. t is the array of time points and data[:,i] is the set of observations for the differential equation system at time point t[i] (or higher dimensional). priors is an array of prior distributions for each parameter, specified via a Distributions.jl type. num_samples is the number of posterior samples. ϵ is the target distance between the data and simulated data. distancefunction is a distance metric specified from the Distances.jl package, the default is euclidean. ABCalgorithm is the ABC algorithm to use, options are ABCSMC or ABCRejection from ApproxBayes.jl, the default is the former which is more efficient. maxiterations is the maximum number of iterations before the algorithm terminates. The extra kwargs are given to the internal differential equation solver."
},

{
    "location": "analysis/parameter_estimation.html#Optimization-Based-ODE-Inference-Examples-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "Optimization-Based ODE Inference Examples",
    "category": "section",
    "text": ""
},

{
    "location": "analysis/parameter_estimation.html#Simple-Local-Optimization-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "Simple Local Optimization",
    "category": "section",
    "text": "We choose to optimize the parameters on the Lotka-Volterra equation. We do so by defining the function as a ParameterizedFunction:function f(du,u,p,t)\n  du[1] = dx = p[1]*u[1] - u[1]*u[2]\n  du[2] = dy = -3*u[2] + u[1]*u[2]\nend\n\nu0 = [1.0;1.0]\ntspan = (0.0,10.0)\np = [1.5]\nprob = ODEProblem(f,u0,tspan,p)We create data using the numerical result with a=1.5:sol = solve(prob,Tsit5())\nt = collect(range(0,stop=10,length=200))\nusing RecursiveArrayTools # for VectorOfArray\nrandomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])\ndata = convert(Array,randomized)Here we used VectorOfArray from RecursiveArrayTools.jl to turn the result of an ODE into a matrix.If we plot the solution with the parameter at a=1.42, we get the following:(Image: Parameter Estimation Not Fit)Notice that after one period this solution begins to drift very far off: this problem is sensitive to the choice of a.To build the objective function for Optim.jl, we simply call the build_loss_objective function:cost_function = build_loss_objective(prob,Tsit5(),L2Loss(t,data),\n                                     maxiters=10000,verbose=false)This objective function internally is calling the ODE solver to get solutions to test against the data. The keyword arguments are passed directly to the solver. Note that we set maxiters in a way that causes the differential equation solvers to error more quickly when in bad regions of the parameter space, speeding up the process. If the integrator stops early (due to divergence), then those parameters are given an infinite loss, and thus this is a quick way to avoid bad parameters. We set verbose=false because this divergence can get noisy.Before optimizing, let\'s visualize our cost function by plotting it for a range of parameter values:vals = 0.0:0.1:10.0\nusing Plots; plotly()\nplot(vals,[cost_function(i) for i in vals],yscale=:log10,\n     xaxis = \"Parameter\", yaxis = \"Cost\", title = \"1-Parameter Cost Function\",\n     lw = 3)(Image: 1 Parameter Likelihood)Here we see that there is a very well-defined minimum in our cost function at the real parameter (because this is where the solution almost exactly fits the dataset).Now this cost function can be used with Optim.jl in order to get the parameters. For example, we can use Brent\'s algorithm to search for the best solution on the interval [0,10] by:using Optim\nresult = optimize(cost_function, 0.0, 10.0)This returns result.minimizer[1]==1.5 as the best parameter to match the data. When we plot the fitted equation on the data, we receive the following:(Image: Parameter Estimation Fit)Thus we see that after fitting, the lines match up with the generated data and receive the right parameter value.We can also use the multivariate optimization functions. For example, we can use the BFGS algorithm to optimize the parameter starting at a=1.42 using:result = optimize(cost_function, [1.42], BFGS())Note that some of the algorithms may be sensitive to the initial condition. For more details on using Optim.jl, see the documentation for Optim.jl. We can improve our solution by noting that the Lotka-Volterra equation requires that the parameters are positive. Thus following the Optim.jl documentation we can add box constraints to ensure the optimizer only checks between 0.0 and 3.0 which improves the efficiency of our algorithm:lower = [0.0]\nupper = [3.0]\nresult = optimize(cost_function, [1.42], lower, upper, Fminbox{BFGS}())Lastly, we can use the same tools to estimate multiple parameters simultaneously. Let\'s use the Lotka-Volterra equation with all parameters free:function f2(du,u,p,t)\n  du[1] = dx = p[1]*u[1] - p[2]*u[1]*u[2]\n  du[2] = dy = -p[3]*u[2] + p[4]*u[1]*u[2]\nend\n\nu0 = [1.0;1.0]\ntspan = (0.0,10.0)\np = [1.5,1.0,3.0,1.0]\nprob = ODEProblem(f2,u0,tspan,p)We can build an objective function and solve the multiple parameter version just as before:cost_function = build_loss_objective(prob,Tsit5(),L2Loss(t,data),\n                                      maxiters=10000,verbose=false)\nresult_bfgs = Optim.optimize(cost_function, [1.3,0.8,2.8,1.2], Optim.BFGS())We can also use First-Differences in L2Loss by passing the kwarg differ_weight which decides the contribution of the differencing loss to the total loss.cost_function = build_loss_objective(prob,Tsit5(),L2Loss(t,data,differ_weight=0.3,data_weight=0.7),\n                                      maxiters=10000,verbose=false)\nresult_bfgs = Optim.optimize(cost_function, [1.3,0.8,2.8,1.2], Optim.BFGS())To solve it using LeastSquaresOptim.jl, we use the build_lsoptim_objective function:cost_function = build_lsoptim_objective(prob1,t,data,Tsit5())The result is a cost function which can be used with LeastSquaresOptim. For more details, consult the documentation for LeastSquaresOptim.jl:x = [1.3,0.8,2.8,1.2]\nres = optimize!(LeastSquaresProblem(x = x, f! = cost_function,\n                output_length = length(t)*length(prob.u0)),\n                LeastSquaresOptim.Dogleg(),LeastSquaresOptim.LSMR())We can see the results are:println(res.minimizer)\n\nResults of Optimization Algorithm\n * Algorithm: Dogleg\n * Minimizer: [1.4995074428834114,0.9996531871795851,3.001556360700904,1.0006272074128821]\n * Sum of squares at Minimum: 0.035730\n * Iterations: 63\n * Convergence: true\n * |x - x\'| < 1.0e-15: true\n * |f(x) - f(x\')| / |f(x)| < 1.0e-14: false\n * |g(x)| < 1.0e-14: false\n * Function Calls: 64\n * Gradient Calls: 9\n * Multiplication Calls: 135and thus this algorithm was able to correctly identify all four parameters.We can also use Multiple Shooting method by creating a multiple_shooting_objectivefunction ms_f(du,u,p,t)\n  dx = p[1]*u[1] - p[2]*u[1]*u[2]\n  dy = -3*u[2] + u[1]*u[2]\nend\nms_u0 = [1.0;1.0]\ntspan = (0.0,10.0)\nms_p = [1.5,1.0]\nms_prob = ODEProblem(ms_f,ms_u0,tspan,ms_p)\nt = collect(range(0,stop=10,length=200))\ndata = Array(solve(ms_prob,Tsit5(),saveat=t,abstol=1e-12,reltol=1e-12))\nbound = Tuple{Float64, Float64}[(0, 10),(0, 10),(0, 10),(0, 10),\n                                (0, 10),(0, 10),(0, 10),(0, 10),\n                                (0, 10),(0, 10),(0, 10),(0, 10),\n                                (0, 10),(0, 10),(0, 10),(0, 10),(0, 10),(0, 10)]\n\n\nms_obj = multiple_shooting_objective(ms_prob,Tsit5(),L2Loss(t,data);discontinuity_weight=1.0,abstol=1e-12,reltol=1e-12)This creates the objective function that can be passed to an optimizer from which we can then get the parameter values and the initial values of the short time periods keeping in mind the indexing.result = bboptimize(ms_obj;SearchRange = bound, MaxSteps = 21e3)\nresult.archive_output.best_candidate[end-1:end]Giving us the results asStarting optimization with optimizer BlackBoxOptim.DiffEvoOpt{BlackBoxOptim.FitPopulation{Float64},BlackBoxOptim.RadiusLimitedSelector,BlackBoxOptim.AdaptiveDiffEvoRandBin{3},BlackBoxOptim.RandomBound{BlackBoxOptim.RangePerDimSearchSpace}}\n\nOptimization stopped after 21001 steps and 136.60030698776245 seconds\nTermination reason: Max number of steps (21000) reached\nSteps per second = 153.7405036862868\nFunction evals per second = 154.43596332393247\nImprovements/step = 0.17552380952380953\nTotal function evaluations = 21096\n\n\nBest candidate found: [0.997396, 1.04664, 3.77834, 0.275823, 2.14966, 4.33106, 1.43777, 0.468442, 6.22221, 0.673358, 0.970036, 2.05182, 2.4216, 0.274394, 5.64131, 3.38132, 1.52826, 1.01721]\n\nFitness: 0.126884213\n\nOut[4]:2-element Array{Float64,1}:\n        1.52826\n        1.01721Here as our model had 2 parameters, we look at the last two indexes of result to get our parameter values and the rest of the values are the initial values of the shorter timespans as described in the reference section.The objective function for Two Stage method can be created and passed to an optimizer astwo_stage_obj = two_stage_method(ms_prob,t,data)\nresult = Optim.optimize(two_stage_obj, [1.3,0.8,2.8,1.2], Optim.BFGS()\n)\nResults of Optimization Algorithm\n * Algorithm: BFGS\n * Starting Point: [1.3,0.8,2.8,1.2]\n * Minimizer: [1.5035938533664717,0.9925731153746833, ...]\n * Minimum: 1.513400e+00\n * Iterations: 9\n * Convergence: true\n   * |x - x\'| ≤ 0.0e+00: false\n     |x - x\'| = 4.58e-10\n   * |f(x) - f(x\')| ≤ 0.0e+00 |f(x)|: false\n     |f(x) - f(x\')| = 5.87e-16 |f(x)|\n   * |g(x)| ≤ 1.0e-08: true\n     |g(x)| = 7.32e-11\n   * Stopped by an increasing objective: false\n   * Reached Maximum Number of Iterations: false\n * Objective Calls: 31\n * Gradient Calls: 31The default kernel used in the method is Epanechnikov others that are available are Uniform,  Triangular, Quartic, Triweight, Tricube, Gaussian, Cosine, Logistic and Sigmoid, this can be passed by the kernel keyword argument. loss_func keyword argument can be used to pass the loss function (cost function) you want  to use and mpg_autodiff enables Auto Differentiation."
},

{
    "location": "analysis/parameter_estimation.html#More-Algorithms-(Global-Optimization)-via-MathProgBase-Solvers-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "More Algorithms (Global Optimization) via MathProgBase Solvers",
    "category": "section",
    "text": "The build_loss_objective function builds an objective function which is able to be used with MathProgBase-associated solvers. This includes packages like IPOPT, NLopt, MOSEK, etc. Building off of the previous example, we can build a cost function for the single parameter optimization problem like:function f(du,u,p,t)\n  dx = p[1]*u[1] - u[1]*u[2]\n  dy = -3*u[2] + u[1]*u[2]\nend\n\nu0 = [1.0;1.0]\ntspan = (0.0,10.0)\np = [1.5]\nprob = ODEProblem(f,u0,tspan,p)\nsol = solve(prob,Tsit5())\n\nt = collect(range(0,stop=10,length=200))\nrandomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])\ndata = convert(Array,randomized)\n\nobj = build_loss_objective(prob,Tsit5(),L2Loss(t,data),maxiters=10000)We can now use this obj as the objective function with MathProgBase solvers. For our example, we will use NLopt. To use the local derivative-free Constrained Optimization BY Linear Approximations algorithm, we can simply do:using NLopt\nopt = Opt(:LN_COBYLA, 1)\nmin_objective!(opt, obj)\n(minf,minx,ret) = NLopt.optimize(opt,[1.3])This finds a minimum at [1.49997]. For a modified evolutionary algorithm, we can use:opt = Opt(:GN_ESCH, 1)\nmin_objective!(opt, obj)\nlower_bounds!(opt,[0.0])\nupper_bounds!(opt,[5.0])\nxtol_rel!(opt,1e-3)\nmaxeval!(opt, 100000)\n(minf,minx,ret) = NLopt.optimize(opt,[1.3])We can even use things like the Improved Stochastic Ranking Evolution Strategy (and add constraints if needed). This is done via:opt = Opt(:GN_ISRES, 1)\nmin_objective!(opt, obj.cost_function2)\nlower_bounds!(opt,[-1.0])\nupper_bounds!(opt,[5.0])\nxtol_rel!(opt,1e-3)\nmaxeval!(opt, 100000)\n(minf,minx,ret) = NLopt.optimize(opt,[0.2])which is very robust to the initial condition. The fastest result comes from the following:using NLopt\nopt = Opt(:LN_BOBYQA, 1)\nmin_objective!(opt, obj)\n(minf,minx,ret) = NLopt.optimize(opt,[1.3])For more information, see the NLopt documentation for more details. And give IPOPT or MOSEK a try!"
},

{
    "location": "analysis/parameter_estimation.html#Using-JuMP-with-DiffEqParamEstim-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "Using JuMP with DiffEqParamEstim",
    "category": "section",
    "text": "JuMP is a domain-specific modeling language  for mathematical optimization embedded in Julia.using OrdinaryDiffEq, DiffEqParamEstim, JuMP, NLopt, PlotsLet\'s define the Lorenz equation to use as our examplefunction g(du,u,p,t)\n  σ,ρ,β = p\n  x,y,z = u\n  du[1] = dx = σ*(y-x)\n  du[2] = dy = x*(ρ-z) - y\n  du[3] = dz = x*y - β*z\nendLet\'s get a solution of the system with parameter values σ=10.0 ρ=28.0 β=8/3 to use as our data. We define some convenience functions model_ode (to create an ODEProblem) and solve_model(to obtain  solution of the ODEProblem) to use in a custom objective function later.u0 = [1.0;0.0;0.0]\nt = 0.0:0.01:1.0\ntspan = (0.0,1.0)\nmodel_ode(p_) = ODEProblem(g, u0, tspan,p_) \nsolve_model(mp_) = OrdinaryDiffEq.solve(model_ode(mp_), Tsit5(),saveat=0.01)\nmock_data = Array(solve_model([10.0,28.0,8/3]))Now we define a custom objective function to pass for optimization to JuMP using the build_loss_objective described above provided by DiffEqParamEstim that defines an objective  function for the parameter estimation problem.loss_objective(mp_, dat) = build_loss_objective(model_ode(mp_), Tsit5(), L2Loss(t,dat))We create a JuMP model, variables, set the objective function and the choice of  optimization algorithm to be used in the JuMP syntax. You can read more about this in  JuMP\'s documentation.juobj(args...) = loss_objective(args, mock_data)(args)\njumodel = Model()\nJuMP.register(jumodel, :juobj, 3, juobj, autodiff=true)\n@variables jumodel begin\n    σ,(start=8)\n    ρ,(start=25.0)\n    β,(start=10/3)\nend\n@NLobjective(jumodel, Min, juobj(σ, ρ, β))\nsetsolver(jumodel, NLoptSolver(algorithm=:LD_MMA))Let\'s call the optimizer to obtain the fitted parameter values.sol = JuMP.solve(jumodel)\nbest_mp = getvalue.(getindex.((jumodel,), Symbol.(jumodel.colNames)))Let\'s compare the solution at the obtained parameter values and our data. sol = OrdinaryDiffEq.solve(best_mp |> model_ode, Tsit5())\nplot(getindex.(sol.(t),1))\nscatter!(mock_data, markersize=2)(Image: jumpestimationplot)"
},

{
    "location": "analysis/parameter_estimation.html#Generalized-Likelihood-Example-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "Generalized Likelihood Example",
    "category": "section",
    "text": "In this example we will demo the likelihood-based approach to parameter fitting. First let\'s generate a dataset to fit. We will re-use the Lotka-Volterra equation but in this case fit just two parameters.f1 = function (du,u,p,t)\n  du[1] = p[1] * u[1] - p[2] * u[1]*u[2]\n  du[2] = -3.0 * u[2] + u[1]*u[2]\nend\np = [1.5,1.0]\nu0 = [1.0;1.0]\ntspan = (0.0,10.0)\nprob1 = ODEProblem(f1,u0,tspan,p)\nsol = solve(prob1,Tsit5())This is a function with two parameters, [1.5,1.0] which generates the same ODE solution as before. This time, let\'s generate 100 datasets where at each point adds a little bit of randomness:using RecursiveArrayTools # for VectorOfArray\nt = collect(range(0,stop=10,length=200))\nfunction generate_data(sol,t)\n  randomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])\n  data = convert(Array,randomized)\nend\naggregate_data = convert(Array,VectorOfArray([generate_data(sol,t) for i in 1:100]))here with t we measure the solution at 200 evenly spaced points. Thus aggregate_data is a 2x200x100 matrix where aggregate_data[i,j,k] is the ith component at time j of the kth dataset. What we first want to do is get a matrix of distributions where distributions[i,j] is the likelihood of component i at take j. We can do this via fit_mle on a chosen distributional form. For simplicity we choose the Normal distribution. aggregate_data[i,j,:] is the array of points at the given component and time, and thus we find the distribution parameters which fits best at each time point via:using Distributions\ndistributions = [fit_mle(Normal,aggregate_data[i,j,:]) for i in 1:2, j in 1:200]Notice for example that we have:julia> distributions[1,1]\nDistributions.Normal{Float64}(μ=1.0022440583676806, σ=0.009851964521952437)that is, it fit the distribution to have its mean just about where our original solution was and the variance is about how much noise we added to the dataset. This this is a good check to see that the distributions we are trying to fit our parameters to makes sense.Note that in this case the Normal distribution was a good choice, and in many cases it\'s a nice go-to choice, but one should experiment with other choices of distributions as well. For example, a TDist can be an interesting way to incorporate robustness to outliers since low degrees of free T-distributions act like Normal distributions but with longer tails (though fit_mle does not work with a T-distribution, you can get the means/variances and build appropriate distribution objects yourself).Once we have the matrix of distributions, we can build the objective function corresponding to that distribution fit:using DiffEqParamEstim\nobj = build_loss_objective(prob1,Tsit5(),LogLikeLoss(t,distributions),\n                                     maxiters=10000,verbose=false)First let\'s use the objective function to plot the likelihood landscape:using Plots; plotly()\nrange = 0.5:0.1:5.0\nheatmap(range,range,[obj([j,i]) for i in range, j in range],\n        yscale=:log10,xlabel=\"Parameter 1\",ylabel=\"Parameter 2\",\n        title=\"Likelihood Landscape\")(Image: 2 Parameter Likelihood)Recall that this is the negative loglikelihood and thus the minimum is the maximum of the likelihood. There is a clear valley where the second parameter is 1.5, while the first parameter\'s likelihood is more muddled. By taking a one-dimensional slice:plot(range,[obj([i,1.0]) for i in range],lw=3,\n     title=\"Parameter 1 Likelihood (Parameter 2 = 1.5)\",\n     xlabel = \"Parameter 1\", ylabel = \"Objective Function Value\")(Image: 1 Parameter Likelihood)we can see that there\'s still a clear minimum at the true value. Thus we will use the global optimizers from BlackBoxOptim.jl to find the values. We set our search range to be from 0.5 to 5.0 for both of the parameters and let it optimize:using BlackBoxOptim\nbound1 = Tuple{Float64, Float64}[(0.5, 5),(0.5, 5)]\nresult = bboptimize(obj;SearchRange = bound1, MaxSteps = 11e3)\n\nStarting optimization with optimizer BlackBoxOptim.DiffEvoOpt{BlackBoxOptim.FitPopulation{Float64},B\nlackBoxOptim.RadiusLimitedSelector,BlackBoxOptim.AdaptiveDiffEvoRandBin{3},BlackBoxOptim.RandomBound\n{BlackBoxOptim.RangePerDimSearchSpace}}\n0.00 secs, 0 evals, 0 steps\n0.50 secs, 1972 evals, 1865 steps, improv/step: 0.266 (last = 0.2665), fitness=-737.311433781\n1.00 secs, 3859 evals, 3753 steps, improv/step: 0.279 (last = 0.2913), fitness=-739.658421879\n1.50 secs, 5904 evals, 5799 steps, improv/step: 0.280 (last = 0.2830), fitness=-739.658433715\n2.00 secs, 7916 evals, 7811 steps, improv/step: 0.225 (last = 0.0646), fitness=-739.658433715\n2.50 secs, 9966 evals, 9861 steps, improv/step: 0.183 (last = 0.0220), fitness=-739.658433715\n\nOptimization stopped after 11001 steps and 2.7839999198913574 seconds\nTermination reason: Max number of steps (11000) reached\nSteps per second = 3951.50873439296\nFunction evals per second = 3989.2242527195904\nImprovements/step = 0.165\nTotal function evaluations = 11106\n\n\nBest candidate found: [1.50001, 1.00001]\n\nFitness: -739.658433715This shows that it found the true parameters as the best fit to the likelihood."
},

{
    "location": "analysis/parameter_estimation.html#Parameter-Estimation-for-Stochastic-Differential-Equations-and-Monte-Carlo-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "Parameter Estimation for Stochastic Differential Equations and Monte Carlo",
    "category": "section",
    "text": "We can use any DEProblem, which not only includes DAEProblem and DDEProblems, but also stochastic problems. In this case, let\'s use the generalized maximum likelihood to fit the parameters of an SDE\'s Monte Carlo evaluation.Let\'s use the same Lotka-Volterra equation as before, but this time add noise:pf_func = function (du,u,p,t)\n  du[1] = p[1] * u[1] - p[2] * u[1]*u[2]\n  du[2] = -3 * u[2] + u[1]*u[2]\nend\n\nu0 = [1.0;1.0]\ntspan = (0.0,10.0)\np = [1.5,1.0]\npg_func = function (du,u,p,t)\n  du[1] = 1e-6u[1]\n  du[2] = 1e-6u[2]\nend\nprob = SDEProblem(pf_func,pg_func,u0,tspan,p)\nsol = solve(prob,SRIW1())Now lets generate a dataset from 10,000 solutions of the SDEusing RecursiveArrayTools # for VectorOfArray\nt = collect(range(0, stop=10, length=200))\nfunction generate_data(t)\n  sol = solve(prob,SRIW1())\n  randomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])\n  data = convert(Array,randomized)\nend\naggregate_data = convert(Array,VectorOfArray([generate_data(t) for i in 1:10000]))Now let\'s estimate the parameters. Instead of using single runs from the SDE, we will use a MonteCarloProblem. This means that it will solve the SDE N times to come up with an approximate probability distribution at each time point and use that in the likelihood estimate.monte_prob = MonteCarloProblem(prob)We use Optim.jl for optimization belowobj = build_loss_objective(monte_prob,SOSRI(),L2Loss(t,aggregate_data),\n                                     maxiters=10000,verbose=false,num_monte = 1000,\n                                     parallel_type = :threads)\nresult = Optim.optimize(obj, [1.0,0.5], Optim.BFGS())Parameter Estimation in case of SDE\'s with a regular L2Loss can have poor accuracy due to only fitting against the mean properties as mentioned in First Differencing.Results of Optimization Algorithm\n * Algorithm: BFGS\n * Starting Point: [1.0,0.5]\n * Minimizer: [6.070728870478734,5.113357737345448]\n * Minimum: 1.700440e+03\n * Iterations: 14\n * Convergence: false\n   * |x - x\'| ≤ 0.0e+00: false\n     |x - x\'| = 1.00e-03\n   * |f(x) - f(x\')| ≤ 0.0e+00 |f(x)|: false\n     |f(x) - f(x\')| = 1.81e-07 |f(x)|\n   * |g(x)| ≤ 1.0e-08: false\n     |g(x)| = 2.34e+00\n   * Stopped by an increasing objective: true\n   * Reached Maximum Number of Iterations: false\n * Objective Calls: 61\n * Gradient Calls: 61Instead when we use L2Loss with first differencing enabled we get much more accurate estimates. obj = build_loss_objective(monte_prob,SRIW1(),L2Loss(t,data,differ_weight=1.0,data_weight=0.5),maxiters=1000,\n                                  verbose=false,verbose_opt=false,verbose_steps=1,num_monte=50)\nresult = Optim.optimize(obj, [1.0,0.5], Optim.BFGS())\nResults of Optimization Algorithm\n * Algorithm: BFGS\n * Starting Point: [1.0,0.5]\n * Minimizer: [1.5010687426045128,1.0023453619050238]\n * Minimum: 1.166650e-01\n * Iterations: 16\n * Convergence: false\n   * |x - x\'| ≤ 0.0e+00: false\n     |x - x\'| = 6.84e-09\n   * |f(x) - f(x\')| ≤ 0.0e+00 |f(x)|: false\n     |f(x) - f(x\')| = 5.85e-06 |f(x)|\n   * |g(x)| ≤ 1.0e-08: false\n     |g(x)| = 1.81e-01\n   * Stopped by an increasing objective: true\n   * Reached Maximum Number of Iterations: false\n * Objective Calls: 118\n * Gradient Calls: 118Here we see that we successfully recovered the drift parameter, and got close to the original noise parameter after searching a two orders of magnitude range."
},

{
    "location": "analysis/parameter_estimation.html#Bayesian-Inference-Examples-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "Bayesian Inference Examples",
    "category": "section",
    "text": ""
},

{
    "location": "analysis/parameter_estimation.html#Stan-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "Stan",
    "category": "section",
    "text": "Like in the previous examples, we setup the Lotka-Volterra system and generate data. Note that using @ode_def here is required.f1 = @ode_def begin\n  dx = a*x - b*x*y\n  dy = -c*y + d*x*y\nend a b c d\np = [1.5,1.0,3.0,1.0]\nu0 = [1.0,1.0]\ntspan = (0.0,10.0)\nprob1 = ODEProblem(f1,u0,tspan,p)\nsol = solve(prob1,Tsit5())\nt = collect(range(1,stop=10,length=10))\nrandomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])\ndata = convert(Array,randomized)Here we now give Stan an array of prior distributions for our parameters. Since the parameters of our differential equation must be positive, we utilize truncated Normal distributions to make sure that is satisfied in the result:priors = [Truncated(Normal(1.5,0.1),0,2),Truncated(Normal(1.0,0.1),0,1.5),\n          Truncated(Normal(3.0,0.1),0,4),Truncated(Normal(1.0,0.1),0,2)]We then give these to the inference function.bayesian_result = stan_inference(prob1,t,data,priors;\n                                 num_samples=100,num_warmup=500,\n                                 vars = (StanODEData(),InverseGamma(4,1)))InverseGamma(4,1) is our starting estimation for the variance hyperparameter of the default Normal distribution. The result is a Mamba.jl chain object. We can pull out the parameter values via:theta1 = bayesian_result.chain_results[:,[\"theta.1\"],:]\ntheta2 = bayesian_result.chain_results[:,[\"theta.2\"],:]\ntheta3 = bayesian_result.chain_results[:,[\"theta.3\"],:]\ntheta4 = bayesian_result.chain_results[:,[\"theta.4\"],:]From these chains we can get our estimate for the parameters via:mean(theta1.value[:,:,1])We can get more of a description via:Mamba.describe(bayesian_result.chain_results)\n\n# Result\n\nIterations = 1:100\nThinning interval = 1\nChains = 1,2,3,4\nSamples per chain = 100\n\nEmpirical Posterior Estimates:\n                  Mean         SD        Naive SE        MCSE         ESS    \n         lp__ -6.15472697 1.657551334 0.08287756670 0.18425029767  80.9314979\naccept_stat__  0.90165904 0.125913744 0.00629568721 0.02781181930  20.4968668\n   stepsize__  0.68014975 0.112183047 0.00560915237 0.06468790087   3.0075188\n  treedepth__  2.68750000 0.524911975 0.02624559875 0.10711170182  24.0159141\n n_leapfrog__  6.77000000 4.121841086 0.20609205428 0.18645821695 100.0000000\n  divergent__  0.00000000 0.000000000 0.00000000000 0.00000000000         NaN\n     energy__  9.12245750 2.518330231 0.12591651153 0.32894488320  58.6109941\n     sigma1.1  0.57164997 0.128579363 0.00642896816 0.00444242658 100.0000000\n     sigma1.2  0.58981422 0.131346442 0.00656732209 0.00397310122 100.0000000\n       theta1  1.50237077 0.008234095 0.00041170473 0.00025803930 100.0000000\n       theta2  0.99778276 0.009752574 0.00048762870 0.00009717115 100.0000000\n       theta3  3.00087782 0.009619775 0.00048098873 0.00020301023 100.0000000\n       theta4  0.99803569 0.008893244 0.00044466218 0.00040886528 100.0000000\n      theta.1  1.50237077 0.008234095 0.00041170473 0.00025803930 100.0000000\n      theta.2  0.99778276 0.009752574 0.00048762870 0.00009717115 100.0000000\n      theta.3  3.00087782 0.009619775 0.00048098873 0.00020301023 100.0000000\n      theta.4  0.99803569 0.008893244 0.00044466218 0.00040886528 100.0000000\n\nQuantiles:\n                  2.5%        25.0%      50.0%      75.0%       97.5%   \n         lp__ -10.11994750 -7.0569000 -5.8086150 -4.96936500 -3.81514375\naccept_stat__   0.54808912  0.8624483  0.9472840  0.98695850  1.00000000\n   stepsize__   0.57975100  0.5813920  0.6440120  0.74276975  0.85282400\n  treedepth__   2.00000000  2.0000000  3.0000000  3.00000000  3.00000000\n n_leapfrog__   3.00000000  7.0000000  7.0000000  7.00000000 15.00000000\n  divergent__   0.00000000  0.0000000  0.0000000  0.00000000  0.00000000\n     energy__   5.54070300  7.2602200  8.7707000 10.74517500 14.91849500\n     sigma1.1   0.38135240  0.4740865  0.5533195  0.64092575  0.89713635\n     sigma1.2   0.39674703  0.4982615  0.5613655  0.66973025  0.88361407\n       theta1   1.48728600  1.4967650  1.5022750  1.50805500  1.51931475\n       theta2   0.97685115  0.9914630  0.9971435  1.00394250  1.01765575\n       theta3   2.98354100  2.9937575  3.0001450  3.00819000  3.02065950\n       theta4   0.97934128  0.9918495  0.9977415  1.00430750  1.01442975\n      theta.1   1.48728600  1.4967650  1.5022750  1.50805500  1.51931475\n      theta.2   0.97685115  0.9914630  0.9971435  1.00394250  1.01765575\n      theta.3   2.98354100  2.9937575  3.0001450  3.00819000  3.02065950\n      theta.4   0.97934128  0.9918495  0.9977415  1.00430750  1.01442975More extensive information about the distributions is given by the plots:plot_chain(bayesian_result)"
},

{
    "location": "analysis/parameter_estimation.html#Turing-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "Turing",
    "category": "section",
    "text": "This case we will build off of the Stan example. Note that turing_inference does not require the use of the @ode_def macro like Stan does, but it will still work with macro-defined functions. Thus, using the same setup as before, we simply give the setup to:bayesian_result = turing_inference(prob,Tsit5(),t,data,priors;num_samples=500)The chain for the ith parameter is then given by:bayesian_result[:theta1]Summary statistics can be also be accessed:Mamba.describe(bayesian_result)The chain can be analysed by the trace plots and other plots obtained by:plot_chain(bayesian_result)"
},

{
    "location": "analysis/parameter_estimation.html#DynamicHMC-1",
    "page": "Parameter Estimation and Bayesian Analysis",
    "title": "DynamicHMC",
    "category": "section",
    "text": "We can use DynamicHMC.jl as the backend for sampling with the dynamic_inference function. It supports any DEProblem, priors can be passed as an array of Distributions.jl distributions, passing initial values is optional and in case where the user has a firm understanding of the domain the parameter values will lie in, transformations can be used to pass an array of constraints for the parameters as an array of Transformations.bayesian_result_hmc = dynamichmc_inference(prob1, Tsit5(), t, data, [Normal(1.5, 1)], [bridge(ℝ, ℝ⁺, )])A tuple with summary statistics and the chain values is returned. The chain for the ith parameter is given by:bayesian_result_hmc[1][i]For accessing the various summary statistics:DynamicHMC.NUTS_statistics(bayesian_result_dynamic[2])Some details about the NUTS sampler can be obtained from:bayesian_result_dynamic[3]In case of dynamic_inference the trace plots for the ith parameter can be obtained by:plot(bayesian_result_hmc[1][i])For a better idea of the summary statistics and plotting you can take a look at the benchmark notebooks"
},

{
    "location": "analysis/bifurcation.html#",
    "page": "Bifurcation Analysis",
    "title": "Bifurcation Analysis",
    "category": "page",
    "text": ""
},

{
    "location": "analysis/bifurcation.html#Bifurcation-Analysis-1",
    "page": "Bifurcation Analysis",
    "title": "Bifurcation Analysis",
    "category": "section",
    "text": "Bifurcation analysis is provided by the wrapper package PyDSTool.jl, which wraps the functionality of PyDSTool. The the package has an interface for directly using PyDSTool itself, included is a higher level interface that makes these tools compatible with more standard JuliaDiffEq types."
},

{
    "location": "analysis/bifurcation.html#Installation-1",
    "page": "Bifurcation Analysis",
    "title": "Installation",
    "category": "section",
    "text": "This functionality does not come standard with DifferentialEquations.jl. To use this functionality, you must install PyDSTool.jl:]add PyDSTool\nusing PyDSTool"
},

{
    "location": "analysis/bifurcation.html#Calcium-Bifurcation-Tutorial-1",
    "page": "Bifurcation Analysis",
    "title": "Calcium Bifurcation Tutorial",
    "category": "section",
    "text": "In this tutorial we will show how to do some simple bifurcation plots. We will follow the PyDSTool tutorial for the calcium channel model and re-create the results using the wrapped functionality."
},

{
    "location": "analysis/bifurcation.html#Specification-of-a-Model-1",
    "page": "Bifurcation Analysis",
    "title": "Specification of a Model",
    "category": "section",
    "text": "We will specify the model using a ParameterizedFunction:using ParameterizedFunctions\nf = @ode_def begin\n  dv = ( i + gl * (vl - v) - gca * 0.5 * (1 + tanh( (v-v1)/v2 )) * (v-vca) )/c\n  dw = v-w\nend vl vca i gl gca c v1 v2(Note that using PyDSTool requires use of the @ode_def macro). Next to build the ODE we need an initial condition and a starting timepoint.u0 = [0;0]\ntspan = [0;30]\np = [-60,120,0.0,2,4,20,-1.2,18]Then we use the following command to build the PyDSTool ODE:dsargs = build_ode(f,u0,tspan,p)Now we need to build the continuation type. Following the setup of PyDSTool\'s tutorial, we need to start near the steady state. The commands translate as:ode = ds[:Generator][:Vode_ODEsystem](dsargs)\node[:set](pars = Dict(\"i\"=>-220))\node[:set](ics  = Dict(\"v\"=>-170))\nPC = ds[:ContClass](ode)Once we have the continuation type, we can call the bifurcation_curve function. Instead of building the args into some object one-by-one, we simply make a function call with keyword arguments. Using the same arguments as the PyDSTool tutorial:bif = bifurcation_curve(PC,\"EP-C\",[\"i\"],\n                        max_num_points=450,\n                        max_stepsize=2,min_stepsize=1e-5,\n                        stepsize=2e-2,loc_bif_points=\"all\",\n                        save_eigen=true,name=\"EQ1\",\n                        print_info=true,calc_stab=true)This returns a BifurcationCurve type. Important fields of this type are:points: the values along the curve\nspecial_points: the values for the bifurcation points\nstab: an array which gives the stability of each point along the curve. \"S\" is for stable, N is for neutral, and U is for unstable.Instead of using the fields directly, we will use the plot recipe. The plot recipe requires you give the x,y coordinates to plot. Here we will plot it in the (i,v) plane:using Plots\nplot(bif,(:i,:v))(Image: bifurcation_plot)"
},

{
    "location": "analysis/bifurcation.html#Bifucation-Curve-Function-Definition-1",
    "page": "Bifurcation Analysis",
    "title": "Bifucation Curve Function Definition",
    "category": "section",
    "text": "function bifurcation_curve(PC,bif_type,freepars;max_num_points=450,\n                          max_stepsize=2,min_stepsize=1e-5,\n                          stepsize=2e-2,loc_bif_points=\"all\",\n                          save_eigen=true,name=\"DefaultName\",\n                          print_info=true,calc_stab=true,\n                          var_tol = 1e-6, func_tol = 1e-6,\n                          test_tol = 1e-4,\n                          initpoint=nothing,solver_sequence=[:forward])"
},

{
    "location": "analysis/sensitivity.html#",
    "page": "Local Sensitivity Analysis (Automatic Differentiation)",
    "title": "Local Sensitivity Analysis (Automatic Differentiation)",
    "category": "page",
    "text": ""
},

{
    "location": "analysis/sensitivity.html#Local-Sensitivity-Analysis-(Automatic-Differentiation)-1",
    "page": "Local Sensitivity Analysis (Automatic Differentiation)",
    "title": "Local Sensitivity Analysis (Automatic Differentiation)",
    "category": "section",
    "text": "Sensitivity analysis for ODE models is provided by the DiffEq suite. The model sensitivities are the derivatives of the solution with respect to the parameters. Specifically, the local sensitivity of the solution to a parameter is defined by how much the solution would change by changes in the parameter, i.e. the sensitivity of the ith independent variable to the jth parameter is fracpartial upartial p_j.Sensitivity analysis serves two major purposes. On one hand, the sensitivities are diagnostics of the model which are useful for understand how it will change in accordance to changes in the parameters. But another use is simply because in many cases these derivatives are useful. Sensitivity analysis provides a cheap way to calculate the gradient of the solution which can be used in parameter estimation and other optimization tasks.There are three types of sensitivity analysis. Local forward sensitivity analysis directly gives the gradient of the solution with respect to each parameter along the time series. The computational cost scales like N*M, where N is the number of states and M is the number of parameters. While this gives all of the information, it can be expensive for models with large numbers of parameters. Local adjoint sensitivity analysis solves directly for the gradient of some functional of the solution, such as a cost function or energy functional, in a manner that is cheaper when the number of parameters is large. Global Sensitivity Analysis methods are meant to be used for exploring the sensitivity over a larger domain without calculating derivatives and are covered on a different page."
},

{
    "location": "analysis/sensitivity.html#Installation-1",
    "page": "Local Sensitivity Analysis (Automatic Differentiation)",
    "title": "Installation",
    "category": "section",
    "text": "This functionality does not come standard with DifferentialEquations.jl. To use this functionality, you must install DiffEqSensitivty.jl:]add DiffEqSensitivity\nusing DiffEqSensitivity"
},

{
    "location": "analysis/sensitivity.html#Efficiency-of-the-Different-Methods-1",
    "page": "Local Sensitivity Analysis (Automatic Differentiation)",
    "title": "Efficiency of the Different Methods",
    "category": "section",
    "text": "For an analysis of which methods will be most efficient for computing the solution derivatives for a given problem, consult our analysis in this arxiv paper. A general rule of thumb is:Discrete Forward Sensitivity Analysis via ForwardDiff.jl is the fastest for ODEs with small numbers of parameters (<100)\nAdjoint senstivity analysis is the fastest when the number of parameters is sufficiently large. There are three configurations of note. Using backsolve is the fastest and uses the least memory, but is not guerenteed to be stable. Checkpointing is the slowest but uses O(1) memory and is stable. Interpolating is the second fastest, is stable, but requires the ability to hold the full forward solution and its interpolation in memory.\nThe methods which use automatic differentiation support the full range of DifferentialEquations.jl features (SDEs, DDEs, events, etc.), but only work on native Julia solvers. The methods which utilize altered ODE systems only work on ODEs (without events), but work on any ODE solver."
},

{
    "location": "analysis/sensitivity.html#Local-Forward-Sensitivity-Analysis-1",
    "page": "Local Sensitivity Analysis (Automatic Differentiation)",
    "title": "Local Forward Sensitivity Analysis",
    "category": "section",
    "text": "Local forward sensitivity analysis gives a solution along with a timeseries of the sensitivities along the solution."
},

{
    "location": "analysis/sensitivity.html#Discrete-Local-Forward-Sensitivity-Analysis-via-ForwardDiff.jl-1",
    "page": "Local Sensitivity Analysis (Automatic Differentiation)",
    "title": "Discrete Local Forward Sensitivity Analysis via ForwardDiff.jl",
    "category": "section",
    "text": "This method is the application of ForwardDiff.jl numbers to the ODE solver. This is done simply by making the u0 state vector a vector of Dual numbers, and multiple dispatch then allows the internals of the solver to propagate the derivatives along the solution."
},

{
    "location": "analysis/sensitivity.html#Examples-using-ForwardDiff.jl-1",
    "page": "Local Sensitivity Analysis (Automatic Differentiation)",
    "title": "Examples using ForwardDiff.jl",
    "category": "section",
    "text": "The easiest way to use ForwardDiff.jl for local forward sensitivity analysis is to simply put the ODE solve inside of a function which you would like to differentiate. For example, let\'s define the ODE system for the Lotka-Volterra equations:function f(du,u,p,t)\n  du[1] = dx = p[1]*u[1] - p[2]*u[1]*u[2]\n  du[2] = dy = -p[3]*u[2] + u[1]*u[2]\nend\n\np = [1.5,1.0,3.0,1.0]\nu0 = [1.0;1.0]\nprob = ODEProblem(f,u0,(0.0,10.0),p)Let\'s say we wanted to get the derivative of the final value w.r.t. each of the parameters. We can define the following function:function test_f(p)\n  _prob = remake(prob;u0=convert.(eltype(p),prob.u0),p=p)\n  solve(_prob,Vern9(),save_everystep=false)[end]\nendWhat this function does is use the remake function from the Problem Interface page to generate a new ODE problem with the new parameters, solves it, and returns the solution at the final time point. Notice that it takes care to make sure that the type of u0 matches the type of p. This is because ForwardDiff.jl will want to use Dual numbers, and thus to propagate the Duals throughout the solver we need to make sure the initial condition is also of the type of Dual number. On this function we can call ForwardDiff.jl and it will return the derivatives we wish to calculate:using ForwardDiff\nfd_res = ForwardDiff.jacobian(test_f,p)If we would like to get the solution and the value at the time point at the same time, we can use DiffResults.jl. For example, the following uses a single ODE solution to calculate the value at the end point and its parameter Jacobian:using DiffResults\nres = DiffResults.JacobianResult(u0,p) # Build the results object\nDiffResults.jacobian!(res,p) # Populate it with the results\nval = DiffResults.value(res) # This is the sol[end]\njac = DiffResults.jacobian(res) # This is dsol/dpIf we would like to get the time series, we can do so by seeding the dual numbers directly. To do this, we use the Dual constructor. The first value is the value of the parameter. The second is a tuple of the derivatives. For each value we want to take the derivative by, we seed a derivative with a 1 in a unique index. For example, we can build our parameter vector like:using ForwardDiff: Dual\nstruct MyTag end\np1dual = Dual{MyTag}(1.5, (1.0, 0.0, 0.0, 0.0))\np2dual = Dual{MyTag}(1.0, (0.0, 1.0, 0.0, 0.0))\np3dual = Dual{MyTag}(3.0, (0.0, 0.0, 1.0, 0.0))\np4dual = Dual{MyTag}(3.0, (0.0, 0.0, 0.0, 1.0))\npdual = [p1dual, p2dual, p3dual, p4dual]or equivalently using the seed_duals convenience function:function seed_duals(x::AbstractArray{V},::Type{T},\n                    ::ForwardDiff.Chunk{N} = ForwardDiff.Chunk(x)) where {V,T,N}\n  seeds = ForwardDiff.construct_seeds(ForwardDiff.Partials{N,V})\n  duals = [Dual{T}(x[i],seeds[i]) for i in eachindex(x)]\nend\npdual = seed_duals(p,MyTag)Next we need to make our initial condition Dual numbers so that these propogate through the solution. We can do this manually like:u0dual = [Dual{MyTag}(1.0, (0.0, 0.0)),Dual{MyTag}(1.0, (0.0, 0.0))]or use the same shorthand from before:u0dual = convert.(eltype(pdual),u0)Now we just use these Dual numbers to solve:prob_dual = ODEProblem(f,u0,tspan,pdual)\nsol_dual = solve(prob_dual,Tsit5(), saveat=0.2)The solution is now in terms of Dual numbers. We can extract the derivatives by looking at the partials of the duals in the solution. For example, sol[1,end] is the Dual number for the x component at the end of the integration, and so sol[1,end].partial[i] is dx(t_end)/dp_i."
},

{
    "location": "analysis/sensitivity.html#Local-Forward-Sensitivity-Analysis-via-ODELocalSensitivityProblem-1",
    "page": "Local Sensitivity Analysis (Automatic Differentiation)",
    "title": "Local Forward Sensitivity Analysis via ODELocalSensitivityProblem",
    "category": "section",
    "text": "For this method local sensitivity is computed using the sensitivity ODE:fracddtfracpartial upartial p_j=fracpartial fpartial ufracpartial upartial p_j+fracpartial fpartial p_j=Jcdot S_j+F_jwhereJ=left(beginarraycccc\nfracpartial f_1partial u_1  fracpartial f_1partial u_2  cdots  fracpartial f_1partial u_k\nfracpartial f_2partial u_1  fracpartial f_2partial u_2  cdots  fracpartial f_2partial u_k\ncdots  cdots  cdots  cdots\nfracpartial f_kpartial u_1  fracpartial f_kpartial u_2  cdots  fracpartial f_kpartial u_k\nendarrayright)is the Jacobian of the system,F_j=left(beginarrayc\nfracpartial f_1partial p_j\nfracpartial f_2partial p_j\nvdots\nfracpartial f_kpartial p_j\nendarrayright)are the parameter derivatives, andS_j=left(beginarrayc\nfracpartial y_1partial p_j\nfracpartial y_2partial p_j\nvdots\nfracpartial y_kpartial p_j\nendarrayright)is the vector of sensitivities. Since this ODE is dependent on the values of the independent variables themselves, this ODE is computed simultaneously with the actual ODE system."
},

{
    "location": "analysis/sensitivity.html#Example-solving-an-ODELocalSensitivityProblem-1",
    "page": "Local Sensitivity Analysis (Automatic Differentiation)",
    "title": "Example solving an ODELocalSensitivityProblem",
    "category": "section",
    "text": "To define a sensitivity problem, simply use the ODELocalSensitivityProblem type instead of an ODE type. For example, we generate an ODE with the sensitivity equations attached for the Lotka-Volterra equations by:function f(du,u,p,t)\n  du[1] = dx = p[1]*u[1] - p[2]*u[1]*u[2]\n  du[2] = dy = -p[3]*u[2] + u[1]*u[2]\nend\n\np = [1.5,1.0,3.0]\nprob = ODELocalSensitivityProblem(f,[1.0;1.0],(0.0,10.0),p)This generates a problem which the ODE solvers can solve:sol = solve(prob,DP8())Note that the solution is the standard ODE system and the sensitivity system combined. We can use the following helper functions to extract the sensitivity information:x,dp = extract_local_sensitivities(sol)\nx,dp = extract_local_sensitivities(sol,i)\nx,dp = extract_local_sensitivities(sol,t)In each case, x is the ODE values and dp is the matrix of sensitivities where dp[i] is the gradient of component i by the parameters. The first gives the full timeseries of values. The second returns the ith values, while the third interpolates to calculate the sensitivities at time t. For example, if we do:x,dp = extract_local_sensitivities(sol)\nda = dp[1]then da is the timeseries for fracpartial u(t)partial p. We can plot thisplot(sol.t,da\',lw=3)transposing so that the rows (the timeseries) is plotted.(Image: Local Sensitivity Solution)Here we see that there is a periodicity to the sensitivity which matches the periodicity of the Lotka-Volterra solutions. However, as time goes on the sensitivity increases. This matches the analysis of Wilkins in Sensitivity Analysis for Oscillating Dynamical Systems.We can also quickly see that these values are equivalent to those given by autodifferentiation and numerical differentiation through the ODE solver:using ForwardDiff, Calculus\nfunction test_f(p)\n  prob = ODEProblem(f,eltype(p).([1.0,1.0]),eltype(p).((0.0,10.0)),p)\n  solve(prob,Vern9(),abstol=1e-14,reltol=1e-14,save_everystep=false)[end]\nend\n\np = [1.5,1.0,3.0]\nfd_res = ForwardDiff.jacobian(test_f,p)\ncalc_res = Calculus.finite_difference_jacobian(test_f,p)Here we just checked the derivative at the end point."
},

{
    "location": "analysis/sensitivity.html#Internal-representation-1",
    "page": "Local Sensitivity Analysis (Automatic Differentiation)",
    "title": "Internal representation",
    "category": "section",
    "text": "For completeness, we detail the internal representation. Therefore, the solution to the ODE are the first n components of the solution. This means we can grab the matrix of solution values like:x = sol[1:sol.prob.indvars,:]Since each sensitivity is a vector of derivatives for each function, the sensitivities are each of size sol.prob.indvars. We can pull out the parameter sensitivities from the solution as follows:da = sol[sol.prob.indvars+1:sol.prob.indvars*2,:]\ndb = sol[sol.prob.indvars*2+1:sol.prob.indvars*3,:]\ndc = sol[sol.prob.indvars*3+1:sol.prob.indvars*4,:]This means that da[1,i] is the derivative of the x(t) by the parameter a at time sol.t[i]. Note that all of the functionality available to ODE solutions is available in this case, including interpolations and plot recipes (the recipes will plot the expanded system)."
},

{
    "location": "analysis/sensitivity.html#Adjoint-Sensitivity-Analysis-1",
    "page": "Local Sensitivity Analysis (Automatic Differentiation)",
    "title": "Adjoint Sensitivity Analysis",
    "category": "section",
    "text": "Adjoint sensitivity analysis is used to find the gradient of the solution with respect to some functional of the solution. In many cases this is used in an optimization problem to return the gradient with respect to some cost function. It is equivalent to \"backpropogation\" or reverse-mode automatic differentiation of a differential equation."
},

{
    "location": "analysis/sensitivity.html#Adjoint-Sensitivity-Analysis-via-adjoint_sensitivities-1",
    "page": "Local Sensitivity Analysis (Automatic Differentiation)",
    "title": "Adjoint Sensitivity Analysis via adjoint_sensitivities",
    "category": "section",
    "text": "This adjoint requires the definition of some scalar functional g(upt) where u is the (numerical) solution to the differential equation. Adjoint sensitivity analysis finds the gradient ofG(up)=G(u(p))=int_t_0^Tg(u(tp))dtsome integral of the solution. It does so by solving the adjoint problemfracdlambda^stardt=g_u(t)-lambda^star(t)f_u(t)thinspacethinspacethinspacelambda^star(T)=0where f_u is the Jacobian of the system with respect to the state u while f_p is the Jacobian with respect to the parameters. The adjoint problem\'s solution gives the sensitivities through the integral:fracdGdp=int_t_0^Tlambda^star(t)f_p(t)+g_p(t)dt+lambda^star(t_0)u_p(t_0)Notice that since the adjoints require the Jacobian of the system at the state, it requires the ability to evaluate the state at any point in time. Thus it requires the continuous forward solution in order to solve the adjoint solution, and the adjoint solution is required to be continuous in order to calculate the resulting integral.There is one extra detail to consider. In many cases we would like to calculate the adjoint sensitivity of some discontinuous functional of the solution. One canonical function is the L2 loss against some data points, that is:L(upt)=sum_i=1^nVerttildeu(t_i)-u(t_ip)Vert^2In this case, we can reinterpret our summation as the distribution integral:G(up)=int_0^Tsum_i=1^nVerttildeu(t_i)-u(t_ip)Vert^2delta(t_i-t)dtwhere δ is the Dirac distribution. In this case, the integral is continuous except at finitely many points. Thus it can be calculated between each t_i. At a given t_i, given that the t_i are unique, we have thatg_y(t_i)=2left(tildeu(t_i)-u(t_ip)right)Thus the adjoint solution is given by integrating between the integrals and applying the jump function g_y at every data point."
},

{
    "location": "analysis/sensitivity.html#Syntax-1",
    "page": "Local Sensitivity Analysis (Automatic Differentiation)",
    "title": "Syntax",
    "category": "section",
    "text": "There are two forms. For discrete adjoints, the form is:s = adjoint_sensitivities(sol,alg,dg,ts;kwargs...)where alg is the ODE algorithm to solve the adjoint problem, dg is the jump function, and ts is the time points for data. dg is given by:dg(out,u,p,t,i)which is the in-place gradient of the cost functional g at time point ts[i] with u=u(t).For continuous functionals, the form is:s = adjoint_sensitivities(sol,alg,g,nothing,dg;kwargs...)for the cost functionalg(u,p,t)with in-place gradientdg(out,u,p,t)If the gradient is omitted, i.e.s = adjoint_sensitivities(sol,alg,g,nothing;kwargs...)then it will be computed automatically using ForwardDiff or finite differencing, depending on the autodiff setting in the SensitivityAlg. Note that the keyword arguments are passed to the internal ODE solver for solving the adjoint problem. Two special keyword arguments are iabstol and ireltol which are the tolerances for the internal quadrature via QuadGK for the resulting functional."
},

{
    "location": "analysis/sensitivity.html#Options-1",
    "page": "Local Sensitivity Analysis (Automatic Differentiation)",
    "title": "Options",
    "category": "section",
    "text": "Options for handling the adjoint computation are set by passing a SensitivityAlg type, e.g. SensitivityAlg(backsolve=true). Additionally, if Gauss-Kronrod quadrature is used, the options ireltol and iabstol into adjoint_sensitivities controls the behavior of the quadrature. Example calls:res = adjoint_sensitivities(sol,Rodas4(),dg,t,ireltol=1e-8)\n\nres = adjoint_sensitivities(sol,Vern9(),dg,t,reltol=1e-8,\n                            sensealg=SensitivityAlg(backsolve=true))checkpointing: When enabled, the adjoint solutions compute the Jacobians by starting from the nearest saved value in sol and computing forward. By default, this is false if sol.dense==true, i.e. if sol has its higher order interpolation then this is by default disabled.\nquad: Use Gauss-Kronrod quadrature to integrate the adjoint sensitivity integral. Disabling it can decrease memory usage but increase computation time. Default is true.\nbacksolve: Solve the differential equation backward to get the past values. Note that for chaotic or non-reversible systems, enabling this option can lead to wildly incorrect results. Enabling it can decrease memory usage but increase computation time. When it is set to true, quad will be automatically set to false. Default is false.\nautodiff: Use automatic differentiation in the internal sensitivity algorithm computations. Default is true.\nchunk_size: Chunk size for forward mode differentiation. Default is 0 for automatic.\nautojacvec: Calculate Jacobian-vector (local sensitivity analysis) or vector-Jacobian (adjoint sensitivity analysis) product via automatic differentiation with special seeding. Default is true if autodiff is true."
},

{
    "location": "analysis/sensitivity.html#Example-discrete-adjoints-on-a-cost-function-1",
    "page": "Local Sensitivity Analysis (Automatic Differentiation)",
    "title": "Example discrete adjoints on a cost function",
    "category": "section",
    "text": "In this example we will show solving for the adjoint sensitivities of a discrete cost functional. First let\'s solve the ODE and get a high quality continuous solution:function f(du,u,p,t)\n  du[1] = dx = p[1]*u[1] - p[2]*u[1]*u[2]\n  du[2] = dy = -p[3]*u[2] + u[1]*u[2]\nend\n\np = [1.5,1.0,3.0]\nprob = ODEProblem(f,[1.0;1.0],(0.0,10.0),p)\nsol = solve(prob,Vern9(),abstol=1e-10,reltol=1e-10)Now let\'s calculate the sensitivity of the L2 error against 1 at evenly spaced points in time, that is:L(upt)=sum_i=1^nfracVert1-u(t_ip)Vert^22for t_i = 05i. This is the assumption that the data is data[i]=1.0. For this function, notice we have that:beginalign\ndg_1=1-u_1 \ndg_2=1-u_2\nendalignand thus:dg(out,u,i) = (out.=1.0.-u)If we had data, we\'d just replace 1.0 with data[i]. To get the adjoint sensitivities, call:res = adjoint_sensitivities(sol,Vern9(),dg,t,abstol=1e-14,\n                            reltol=1e-14,iabstol=1e-14,ireltol=1e-12)This is super high accuracy. As always, there\'s a tradeoff between accuracy and computation time. We can check this almost exactly matches the autodifferentiation and numerical differentiation results:using ForwardDiff,Calculus\nfunction G(p)\n  tmp_prob = remake(prob,u0=convert.(eltype(p),prob.u0),p=p)\n  sol = solve(tmp_prob,Vern9(),abstol=1e-14,reltol=1e-14,saveat=t)\n  A = convert(Array,sol)\n  sum(((1-A).^2)./2)\nend\nG([1.5,1.0,3.0])\nres2 = ForwardDiff.gradient(G,[1.5,1.0,3.0])\nres3 = Calculus.gradient(G,[1.5,1.0,3.0])\nres4 = Flux.Tracker.gradient(G,[1.5,1.0,3.0])\nres5 = ReverseDiff.gradient(G,[1.5,1.0,3.0])and see this gives the same values."
},

{
    "location": "analysis/sensitivity.html#Example-controlling-adjoint-method-choices-and-checkpointing-1",
    "page": "Local Sensitivity Analysis (Automatic Differentiation)",
    "title": "Example controlling adjoint method choices and checkpointing",
    "category": "section",
    "text": "In the previous examples, all calculations were done using the interpolating method. This maximizes speed but at a cost of requiring a dense sol. If it is not possible to hold a dense forward solution in memory, then one can use checkpointing. This is enabled by default if sol is not dense, so for examplesol = solve(prob,Vern9(),saveat=[0.0,0.2,0.5,0.7])Creates a non-dense solution with checkpoints at [0.0,0.2,0.5,0.7]. Now we can dores = adjoint_sensitivities(sol,Vern9(),dg,t)When grabbing a Jacobian value during the backwards solution, it will no longer interpolate to get the value. Instead, it will start a forward solution at the nearest checkpoint and solve until the necessary time.To eliminate the extra forward solutions, one can instead pass the SensitivityAlg with the backsolve=true option:sol = solve(prob,Vern9(),save_everystep=false,save_start=false)\nres = adjoint_sensitivities(sol,Vern9(),dg,t,sensealg=SensitivityAlg(backsolve=true))When this is done, the values for the Jacobian will be computing the original ODE in reverse. Note that this only requires the final value of the solution."
},

{
    "location": "analysis/sensitivity.html#Applicability-of-Backsolve-and-Caution-1",
    "page": "Local Sensitivity Analysis (Automatic Differentiation)",
    "title": "Applicability of Backsolve and Caution",
    "category": "section",
    "text": "When backsolve is applicable it is the fastest method and requires the least memory. However, one must be cautious because not all ODEs are stable under backwards integration by the majority of ODE solvers. An example of such an equation is the Lorenz equation. Notice that if one solves the Lorenz equation forward and then in reverse with any adaptive time step and non-reversible integrator, then the backwards solution diverges from the forward solution. As a quick demonstration:using Sundials, DiffEqBase\nfunction lorenz(du,u,p,t)\n du[1] = 10.0*(u[2]-u[1])\n du[2] = u[1]*(28.0-u[3]) - u[2]\n du[3] = u[1]*u[2] - (8/3)*u[3]\nend\nu0 = [1.0;0.0;0.0]\ntspan = (0.0,100.0)\nprob = ODEProblem(lorenz,u0,tspan)\nsol = solve(prob,Tsit5(),reltol=1e-12,abstol=1e-12)\nprob2 = ODEProblem(lorenz,sol[end],(100.0,0.0))\nsol = solve(prob,Tsit5(),reltol=1e-12,abstol=1e-12)\n@show sol[end]-u0 #[-3.22091, -1.49394, 21.3435]Thus one should check the stability of the backsolve on their type of problem before enabling this method."
},

{
    "location": "analysis/sensitivity.html#Example-continuous-adjoints-on-an-energy-functional-1",
    "page": "Local Sensitivity Analysis (Automatic Differentiation)",
    "title": "Example continuous adjoints on an energy functional",
    "category": "section",
    "text": "In this case we\'d like to calculate the adjoint sensitivity of the scalar energy functionalG(up)=int_0^Tfracsum_i=1^nu_i^2(t)2dtwhich isg(u,p,t) = (sum(u).^2) ./ 2Notice that the gradient of this function with respect to the state u is:function dg(out,u,p,t)\n  out[1]= u[1] + u[2]\n  out[2]= u[1] + u[2]\nendTo get the adjoint sensitivities, we call:res = adjoint_sensitivities(sol,Vern9(),g,nothing,dg,abstol=1e-8,\n                                 reltol=1e-8,iabstol=1e-8,ireltol=1e-8)Notice that we can check this against autodifferentiation and numerical differentiation as follows:function G(p)\n  tmp_prob = remake(prob,p=p)\n  sol = solve(tmp_prob,Vern9(),abstol=1e-14,reltol=1e-14)\n  res,err = quadgk((t)-> (sum(sol(t)).^2)./2,0.0,10.0,abstol=1e-14,reltol=1e-10)\n  res\nend\nres2 = ForwardDiff.gradient(G,[1.5,1.0,3.0])\nres3 = Calculus.gradient(G,[1.5,1.0,3.0])"
},

{
    "location": "analysis/global_sensitivity.html#",
    "page": "Global Sensitivity Analysis",
    "title": "Global Sensitivity Analysis",
    "category": "page",
    "text": ""
},

{
    "location": "analysis/global_sensitivity.html#Global-Sensitivity-Analysis-1",
    "page": "Global Sensitivity Analysis",
    "title": "Global Sensitivity Analysis",
    "category": "section",
    "text": "Global Sensitivity Analysis (GSA) methods are used to quantify the uncertainty in output of a model w.r.t. the parameters, their individual contributions, or the contribution of their interactions. The type of GSA method to use depends on the interest of the user, below we describe the methods available in the suite at the moment (some more are already in development) and explain what is the output of each of the methods and what it represents."
},

{
    "location": "analysis/global_sensitivity.html#Installation-1",
    "page": "Global Sensitivity Analysis",
    "title": "Installation",
    "category": "section",
    "text": "This functionality does not come standard with DifferentialEquations.jl. To use this functionality, you must install DiffEqSensitivty.jl:]add DiffEqSensitivity\nusing DiffEqSensitivty"
},

{
    "location": "analysis/global_sensitivity.html#Morris-Method-1",
    "page": "Global Sensitivity Analysis",
    "title": "Morris Method",
    "category": "section",
    "text": "The Morris method also known as Morris’s OAT method where OAT stands for One At a Time can be described in the following steps:We calculate local sensitivity measures known as “elementary effects”, which are calculated by measuring the perturbation in the output of the model on changing one parameter.EE_i = fracf(x_1x_2x_i+ Deltax_k) - yDeltaThese are evaluated at various points in the input chosen such that a wide “spread” of the parameter space is explored and considered in the analysis, to provide an approximate global importance measure. The mean and variance of these elementary effects is computed. A high value of the mean implies that a parameter is important, a high variance implies that its effects are non-linear or the result of interactions with other inputs. This method does not evaluate separately the contribution from the interaction and the contribution of the parameters individually and gives the effects for each parameter which takes into cpnsideration all the interactions and its individual contribution.morris_effects = morris_sensitivity(f,param_range,param_steps;relative_scale=false,kwargs...)morris_effects = morris_sensitivity(prob::DiffEqBase.DEProblem,alg,t,param_range,param_steps;kwargs...)Here, f is just the model (as a julia function  f(input_vector) -> output_vector or a DEProblem) you want to run the analysis on, param_range requires an array of 2-tuples with the lower bound and the upper bound, param_steps decides the value of Delta in the equation above and relative_scale, the above equation takes the assumption that the parameters lie in the range [0,1] but as this is not always the case scaling is used to get more informative, scaled effects."
},

{
    "location": "analysis/global_sensitivity.html#Sobol-Method-1",
    "page": "Global Sensitivity Analysis",
    "title": "Sobol Method",
    "category": "section",
    "text": "Sobol is a variance-based method and it decomposes the variance of the output of the model or system into fractions which can be attributed to inputs or sets of inputs. This helps to get not just the individual parameter\'s sensitivities but also gives a way to quantify the affect and sensitivity from the interaction between the parameters. Y = f_0+ sum_i=1^d f_i(X_i)+ sum_i  j^d f_ij(X_iX_j)  + f_12d(X_1X_2X_d) Var(Y) = sum_i=1^d V_i + sum_i  j^d V_ij +  + V_12dThe Sobol Indices are \"order\"ed, the first order indices given by S_i = fracV_iVar(Y) the contribution to the output variance of the main effect of X_i, therefore it measures the effect of varying X_i alone, but averaged over variations in other input parameters. It is standardised by the total variance to provide a fractional contribution. Higher-order interaction indices S_ij S_ijk and so on can be formed by dividing other terms in the variance decomposition by Var(Y).sobol_second_order = sobol_sensitivity(f,param_range,N,order=2)sobol_second_order = sobol_sensitivity(prob::DiffEqBase.DEProblem,alg,t,param_range,N,order=2)Here f and param_range are the same as Morris\'s, providing a uniform interface."
},

{
    "location": "analysis/global_sensitivity.html#Regression-Method-1",
    "page": "Global Sensitivity Analysis",
    "title": "Regression Method",
    "category": "section",
    "text": "If a sample of inputs and outputs (X^n Y^n) = (X^i_1     X^i_d Y_i)_i=1n􏰁 is available, it is possible to fit a linear model explaining the behavior of Y given the values of X, provided that the sample size n is sufficiently large (at least n > d).The measures provided for this analysis by us in DiffEqSensitivity.jl area) Pearson Correlation Coefficient:r = fracsum_i=1^n (x_i - overlinex)(y_i - overliney)sqrtsum_i=1^n (x_i - overlinex)^2(y_i - overliney)^2b) Standard Regression Coefficient (SRC):SRC_j = beta_j sqrtfracVar(X_j)Var(Y)where beta_j is the linear regression coefficient associated to X_j.c) Partial Correlation Coefficient (PCC):PCC_j = rho(X_j - hatX_-jY_j - hatY_-j)where hatX_-j􏰈 is the prediction of the linear model, expressing X_j   with respect to the other inputs and hatY_-j is the prediction of the   linear model where X_j is absent. PCC measures the sensitivity of Y to   X_j when the effects of the other inputs have been canceled.regre_sensitivity = regression_sensitivity(f,param_range,param_fixed,n;coeffs=:rank)regre_sensitivity = regression_sensitivity(prob::DiffEqBase.DEProblem,alg,t,param_range,param_fixed,n;coeffs=:rank)Again, f and param_range are the same as above. An array of the true parameter values that lie within the param_range bounds are passed through the param_fixed argument. n determines the number of simulations of the model run to generate the data points of the solution and parameter values and the coeffs kwarg lets you decide the coefficients you want."
},

{
    "location": "analysis/global_sensitivity.html#GSA-example-1",
    "page": "Global Sensitivity Analysis",
    "title": "GSA example",
    "category": "section",
    "text": "Let\'s create the ODE problem to run our GSA on.function f(du,u,p,t)\n  du[1] = p[1]*u[1] - p[2]*u[1]*u[2]\n  du[2] = -3*u[2] + u[1]*u[2]\nend\nu0 = [1.0;1.0]\ntspan = (0.0,10.0)\np = [1.5,1.0]\nprob = ODEProblem(f,u0,tspan,p)\nt = collect(range(0, stop=10, length=200))For Morris Methodm = DiffEqSensitivity.morris_sensitivity(prob,Tsit5(),t,[[1,5],[0.5,5]],[10,10],len_trajectory=1500,total_num_trajectory=1000,num_trajectory=150)Let\'s get the means and variances from the MorrisSensitivity struct.m.means\n\nOut[9]: 2-element Array{Array{Float64,2},1}:\n [0.0 0.0513678 … 7.91336 7.93783; 0.0 0.00115769 … 3.66156 3.67284]\n [0.0 0.0488899 … 2.50728 2.359; 0.0 0.00112006 … 2.23431 2.44946]\n\nm.variances\n\nOut[10]: 2-element Array{Array{Float64,2},1}:\n [0.0 1.94672e-5 … 26.4223 24.8513; 0.0 4.81347e-9 … 37.4061 30.3068]\n [0.0 1.77615e-5 … 17.9555 14.9231; 0.0 4.47931e-9 … 48.074 51.9312]This gives the means of the effects and it\'s variances over the entire timespan and thus we get 200-length arrays for each paramter and dependent variable pair.We can plot the trajectory of the sensitivity with the standard deviation bars.# For the first parameter (a)\nstdv1 = sqrt.(m.variances[1])\np = plot(m.means[1]\', yerror=stdv1)(Image: morrisparameter1)# For the second parameter (b)\nstdv2 = sqrt.(m.variances[2])\np = plot(m.means[2]\', yerror=stdv2)(Image: morrisparameter2)For Sobol Method\ns0 = sobol_sensitivity(prob,Tsit5(),t,[[1,5],[0.5,5]],N,0)\nOut[8]: 2-element Array{Array{Float64,2},1}:\n [NaN 0.507831 … 1.00731 1.00436; NaN 1.92336 … 0.732384 0.730945]\n [NaN 0.47214 … 0.676224 0.681525; NaN -1.68656 … 0.879557 0.877603]\n\ns1 = sobol_sensitivity(prob,Tsit5(),t,[[1,5],[0.5,5]],N,1)\nOut[9]: 2-element Array{Array{Float64,2},1}:\n [NaN 0.39537 … 0.341697 0.343645; NaN -2.06101 … 0.10922 0.106976]\n [NaN 0.652815 … 0.00910675 0.00815206; NaN 5.24832 … 0.296978 0.296639]\n\ns2 = sobol_sensitivity(prob,Tsit5(),t,[[1,5],[0.5,5]],N,2)\nOut[10]: 1-element Array{Array{Float64,2},1}:\n [NaN -0.0596478 … 0.652303 0.657847; NaN -1.84504 … 0.645139 0.620036]We can decide which order of Sobol Indices we are interested in my passing an argument for it, by default it gives the second order indices. Again the result is obtained over the entire timespanWe plot the first order and total order Sobol Indices for some timepoints for each of the parameters (a and b).\np1 = bar([\"a\",\"b\"],[s0[1][end-2],s0[2][end-2]],color=[:red,:blue],title=\"Total Order Indices at t=9.949748743718592\",legend=false)\np2 = bar([\"a\",\"b\"],[s1[1][end-2],s1[2][end-2]],color=[:red,:blue],title=\"First Order Indices at t=9.949748743718592\",legend=false)\np3 = bar([\"a\",\"b\"],[s0[1][3],s0[2][3]],color=[:red,:blue],title=\"Total Order Indices at t=0.05025125628140704\",legend=false)\np4 = bar([\"a\",\"b\"],[s1[1][3],s1[2][3]],color=[:red,:blue],title=\"First Order Indices at t=0.05025125628140704\",legend=false)\nplo = plot(p1,p2,p3,p4,layout=(4,1),size=(600,500))\n(Image: sobolplot)Here we plot the Sobol indices of first order and the total Sobol indices for the parameters a and b. The plots are obtained by getting the Sobol Indices at the t = 9.949748743718592 and the t = 0.05025125628140704 time point of the first dependent variable x(t) from the 200-length sensitivities over the entire time span. The length of the bar represents the quantification of the sensitivity of the output to that parameter and hence for the 199th time point you can say that x(t) is more sensitive to b, also you can observe how the relative difference between a and b is larger in the first order than the total order indices, this tells us that most of the contribution of a to x(t) arises from interactions and it\'s individual non-interaction contribution is significantly lesser than b and vice-versa for b as it\'s first order plot indicates quite high value."
},

{
    "location": "analysis/uncertainty_quantification.html#",
    "page": "Uncertainty Quantification",
    "title": "Uncertainty Quantification",
    "category": "page",
    "text": ""
},

{
    "location": "analysis/uncertainty_quantification.html#Uncertainty-Quantification-1",
    "page": "Uncertainty Quantification",
    "title": "Uncertainty Quantification",
    "category": "section",
    "text": "Uncertainty quantification allows a user to identify the uncertainty associated with the numerical approximation given by DifferentialEquations.jl. This page describes the different methods available for quantifying such uncertainties. Note that this requires one of the native Julia solvers like OrdinaryDiffEq.jl, StochasticDiffEq.jl, or DelayDiffEq.jl."
},

{
    "location": "analysis/uncertainty_quantification.html#Installation-1",
    "page": "Uncertainty Quantification",
    "title": "Installation",
    "category": "section",
    "text": "This functionality does not come standard with DifferentialEquations.jl. To use this functionality, you must install DiffEqUncertainty.jl:]add DiffEqUncertainty\nusing DiffEqUncertainty"
},

{
    "location": "analysis/uncertainty_quantification.html#ProbInts-1",
    "page": "Uncertainty Quantification",
    "title": "ProbInts",
    "category": "section",
    "text": "The ProbInts method for uncertainty quantification involves the transformation of an ODE into an associated SDE where the noise is related to the timesteps and the order of the algorithm. This is implemented into the DiffEq system via a callback function. The first form is:ProbIntsUncertainty(σ,order,save=true)σ is the noise scaling factor and order is the order of the algorithm. save is for choosing whether this callback should control the saving behavior. Generally this is true unless one is stacking callbacks in a CallbackSet. It is recommended that σ is representative of the size of the errors in a single step of the equation.If you are using an adaptive algorithm, the callbackAdaptiveProbIntsUncertainty(order,save=true)determines the noise scaling automatically using an internal error estimate."
},

{
    "location": "analysis/uncertainty_quantification.html#Example-1:-FitzHugh-Nagumo-1",
    "page": "Uncertainty Quantification",
    "title": "Example 1: FitzHugh-Nagumo",
    "category": "section",
    "text": "In this example we will determine our uncertainty when solving the FitzHugh-Nagumo model with the Euler() method. We define the FitzHugh-Nagumo model:function fitz(du,u,p,t)\n  V,R = u\n  a,b,c = p\n  du[1] = c*(V - V^3/3 + R)\n  du[2] = -(1/c)*(V -  a - b*R)\nend\nu0 = [-1.0;1.0]\ntspan = (0.0,20.0)\np = (0.2,0.2,3.0)\nprob = ODEProblem(fitz,u0,tspan,p)Now we define the ProbInts callback. In this case, our method is the Euler method and thus it is order 1. For the noise scaling, we will try a few different values and see how it changes. For σ=0.2, we define the callback as:cb = ProbIntsUncertainty(0.2,1)This is akin to having an error of approximately 0.2 at each step. We now build and solve a EnsembleProblem for 100 trajectories:ensemble_prob = EnsembleProblem(prob)\nsim = solve(ensemble_prob,Euler(),trajectories=100,callback=cb,dt=1/10)Now we can plot the resulting Monte Carlo solution:using Plots; plotly(); plot(sim,vars=(0,1),linealpha=0.4)(Image: uncertainty_02)If we increase the amount of error, we see that some parts of the equation have less uncertainty than others. For example, at σ=0.5:cb = ProbIntsUncertainty(0.5,1)\nensemble_prob = EnsembleProblem(prob)\nsim = solve(ensemble_prob,Euler(),trajectories=100,callback=cb,dt=1/10)\nusing Plots; plotly(); plot(sim,vars=(0,1),linealpha=0.4)(Image: uncertainty_05)But at this amount of noise, we can see how we contract to the true solution by decreasing dt:cb = ProbIntsUncertainty(0.5,1)\nensemble_prob = EnsembleProblem(prob)\nsim = solve(ensemble_prob,Euler(),trajectories=100,callback=cb,dt=1/100)\nusing Plots; plotly(); plot(sim,vars=(0,1),linealpha=0.4)(Image: uncertainty_lowh)"
},

{
    "location": "analysis/uncertainty_quantification.html#Example-2:-Adaptive-ProbInts-on-FitzHugh-Nagumo-1",
    "page": "Uncertainty Quantification",
    "title": "Example 2: Adaptive ProbInts on FitzHugh-Nagumo",
    "category": "section",
    "text": "While the first example is academic and shows how the ProbInts method scales, the fact that one should have some idea of the error in order to calibrate σ can lead to complications. Thus the more useful method in many cases is the AdaptiveProbIntsUncertainty version. In this version, no σ is required since this is calculated using an internal error estimate. Thus this gives an accurate representation of the possible error without user input.Let\'s try this with the order 5 Tsit5() method on the same problem as before:cb = AdaptiveProbIntsUncertainty(5)\nsol = solve(prob,Tsit5())\nensemble_prob = EnsembleProblem(prob)\nsim = solve(ensemble_prob,Tsit5(),trajectories=100,callback=cb)\nusing Plots; plotly(); plot(sim,vars=(0,1),linealpha=0.4)(Image: uncertainty_adaptive_default)In this case, we see that the default tolerances give us a very good solution. However, if we increase the tolerance a lot:cb = AdaptiveProbIntsUncertainty(5)\nsol = solve(prob,Tsit5())\nensemble_prob = EnsembleProblem(prob)\nsim = solve(ensemble_prob,Tsit5(),trajectories=100,callback=cb,abstol=1e-3,reltol=1e-1)\nusing Plots; plotly(); plot(sim,vars=(0,1),linealpha=0.4)(Image: uncertainty_adaptive_default)we can see that the moments just after the rise can be uncertain."
},

{
    "location": "analysis/uncertainty_quantification.html#Example-3:-Adaptive-ProbInts-on-the-Lorenz-Attractor-1",
    "page": "Uncertainty Quantification",
    "title": "Example 3: Adaptive ProbInts on the Lorenz Attractor",
    "category": "section",
    "text": "One very good use of uncertainty quantification is on chaotic models. Chaotic equations diverge from the true solution according to the error exponentially. This means that as time goes on, you get further and further from the solution. The ProbInts method can help diagnose how much of the timeseries is reliable.As in the previous example, we first define the model:function g(du,u,p,t)\n du[1] = p[1]*(u[2]-u[1])\n du[2] = u[1]*(p[2]-u[3]) - u[2]\n du[3] = u[1]*u[2] - p[3]*u[3]\nend\nu0 = [1.0;0.0;0.0]\ntspan = (0.0,30.0)\np = [10.0,28.0,8/3]\nprob = ODEProblem(g,u0,tspan,p)and then we build the ProbInts type. Let\'s use the order 5 Tsit5 again.cb = AdaptiveProbIntsUncertainty(5)Then we solve the MonteCarloProblemensemble_prob = EnsembleProblem(prob)\nsim = solve(ensemble_prob,Tsit5(),trajectories=100,callback=cb)\nusing Plots; plotly(); plot(sim,vars=(0,1),linealpha=0.4)(Image: uncertainty_chaos)Here we see that by t about 22 we start to receive strong deviations from the \"true\" solution. We can increase the amount of time before error explosion by using a higher order method with stricter tolerances:tspan = (0.0,40.0)\nprob = ODEProblem(g,u0,tspan,p)\ncb = AdaptiveProbIntsUncertainty(7)\nensemble_prob = EnsembleProblem(prob)\nsim = solve(ensemble_prob,Vern7(),trajectories=100,callback=cb,reltol=1e-6)\nusing Plots; plotly(); plot(sim,vars=(0,1),linealpha=0.4)(Image: uncertainty_high_order)we see that we can extend the amount of time until we deviate strongly from the \"true\" solution. Of course, for a chaotic system like the Lorenz one presented here, it is impossible to follow the true solution for long times, due to the fact that the system is chaotic and unavoidable deviations due to the numerical precision of a cumputer get amplified exponentially.However, not all hope is lost. The shadowing theorem is a strong statement for having confidence in numerical evolution of chaotic systems:Although a numerically computed chaotic trajectory diverges exponentially from the true trajectory with the same initial coordinates, there exists an errorless trajectory with a slightly different initial condition that stays near (\"shadows\") the numerically computed one.For more info on the shadowing theorem, please see the book Chaos in Dynamical Systems by E. Ott."
},

{
    "location": "analysis/neural_networks.html#",
    "page": "Neural Networks",
    "title": "Neural Networks",
    "category": "page",
    "text": ""
},

{
    "location": "analysis/neural_networks.html#Neural-Networks-1",
    "page": "Neural Networks",
    "title": "Neural Networks",
    "category": "section",
    "text": "To use DifferentialEquations.jl with the Flux.jl neural network package, consult the documentation at DiffEqFlux.jl."
},

{
    "location": "analysis/dev_and_test.html#",
    "page": "Algorithm Development and Testing",
    "title": "Algorithm Development and Testing",
    "category": "page",
    "text": ""
},

{
    "location": "analysis/dev_and_test.html#Algorithm-Development-and-Testing-1",
    "page": "Algorithm Development and Testing",
    "title": "Algorithm Development and Testing",
    "category": "section",
    "text": "Algorithm developing and testing tools are provided by DiffEqDevTools.jl and are documented in the developer documentation."
},

{
    "location": "analysis/dev_and_test.html#Installation-1",
    "page": "Algorithm Development and Testing",
    "title": "Installation",
    "category": "section",
    "text": "This functionality does not come standard with DifferentialEquations.jl. To use this functionality, you must install DiffEqDevTools.jl:]add DiffEqDevTools\nusing DiffEqDevTools"
},

{
    "location": "models/multiscale.html#",
    "page": "Multi-Scale Models",
    "title": "Multi-Scale Models",
    "category": "page",
    "text": ""
},

{
    "location": "models/multiscale.html#Multi-Scale-Models-1",
    "page": "Multi-Scale Models",
    "title": "Multi-Scale Models",
    "category": "section",
    "text": "The multi-scale modeling functionality is provided by MultiScaleArrays.jl. It allows for designing a multi-scale model as an extension of an array, which in turn can be directly used in the native Julia solvers of DifferentialEquations.jl."
},

{
    "location": "models/multiscale.html#More-Information-1",
    "page": "Multi-Scale Models",
    "title": "More Information",
    "category": "section",
    "text": "For more information, please see the MultiScaleArrays.jl README."
},

{
    "location": "models/physical.html#",
    "page": "Physical Models",
    "title": "Physical Models",
    "category": "page",
    "text": ""
},

{
    "location": "models/physical.html#Physical-Models-1",
    "page": "Physical Models",
    "title": "Physical Models",
    "category": "section",
    "text": "The physical modeling functionality is provided by DiffEqPhysics.jl and helps the user build and solve the differential equation based physical models."
},

{
    "location": "models/physical.html#Hamiltonian-Problems-1",
    "page": "Physical Models",
    "title": "Hamiltonian Problems",
    "category": "section",
    "text": "ODEs defined by Hamiltonians is described in the Dynamical ODEs section."
},

{
    "location": "models/physical.html#N-Body-Problems-1",
    "page": "Physical Models",
    "title": "N-Body Problems",
    "category": "section",
    "text": "nprob = NBodyProblem(f, mass, vel, pos, tspan)where f is the potential function, mass is the mass matrix, pos and vel are ArrayPartitions for the intial positions and velocity, and tspan is the timespan to solve on."
},

{
    "location": "models/physical.html#Example-1",
    "page": "Physical Models",
    "title": "Example",
    "category": "section",
    "text": "In this example we will model the outer solar system planets.G = 2.95912208286e-4\nM = [1.00000597682, 0.000954786104043, 0.000285583733151, 0.0000437273164546, 0.0000517759138449, 1/1.3e8]\ninvM = inv.(M)\nplanets = [\"Sun\", \"Jupiter\", \"Saturn\", \"Uranus\", \"Neptune\", \"Pluto\"]\n\npos_x = [0.0,-3.5023653,9.0755314,8.3101420,11.4707666,-15.5387357]\npos_y = [0.0,-3.8169847,-3.0458353,-16.2901086,-25.7294829,-25.2225594]\npos_z = [0.0,-1.5507963,-1.6483708,-7.2521278,-10.8169456,-3.1902382]\npos = ArrayPartition(pos_x,pos_y,pos_z)\n\nvel_x = [0.0,0.00565429,0.00168318,0.00354178,0.00288930,0.00276725]\nvel_y = [0.0,-0.00412490,0.00483525,0.00137102,0.00114527,-0.00170702]\nvel_z = [0.0,-0.00190589,0.00192462,0.00055029,0.00039677,-0.00136504]\nvel = ArrayPartition(vel_x,vel_y,vel_z)\n\ntspan = (0.,200_000)\n\nconst ∑ = sum\nconst N = 6\npotential(p, t, x, y, z, M) = -G*∑(i->∑(j->(M[i]*M[j])/sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2 + (z[i]-z[j])^2), 1:i-1), 2:N)\nnprob = NBodyProblem(potential, M, vel, pos, tspan)\nsol = solve(nprob,Yoshida6(), dt=100)"
},

{
    "location": "models/financial.html#",
    "page": "Financial Models",
    "title": "Financial Models",
    "category": "page",
    "text": ""
},

{
    "location": "models/financial.html#Financial-Models-1",
    "page": "Financial Models",
    "title": "Financial Models",
    "category": "section",
    "text": "The financial models functionality is provided by DiffEqFinancial.jl and helps the user build and solve the differential equation based financial models."
},

{
    "location": "models/financial.html#SDE-Model-Library-1",
    "page": "Financial Models",
    "title": "SDE Model Library",
    "category": "section",
    "text": "The following constructors create SDEProblem types which can be solved using the stochastic differential equation solvers."
},

{
    "location": "models/financial.html#HestonProblem-1",
    "page": "Financial Models",
    "title": "HestonProblem",
    "category": "section",
    "text": "dS = μSdt + sqrtvSdW_1 \ndv = κ(Θ-v)dt + σsqrtvdW_2 \ndW_1 dW_2 = ρ dtConstructor:HestonProblem(μ,κ,Θ,σ,ρ,u0,tspan)"
},

{
    "location": "models/financial.html#GeneralizedBlackScholesProblem-1",
    "page": "Financial Models",
    "title": "GeneralizedBlackScholesProblem",
    "category": "section",
    "text": "d ln S(t) = (r(t) - q(t) - fracΘ(tS)^22)dt + σ dW_tSolves for log S(t). Constructor:GeneralizedBlackScholesProblem(r,q,Θ,σ,u0,tspan)"
},

{
    "location": "models/financial.html#BlackScholesProblem-1",
    "page": "Financial Models",
    "title": "BlackScholesProblem",
    "category": "section",
    "text": "d ln S(t) = (r(t) - fracΘ(tS)^22)dt + σ dW_tSolves for log S(t). Constructor:BlackScholesProblem(r,Θ,σ,u0,tspan)"
},

{
    "location": "models/financial.html#ExtendedOrnsteinUhlenbeckProblem-1",
    "page": "Financial Models",
    "title": "ExtendedOrnsteinUhlenbeckProblem",
    "category": "section",
    "text": "dx = a(b(t)-x)dt + σ dW_tConstructor:ExtendedOrnsteinUhlenbeckProblem(a,b,σ,u0,tspan)"
},

{
    "location": "models/financial.html#OrnsteinUhlenbeckProblem-1",
    "page": "Financial Models",
    "title": "OrnsteinUhlenbeckProblem",
    "category": "section",
    "text": "dx = a(r-x)dt + σ dW_tConstructor:OrnsteinUhlenbeckProblem(a,r,σ,u0,tspan)"
},

{
    "location": "models/biological.html#",
    "page": "Chemical Reaction Models",
    "title": "Chemical Reaction Models",
    "category": "page",
    "text": ""
},

{
    "location": "models/biological.html#Chemical-Reaction-Models-1",
    "page": "Chemical Reaction Models",
    "title": "Chemical Reaction Models",
    "category": "section",
    "text": "The biological models functionality is provided by DiffEqBiological.jl and helps the user to build discrete stochastic and differential equation based systems biological models. These tools allow one to define the models at a high level by specifying reactions and rate constants, and the creation of the actual problems is then handled by the modelling package."
},

{
    "location": "models/biological.html#Installation-1",
    "page": "Chemical Reaction Models",
    "title": "Installation",
    "category": "section",
    "text": "This functionality does not come standard with DifferentialEquations.jl. To use this functionality, you must install DiffEqBiological.jl:]add DiffEqBiological\nusing DiffEqBiological"
},

{
    "location": "models/biological.html#The-Reaction-DSL-Basic-1",
    "page": "Chemical Reaction Models",
    "title": "The Reaction DSL - Basic",
    "category": "section",
    "text": "This section covers some of the basic syntax for building chemical reaction network models. Examples showing how to both construct and solve network models are provided in Chemical Reaction Network Examples."
},

{
    "location": "models/biological.html#Basic-syntax-1",
    "page": "Chemical Reaction Models",
    "title": "Basic syntax",
    "category": "section",
    "text": "The @reaction_network macro allows the (symbolic) specification of reaction networks with a simple format. Its input is a set of chemical reactions, and from them it generates a reaction network object which can be used as input to ODEProblem, SteadyStateProblem, SDEProblem and JumpProblem constructors. The @min_reaction_network macro constructs a more simplified reaction network, deferring construction of all the various functions needed for each of these problem types. It can then be incrementally filled in for specific problem types as needed, which reduces network construction time for very large networks (see The Min Reaction Network Object for a detailed description).The basic syntax is:rn = @reaction_network begin\n  2.0, X + Y --> XY               \n  1.0, XY --> Z1 + Z2            \nendwhere each line corresponds to a chemical reaction. Each reaction consists of a reaction rate (the expression on the left hand side of  ,), a set of substrates (the expression in-between , and -->), and a set of products (the expression on the right hand side of -->). The substrates and the products may contain one or more reactants, separated by +.  The naming convention for these are the same as for normal variables in Julia.The chemical reaction model is generated by the @reaction_network macro and stored in the rn variable (a normal variable, do not need to be called rn). The macro generates a differential equation model according to the law of mass action, in the above example the ODEs become:fracdXdt = -2 X Y\nfracdYdt = -2 X Y\nfracdXYdt = 2 X Y - XY\nfracdZ1dt= XY\nfracdZ2dt = XY"
},

{
    "location": "models/biological.html#Arrow-variants-1",
    "page": "Chemical Reaction Models",
    "title": "Arrow variants",
    "category": "section",
    "text": "Several types of arrows are accepted by the DSL and works instead of -->. All of these works:  >, → ↣, ↦, ⇾, ⟶, ⟼, ⥟, ⥟, ⇀, ⇁. Backwards arrows can also be used to write the reaction in the opposite direction. Hence all of these three reactions are equivalent:rn = @reaction_network begin\n  1.0, X + Y --> XY               \n  1.0, X + Y → XY      \n  1.0, XY ← X + Y      \nend(note that due to technical reasons <-- cannot be used)"
},

{
    "location": "models/biological.html#Using-bi-directional-arrows-1",
    "page": "Chemical Reaction Models",
    "title": "Using bi-directional arrows",
    "category": "section",
    "text": "Bi-directional arrows can be used to designate a reaction that goes two ways. These two models are equivalent:rn = @reaction_network begin\n  2.0, X + Y → XY             \n  2.0, X + Y ← XY          \nend\nrn = @reaction_network begin\n  2.0, X + Y ↔ XY               \nendIf the reaction rate in the backwards and forwards directions are different they can be designated in the following way:rn = @reaction_network begin\n  (2.0,1.0) X + Y ↔ XY               \nendwhich is identical torn = @reaction_network begin\n  2.0, X + Y → XY             \n  1.0, X + Y ← XY          \nend"
},

{
    "location": "models/biological.html#Combining-several-reactions-in-one-line-1",
    "page": "Chemical Reaction Models",
    "title": "Combining several reactions in one line",
    "category": "section",
    "text": "Several similar reactions can be combined in one line by providing a tuple of reaction rates and/or substrates and/or products. If several tuples are provided they much all be of identical length. These pairs of reaction networks are all identical:rn1 = @reaction_network begin\n  1.0, S → (P1,P2)               \nend\nrn2 = @reaction_network begin\n  1.0, S → P1     \n  1.0, S → P2\nendrn1 = @reaction_network begin\n  (1.0,2.0), (S1,S2) → P             \nend\nrn2 = @reaction_network begin\n  1.0, S1 → P     \n  2.0, S2 → P\nendrn1 = @reaction_network begin\n  (1.0,2.0,3.0), (S1,S2,S3) → (P1,P2,P3)        \nend\nrn2 = @reaction_network begin\n  1.0, S1 → P1\n  2.0, S2 → P2   \n  3.0, S3 → P3  \nendThis can also be combined with bi-directional arrows in which case separate tuples can be provided for the backward and forward reaction rates separately. These reaction networks are identicalrn1 = @reaction_network begin\n (1.0,(1.0,2.0)), S ↔ (P1,P2)  \nend\nrn2 = @reaction_network begin\n  1.0, S → P1\n  1.0, S → P2\n  1.0, P1 → S   \n  2.0, P2 → S\nend"
},

{
    "location": "models/biological.html#Production-and-Destruction-and-Stoichiometry-1",
    "page": "Chemical Reaction Models",
    "title": "Production and Destruction and Stoichiometry",
    "category": "section",
    "text": "Sometimes reactants are produced/destroyed from/to nothing. This can be designated using either 0 or ∅:rn = @reaction_network begin\n  2.0, 0 → X\n  1.0, X → ∅\nendSometimes several molecules of the same reactant is involved in a reaction, the stoichiometry of a reactant in a reaction can be set using a number. Here two species of X forms the dimer X2:rn = @reaction_network begin\n  1.0, 2X → X2\nendthis corresponds to the differential equation:fracdXdt = -X^2\nfracdX2dt = frac12 X^2Other numbers than 2 can be used and parenthesises can be used to use the same stoichiometry for several reactants:rn = @reaction_network begin\n  1.0, X + 2(Y + Z) → XY2Z2\nend"
},

{
    "location": "models/biological.html#Variable-reaction-rates-1",
    "page": "Chemical Reaction Models",
    "title": "Variable reaction rates",
    "category": "section",
    "text": "Reaction rates do not need to be constant, but can also depend on the current concentration of the various reactants (when e.g. one reactant activate the production of another one). E.g. this is a valid notation:rn = @reaction_network begin\n  X, Y → ∅\nendand will have Y degraded at ratefracdYdt = -XYNote that this is actually equivalent to the reactionrn = @reaction_network begin\n  1.0, X + Y → X\nendMost expressions and functions are valid reaction rates, e.g:rn = @reaction_network begin\n  2.0*X^2, 0 → X + Y\n  gamma(Y)/5, X → ∅\n  pi*X/Y, Y → ∅\nendplease note that user defined functions cannot be used directly (see later section \"User defined functions in reaction rates\")."
},

{
    "location": "models/biological.html#Defining-parameters-1",
    "page": "Chemical Reaction Models",
    "title": "Defining parameters",
    "category": "section",
    "text": "Just as when defining normal differential equations using DifferentialEquations parameter values does not need to be set when the model is created. Components can be designated as parameters by declaring them at the end:rn = @reaction_network begin\n  p, ∅ → X\n  d, X → ∅\nend p dParameters can only exist in the reaction rates (where they can be mixed with reactants). All variables not declared at the end will be considered a reactant."
},

{
    "location": "models/biological.html#Pre-defined-functions-1",
    "page": "Chemical Reaction Models",
    "title": "Pre-defined functions",
    "category": "section",
    "text": "Hill functions and a Michaelis-Menten function are pre-defined and can be used as rate laws. Below, the pair of reactions within rn1 are equivalent, as are the pair of reactions within rn2:rn1 = @reaction_network begin\n  hill(X,v,K,n), ∅ → X\n  v*X^n/(X^n+K^n), ∅ → X\nend v K n\nrn2 = @reaction_network begin\n  mm(X,v,K), ∅ → X\n  v*X/(X+K), ∅ → X\nend v KRepressor Hill (hillr) and Michaelis-Menten (mmr) functions are also provided:rn1 = @reaction_network begin\n  hillr(X,v,K,n), ∅ → X\n  v*K^n/(X^n+K^n), ∅ → X\nend v K n\nrn2 = @reaction_network begin\n  mmr(X,v,K), ∅ → X\n  v*K/(X+K), ∅ → X\nend v K"
},

{
    "location": "models/biological.html#Model-Simulation-1",
    "page": "Chemical Reaction Models",
    "title": "Model Simulation",
    "category": "section",
    "text": "Once created, a reaction network can be used as input to various problem types which can be solved by DifferentialEquations.jl."
},

{
    "location": "models/biological.html#Deterministic-simulations-using-ODEs-1",
    "page": "Chemical Reaction Models",
    "title": "Deterministic simulations using ODEs",
    "category": "section",
    "text": "A reaction network can be used as input to an ODEProblem instead of a function, using probODE = ODEProblem(rn, args...; kwargs...) E.g. a model can be created and simulated using:rn = @reaction_network begin\n  p, ∅ → X\n  d, X → ∅\nend p d\np = [1.0,2.0]\nu0 = [0.1]\ntspan = (0.,1.)\nprob = ODEProblem(rn,u0,tspan,p)\nsol = solve(prob)(if no parameters are given p does not need to be provided)To solve for a steady-state starting from the guess u0, one can useprob = SteadyStateProblem(rn,u0,p)\nsol = solve(prob, SSRootfind())orprob = SteadyStateProblem(rn,u0,p)\nsol = solve(prob, DynamicSS(Tsit5()))"
},

{
    "location": "models/biological.html#Stochastic-simulations-using-SDEs-1",
    "page": "Chemical Reaction Models",
    "title": "Stochastic simulations using SDEs",
    "category": "section",
    "text": "In a similar way a SDE can be created using probSDE = SDEProblem(rn, args...; kwargs...). In this case the chemical Langevin equations (as derived in Gillespie 2000) will be used to generate stochastic differential equations."
},

{
    "location": "models/biological.html#Stochastic-simulations-using-discrete-stochastic-simulation-algorithms-1",
    "page": "Chemical Reaction Models",
    "title": "Stochastic simulations using discrete stochastic simulation algorithms",
    "category": "section",
    "text": "Instead of solving SDEs one can create a stochastic jump process model using integer copy numbers and a discrete stochastic simulation algorithm. This can be done using:rn = @reaction_network begin\n  p, ∅ → X\n  d, X → ∅\nend p d\np = [1.0,2.0]\nu0 = [10]\ntspan = (0.,1.)\ndiscrete_prob = DiscreteProblem(rn, u0,tspan,p)\njump_prob = JumpProblem(discrete_prob,Direct(),rn)\nsol = solve(jump_prob,SSAStepper())Here we used Gillespie\'s Direct method as the underlying stochastic simulation algorithm."
},

{
    "location": "models/biological.html#Reaction-rate-laws-used-in-simulations-1",
    "page": "Chemical Reaction Models",
    "title": "Reaction rate laws used in simulations",
    "category": "section",
    "text": "In generating mathematical models from a reaction_network, reaction rates are treated as microscopic rates. That is, for a general mass action reaction of the form n_1 S_1 + n_2 S_2 + dots n_M S_M to dots with stoichiometric substrate coefficients n_i_i=1^M and rate constant k, the corresponding ODE rate law is taken to bek prod_i=1^M frac(S_i)^n_in_iwhile the jump process transition rate (i.e. propensity function) isk prod_i=1^M fracS_i (S_i-1) dots (S_i-n_i+1)n_iFor example, the ODE model of the reaction 2X + 3Y to Z with rate constant k would befracdXdt =  -2 k fracX^22 fracY^33 = -k fracX^2 Y^33 \nfracdYdt =  -3 k fracX^22 fracY^33 = -k fracX^2 Y^34 \nfracdZdt = k fracX^22 fracY^33"
},

{
    "location": "models/biological.html#The-Reaction-DSL-Advanced-1",
    "page": "Chemical Reaction Models",
    "title": "The Reaction DSL - Advanced",
    "category": "section",
    "text": "This section covers some of the more advanced syntax for building chemical reaction network models (still not very complicated!)."
},

{
    "location": "models/biological.html#User-defined-functions-in-reaction-rates-1",
    "page": "Chemical Reaction Models",
    "title": "User defined functions in reaction rates",
    "category": "section",
    "text": "The reaction network DSL cannot \"see\" user defined functions. E.g. this is not correct syntax:myHill(x) = 2.0*x^3/(x^3+1.5^3)\nrn = @reaction_network begin\n  myHill(X), ∅ → X\nendHowever, it is possible to define functions in such a way that the DSL can see them using the @reaction_func macro:@reaction_func myHill(x) = 2.0*x^3/(x^3+1.5^3)\nrn = @reaction_network begin\n  myHill(X), ∅ → X\nend"
},

{
    "location": "models/biological.html#Defining-a-custom-reaction-network-type-1",
    "page": "Chemical Reaction Models",
    "title": "Defining a custom reaction network type",
    "category": "section",
    "text": "While the default type of a reaction network is reaction_network (which inherits from AbstractReactionNetwork) it is possible to define a custom type (which also will inherit from AbstractReactionNetwork) by adding the type name as a first argument to the @reaction_network macro:rn = @reaction_network my_custom_type begin\n  1.0, ∅ → X\nend"
},

{
    "location": "models/biological.html#Scaling-noise-in-the-chemical-Langevin-equations-1",
    "page": "Chemical Reaction Models",
    "title": "Scaling noise in the chemical Langevin equations",
    "category": "section",
    "text": "When making stochastic simulations using SDEs it is possible to scale the amount of noise in the simulations by declaring a noise scaling parameter. This parameter is declared as a second argument to the @reaction_network macro (when scaling the noise one have to declare a custom type).rn = @reaction_network my_custom_type ns begin\n  1.0, ∅ → X\nendThe noise scaling parameter is automatically added as a last argument to the parameter array (even if not declared at the end). E.g. this is correct syntax:rn = @reaction_network my_custom_type ns begin\n  1.0, ∅ → X\nend\np = [0.1,]\nu0 = [0.1]\ntspan = (0.,1.)\nprob = SDEProblem(rn,u0,tspan,p)\nsol = solve(prob)Here the amount of noise in the stochastic simulation will be reduced by a factor 10."
},

{
    "location": "models/biological.html#Ignoring-mass-kinetics-1",
    "page": "Chemical Reaction Models",
    "title": "Ignoring mass kinetics",
    "category": "section",
    "text": "While one in almost all cases want the reaction rate to take the law of mass action into account, so the reactionrn = @reaction_network my_custom_type ns begin\n  k, X → ∅\nend koccur at the rate dXdt = -kX, it is possible to ignore this by using any of the following non-filled arrows when declaring the reaction: ⇐, ⟽, ⇒, ⟾, ⇔, ⟺. This means that the reactionrn = @reaction_network my_custom_type ns begin\n  k, X ⇒ ∅\nend kwill occur at rate dXdt = -k (which might become a problem since X will be degraded at a constant rate even when very small or equal to 0."
},

{
    "location": "models/biological.html#The-Reaction-Network-Object-1",
    "page": "Chemical Reaction Models",
    "title": "The Reaction Network Object",
    "category": "section",
    "text": "The @reaction_network macro generates a reaction_network object, which has a number of fields  which can be accessed.rn.f is a function encoding the right hand side of the ODEs (i.e. the time derivatives of the chemical species).\nrn.f_func is a vector of expressions corresponding to the time derivatives of the chemical species.\nrn.f_symfuncs is a vector of SymEngine expressions corresponding to the time derivatives of the chemical species.\nrn.g is a function encoding the noise terms for the SDEs (see rn.g_func for details).\nrn.g_func is a vector containing expressions corresponding to the noise terms used when creating the SDEs (n*m elements when there are n reactants and m reactions. The first m elements correspond to the noise terms for the first reactant and each reaction, the next m elements for the second reactant and all reactions, and so on).\nrn.jac is a function that evaluates the Jacobian of rn.f in place. i.e. has the form rn.jac(dJ,u,p,t), for pre-allocated Jacobian matrix dJ.\nrn.jump_affect_expr is a vector of expressions for how each reaction causes the species populations to change.\nrn.jump_rate_expr is a vector of expressions for how the transition rate (i.e. propensity) of each reaction is calculated from the species populations.\nrn.jumps is a vector storing a jump corresponding to each reaction (i.e. ConstantRateJump, VariableRateJump, etc...)\nrn.odefun stores an ODEFunction that can be used to create an ODEProblem corresponding to the reaction network.\nrn.p_matrix is a prototype matrix with the same size as the noise term.\nrn.paramjac is a function that evaluates the Jacobian of rn.f with respect to the parameters, p, in-place. It has the form rn.paramjac(dpJ,u,p,t) for pre-allocated parameter Jacobian matrix dpJ.\nrn.params is a vector containing symbols corresponding to all the parameters of the network.\nrn.params_to_ints provides a mapping from parameter symbol to the integer id of the parameter (i.e. where it is stored in the parameter vector passed to ODEProblem, SDEProblem, etc...)\nrn.reactions stores a vector of DiffEqBiological.ReactionStructs, which collect info for their corresponding reaction (such as stoichiometric coefficients).\nrn.regular_jumps stores a RegularJump representation of the network, for use in tau-leaping methods.\nrn.scale_noise is the noise scaling parameter symbol (if provided).\nrn.sdefun is a SDEFunction that can be used to create an SDEProblem corresponding to the reaction network.\nrn.symjac is the symbolically calculated Jacobian of the ODEs corresponding to the model.\nrn.syms is a vector containing symbols for all species of the network.\nrn.syms_to_ints is a map from the symbol of a species to its integer index within the solution vector."
},

{
    "location": "models/biological.html#The-Min-Reaction-Network-Object-1",
    "page": "Chemical Reaction Models",
    "title": "The Min Reaction Network Object",
    "category": "section",
    "text": "The @min_reaction_network macro works similarly to the @reaction_network macro, but initially only fills in fields corresponding to basic reaction network properties (i.e. rn.params, rn.params_to_ints, rn.scale_noise, rn.reactions, rn.syms, and rn.syms_to_ints). To fill in the remaining fields call (in the following [val] denotes the default value of a keyword argument):addodes!(rn) to complete ODE-related fields, optional keyword arguments include:\nbuild_jac=[true], is true if rn.jac and rn.symjac should be constructed. (Currently these build a dense Jacobian, so should be set to false for sufficiently large systems.)\nbuild_symfuncs=[true], is true if symbolic functions should be constructed for each ODE rhs. It is recommended to disable this for larger systems to reduce memory usage and speedup network construction.\naddsdes!(rn) to complete SDE-related fields. \naddjumps!(rn) to complete jump-related fields. addjumps! accepts several keyword arguments to control which jumps get created: \nbuild_jumps=[true] is true if rn.jumps should be constructed. This can be set to false for regular jump problems, where only rn.regular_jumps is needed.\nbuild_regular_jumps=[true] is true if rn.regular_jumps should be constructed. This can be set to false for Gillespie-type jump problems, where regular_jumps are not used.\nminimal_jumps=[false] is false if rn.jumps should contain a jump for each possible reaction. If set to true jumps are only added to rn.jumps for non-mass action jumps. (Note, mass action jumps are still resolved within any jump simulation. This option simply speeds up the construction of the jump problem since entries in rn.jumps that correspond to mass action jumps are never directly called within jump simulations.)For example, to simulate a jump process (i.e. Gillespie) simulation without constructing any RegularJumps, and only constructing a minimal set of jumps:rs = @min_reaction_network begin\n  c1, X --> 2X\n  c2, X --> 0\n  c3, 0 --> X\nend c1 c2 c3\np = (2.0,1.0,0.5)\naddjumps!(rs; build_regular_jumps=false, minimal_jumps=true)\nprob = DiscreteProblem(rs, [5], (0.0, 4.0), p)\njump_prob = JumpProblem(prob, Direct(), rs)\nsol = solve(jump_prob, SSAStepper())"
},

{
    "location": "models/biological.html#Chemical-Reaction-Network-Examples-1",
    "page": "Chemical Reaction Models",
    "title": "Chemical Reaction Network Examples",
    "category": "section",
    "text": ""
},

{
    "location": "models/biological.html#Example:-Birth-Death-Process-1",
    "page": "Chemical Reaction Models",
    "title": "Example: Birth-Death Process",
    "category": "section",
    "text": "rs = @reaction_network begin\n  c1, X --> 2X\n  c2, X --> 0\n  c3, 0 --> X\nend c1 c2 c3\np = (1.0,2.0,50.)\ntspan = (0.,4.)\nu0 = [5.]\n\n# solve ODEs\noprob = ODEProblem(rs, u0, tspan, p)\nosol  = solve(oprob, Tsit5())\n\n# solve for Steady-States\nssprob = SteadyStateProblem(rs, u0, p)\nsssol  = solve(ssprob, SSRootfind())\n\n# solve SDEs\nsprob = SDEProblem(rs, u0, tspan, p)\nssol  = solve(sprob, EM(), dt=.01)\n\n# solve JumpProblem\nu0 = [5]\ndprob = DiscreteProblem(rs, u0, tspan, p)\njprob = JumpProblem(dprob, Direct(), rs)\njsol = solve(jprob, SSAStepper())"
},

{
    "location": "models/biological.html#Example:-Michaelis-Menten-Enzyme-Kinetics-1",
    "page": "Chemical Reaction Models",
    "title": "Example: Michaelis-Menten Enzyme Kinetics",
    "category": "section",
    "text": "rs = @reaction_network begin\n  c1, S + E --> SE\n  c2, SE --> S + E\n  c3, SE --> P + E\nend c1 c2 c3\np = (0.00166,0.0001,0.1)\ntspan = (0., 100.)\nu0 = [301., 100., 0., 0.]  # S = 301, E = 100, SE = 0, P = 0\n\n# solve ODEs\noprob = ODEProblem(rs, u0, tspan, p)\nosol  = solve(oprob, Tsit5())\n\n# solve JumpProblem\nu0 = [301, 100, 0, 0] \ndprob = DiscreteProblem(rs, u0, tspan, p)\njprob = JumpProblem(dprob, Direct(), rs)\njsol = solve(jprob, SSAStepper())"
},

{
    "location": "models/external_modeling.html#",
    "page": "External Modeling Packages",
    "title": "External Modeling Packages",
    "category": "page",
    "text": ""
},

{
    "location": "models/external_modeling.html#External-Modeling-Packages-1",
    "page": "External Modeling Packages",
    "title": "External Modeling Packages",
    "category": "section",
    "text": "This is a list of modeling packages built upon the JuliaDiffEq ecosystem."
},

{
    "location": "models/external_modeling.html#DynamicalSystems.jl-1",
    "page": "External Modeling Packages",
    "title": "DynamicalSystems.jl",
    "category": "section",
    "text": "DynamicalSystems.jl is a package for the exploration of continuous and discrete dynamical systems, with focus on nonlinear dynamics and chaos.It uses DifferentialEquations.jl for all evolution regarding continuous systems while still retaining unified interface for discrete systems.A quick summary of features: Lyapunov exponents, generalized entropies (Renyi entropy), generalized & fractal dimensions, delay coordinates embedding (reconstruction), chaos detection, Lyapunov exponents of a numerical timeseries, finding periodic orbits of any order for maps."
},

{
    "location": "models/external_modeling.html#BioEnergeticFoodWebs.jl-1",
    "page": "External Modeling Packages",
    "title": "BioEnergeticFoodWebs.jl",
    "category": "section",
    "text": "BioEnergeticFoodWebs.jl is a package for simulations of biomass flows in food webs."
},

{
    "location": "models/external_modeling.html#SwitchTimeOpt.jl-1",
    "page": "External Modeling Packages",
    "title": "SwitchTimeOpt.jl",
    "category": "section",
    "text": "SwitchTimeOpt.jl is a Julia package to easily define and efficiently solve switching time optimization (STO) problems for linear and nonlinear systems."
},

{
    "location": "models/external_modeling.html#VehicleModels.jl-1",
    "page": "External Modeling Packages",
    "title": "VehicleModels.jl",
    "category": "section",
    "text": "VehicleModels.jl is a package for simulating vehicle models."
},

{
    "location": "models/external_modeling.html#MADS.jl-1",
    "page": "External Modeling Packages",
    "title": "MADS.jl",
    "category": "section",
    "text": "MADS.jl is a package data and model analysis. It adds many sensitivity analysis, uncertainty quantification, and model selection routines."
},

{
    "location": "models/external_modeling.html#QuantumOptics.jl-1",
    "page": "External Modeling Packages",
    "title": "QuantumOptics.jl",
    "category": "section",
    "text": "QunatumOptics.jl is a package for simulation of quantum systems."
},

{
    "location": "apis/diffeqbio.html#",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.jl API",
    "category": "page",
    "text": ""
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.jl-API-1",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.jl API",
    "category": "section",
    "text": "CurrentModule = DiffEqBiological"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.@reaction_network",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.@reaction_network",
    "category": "macro",
    "text": "@reaction_network\n\nGenerates a subtype of an AbstractReactionNetwork that encodes a chemical reaction network, and complete ODE, SDE and jump representations of the system. See the Chemical Reaction Model docs for details on parameters to the macro.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.@min_reaction_network",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.@min_reaction_network",
    "category": "macro",
    "text": "@min_reaction_network\n\nGenerates a subtype of an AbstractReactionNetwork that only encodes a chemical reaction network. Use addodes!, addsdes! or addjumps! to complete the network for specific problem types. It accepts the same arguments as @reaction_network.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.@empty_reaction_network",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.@empty_reaction_network",
    "category": "macro",
    "text": "@empty_reaction_network networktype\n\nGenerates a subtype of an AbstractReactionNetwork that encodes an empty chemical reaction network. networktype is an optional parameter that specifies the type of the generated network. Use addspecies!, addparam! and addreaction! to extend the network.  Use addodes!, addsdes! or addjumps! to complete the network for specific problem types.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#Reaction-Network-Generation-Macros-1",
    "page": "DiffEqBiological.jl API",
    "title": "Reaction Network Generation Macros",
    "category": "section",
    "text": "DiffEqBiological has three macros for generating reaction networks. @reaction_network generates a complete network, including everything needed to construct ODE, SDE and Jump problems directly from the network. For small systems it is the recommended form to use.@min_reaction_network constructs a network that just stores the basic information needed to represent the species, parameters and chemical reactions. This is sufficient for network analysis, such as calculating dependency graphs, but means the network must be extended to build mathematical models, see addodes!, addsdes!, and addjumps!. @empty_reaction_network constructs an empty network. Both min_reaction_networks and empty_reaction_networks can be enlarged using addspecies!, addparam!, and addreaction!. Once the final chemistry for the network is set, addodes!, addsdes!, and addjumps! can be called to build corresponding mathematical models.It is important to note for @reaction_network and @min_reaction_network that species which are used within the macro as part of a rate expression, but not as a substrate or product of some reaction, are not recognized as either a species or parameter. i.e. avoidrn = @reaction_network begin\n    k*X, Y --> W\nend kas here X is never defined as either a species or parameter. This leads to internal problems in the representation of reactions that can not be corrected by subsequently calling addspecies!. @reaction_network\n@min_reaction_network\n@empty_reaction_network"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.species",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.species",
    "category": "function",
    "text": "species(network)\n\nGiven an AbstractReactionNetwork, return a vector of species symbols.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.speciesmap",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.speciesmap",
    "category": "function",
    "text": "speciesmap(network)\n\nGiven an AbstractReactionNetwork, return a Dictionary mapping from species symbol to species index. \n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.params",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.params",
    "category": "function",
    "text": "params(network)\n\nGiven an AbstractReactionNetwork, return a vector of parameter symbols.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.paramsmap",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.paramsmap",
    "category": "function",
    "text": "paramsmap(network)\n\nGiven an AbstractReactionNetwork, return a Dictionary mapping from parameter symbol to parameter index.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.numspecies",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.numspecies",
    "category": "function",
    "text": "numspecies(network)\n\nReturn the number of species within the given AbstractReactionNetwork.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.numparams",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.numparams",
    "category": "function",
    "text": "numparams(network)\n\nReturn the number of parameters within the given AbstractReactionNetwork.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.numreactions",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.numreactions",
    "category": "function",
    "text": "numreactions(network)\n\nReturn the number of reactions within the given AbstractReactionNetwork.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#Basic-properties-1",
    "page": "DiffEqBiological.jl API",
    "title": "Basic properties",
    "category": "section",
    "text": "species\nspeciesmap\nparams\nparamsmap \nnumspecies \nnumparams\nnumreactions"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.substrates",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.substrates",
    "category": "function",
    "text": "substrates(network, rxidx)\n\nGiven an AbstractReactionNetwork and a reaction index, rxidx, return a vector of symbols of species that correspond to substrates in the reaction.  i.e. for\n\nk*W, X + 3Y --> X + W\n\nthe returned vector would be [:X,:Y].\n\nAllocates a new vector to store the symbols.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.products",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.products",
    "category": "function",
    "text": "products(network, rxidx)\n\nGiven an AbstractReactionNetwork and a reaction index, rxidx, return a vector of symbols of species that correspond to products in the reaction.  i.e. for\n\nk*W, X + 3Y --> X + W\n\nthe returned vector would be [:X,:W].\n\nAllocates a new vector to store the symbols.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.dependents",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.dependents",
    "category": "function",
    "text": "dependents(network, rxidx)\n\nGiven an AbstractReactionNetwork and a reaction index, rxidx, return a vector of symbols of species the reaction rate law depends on. i.e. for\n\nk*W, 2X + 3Y --> 5Z + W\n\nthe returned vector would be [:W,:X,:Y].\n\nNon-allocating, returns underlying field within the reaction_network.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.dependants",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.dependants",
    "category": "function",
    "text": "dependants(network, rxidx)\n\nSee documentation for dependents(network, rxidx).\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.ismassaction",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.ismassaction",
    "category": "function",
    "text": "ismassaction(network, rxidx)\n\nGiven an AbstractReactionNetwork and a reaction index, rxidx, return a boolean indicating whether the given reaction is of mass action form. For example, the reaction\n\n2*k, 2X + 3Y --> 5Z + W\n\nwould return true, while reactions with state-dependent rates like\n\nk*X, X + Y --> Z\n\nwould return false.\n\nNon-allocating, returns underlying field within the reaction_network.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.substratestoich",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.substratestoich",
    "category": "function",
    "text": "substratestoich(network, rxidx)\n\nGiven an AbstractReactionNetwork and a reaction index, rxidx, return a vector of pairs, mapping ids of species that serve as substrates in the reaction to the corresponding stoichiometric coefficient as a substrate. \n\nAllocates a new vector to store the pairs.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.substratesymstoich",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.substratesymstoich",
    "category": "function",
    "text": "substratesymstoich(network, rxidx)\n\nGiven an AbstractReactionNetwork and a reaction index, rxidx, return a Vector of ReactantStructs, mapping the symbols of species that serve as substrates in the reaction to the corresponding stoichiometric coefficient as a substrate. \n\nNon-allocating, returns underlying field within the reaction_network.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.productstoich",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.productstoich",
    "category": "function",
    "text": "productstoich(network, rxidx)\n\nGiven an AbstractReactionNetwork and a reaction index, rxidx, return a vector of pairs, mapping ids of species that are products in the reaction to the corresponding stoichiometric coefficient as a product.\n\nAllocates a new vector to store the pairs.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.productsymstoich",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.productsymstoich",
    "category": "function",
    "text": "productsymstoich(network, rxidx)\n\nGiven an AbstractReactionNetwork and a reaction index, rxidx, return a Vector of ReactantStructs, mapping the symbols of species that are products in the reaction to the corresponding stoichiometric coefficient as a product. \n\nNon-allocating, returns underlying field within the reaction_network.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.netstoich",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.netstoich",
    "category": "function",
    "text": "netstoich(network, rxidx)\n\nGiven an AbstractReactionNetwork and a reaction index, rxidx, return a vector of pairs, mapping ids of species that change numbers due to the reaction to the net stoichiometric coefficient of the species (i.e. net change in the species due to the reaction).\n\nAllocates a new vector to store the pairs.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#Reaction-Properties-1",
    "page": "DiffEqBiological.jl API",
    "title": "Reaction Properties",
    "category": "section",
    "text": "substrates\nproducts\ndependents\ndependants\nismassaction\nsubstratestoich\nsubstratesymstoich\nproductstoich\nproductsymstoich\nnetstoich"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.addspecies!",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.addspecies!",
    "category": "function",
    "text": "addspecies!(network, speciessym::Symbol)\n\nGiven an AbstractReaction network, add the species corresponding to the passed in symbol to the network (if it is not already defined).\n\n\n\n\n\naddspecies!(network, speciesname::String)\n\nGiven an AbstractReaction network, add the species with name given by the passed in string to the network (if it is not already defined.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.addparam!",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.addparam!",
    "category": "function",
    "text": "addparam!(network, param::Symbol)\n\nGiven an AbstractReaction network, add the parameter corresponding to the passed in symbol to the network (if it is not already defined).\n\n\n\n\n\naddparam!(network, paramname::String)\n\nGiven an AbstractReaction network, add the parameter with name given by the passed in string to the network (if it is not already defined).\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.addreaction!",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.addreaction!",
    "category": "function",
    "text": "addreaction!(network, rateex::Union{Expr,Symbol,Int,Float64}, rxexpr::Expr)\n\nGiven an AbstractReaction network, add a reaction with the passed in rate and reaction expressions. i.e. a reaction of the form\n\nk*X, 2X + Y --> 2W\n\nwould have rateex=:(k*X) and rxexpr=:(2X + Y --> W),\n\n10.5, 0 --> X\n\nwould have rateex=10.5 and rxexpr=:(0 --> X), and\n\nk, X+X --> Z\n\nwould have rateex=:k and rxexpr=:(X+X --> Z). All normal DSL reaction definition notation should be supported.\n\n\n\n\n\naddreaction!(network, rateex::Union{Expr,Symbol,Int,Float64}, substrates, products)\n\nGiven an AbstractReaction network, add a reaction with the passed in rate, rateex, substrate stoichiometry, and product stoichiometry. Stoichiometries are represented as tuples of Pair{Symbol,Int}. i.e. a reaction of the form\n\nk*X, 2X + Y --> 2W\n\nwould have rateex=:(k*X), substrates=(:X=>2, :Y=>2)andproducts=(W=>2,)`,\n\n10.5, 0 --> X\n\nwould have rateex=10.5, substrates=() and products=(:X=>1,), and\n\nk, X+X --> Z\n\nwould have rateex=:k, substrates=(:X=>2,) and products=(:Z=>2,). All normal DSL reaction definition notation should be supported for the rateex.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.add_scale_noise_param!",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.add_scale_noise_param!",
    "category": "function",
    "text": "add_scale_noise_param!(network, scale_noise::Symbol)\n\nGiven an AbstractReaction network, add the parameter corresponding to the passed in symbol to the network (if it is not already defined), and register it as the noise scaling coefficient.\n\n\n\n\n\nadd_scale_noise_param!(network, scale_noise_name::String)\n\nGiven an AbstractReaction network, add the parameter with the passed in string as its name to the network (if it is not already defined), and register it as the noise scaling coefficient.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#Functions-to-Add-Species,-Parameters-and-Reactions-to-a-Network-1",
    "page": "DiffEqBiological.jl API",
    "title": "Functions to Add Species, Parameters and Reactions to a Network",
    "category": "section",
    "text": "Both @min_reaction_network and @empty_reaction_network can be extended with additional species, parameters, and reactions. Note, always add all species and parameter definitions before adding any reaction definitions. Other orderings may result in incorrect information stored within the generated network.addspecies!\naddparam!\naddreaction!\nadd_scale_noise_param!"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.addodes!",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.addodes!",
    "category": "function",
    "text": "addodes!(network; build_jac=true, build_symfuncs=true)\n\nExtend an AbstractReactionNetwork generated with the @min_reaction_network or @empty_reaction_network macros with everything needed to use ODE solvers.\n\nOptional kwargs can be used to disable the construction of additional ODE solver components.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.addsdes!",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.addsdes!",
    "category": "function",
    "text": "addsdes!(network; build_jac=true, build_symfuncs=true)\n\nExtend an AbstractReactionNetwork generated with the @min_reaction_network or @empty_reaction_network macros with everything needed to use SDE solvers.\n\nOptional kwargs can be used to disable the construction of additional SDE solver components.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.addjumps!",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.addjumps!",
    "category": "function",
    "text": "addjumps!(network; build_jumps=true, build_regular_jumps=true, minimal_jumps=false)\n\nExtend an AbstractReactionNetwork generated with the @min_reaction_network or @empty_reaction_network macros with everything needed to use jump SSA solvers.\n\nOptional kwargs can be used to disable the construction of additional jump solver components.\n\nKeyword arguments:\n\nbuild_jumps: if true jump rates and affects will be calculated for use in DiffEqJump SSAs.\nbuild_regular_jumps: if true a RegularJump representation of the stochastic chemical kinetics model will be calculated for use in τ-leaping methods.\nminimal_jumps: if true ConstantRate jumps are only constructed for non-mass action jumps. (Note, mass action jumps are still resolved within any jump simulation. This option simply speeds up the construction of the jump problem since it avoids building redundant ConstantRate jumps that encode MassActionJumps, which are subsequently ignored within jump simulations.)\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#Functions-to-Add-ODEs,-SDEs-or-Jumps-to-a-Network-1",
    "page": "DiffEqBiological.jl API",
    "title": "Functions to Add ODEs, SDEs or Jumps to a Network",
    "category": "section",
    "text": "addodes!\naddsdes!\naddjumps!"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.oderhsfun",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.oderhsfun",
    "category": "function",
    "text": "oderhsfun(network)\n\nGiven an AbstractReactionNetwork, return a function, f!(du,u,p,t), that evaluates the current value of the ODE model derivative functions, dudt = f(ut),  within du.\n\nNote, for a network generated with the @min_reaction_network macro addodes! must be called first.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.jacfun",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.jacfun",
    "category": "function",
    "text": "jacfun(network)\n\nGiven an AbstractReactionNetwork, return a function, jac!(J,u,p,t), that evaluates the current Jacobian matrix, J, of the ODE model, dudt = f(ut).  The Jacobian matrix has entries \n\nJ_ij = partial f_i(ut)  partial u_j.\n\nNote, for a network generated with the @min_reaction_network macro addodes!  must be called first.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.paramjacfun",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.paramjacfun",
    "category": "function",
    "text": "paramjacfun(network)\n\nGiven an AbstractReactionNetwork, return a function, pjac(pJ,u,p,t), that evaluates the current parameter Jacobian matrix, pJ, of the ODE model, dudt = f(ut).  The parameter Jacobian matrix has entries \n\npJ_ij = partial f_i(ut)  partial p_j.\n\nNote, for a network generated with the @min_reaction_network macro addodes!  must be called first.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.odefun",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.odefun",
    "category": "function",
    "text": "odefun(network)\n\nGiven an AbstractReactionNetwork, return a DiffEqBase.ODEFunction encoding an ODE model for the reaction network. \n\nNote, for a network generated with the @min_reaction_network macro addodes!  must be called first.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.noisefun",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.noisefun",
    "category": "function",
    "text": "noisefun(network)\n\nGiven an AbstractReactionNetwork, return a function, g(η,u,p,t), that evaluates the current noise coefficients for each reaction in the Chemical Langevin Equation representation within η.\n\nNote, for a network generated with the @min_reaction_network macro addsdes!  must be called first.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.sdefun",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.sdefun",
    "category": "function",
    "text": "sdefun(network)\n\nGiven an AbstractReactionNetwork, return a DiffEqBase.SDEFunction encoding a Chemical Langevin Equation SDE model for the reaction network. \n\nNote, for a network generated with the @min_reaction_network macro addsdes!  must be called first.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.jumps",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.jumps",
    "category": "function",
    "text": "jumps(network)\n\nGiven an AbstractReactionNetwork, return a tuple of AbstractJumps encoding a stochastical chemical kinetics representation for the reaction network.\n\nNote, for a network generated with the @min_reaction_network macro addjumps!  must be called first.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.regularjumps",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.regularjumps",
    "category": "function",
    "text": "regularjumps(network)\n\nGiven an AbstractReactionNetwork, return a RegularJump encoding a stochastical chemical kinetics representation of the reaction network for use in tau-leaping approximations.\n\nNote, for a network generated with the @min_reaction_network macro addjumps!  must be called first.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#Generated-Functions-for-Models-1",
    "page": "DiffEqBiological.jl API",
    "title": "Generated Functions for Models",
    "category": "section",
    "text": "oderhsfun\njacfun\nparamjacfun\nodefun\nnoisefun\nsdefun\njumps\nregularjumps"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.odeexprs",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.odeexprs",
    "category": "function",
    "text": "odeexprs(network)\n\nGiven an AbstractReactionNetwork, return a vector of the ODE expressions.\n\nNote, for a network generated with the @min_reaction_network macro addodes!  must be called first.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.jacobianexprs",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.jacobianexprs",
    "category": "function",
    "text": "jacobianexprs(network)\n\nGiven an AbstractReactionNetwork, return a matrix with the ODE Jacobian expressions.\n\nNote, for a network generated with the @min_reaction_network macro addodes!  must be called first.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.noiseexprs",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.noiseexprs",
    "category": "function",
    "text": "noiseexprs(network)\n\nGiven an AbstractReactionNetwork, return a vector of the SDE noise expressions for each reaction.\n\nNote, for a network generated with the @min_reaction_network macro addsdes!  must be called first.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.jumpexprs",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.jumpexprs",
    "category": "function",
    "text": "jumpexprs(network)\n\nGiven an AbstractReactionNetwork, return a tuple of the jump rates and affects expressions.\n\nNote, for a network generated with the @min_reaction_network macro addjumps!  must be called first.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.rateexpr",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.rateexpr",
    "category": "function",
    "text": "rateexpr(network, rxidx)\n\nGiven an AbstractReactionNetwork, return the reaction rate expression for the reaction with index rxidx. Note, for a reaction defined by\n\nk*X*Y, X+Z --> 2X + Y\n\nthe expression that is returned will be :(k*X*Y), while the rate law used in ODEs and SDEs would be k*X^2*Y*Z.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.oderatelawexpr",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.oderatelawexpr",
    "category": "function",
    "text": "oderatelawexpr(network, rxidx)\n\nGiven an AbstractReactionNetwork, return the reaction rate law expression used in generated ODEs for the reaction with index rxidx. Note, for a reaction defined by\n\nk*X*Y, X+Z --> 2X + Y\n\nthe expression that is returned will be :(k*X^2*Y*Z). For a reaction of the form \n\nk, 2X+3Y --> Z\n\nthe expression that is returned will be :(k * (X^2/2) * (Y^3/6)).\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.ssaratelawexpr",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.ssaratelawexpr",
    "category": "function",
    "text": "ssaratelawexpr(network, rxidx)\n\nGiven an AbstractReactionNetwork, return the reaction rate law expression used in generated stochastic chemical kinetic model SSAs for the reaction with index rxidx. Note, for a reaction defined by\n\nk*X*Y, X+Z --> 2X + Y\n\nthe expression that is returned will be :(k*X^2*Y*Z). For a reaction of the form \n\nk, 2X+3Y --> Z\n\nthe expression that is returned will be :(k * binomial(X,2) * binomial(Y,3)).\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#Generated-Expressions-1",
    "page": "DiffEqBiological.jl API",
    "title": "Generated Expressions",
    "category": "section",
    "text": "odeexprs\njacobianexprs\nnoiseexprs\njumpexprs\nrateexpr\noderatelawexpr\nssaratelawexpr"
},

{
    "location": "apis/diffeqbio.html#Base.:==-Tuple{DiffEqBase.AbstractReactionNetwork,DiffEqBase.AbstractReactionNetwork}",
    "page": "DiffEqBiological.jl API",
    "title": "Base.:==",
    "category": "method",
    "text": "==(rn1::DiffEqBase.AbstractReactionNetwork, rn2::DiffEqBase.AbstractReactionNetwork)\n\nTests whether the underlying species symbols, parameter symbols and reactions are the same in the two networks. Ignores order network components were defined, so the integer id of any individual species/parameters/reactions may be different between the two networks. Does not currently account for different reaction definitions, so \"k, X+Y –> Y + Z\" will not be the same as \"k, Y+X –> Y + Z\"\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#Network-Comparison-Functions-1",
    "page": "DiffEqBiological.jl API",
    "title": "Network Comparison Functions",
    "category": "section",
    "text": "==(rn1::DiffEqBase.AbstractReactionNetwork, rn2::DiffEqBase.AbstractReactionNetwork)"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.rxtospecies_depgraph",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.rxtospecies_depgraph",
    "category": "function",
    "text": "rxtospecies_depgraph(network)\n\nGiven an AbstractReactionNetwork, returns a Vector{Vector{Int}} mapping each reaction index to the indices of species that depend on it.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.speciestorx_depgraph",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.speciestorx_depgraph",
    "category": "function",
    "text": "speciestorx_depgraph(network)\n\nGiven an AbstractReactionNetwork, returns a Vector{Vector{Int}} mapping each species index to the indices of reactions that depend on it.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#DiffEqBiological.rxtorx_depgraph",
    "page": "DiffEqBiological.jl API",
    "title": "DiffEqBiological.rxtorx_depgraph",
    "category": "function",
    "text": "rxtorx_depgraph(network)\n\nGiven an AbstractReactionNetwork, returns a Vector{Vector{Int}} mapping each reaction index to the indices of reactions that depend on it.\n\n\n\n\n\n"
},

{
    "location": "apis/diffeqbio.html#Dependency-Graphs-1",
    "page": "DiffEqBiological.jl API",
    "title": "Dependency Graphs",
    "category": "section",
    "text": "rxtospecies_depgraph\nspeciestorx_depgraph\nrxtorx_depgraph"
},

{
    "location": "extras/timestepping.html#",
    "page": "Timestepping Method Descriptions",
    "title": "Timestepping Method Descriptions",
    "category": "page",
    "text": ""
},

{
    "location": "extras/timestepping.html#Timestepping-Method-Descriptions-1",
    "page": "Timestepping Method Descriptions",
    "title": "Timestepping Method Descriptions",
    "category": "section",
    "text": ""
},

{
    "location": "extras/timestepping.html#Common-Setup-1",
    "page": "Timestepping Method Descriptions",
    "title": "Common Setup",
    "category": "section",
    "text": "All methods start by calculating a scaled error estimate on each scalar component of u:err^scaled_i = norm(err_i(abstol_i + max(uprev_iu_i)reltol_i))On this scaled error estimate, we calculate the norm. This norm is usually the Hairer semi-norm:norm(x) = sqrt(sum(x^2)length(x))This norm works well because it does not change if we add new pieces to the differential equation: it scales our error by the number of equations so that independent equations will not step differently than a single solve.In all cases, the step is rejected if err^scaled1 since that means the error is larger than the tolerances, and the step is accepted if err^scaled1."
},

{
    "location": "extras/timestepping.html#Proportional-Control-(Standard-Control)-1",
    "page": "Timestepping Method Descriptions",
    "title": "Proportional Control (Standard Control)",
    "category": "section",
    "text": "The proportional control algorithm is the \"standard algorithm\" for adaptive timestepping. Note that it is not the default in DifferentialEquations.jl because it is usually awful for performance, but it is explained first because it is the most widely taught algorithm and others build off of its techniques.The control simply changes dt proportional to the error. There is an exponentiation based on the order of the algorithm which goes back to a result by Cechino for the optimal stepsize to reduce the error. The algorithm is:qtmp = integrator.EEst^(1/(alg_adaptive_order(integrator.alg)+1))/integrator.opts.gamma\n@fastmath q = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),qtmp))\nintegrator.dtnew = integrator.dt/qThus q is the scaling factor for dt, and it must be between qmin and qmax. gamma is the safety factor, 0.9, for how much dt is decreased below the theoretical \"optimal\" value.Since proportional control is \"jagged\", i.e. can cause large changes between one step to the next, it can effect the stability of explicit methods. Thus it\'s only applied by default to low order implicit solvers."
},

{
    "location": "extras/timestepping.html#Proportional-Integral-Control-(PI-Control)-1",
    "page": "Timestepping Method Descriptions",
    "title": "Proportional-Integral Control (PI-Control)",
    "category": "section",
    "text": "The proportional-integral control algorithm is a standard control algorithm from control theory. It mixes proportional control with memory in order to make the timesteps more stable, which actually increases the adaptive stability region of the algorithm. This stability property means that it\'s well-suited for explicit solvers, and it\'s applied by default to the Rosenbrock methods as well. The form for the updates is:EEst,beta1,q11,qold,beta2 = integrator.EEst, integrator.opts.beta1, integrator.q11,integrator.qold,integrator.opts.beta2\n@fastmath q11 = EEst^beta1\n@fastmath q = q11/(qold^beta2)\nintegrator.q11 = q11\n@fastmath q = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),q/integrator.opts.gamma))\nif q <= integrator.opts.qsteady_max && q >= integrator.opts.qsteady_min\n  q = one(q)\nend\nqbeta1 is the gain on the proportional part, and beta2 is the gain for the history portion. qoldinit is the initialized value for the gain history."
},

{
    "location": "extras/timestepping.html#Gustafsson-Acceleration-1",
    "page": "Timestepping Method Descriptions",
    "title": "Gustafsson Acceleration",
    "category": "section",
    "text": "The Gustafsson acceleration algorithm accelerates changes so that way algorithms can more swiftly change to handle quick transients. This algorithm is thus well-suited for stiff solvers where this can be expected, and is the default for algorithms like the (E)SDIRK methods.gamma = integrator.opts.gamma\nniters = integrator.cache.newton_iters\nfac = min(gamma,(1+2*integrator.alg.max_newton_iter)*gamma/(niters+2*integrator.alg.max_newton_iter))\nexpo = 1/(alg_order(integrator.alg)+1)\nqtmp = (integrator.EEst^expo)/fac\n@fastmath q = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),qtmp))\nif q <= integrator.opts.qsteady_max && q >= integrator.opts.qsteady_min\n  q = one(q)\nend\nintegrator.qold = q\nqIn this case, niters is the number of Newton iterations which was required in the most recent step of the algorithm. Note that these values are used differently depending on acceptance and rejectance. When the step is accepted, the following logic is applied:if integrator.success_iter > 0\n  expo = 1/(alg_adaptive_order(integrator.alg)+1)\n  qgus=(integrator.dtacc/integrator.dt)*(((integrator.EEst^2)/integrator.erracc)^expo)\n  qgus = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),qgus/integrator.opts.gamma))\n  qacc=max(q,qgus)\nelse\n  qacc = q\nend\nintegrator.dtacc = integrator.dt\nintegrator.erracc = max(1e-2,integrator.EEst)\nintegrator.dt/qaccWhen it rejects, its the same as the proportional control:if integrator.success_iter == 0\n  integrator.dt *= 0.1\nelse\n  integrator.dt = integrator.dt/integrator.qold\nend"
},

]}
