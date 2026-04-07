Here are some notes:

1. Important for Newton-Krylov methods that they converged on the first sub-step; otherwise, the displacement will remain zero.
2. Newton-Krylov methods converge fast, but require lots of sub-steps. Additionally, they aren't exact, which makes them especially unhelpful for contact problems.
3. For all modified Newton-based methods, preferably apply a linear subloading scheme.
4.   
