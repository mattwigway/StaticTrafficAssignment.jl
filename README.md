# StaticTrafficAssignment.jl

This is code from [my dissertation](https://hdl.handle.net/2286/R.2.N.161665) for doing Frank-Wolfe static traffic assignment. It's currently a WIP as it's been copied out of individual scripts into package that is not complete yet, so it does not work at this time, but the code that actually implements the algorithm is in `src/assignment.jl` It also needs additional testing to ensure correctness; in particular it uses the `Threads.threadid()` pattern that can cause race conditions in recent versions of Julia (though not with Julia 1.5, which was what was originally used for the dissertation).