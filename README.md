# Well-balanced scheme for Euler equations with gravity

A 2-d well-balanced finite volume code for Euler equations with gravity. The scheme is explained in the paper

Praveen C and C. Klingenberg, "A second order well-balanced finite volume scheme for Euler equations with gravity", SIAM J. Sci. Comp., vol. 37, Issue 3, pp. 382-402.

Initial conditions and problem definition are in separate files called ```init_cond_xyz.f95```; uncomment the required case in ```init_cond.f95```.

The boundaries can be treated as solid walls. If you do not want solid wall, set

```
wall = 0
```

in ```solveFVM.f95```, otherwise set it to 1.

To switch off gravity, return 1.0 in function ```potential```.
