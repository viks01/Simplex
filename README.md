# Simplex Algorithm 
- A programming assignment in optimization methods.
- Uses simplex algorithm to solve linear programming problems.
- Optimization problem of the form 
```math
min c^{T}x <br>
subject to Ax = b <br> 
x \geq 0
```
- Input from stdin.
- Output to stdout.
- Uses two phase method for initialization of simplex algorithm if initial basic feasible solution is not easy to find. 
- Uses Bland’s rule only if cycling is detected.
- Handles unbounded solution and no solution cases.
- Created in python.

### Input format
- The first line contains three integers: n, u, and v.
- The second line contains n numbers, denoting the vector c.
- The next u + v lines contain n numbers each.
- The first u lines denote the coefficients of the decision variables in inequalities of the form ≤.
- The next v lines denote the coefficients of the decision variables in inequalities of the form ≥.
- The final line contains u + v numbers, which denote the R.H.S values of the u + v constraints given, or in other terms, the vector b.

### Output format
- If the optimal solution exists and no cycling is detected: Print the optimal value on line 1. Print the vector of optimal values of the variables on line 2.
- If the optimal solution exists and cycling is detected: Print “Cycling detected” on line 1. Print the optimal value on line 2. Print the vector of optimal values of the variables on line 3.
- If the problem is unbounded: Print “Unbounded”.
- If the problem is infeasible: Print “Infeasible”.
