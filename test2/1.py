import numpy as np
from simplex import *

# Get the solution vector after simplex is run
def get_solution(table, basic, size):
    solution = []
    for i in range(size):
        if i not in basic:
            solution.append(0)
        else:
            idx = 0
            for j in range(len(basic)):
                if basic[j] == i:
                    idx = j
                    break
            solution.append(table[idx + 1][-1])
    return solution

# Format and display output based on status
def display_result(table, basic, status):
    if status == "Success" or status == "Cycling detected":
        solution = get_solution(table, basic, len(table[0]) - 1)

        if status == "Cycling detected":
            print(status)

        print(table[0][-1])
        print(*solution)
    else:
        print(status)

################################################################################################################################

# Input
N, u, v = list(map(int, input().strip().split()))                        # N - dimension of c
c = list(map(float, input().strip().split()))                            # objective coefficients
U = [[float(i) for i in input().strip().split()] for x in range(u)]      # <= inequalities
V = [[float(i) for i in input().strip().split()] for x in range(v)]      # >= inequalities
b = list(map(float, input().strip().split()))                            # RHS of inequalities

# A --> m x n
m = u + v      # total number of slack variables
n = N + m      # total number of variables (objective + slack)

# Augment objective vector by m values
for i in range(m):
    c.append(0)

# Augment each inequality by m values for slack variables, converting the inequalities to equalities
for i in range(u):
    for j in range(m):
        val = 1 if i == j else 0
        U[i].append(val)
for i in range(v):
    for j in range(m):
        val = -1 if i + u == j else 0
        V[i].append(val)

# Combine U and V matrices into one A matrix
A = []
for i in range(u):
    A.append(U[i])
for i in range(v):
    A.append(V[i])

################################################################################################################################

# Check for two phased method
twoPhased = False
# Case 1: <= inequalities
for i in range(u):
    if b[i] < 0:
        twoPhased = True
# Case 2: >= inequalities
for i in range(u, m):
    if b[i] > 0:
        twoPhased = True

# Make b >= 0 without changing the equalities in Ax = b
for i in range(m):
    if b[i] < 0:
        b[i] *= -1
        for j in range(n):
            A[i][j] *= -1

################################################################################################################################

x = []
basic = []

if twoPhased:
    print('Two-phased')
    print()

    # Phase 1
    
    # Check for identity columns (basic variables) already existing among the slack variables
    basic = []
    rowIdxs = []                                # Keep track of rows of A to skip adding ones

    # print(*A, sep='\n')
    # print()

    for i in range(m):
        j = i + N
        if A[i][j] == 1:
            basic.append(j)
            rowIdxs.append(i)

    # print(basic)
    # print(rowIdxs)
    # print()
    
    # Find number of artificial variables to include
    num_artificial = m - len(rowIdxs)

    # Make a copy of A
    A1 = [row[:] for row in A]
    
    # Append identity columns to A corresponding to artificial variables
    for j in range(num_artificial):
        flag = True
        for i in range(m):
            val = 0
            if flag and i not in rowIdxs:
                val = 1
                rowIdxs.append(i)
                flag = False
            A1[i].append(val)

    # print(*A1, sep='\n')
    # print()

    # Update the basic variables list to include artificial variables
    for i in range(n, n + num_artificial):
        basic.append(i)

    # Order the basic variables by position of 1 in identity column (based on row)
    ordered_basic = []
    for i in range(m):
        for j in basic:
            if A1[i][j] == 1:
                ordered_basic.append(j)

    # print(ordered_basic)
    # print(basic)
    # print()

    # Initialize the objective vector
    e = []
    for i in range(n):
        e.append(0) 
    for i in range(n, n + num_artificial):
        e.append(1)

    # Initialize vector of decision variables (non-basic variables are 0 and basic variables are elements of b)
    x = []
    i = 0
    for j in range(n + num_artificial):
        val = 0
        if j in ordered_basic:
            val = b[i]
            i += 1
        x.append(val)

    # Initialize tableau form
    # First row (augmented with objective value in last column)
    table = [[-i for i in e]]
    table[0].append(0)

    # Rows corresponding to slack variables
    for i in range(m):
        table.append(A1[i])

    # Last column or RHS
    for i in range(1, m+1):
        table[i].append(b[i-1])

    print(*table, sep="\n")
    print()

    # Make basic variables zero in 1st row
    for j in range(n, n + num_artificial):
        pivot = table[0][j]
        if pivot < 0:
            idx = 0
            for i in range(1, m+1):
                if table[i][j] == 1:
                    idx = i
                    break
            for i in range(n + num_artificial + 1):
                table[0][i] -= pivot * table[idx][i]        # Equivalent to R_0 --> R_0 + R_idx
    
    print(*table, sep="\n")
    print()

    # Run simplex algorithm
    table, basic, status = simplex(table, ordered_basic)
    solution = get_solution(table, basic, len(table[0]) - 1)

    print(*table, sep="\n")
    print()
    print(basic)
    print(solution)
    print()
    
    # Check for infeasibility
    infeasible = False
    for i in range(n, n + num_artificial):
        if solution[i] != 0:
            infeasible = True
            status = "Infeasible"
            break

    if table[0][-1] != 0:
        infeasible = True
        status = "Infeasible"

    if infeasible:
        display_result(table, basic, status)
    else:
        # Phase 2

        # Form a new table with lesser columns without the artificial variables
        new_table = [[-i for i in c]]
        new_table[0].append(0)
        for i in range(1, len(table)):
            row = []
            for j in range(n):
                row.append(table[i][j])
            row.append(table[i][-1])
            new_table.append(row)

        # Make basic variables zero in 1st row
        for j in basic:
            pivot = new_table[0][j]
            if pivot != 0:
                idx = 0
                for i in range(1, m+1):
                    if new_table[i][j] == 1:
                        idx = i
                        break
                for i in range(n+1):
                    new_table[0][i] -= pivot * new_table[idx][i]

        # Run simplex algorithm and show output
        new_table, basic, status = simplex(new_table, basic)
        print(*new_table, sep='\n')
        print(basic)
        display_result(new_table, basic, status)
                
else:
    # Initialize vector of decision variables (non-basic variables are 0 and basic variables are elements of b)
    x = []
    for i in range(N):
        x.append(0)
    for i in range(m):
        x.append(b[i])

    # Initialize basic variables index set
    basic = [i for i in range(N, n)]

    # Initialize tableau form
    # First row (augmented with objective value in last column)
    table = [[-i for i in c]]
    table[0].append(0)

    # Rows corresponding to slack variables
    for i in range(m):
        table.append(A[i])

    # Last column or RHS
    for i in range(1, m+1):
        table[i].append(b[i-1])

    # print(*table, sep="\n")
    # print()

    # Run simplex algorithm and show output
    table, basic, status = simplex(table, basic)
    display_result(table, basic, status)

