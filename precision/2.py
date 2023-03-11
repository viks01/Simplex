import math

def simplex(table, basic):
    status = "Success"

    # Number of basic variables
    m = len(table) - 1

    # Initial check of first row: (c_B)(A_B^-1)(a_j) - (c_j) > 0
    maxVal = float("-inf")
    maxIdx = -1                                  # index of entering variable
    for i in range(len(table[0]) - 1):
        val = table[0][i]
        if val > maxVal:                         # entering variable has largest positive reduced cost
            maxVal = val
            maxIdx = i
    
    # If maxVal <= 0, then the while loop is skipped and the initial basis selected is optimal
    while maxVal > 0:

        # Check for leaving variable along the column of entering variable (maxIdx)
        minVal = float("inf")
        minIdx = -1                              # index of leaving variable
        for i in range(1, len(table)):
            y = table[i][maxIdx]
            if y > 0:
                val = table[i][-1] / y
                if val < minVal:                 # leaving variable is from lowest index row 
                    minVal = val
                    minIdx = i
                if val == 0:                     # Check for cycling
                    status = "Cycling detected"
                    break

        # Check for unboundedness
        if minIdx == -1:
            status = "Unbounded"
            break

        # Check for cycling
        if status == "Cycling detected":
            break
        
        # Update basic variables set. Here, table basic variables are 1-indexed (because we skip 1st row) and basic variables index set is 0-indexed
        basic[minIdx - 1] = maxIdx

        # Divide minIdx row by corresponding y. (minIdx = r, y = y_r)
        y = table[minIdx][maxIdx]
        for i in range(len(table[0])):
            table[minIdx][i] /= y

        # Update other rows
        for i in range(len(table)):
            if i != minIdx:
                pivot = table[i][maxIdx]
                for j in range(len(table[0])):
                    table[i][j] -= pivot * table[minIdx][j]
        
        # Initial check of first row: (c_B)(A_B^-1)(a_j) - (c_j) > 0
        maxVal = float("-inf")
        maxIdx = -1                                  # index of entering variable
        for i in range(len(table[0]) - 1):
            val = table[0][i]
            if val > maxVal:                         # entering variable has largest positive reduced cost
                maxVal = val
                maxIdx = i

    # Bland's rule
    if status == "Cycling detected":

        # Initial check of first row: (c_B)(A_B^-1)(a_j) - (c_j) > 0
        maxVal = float("inf")
        maxIdx = -1                                  # index of entering variable
        for i in range(len(table[0]) - 1):
            val = table[0][i]
            if val > 0 and val < maxVal:             # entering variable has smallest positive reduced cost
                maxVal = val
                maxIdx = i
        
        # If maxIdx < 0, then the while loop is skipped and the optimal solution is reached
        while maxIdx >= 0:

            # Check for leaving variable along the column of entering variable (maxIdx)
            minVal = float("inf")
            minIdx = -1                              # index of leaving variable
            for i in range(1, len(table)):
                y = table[i][maxIdx]
                if y > 0:
                    val = table[i][-1] / y
                    if val < minVal:                 # leaving variable is from lowest index row 
                        minVal = val
                        minIdx = i

            # Check for unboundedness
            if minIdx == -1:
                status = "Unbounded"
                break
            
            # Update basic variables set. Here, table basic variables are 1-indexed (because we skip 1st row) and basic variables index set is 0-indexed
            basic[minIdx - 1] = maxIdx

            # Divide minIdx row by corresponding y. (minIdx = r, y = y_r)
            y = table[minIdx][maxIdx]
            for i in range(len(table[0])):
                table[minIdx][i] /= y

            # Update other rows
            for i in range(len(table)):
                if i != minIdx:
                    pivot = table[i][maxIdx]
                    for j in range(len(table[0])):
                        table[i][j] -= pivot * table[minIdx][j]
            
            # Initial check of first row: (c_B)(A_B^-1)(a_j) - (c_j) > 0
            maxVal = float("inf")
            maxIdx = -1                                  # index of entering variable
            for i in range(len(table[0]) - 1):
                val = table[0][i]
                if val > 0 and val < maxVal:             # entering variable has smallest positive reduced cost
                    maxVal = val
                    maxIdx = i

    return table, basic, status

# Initialize vector of decision variables (non-basic variables are 0 and basic variables are elements of b)
def initial_x0(basic, size, b):
    x = []
    i = 0
    for j in range(size):
        val = 0
        if j in basic:
            val = b[i]
            i += 1
        x.append(val)
    return x

# Get the solution vector after simplex is run
def get_solution(table, basic, size):
    solution = []
    rhs = []
    for i in range(len(table)):
        rhs.append(table[i][-1])

    # size = total number of variables in the solution vector
    for i in range(size):
        idx = -1
        for j in range(len(basic)):
            if basic[j] == i:
                idx = j
                break
        if idx == -1:
            solution.append(0)
        else:
            solution.append(rhs[idx + 1])
    return solution

# Format and display output based on status
def display_result(table, basic, status, N=0):
    if N > 0 and (status == "Success" or status == "Cycling detected"):
        solution = get_solution(table, basic, len(table[0]) - 1)
        res = solution[:N]
        
        if status == "Cycling detected":
            print(status)

        print("%.7f" % table[0][-1])
        for val in res:
            print("%.7f" % val, end=" ")
        print()
    else:
        print(status)

# Does pivoting
def dual_simplex(table, basic):
    status = "Success"

    # Number of basic variables
    m = len(table) - 1

    # Leaving variable is the row index of most negative value in RHS (last) column
    minVal = float("inf")
    minIdx = -1                              # index of leaving variable
    for i in range(1, len(table)):
        val = table[i][-1]
        if val < minVal:
            minVal = val
            minIdx = i

    # Assert minIdx = len(table) - 1
    # if minIdx != len(table) - 1:
    #     print('minIdx error in dual simplex')
    #     status = "Failure"

    # If minVal >= 0, then the current solution is optimal and the while loop ends
    while minVal < 0:

        # Entering variable is from the minimum positive absolute value of first row divided by value in row of leaving variable (minIdx)
        minAbsVal = float("inf")
        minAbsValIdx = -1
        for j in range(len(table[0]) - 1):
            y = table[minIdx][j]
            if y < 0:
                val = abs(table[0][j] / y)
                if val < minAbsVal:
                    minAbsVal = val
                    minAbsValIdx = j

        # Assert minAbsValIdx != -1 (unbounded)
        # if minAbsValIdx == -1:
        #     print('minAbsVal error in dual simplex')
        #     status = "Failure"

        if minAbsValIdx == -1:
            status = "Unbounded"
            break
        
        # Update basic variables set. Here, table basic variables are 1-indexed (because we skip 1st row) and basic variables index set is 0-indexed
        basic[minIdx - 1] = minAbsValIdx

        # Divide minIdx row by corresponding pivot
        pivot = table[minIdx][minAbsValIdx]
        for j in range(len(table[0])):
            table[minIdx][j] /= pivot

        # Update other rows
        for i in range(len(table)):
            if i != minIdx:
                pivot = table[i][minAbsValIdx]
                for j in range(len(table[0])):
                    table[i][j] -= pivot * table[minIdx][j]

        # Leaving variable is the row index of most negative value in RHS (last) column
        minVal = float("inf")
        minIdx = -1                              # index of leaving variable
        for i in range(1, len(table)):
            val = table[i][-1]
            if val < minVal:
                minVal = val
                minIdx = i
    
    return table, basic, status

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

x = []
basic = []
solution = []
table = [[]]
status = None

if twoPhased:
    ########################################### Phase 1 ###########################################
    
    # Check for identity columns (basic variables) already existing among the slack variables
    basic = []
    rowIdxs = []                                # Keep track of rows of A to skip adding ones

    for i in range(m):
        j = i + N
        if A[i][j] == 1:
            basic.append(j)
            rowIdxs.append(i)
    
    # Find number of artificial variables to include
    num_artificial = m - len(rowIdxs)
    
    # Store total number of variables for phase 1 in a separate variable
    t = n + num_artificial

    # Make a copy of A
    A1 = [row[:] for row in A]
    
    # Append identity columns to A1 corresponding to artificial variables
    for j in range(num_artificial):
        flag = True
        for i in range(m):
            val = 0
            if flag and i not in rowIdxs:
                val = 1
                rowIdxs.append(i)
                flag = False
            A1[i].append(val)

    # Update the basic variables list to include artificial variables
    for i in range(n, t):
        basic.append(i)

    # Order the basic variables by position of 1 in identity column (based on row)
    ordered_basic = []
    for i in range(m):
        for j in basic:
            if A1[i][j] == 1:
                ordered_basic.append(j)

    # Initialize the objective vector
    e = [0 for i in range(n)]
    for i in range(n, t):
        e.append(1)

    # Initialize vector of decision variables (non-basic variables are 0 and basic variables are elements of b)
    x = initial_x0(ordered_basic, t, b)

    ################################### Initialize tableau form ###################################
    # First row (augmented with objective value in last column)
    table = [[-i for i in e]]
    table[0].append(0)

    # Rows corresponding to slack variables
    for i in range(m):
        table.append(A1[i])

    # Last column or RHS
    for i in range(1, m+1):
        table[i].append(b[i-1])

    # Make basic variables zero in 1st row
    # Find rows to add to first row based on columns corresponding to artificial variables
    idxs = []                                
    for j in range(n, t):
        for i in range(1, m+1):
            if table[i][j] == 1:
                idxs.append(i)
                
    # Add rows to R_0 to complete the operation R_0 --> R_0 + (sum(R_i) for i in idxs)
    for i in idxs:
        for j in range(t + 1):
            table[0][j] += table[i][j]
    ###############################################################################################

    # Run simplex algorithm
    table, basic, status = simplex(table, ordered_basic)
    solution = get_solution(table, basic, len(table[0]) - 1)
    
    ###############################################################################################
    
    # Check for infeasibility 
    # 1. Artificial variables need to be 0 (account for floating point error)
    for i in range(n, t):
        if abs(solution[i]) > 0.001:
            status = "Infeasible"
            break
            
    # 2. Objective value needs to be 0 (account for floating point error)
    if abs(table[0][-1]) > 0.001:
        status = "Infeasible"

    if status != "Infeasible":
        ########################################### Phase 2 ###########################################
        
        # Form a new table with lesser columns without the artificial variables
        new_table = [[-i for i in c]]
        new_table[0].append(0)
        for i in range(1, len(table)):
            # Don't add rows corresponding to artificial variables (in case of redundant constraint)
            if basic[i-1] < n:
                row = table[i][:n]
                row.append(table[i][-1])
                new_table.append(row)

        # Modify basic to exclude artificial variables (in case of redundant constraint)
        new_basic = []
        for j in basic:
            if j < n:
                new_basic.append(j)
        basic = new_basic

        # Make basic variables zero in 1st row
        for j in range(len(basic)):
            basicVarIdx = basic[j]
            pivot = new_table[0][basicVarIdx]
            if pivot != 0:
                rowIdx = j+1   # Since the 1 in the column of a basic variable is located at the corresponding row of the basic variable
                for i in range(n+1):
                    new_table[0][i] -= pivot * new_table[rowIdx][i]

        # Run simplex algorithm
        new_table, basic, status = simplex(new_table, basic)
        # display_result(new_table, basic, status, N)
        table = new_table
        ###############################################################################################
                
else:
    # Initialize vector of decision variables (non-basic variables are 0 and basic variables are elements of b)
    x = [0 for i in range(N)]
    for i in range(m):
        x.append(b[i])

    # Initialize basic variables index set
    basic = [i for i in range(N, n)]

    ################################### Initialize tableau form ###################################
    # First row (augmented with objective value in last column)
    table = [[-i for i in c]]
    table[0].append(0)

    # Rows corresponding to slack variables
    for i in range(m):
        table.append(A[i])

    # Last column or RHS
    for i in range(1, m+1):
        table[i].append(b[i-1])
    
    # Here, basic variables are already 0 in the first row
    ###############################################################################################

    # Run simplex algorithm
    table, basic, status = simplex(table, basic)
    # display_result(table, basic, status, N)

# Store result of initial LP relaxation in variables
solution = get_solution(table, basic, len(table[0]) - 1)
optimal_point = solution[:N]

# Check for integer solution
cutting_plane = False
for val in optimal_point:
    if abs(val - int(val)) > 0.001:
        cutting_plane = True
        break

while status != "Infeasible" and cutting_plane:
    #################################### Cutting Plane Method ####################################

    # Find first cutting plane by finding rowIdx, the row of maximum fractional value in last/RHS column
    maxFrac = float("-inf")
    rowIdx = -1
    for i in range(1, len(table)):
        if basic[i-1] < N:
            val = table[i][-1]
            frac = abs(val - int(val))
            if frac > maxFrac:
                maxFrac = frac
                rowIdx = i
    
    # Error
    # if rowIdx == -1:
    #     print('Row idx is -1')
    #     break

    cp = table[rowIdx][:]

    # Modify the inequality at rowIdx
    for i in range(len(cp)):
        cp[i] = math.floor(cp[i])

    for i in range(len(cp)):
        cp[i] -= table[rowIdx][i]

    # Add slack variable and extra row (new equality with new slack variable)
    cp.append(cp[-1])
    cp[-2] = 1
    for i in range(len(table)):
        table[i].append(table[i][-1])
        table[i][-2] = 0
    table.append(cp)

    # Add new slack variable to basic variables
    basic.append(len(cp) - 2)

    # Run simplex algorithm
    table, basic, status = dual_simplex(table, basic)

    # Store result in variables
    solution = get_solution(table, basic, len(table[0]) - 1)
    optimal_point = solution[:N]

    # Check for integer solution
    cutting_plane = False
    for val in optimal_point:
        if abs(val - int(val)) > 0.001:
            cutting_plane = True
            break

    ##############################################################################################
    
display_result(table, basic, status, N)
