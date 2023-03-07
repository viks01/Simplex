import numpy as np

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
    
    # Check for leaving variable along the column of entering variable (maxIdx)
    # If maxVal <= 0, then the while loop is skipped and the initial basis selected is optimal
    while maxVal > 0:
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
        table[minIdx] /= y

        # Update other rows
        for i in range(len(table)):
            if i != minIdx:
                pivot = table[i][maxIdx]
                table[i] -= pivot * table[minIdx]
        
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
        print(status)
        
        # Initial check of first row: (c_B)(A_B^-1)(a_j) - (c_j) > 0
        maxVal = float("inf")
        maxIdx = -1                                  # index of entering variable
        for i in range(len(table[0]) - 1):
            val = table[0][i]
            if val > 0 and val < maxVal:             # entering variable has smallest positive reduced cost
                maxVal = val
                maxIdx = i
        
        # Check for leaving variable along the column of entering variable (maxIdx)
        # If maxIdx < 0, then the while loop is skipped and the optimal solution is reached
        while maxIdx >= 0:
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
            table[minIdx] /= y

            # Update other rows
            for i in range(len(table)):
                if i != minIdx:
                    pivot = table[i][maxIdx]
                    table[i] -= pivot * table[minIdx]
            
            # Initial check of first row: (c_B)(A_B^-1)(a_j) - (c_j) > 0
            maxVal = float("inf")
            maxIdx = -1                                  # index of entering variable
            for i in range(len(table[0]) - 1):
                val = table[0][i]
                if val > 0 and val < maxVal:             # entering variable has smallest positive reduced cost
                    maxVal = val
                    maxIdx = i

    return table, basic, status

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
def display_result(table, basic, status, N=0):
    if N > 0 and (status == "Success" or status == "Cycling detected"):
        solution = get_solution(table, basic, len(table[0]) - 1)
        res = solution[:N]
        
        print(table[0][-1])
        print(*res)
    else:
        print(status)

# Input
N, u, v = list(map(int, input().strip().split()))                                          # N - dimension of c
c = np.array(list(map(float, input().strip().split())))                                    # objective coefficients
U = np.array([np.array([float(i) for i in input().strip().split()]) for x in range(u)])    # <= inequalities
V = np.array([np.array([float(i) for i in input().strip().split()]) for x in range(v)])    # >= inequalities
b = np.array(list(map(float, input().strip().split())))                                    # RHS of inequalities

# A --> m x n
m = u + v      # total number of slack variables
n = N + m      # total number of variables (objective + slack)

# Augment objective vector by m values
c = np.append(c, [0] * m)

# Combine U and V matrices into one A matrix
if u <= 0:
    A = np.copy(V)
elif v <= 0:
    A = np.copy(U)
else:
    A = np.copy(U)
    A = np.append(A, V, axis=0)

# Augment each inequality by m values for slack variables, converting the inequalities to equalities
I = np.identity(m)
A = np.append(A, I, axis=1)

# Check for two phased method
twoPhased = False
# Case 1: <= inequalities
if np.any(b[:u] < 0):
    twoPhased = True
# Case 2: >= inequalities
if np.any(b[u:m] > 0):
    twoPhased = True

# Make b >= 0 without changing the equalities in Ax = b
for i in range(m):
    if b[i] < 0:
        b[i] *= -1
        A[i] *= -1

x = []
basic = []

if twoPhased:
    ########################################### Phase 1 ###########################################
    
    # Check for identity columns (basic variables) already existing among the slack variables
    basic = np.array([])
    rowIdxs = np.array([])                                # Keep track of rows of A to skip adding ones

    for i in range(m):
        j = i + N
        if A[i][j] == 1:
            basic = np.append(basic, [j])
            rowIdxs = np.append(rowIdxs, [i])

    # Convert basic and rowIdxs to integer arrays
    basic = basic.astype('i')
    rowIdxs = rowIdxs.astype('i')
    
    # Find number of artificial variables to include
    num_artificial = m - len(rowIdxs)
    
    # Store total number of variables for phase 1 in a separate variable
    t = n + num_artificial

    # Make a copy of A
    A1 = np.copy(A)
    
    # Append identity columns to A1 corresponding to artificial variables
    for j in range(num_artificial):
        flag = True
        for i in range(m):
            val = 0
            if flag and i not in rowIdxs:
                val = 1
                rowIdxs = np.append(rowIdxs, [i])
                flag = False
            A1[i] = np.append(A1[i], [val])

    # Update the basic variables list to include artificial variables
    for i in range(n, t):
        basic = np.append(basic, [i])

    # Order the basic variables by position of 1 in identity column (based on row)
    ordered_basic = np.array([])
    for i in range(m):
        for j in basic:
            if A1[i, j] == 1:
                ordered_basic = np.append(ordered_basic, [j])

    # Initialize the objective vector
    e = np.array([0] * n)
    ones = [1] * num_artificial
    e = np.append(e, ones)

    # Initialize vector of decision variables (non-basic variables are 0 and basic variables are elements of b)
    x = np.array([])
    i = 0
    for j in range(t):
        val = 0
        if j in ordered_basic:
            val = b[i]
            i += 1
        x = np.append(x, [val])

    ################################### Initialize tableau form ###################################
    # First row 
    table = np.array([[-i for i in e]])

    # Rows corresponding to slack variables
    table = np.append(table, A1, axis=0)

    # Last column or RHS (augmented with objective value in 1st row)
    col = np.array([0])
    col = np.append(col, b)
    col = [[i] for i in col]
    table = np.append(table, col, axis=1)

    # Make basic variables zero in 1st row
    # Find rows to add to first row based on columns corresponding to artificial variables
    idxs = []                                
    for j in range(n, t):
        for i in range(1, m+1):
            if table[i][j] == 1:
                idxs.append(i)
                
    # Add rows to R_0 to complete the operation R_0 --> R_0 + (sum(R_i) for i in idxs)
    for i in idxs:
        table[0] += table[i]
            
    # Directly makes artificial variable coefficients in 1st row 0, irrespective of objective vector e 
    # for j in range(n, t)
    #     pivot = table[0][j]
    #     if pivot != 0:
    #         idx = 0
    #         for i in range(1, m+1):
    #             if table[i][j] == 1:
    #                 idx = i
    #                 break
    #         for i in range(t + 1):
    #             table[0][i] -= pivot * table[idx][i]
    ###############################################################################################

    # Run simplex algorithm
    table, basic, status = simplex(table, ordered_basic)
    solution = get_solution(table, basic, len(table[0]) - 1)

    # Convert basic to integer array
    basic = basic.astype('i')
    
    ###############################################################################################
    
    # Check for infeasibility 
    # 1. Artificial variables need to be 0 (account for floating point error)
    infeasible = False
    for i in range(n, t):
        if abs(solution[i]) > 0.001:
            infeasible = True
            status = "Infeasible"
            break
            
    # 2. Objective value needs to be 0 (account for floating point error)
    if abs(table[0][-1]) > 0.001:
        infeasible = True
        status = "Infeasible"

    if infeasible:
        display_result(table, basic, status)
    else:
        ########################################### Phase 2 ###########################################

        # Form a new table with lesser columns without the artificial variables
        first_row = np.array([-i for i in c])
        first_row = np.append(first_row, [0])
        new_table = np.array([first_row])
        for i in range(1, len(table)):
            row = table[i][:n]
            row = np.append(row, [table[i][-1]])
            row = np.array([row])
            new_table = np.append(new_table, row, axis=0)

        # Make basic variables zero in 1st row
        for j in basic:
            pivot = new_table[0][j]
            if pivot != 0:
                idx = 0
                for i in range(1, m+1):
                    if new_table[i][j] == 1:
                        idx = i
                        break
                new_table[0] -= pivot * new_table[idx]

        # Run simplex algorithm and show output
        new_table, basic, status = simplex(new_table, basic)
        display_result(new_table, basic, status, N)
        ###############################################################################################
                
else:
    # Initialize vector of decision variables (non-basic variables are 0 and basic variables are elements of b)
    x = np.array([0] * N)
    x = np.append(x, b)

    # Initialize basic variables index set
    basic = np.array([i for i in range(N, n)])

    ################################### Initialize tableau form ###################################
    # First row
    table = np.array([[-i for i in c]])

    # Rows corresponding to slack variables
    table = np.append(table, A, axis=0)

    # Last column or RHS (augmented with objective value in 1st row)
    col = np.array([0])
    col = np.append(col, b)
    col = [[i] for i in col]
    table = np.append(table, col, axis=1)
    ###############################################################################################

    # Run simplex algorithm and show output
    table, basic, status = simplex(table, basic)
    display_result(table, basic, status, N)
