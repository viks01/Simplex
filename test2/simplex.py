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
        for i in range(len(table[0])):
            table[minIdx][i] /= y

        # Update other rows
        for i in range(len(table)):
            if i != minIdx:
                pivot = table[i][maxIdx]
                for j in range(len(table[0])):
                    table[i][j] -= pivot * table[minIdx][j]

        # print(*table, sep='\n')
        # print()
        
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
            for i in range(len(table[0])):
                table[minIdx][i] /= y

            # Update other rows
            for i in range(len(table)):
                if i != minIdx:
                    pivot = table[i][maxIdx]
                    for j in range(len(table[0])):
                        table[i][j] -= pivot * table[minIdx][j]

            # print(*table, sep='\n')
            # print()
            
            # Initial check of first row: (c_B)(A_B^-1)(a_j) - (c_j) > 0
            maxVal = float("inf")
            maxIdx = -1                                  # index of entering variable
            for i in range(len(table[0]) - 1):
                val = table[0][i]
                if val > 0 and val < maxVal:             # entering variable has smallest positive reduced cost
                    maxVal = val
                    maxIdx = i

    return table, basic, status