# Check for 2-phased method
# flag = False
# for i in range(m-1, u-1, -1):
#     if b[i] > 0:
#         flag = True
#         break

# Convert >= inequalities to <= inequalities
# for i in range(v):
#     for j in range(N):
#         V[i][j] *= -1
# for i in range(m-1, u-1, -1):
#     b[i] *= -1

    # for i in range(N):
    #     row = []
    #     for j in range(n):
    #         val = 1 if i == j else 0
    #         row.append(val)
    #     table.append(row)


# Check for two phased method
# flag = False
# for i in range(m):
#     if A[i][i + N] < 0:                         # Check if slack variables can be basic variables
#         flag = True
#         break
# for j in range(n):
#     check = True
#     for i in range(m):
#         if A[i][j] == 1 and check:
#             check = False
#         elif A[i][j] != 0:
#             flag = True
#             break

# Convert some >= inequalities into <= inequalities in matrix A while keeping b >= 0
# for i in range(u, m):
#     if b[i] < 0:
#         b[i] *= -1
#         for j in range(n):
#             A[i][j] *= -1


    # rowIdx = -1
    # for j in range(N, n):
    #     zeroes = 0
    #     ones = 0
    #     rowIdx = -1
    #     for i in range(m):
    #         val = A[i][j]
    #         if val == 0:
    #             zeroes += 1
    #         elif val == 1:
    #             ones += 1
    #             rowIdx = i
    #     if zeroes == m-1 and ones == 1:
    #         basic.append(j)
    #         rowIdxs.append(rowIdx)
    #     if len(basic) == m:
    #         break

        # for j in basic:
        #     pivot = new_table[0][j]
        #     if pivot != 0:
        #         idx = 0
        #         for i in range(1, m+1):
        #             if new_table[i][j] == 1:
        #                 idx = i
        #                 break
        #         for i in range(n+1):
        #             new_table[0][i] -= pivot * new_table[idx][i]

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