import numpy as np

W1 = np.array([8, 7, 5])
P1 = np.array([0, 1, 5])

W2 = np.array([8, 8, 5])
P2 = np.array([0, 3, 6])

W3 = np.array([9, 9, 5])
P3 = np.array([0, 5, 2.25])

W4 = np.array([9, 9, 5, 7, 7, 5, 3])
P4 = np.array([0, -1, -3, -3, -5, -7, -9])

build_corr = np.array([0, 0, 0, 1])
stable_corr = np.array([0, 0, 1, 1])

def buildablegen(W, P):
    B = True
    for i in range(1, len(P)):
        P_stack = P[1:i+1]
        W_stack = W[1:i+1]

        stack_center = np.sum(W_stack * P_stack) / np.sum(W_stack)

        bottom_l_edge = P[0] - W[0] / 2
        bottom_r_edge = P[0] + W[0] / 2

        below_l_edge = P[i-1] - W[i-1] / 2
        below_r_edge = P[i-1] + W[i-1] / 2

        if stack_center > bottom_r_edge or stack_center < bottom_l_edge:
            B = False
            break
        elif P[i] > below_r_edge or P[i] < below_l_edge:
            B = False
            break
    
    return B

def stablegen(W, P):
    S = True
    for i in range(1, len(P)):
        P_stack = P[i:]
        W_stack = W[i:]

        stack_center = np.sum(W_stack * P_stack) / np.sum(W_stack)
        l_edge = P[i-1] - W[i-1] / 2
        r_edge = P[i-1] + W[i-1] / 2

        if stack_center > r_edge or stack_center < l_edge:
            S = False
            break
    
    return S

test_build = np.array([buildablegen(W1, P1), buildablegen(W2, P2), buildablegen(W3, P3), buildablegen(W4, P4)])
test_stabl = np.array([stablegen(W1, P1), stablegen(W2, P2), stablegen(W3, P3), stablegen(W4, P4)])

print(test_build, test_stabl)

if np.array_equal(build_corr, test_build):
    print('build test passed')

if np.array_equal(stable_corr, test_stabl):
    print('stable test passed')