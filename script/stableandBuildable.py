def stable(W, P):
    S = True
    for i in range(1, len(P)):
        P_stack = P[i:]
        W_stack = W[i:]
        
        stack_center = sum([w * p for w, p in zip(W_stack, P_stack)]) / sum(W_stack)
        l_edge = P[i-1] - W[i-1] / 2
        r_edge = P[i-1] + W[i-1] / 2
        if stack_center > r_edge or stack_center < l_edge:
            S = False
            break
    return S


def buildable(W, P):
    B = True
    for i in range(1, len(P)):
        P_stack = P[1:i+1]
        W_stack = W[1:i+1]
        
        stack_center = sum([w * p for w, p in zip(W_stack, P_stack)]) / sum(W_stack)
        l_edge = P[i-1] - W[i-1] / 2
        r_edge = P[i-1] + W[i-1] / 2
        l_edge_bottom = P[0] - W[0] / 2
        r_edge_bottom = P[0] + W[0] / 2
        if P[i] > r_edge or P[i] < l_edge:
            B = False
            break
        elif stack_center > r_edge_bottom or stack_center < l_edge_bottom:
            B = False
            break
    return B


# Example usage:
W=[9, 9, 5, 7, 7, 5, 3]
P=[0, -1, -3, -3, -5, -7, -9]

# Call the functions
is_stable = stable(W, P)
is_buildable = buildable(W, P)

print("Stable:", is_stable)
print("Buildable:", is_buildable)
