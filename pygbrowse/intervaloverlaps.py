def compute_interval_overlaps(A, B, min_overlap=0, is_sorted=False):
    """
    Given two lists, A and B, of interval tuples in the form (name, start, end)
    return a list of tuples of the form:

    (A_name, B_name, overlap_start, overlap_end)

    for every pair of intervals that overlaps by at least <min_overlap>
    """
    overlaps = []

    if not is_sorted:
        A = sorted(A, key=lambda x: x[1])
        B = sorted(B, key=lambda x: x[1])

    A_ptr = 0
    B_ptr = 0

    # Two conditions must hold to have an overlap:
    # B start <= A end
    # A start <= B end

    # initialize the loop by checking for overlaps at the start (all B intervals that overlap with the first interval in A)
    if B[B_ptr][1] <= A[A_ptr][2]:
        while B[B_ptr][1] <= A[A_ptr][2]:
            overlap_start = max(A[A_ptr][1], B[B_ptr][1])
            overlap_end = min(A[A_ptr][2], B[B_ptr][2])
            if overlap_end - overlap_start >= min_overlap:
                overlaps.append((A[A_ptr][0], B[B_ptr][0], overlap_start, overlap_end))
            if B_ptr < len(B) - 1:
                B_ptr += 1
            else:
                break
        B_ptr -= 1
    #     print('initialized: {}, {}, {}, {} '.format(A_ptr, B_ptr, A[A_ptr], B[B_ptr]))
    # advance the A pointer until B start is upstream of A end
    while True:
        #         print('\tnew run: {}, {}, {}, {} '.format(A_ptr, B_ptr, A[A_ptr], B[B_ptr]))
        if A_ptr < len(A) - 1:
            A_ptr += 1
        #             print('\tadvanced A: {}, {}, {}, {} '.format(A_ptr, B_ptr, A[A_ptr], B[B_ptr]))
        else:
            break
        # advance the B pointer until A start is upstream of B end
        while A[A_ptr][1] > B[B_ptr][2]:
            if B_ptr < len(B) - 1:
                B_ptr += 1
            #                 print('\tadvanced B: {}, {}, {}, {} '.format(A_ptr, B_ptr, A[A_ptr], B[B_ptr]))
            else:
                break
        #         print('aligned: {}, {}, {}, {} '.format(A_ptr, B_ptr, A[A_ptr], B[B_ptr]))
        # capture the overlaps in B until B start is no longer upstream of A end
        if B[B_ptr][1] <= A[A_ptr][2]:
            while B[B_ptr][1] <= A[A_ptr][2]:
                overlap_start = max(A[A_ptr][1], B[B_ptr][1])
                overlap_end = min(A[A_ptr][2], B[B_ptr][2])
                #                 print('grabbing: {}, {}, {}, {} '.format(A_ptr, B_ptr, A[A_ptr], B[B_ptr]))
                if overlap_end - overlap_start >= min_overlap:
                    overlaps.append((A[A_ptr][0], B[B_ptr][0],
                                     overlap_start, overlap_end))
                if B_ptr < len(B) - 1:
                    B_ptr += 1
                #                     print('\tadvanced B: {}, {}, {}, {} '.format(A_ptr, B_ptr, A[A_ptr], B[B_ptr]))
                else:
                    break
            #             A_ptr += 1
            B_ptr -= 1
    return overlaps


def test():
    test_A = [('A_1', 4, 8), ('A_2', 11, 19)]
    test_B = [('B_1', 1, 6), ('B_2', 9, 10), ('B_3', 13, 15), ('B_4', 17, 20)]
    print(compute_interval_overlaps(test_A, test_B))
