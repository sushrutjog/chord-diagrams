import math
import svgwrite
import os
import random
import json
from collections import deque

#MARK: get support quadrilateral array

def get_supp(pairs_input, moveAbout):
    """
    determine the partners of the neighbours of the move_about pair. get support array of 4 pairs, whose first and third entries will be used to get the index of the conjugate arc which will be inserted back
    """
    k = 2 * len(pairs_input)

    # create a map for partner lookup 
    partnerMap = {p1: p2 for p1, p2 in pairs_input}
    partnerMap.update({p2: p1 for p1, p2 in pairs_input})

    #both these helper functions should be in mod k

    def nextInteger(n):
        return (n % k) + 1

    def prevInteger(n):
        return k if n == 1 else n - 1

    a, b = moveAbout

    #get the four neighbours of the moveAbout pair
    neighbors = [
        nextInteger(a),      # a+1
        nextInteger(b),      # b+1
        prevInteger(a),  # a-1
        prevInteger(b)   # b-1
    ]

    #get partners of neighbours

    partners = sorted([partnerMap[num] for num in neighbors])

    finalPartners = partners

    # MARK: edge cases
    # there are two edge cases:
    # [1, a, b, k] and [1, 2, k-1, k]

    if len(partners) == 4 and partners[0] == 1 and partners[-1] == k:

        # for [1, 2, k-1, k], do nothing
        is_special_case = (partners[1] == 2 and partners[2] == k - 1)

        # for [1, a, b, k], convert to [a, b, k, 1]
        if not is_special_case:
            finalPartners = partners[1:] + [partners[0]]

    # check if array is of the form [a, a+1, b, b+1] in mod k,throw an error if not

    if len(finalPartners) == 4:
        is_first_pair_successor = (nextInteger(finalPartners[0]) == finalPartners[1])
        is_second_pair_successor = (nextInteger(finalPartners[2]) == finalPartners[3])

        # this was used for the second edge case.
        if not (is_first_pair_successor and is_second_pair_successor):
            raise ValueError(f"Support array {finalPartners} is not of the form [a, a+1, b, b+1]. Stopping iterations.")

    return finalPartners

#MARK: reconstruct

def reconstruct(pairsInput, moveAbout, k1, k2):

    # remove the moveAbout pair

    a, b = moveAbout
    min_ab, max_ab = min(a, b), max(a, b)

    remainingPairs = [p for p in pairsInput if p != tuple(sorted(moveAbout)) and p != tuple(sorted(moveAbout))[::-1]]

    def adjustIndexDown(i):
        if i > max_ab:
            return i - 2
        if min_ab < i < max_ab:
            return i - 1
        return i

    pairsSupport = [(adjustIndexDown(p1), adjustIndexDown(p2)) for p1, p2 in remainingPairs]

    # adding the new pair back
    min_k, max_k = min(k1, k2), max(k1, k2)

    def adjustIndexUp(i):
        if i > max_k:
            return i + 2
        if min_k < i <= max_k:
            return i + 1
        return i

    pairsOutput = [(adjustIndexUp(p1), adjustIndexUp(p2)) for p1, p2 in pairsSupport]

    # add the new pair that results from the move
    newPair = (min_k + 1, max_k + 2)
    pairsOutput.append(newPair)

    return pairsOutput, pairsSupport



def decreaseIndexHelper(moveAbout, n):
    # helper function to adjust a single index after a pair is removed.
    a, b = moveAbout
    min_moveAbout, max_moveAbout = min(a, b), max(a, b)

    if n < min_moveAbout:
        return n
    if min_moveAbout < n < max_moveAbout:
        return n - 1
    return n - 2


#MARK: draw svg -- function

def draw_chord_diagram(dwg, pairs, x_offset, highlightPair=None):
    """
    Draws a chord diagram in an SVG group.

    Args:
        dwg: The svgwrite.Drawing object.
        pairs: A list of tuples representing the chords.
        x_offset: The horizontal offset for the diagram.
        highlightPair: A specific pair (a, b) to highlight. If provided, this
                      pair is drawn in red, and chords connected to the
                      neighbors of a and b (a-1, a+1, b-1, b+1) are drawn in green.
    """
    if not pairs:
        print("Empty list!!!!")
        return

    k = 2 * len(pairs)

    # helper functions for working in mod k
    def nextInteger(n):
        return (n % k) + 1

    def prevInteger(n):
        return k if n == 1 else n - 1

    # sort the main highlight pair for easier comparison
    sortedHighlight = tuple(sorted(highlightPair)) if highlightPair else None

    # determining the neighbours of moveAbout = highlight_pair
    neighbor_points = set()
    if highlightPair:
        a, b = highlightPair
        neighbor_points.add(nextInteger(a))
        neighbor_points.add(prevInteger(a))
        neighbor_points.add(nextInteger(b))
        neighbor_points.add(prevInteger(b))

    # MARK: svg core logic

    # map internal coordinates (-1.2 to 1.2) to a 500px box
    scale_factor = 500 / 2.4

    # invert the y-axis
    transform_str = f"translate({x_offset + 250}, 250) scale({scale_factor}, {-scale_factor})"

    main_group = dwg.g(transform=transform_str)
    dwg.add(main_group)

    # draw the circle
    main_group.add(dwg.circle(center=(0, 0), r=1, stroke='black', fill='none', stroke_width=0.015))

    # get the positions of k points on the circumference.
    # '1' is placed at the northernmost point, and points are enumerated counter-clockwise -- CCW
    #MARK: plus = CCW (AFTER inverting the y axis)

    points = []
    for i in range(k):
        angle = (math.pi / 2) + (2 * math.pi * i / k)
        points.append((math.cos(angle), math.sin(angle)))

    # draw the chords
    for i, j in pairs:
        if not (1 <= i <= k and 1 <= j <= k):
            raise ValueError(f"Invalid pair ({i}, {j}). Indices must be between 1 and {k}.")

        # chord colours
        current_pair_normalized = tuple(sorted((i, j)))
        color = 'royalblue'  # default color

        if sortedHighlight:
            if current_pair_normalized == sortedHighlight:
                color = 'red'  # main highlighted pair
            elif i in neighbor_points or j in neighbor_points:
                color = 'orange'  # neighbour chords

        p1 = points[i - 1]
        p2 = points[j - 1]
        x1, y1 = p1
        x2, y2 = p2

        """
        orthogonal circles: if two circles of radii r1 and r2 are d distance apart, then we have that (r1)^2 + (r2)^2 = d^2.

        If two circles having cartesian equaltions

        x^2 + y^2 + 2 g x + 2 f y + c = 0   and
        x^2 + y^2 + 2 g' x + 2 f' y + c' = 0

        are orthogonal, then we have:

        2 g g' + 2 f f' + c + c'
        """

        # check for diametrically opposite points
        if abs(x1 + x2) < 1e-9 and abs(y1 + y2) < 1e-9:
            # draw a straight line if they are
            main_group.add(dwg.line(start=p1, end=p2, stroke=color, stroke_width=0.015))
        else:
            # draw hyperbolic gedesic for others
            denominator = x1 * y2 - x2 * y1
            if abs(denominator) < 1e-9: # Should not happen if points are not opposite
                main_group.add(dwg.line(start=p1, end=p2, stroke=color, stroke_width=0.015))
                continue

            # center and radius for the arc's circle
            cx = (y2 - y1) / denominator
            cy = (x1 - x2) / denominator
            radius = math.sqrt(cx**2 + cy**2 - 1)

            # Determine sweep_flag for correct arc orientation. 
            orientation = (cx - x1) * (y2 - y1) - (cy - y1) * (x2 - x1)
            sweep_flag = 1 if orientation < 0 else 0

            # SVG path data string for the arc
            path_d = f"M {x1} {y1} A {radius} {radius} 0 0 {sweep_flag} {x2} {y2}"
            main_group.add(dwg.path(d=path_d, fill="none", stroke=color, stroke_width=0.015))


#MARK: iterations
# obsolete function

def iterations(initial_pairs_input, num_iterations):

    current_pairs = initial_pairs_input.copy()

    for k in range(1, num_iterations + 1):
        print(f"--- Starting Iteration {k} of {num_iterations} ---")
        try:
            move_about = random.choice(current_pairs)
            print(f"Randomly selected pair to move: {move_about}")

            supp_array = get_supp(current_pairs, move_about)
            k1 = decreaseIndexHelper(move_about, supp_array[0])
            k2 = decreaseIndexHelper(move_about, supp_array[2])
            pairs_output, pairs_support = reconstruct(current_pairs, move_about, k1, k2)

            insert_pair = (min(k1, k2) + 1, max(k1, k2) + 2)
            output_folder = 'coloured_iterations'
            os.makedirs(output_folder, exist_ok=True)
            output_filename = os.path.join(output_folder, f'{k}.svg')
            # dwg = svgwrite.Drawing(output_filename, profile='tiny', size=('1500px', '500px'))
            dwg = svgwrite.Drawing(output_filename, profile='tiny', size=('1000', '500px'))

            draw_chord_diagram(dwg, current_pairs, 0, highlightPair=move_about)

            draw_chord_diagram(dwg, pairs_support, 500)

            draw_chord_diagram(dwg, pairs_output, 1000, highlightPair=insert_pair)

            dwg.save()
            print(f"Success: '{dwg.filename}' saved.\n")

            current_pairs = pairs_output

        except (ValueError, FileNotFoundError, IndexError) as e:
            print(f"An error occurred during iteration {k}: {e}")
            print("Stopping further iterations.")
            break


# MARK: check isometries
# obsolete function

def check_isometry(current_pair, pairs_output):
    """
    checks if two diagrams differ by a rotation. 

    determines if the pairs_output and current_pair differ by a constant c in the range 0 to k-1
    """

    # check if the number of chord are different
    if len(current_pair) != len(pairs_output):
        return False

    # edge case: empty diagrams are always isomorphic
    num_pairs = len(current_pair)
    if num_pairs == 0:
        return True  
    k = 2 * num_pairs

    # create a canonical representation using set -- O(1) average case lookups. 
    # canonical pairing sorts each pair
    canonical_set = {tuple(sorted(p)) for p in current_pair}

    # check for all possible rotational differences
    for c in range(k):
        transformed_pairs = set()

        # apply the current rotation c to the second list
        for p1, p2 in pairs_output:
            # work in mod k
            new_p1 = (p1 - 1 + c) % k + 1
            new_p2 = (p2 - 1 + c) % k + 1

            transformed_pairs.add(tuple(sorted((new_p1, new_p2))))

        # compare the resultant set with the canonical pairing.
        if transformed_pairs == canonical_set:
            return True

    # return false if they don't differ by a rotation
    return False

# MARK: canonical repn

def get_canonical(pairs):
    """
    computes the canonical representation of the chord diagram.

    sort numbers within pairs, and then sort the pairs themselves lexicographically

    Args:
        pairs: A list of tuples representing the chords.

    Returns:
        A tuple of sorted tuples, representing the canonical form.
    """
    if not pairs:
        return tuple()

    k = 2 * len(pairs)
    canonical_form = None

    for c in range(k):
        rotated_pairs = []
        for p1, p2 in pairs:
            # apply rotations
            new_p1 = (p1 - 1 + c) % k + 1
            new_p2 = (p2 - 1 + c) % k + 1
            rotated_pairs.append(tuple(sorted((new_p1, new_p2))))

        # sort the lsit of pairs -- standardize
        rotated_pairs.sort()

        # convert to a tuple
        current_form = tuple(rotated_pairs)

        # if this is the first one generated, keep it
        # if this is 'more' canonical than what we already have, keep it
        if canonical_form is None or current_form < canonical_form:
            canonical_form = current_form

    return canonical_form

#MARK: enumerate diagrams

def enumerate_diagrams(initial_pairs):
    """
    generates ALL possible chord diagrams up to rotations by considering elementary moves about ALL possible chords.

    uses BFS. we maintain the set of the canonical forms of all diagrams seen so far.

    Args:
        initial_pairs: A list of tuples representing the starting chord diagram.
    """
    output_folder = 'enumerate'
    os.makedirs(output_folder, exist_ok=True)

    # queue for BFS. storing tuples of (list, move_about)
    queue = deque()

    # use sets for constant (average case) lookups
    seen_diagrams = set()

    # a list to store all the unique pairings found, stored to a file
    pairings_data = []

    # start with the initial diagram
    initial_canonical = get_canonical(initial_pairs)
    seen_diagrams.add(initial_canonical)
    pairings_data.append(initial_pairs)

    # add all possible chords in the initial diagram to the queue
    for pair in initial_pairs:
        queue.append((initial_pairs, pair))

    svg_counter = 1

    # start BFS
    while queue:
        current_pairs, move_about = queue.popleft()
        try:
            # transform
            supp_array = get_supp(current_pairs, move_about)
            k1 = decreaseIndexHelper(move_about, supp_array[0])
            k2 = decreaseIndexHelper(move_about, supp_array[2])
            output_pairs, pairs_support = reconstruct(current_pairs, move_about, k1, k2)

            # check if the resultant diagram is new
            output_canonical = get_canonical(output_pairs)

            if output_canonical not in seen_diagrams:
                svg_counter += 1
                seen_diagrams.add(output_canonical)
                pairings_data.append(output_pairs)
                print(f"Found a new diagram {svg_counter}")

                # if it is new, create an svg, save sgv.
                output_filename = os.path.join(output_folder, f'{svg_counter}.svg')
                # dwg = svgwrite.Drawing(output_filename, profile='tiny', size=('1500px', '500px'))
                dwg = svgwrite.Drawing(output_filename, profile='tiny', size=('1000px', '500px'))
                insert_pair = (min(k1, k2) + 1, max(k1, k2) + 2)
                draw_chord_diagram(dwg, current_pairs, 0, highlightPair=move_about)
                # draw_chord_diagram(dwg, pairs_support, 500)
                # draw_chord_diagram(dwg, output_pairs, 1000, highlight_pair=insert_pair)
                draw_chord_diagram(dwg, output_pairs, 500, highlightPair=insert_pair)
                dwg.save()

                # add all chords from this new diagram to the queue
                for new_move in output_pairs:
                    queue.append((output_pairs, new_move))
            else:
                print("skipping")

        # when you encounter an error, stop.
        except ValueError as e:
            print(f"Error: {e}")

    # save the unique pairing data to a json
    pairings_filename = 'pairings_data.json'
    with open(pairings_filename, 'w') as f:
        json.dump(pairings_data, f, indent=4)

    print(f"\n Enumeration Complete")
    print(f"Found a total of {len(pairings_data)} unique diagrams.")


#MARK: get only n_g
# the main function used for higher genera

def total_n_g(initial_pairs):
    """
    for more comments, look at the function enumerate_diagrams

    the minimal version of the function enumerate_diagrams. only outputs the total number of diagrams seen up to rotation.

    sends updates every 10,000 new diagrams found.

    peace loving function which does not draw or save anything.

    rest is the same as enumerate_diagrams.

    Args:
        initial_pairs: A list of tuples representing the starting chord diagram.

    Returns:
        An integer representing the total number of unique diagrams found.
    """
    if not initial_pairs:
        return 0
    
    print("Running...")

    # queue for BFS
    queue = deque()
    
    # get a set to store all the canonical pairings seen so far.
    seen_diagrams = set()

    initial_canonical = get_canonical(initial_pairs)
    seen_diagrams.add(initial_canonical)
    for pair in initial_pairs:
        queue.append((initial_pairs, pair))
    
    found_count = 1

    while queue:
        current_pairs, move_about = queue.popleft()

        try:
            supp_array = get_supp(current_pairs, move_about)
            k1 = decreaseIndexHelper(move_about, supp_array[0])
            k2 = decreaseIndexHelper(move_about, supp_array[2])
            
            output_pairs, _ = reconstruct(current_pairs, move_about, k1, k2)
            
            output_canonical = get_canonical(output_pairs)
            
            if output_canonical not in seen_diagrams:
                seen_diagrams.add(output_canonical)
                found_count += 1

                # print progress
                if found_count % 10000 == 0:
                    print(f"{found_count:,} diagrams found")

                for new_move in output_pairs:
                    queue.append((output_pairs, new_move))

        # except (ValueError, IndexError):
        #     continue

        except ValueError as e:
            print(f"Error: {e}")

    # return or print the final output
    # return len(seen_diagrams)
    print(f"Found a total of {len(seen_diagrams)} unique diagrams.")


#MARK: main

if __name__ == '__main__':

    # genus 1:
    # initial_pairs = [(1, 4), (2, 5), (3, 6)]

    # genus 2:
    # initial_pairs = [(2, 5), (18, 15), (17, 8), (16, 7), (9, 6), (1, 10), (4, 13), (3, 12), (11, 14)]

    # genus3:
    initial_pairs = [(1, 16), (2, 5), (3, 28), (4, 29), (6, 15), (7, 10), (8, 12), (9, 13), (11, 14), (17, 26), (18, 21), (19, 23), (20, 24), (22, 25), (27, 30)]

    # genus 4:
    # initial_pairs = [(2, 5), (3, 8), (4, 9), (10, 7), (12, 15), (13, 17), (14, 18), (16, 19), (42, 39), (41, 36), (40, 35), (34, 37), (32, 29), (31, 27), (30, 26), (28, 25), (1, 22), (6, 21), (38, 23), (11, 20), (33, 24)]

    # To run the random iterations:
    # number_of_iterations = 100
    # iterations(initial_pairs, number_of_iterations)
    
    # enumerate_diagrams(initial_pairs)

    total_n_g(initial_pairs)