import math
import svgwrite
import os

genus = int(input("Enter the genus -- SHOULD BE EVEN: "))
k = 12 * genus - 6

pairs = []

k2 = k + 2
h = (k2) / 2
p = (k + 6) / 24
pInt = int(p)

for i in range(pInt):
    if i == pInt - 1:
        pairs.append((2 + 10 * i, 5 + 10 * i))
        pairs.append((3 + 10 * i, 7 + 10 * i))
        pairs.append((4 + 10 * i, 8 + 10 * i))
        pairs.append((6 + 10 * i, 9 + 10 * i))
    else:
        pairs.append((2 + 10 * i, 5 + 10 * i))
        pairs.append((3 + 10 * i, 8 + 10 * i))
        pairs.append((4 + 10 * i, 9 + 10 * i))
        pairs.append((10 * (i + 1), 7 + 10 * i))

for i in range(pInt):
    base = 10 * i
    if i == pInt - 1:
        pairs.append((k2 - (2 + base), k2 - (5 + base)))
        pairs.append((k2 - (3 + base), k2 - (7 + base)))
        pairs.append((k2 - (4 + base), k2 - (8 + base)))
        pairs.append((k2 - (6 + base), k2 - (9 + base)))
    else:
        pairs.append((k2 - (2 + base), k2 - (5 + base)))
        pairs.append((k2 - (3 + base), k2 - (8 + base)))
        pairs.append((k2 - (4 + base), k2 - (9 + base)))
        pairs.append((k2 - (10 * (i + 1)), k2 - (7 + base)))

pairs.append((1, int(h)))

limit = 2 * pInt - 2

for i in range(1, limit + 1):
    x_orig = 5 * i + 1
    y_orig = h - i
    pairs.append((x_orig, int(y_orig)))

    x_transformed = k2 - x_orig
    y_transformed = k2 - y_orig
    pairs.append((int(x_transformed), int(y_transformed)))

# done generating the fan triangulation

def getFan(kVal, pairsVal):
    """
    Draws the fan triangulation of the given even genus.

    Args:
        kVal: number of ends depending on the genus, 12g-6.
        pairsVal: the pairing
    """
    if kVal % 2 != 0:
        raise ValueError("please enter an even genus")

    #svg setup
    dwg = svgwrite.Drawing(profile='tiny', size=('500px', '500px'))
    dwg.viewbox(minx=-1.2, miny=-1.2, width=2.4, height=2.4)

    #invert the y-axis
    main_group = dwg.g(transform="scale(1,-1)")
    dwg.add(main_group)

    # draw the circle
    main_group.add(dwg.circle(center=(0, 0), r=1, stroke='black', fill='none', stroke_width=0.015))

    # get the positions of k points on the circumference.
    # '1' is placed at the northernmost point, and points are enumerated counter-clockwise -- CCW
    #MARK: plus = CCW (AFTER inverting the y axis)

    points = []
    for i in range(kVal):
        angle = (math.pi / 2) + (2 * math.pi * i / kVal)
        points.append((math.cos(angle), math.sin(angle)))

    for i, j in pairsVal:
        if not (1 <= i <= kVal and 1 <= j <= kVal):
            raise ValueError(f"Invalid pair ({i}, {j}). Indices must be between 1 and {kVal}.")

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

        if abs(x1 + x2) < 1e-9 and abs(y1 + y2) < 1e-9:
            main_group.add(dwg.line(start=p1, end=p2, stroke='royalblue', stroke_width=0.015))
        else:
            denominator = x1 * y2 - x2 * y1
            if abs(denominator) < 1e-9:
                main_group.add(dwg.line(start=p1, end=p2, stroke='royalblue', stroke_width=0.015))
                continue

            cx = (y2 - y1) / denominator
            cy = (x1 - x2) / denominator
            radius = math.sqrt(cx**2 + cy**2 - 1)
            orientation = (cx - x1) * (y2 - y1) - (cy - y1) * (x2 - x1)
            flag = 1 if orientation < 0 else 0

            path = dwg.path(d=f"M {x1} {y1} A {radius} {radius} 0 0 {flag} {x2} {y2}",
                            fill="none", stroke="royalblue", stroke_width=0.015)
            main_group.add(path)

    return dwg.tostring()

if __name__ == '__main__':
    try:
        svg_output = getFan(k, pairs)
        output_dir = "fan-triangulations"
        os.makedirs(output_dir, exist_ok=True)
        filename = f"genus-{genus}-fan.svg"
        output_filename = os.path.join(output_dir, filename)

        with open(output_filename, "w") as f:
            f.write(svg_output)
        print("successfully generated the fan triangulation")
        # print(pairs)
    except (ValueError, ImportError) as e:
        if isinstance(e, ImportError):
            print(f"Error: {e}")