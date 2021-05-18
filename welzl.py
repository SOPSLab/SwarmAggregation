# Project:     SwarmAggregation
# Filename:    welzl.py
# Authors:     Joshua J. Daymude (jdaymude@asu.edu) and Noble C. Harasha
#              (nharasha1202@gmail.com).

"""
welzl: An implementation of Welzl's smallest enclosing circle algorithm based
       on that of the Nayuki Project. The original code can be found at:
       https://www.nayuki.io/page/smallest-enclosing-circle
"""

# Copyright (c) 2020 Project Nayuki
# https://www.nayuki.io/page/smallest-enclosing-circle
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program (see COPYING.txt and COPYING.LESSER.txt).
# If not, see <http://www.gnu.org/licenses/>.

from math import hypot
import random


def welzl(points):
    """
    Takes as input an N x 2 array representing (x, y) positional data and
    returns the center and radius of the smallest enclosing circle.
    """
    # Convert to float and randomize order
    shuffled = [(float(x), float(y)) for (x, y) in points]
    random.shuffle(shuffled)

    # Progressively add points to circle or recompute circle
    c = None
    for (i, p) in enumerate(shuffled):
        if c is None or not in_circle(p, c):
            c = circle_from_one(shuffled[:i+1], p)
    return c


def circle_from_one(points, p):
    """
    Finds the smallest enclosing circle when one boundary point is known.
    """
    c = (p[0], p[1], 0.0)
    for (i, q) in enumerate(points):
        if not in_circle(q, c):
            if c[2] == 0.0:
                c = diameter(p, q)
            else:
                c = circle_from_two(points[:i+1], p, q)
    return c


def circle_from_two(points, p, q):
    """
    Finds the smallest enclosing circle when two boundary points are known.
    """
    circ = diameter(p, q)
    left  = None
    right = None
    px, py = p
    qx, qy = q

    # For each point not in the two-point circle, form a circumcircle and
    # classify it on left or right side.
    for r in points:
        if in_circle(r, circ):
            continue

        cross = cross_prod(px, py, qx, qy, r[0], r[1])
        c = circumcircle(p, q, r)
        if c is None:
            continue
        elif cross > 0.0 and (left is None or cross_prod(px, py, qx, qy, c[0], c[1]) > cross_prod(px, py, qx, qy, left[0], left[1])):
            left = c
        elif cross < 0.0 and (right is None or cross_prod(px, py, qx, qy, c[0], c[1]) < cross_prod(px, py, qx, qy, right[0], right[1])):
            right = c

    # Select which circle to return.
    if left is None and right is None:
        return circ
    elif left is None:
        return right
    elif right is None:
        return left
    else:
        return left if (left[2] <= right[2]) else right


def diameter(a, b):
    """
    Finds the midpoint and distance to center between two points.
    """
    cx = (a[0] + b[0]) / 2
    cy = (a[1] + b[1]) / 2
    r0 = hypot(cx - a[0], cy - a[1])
    r1 = hypot(cx - b[0], cy - b[1])
    return (cx, cy, max(r0, r1))


def circumcircle(a, b, c):
    """
    Finds the inscribed circumcircle from three points.
    """
    ox = (min(a[0], b[0], c[0]) + max(a[0], b[0], c[0])) / 2
    oy = (min(a[1], b[1], c[1]) + max(a[1], b[1], c[1])) / 2
    ax = a[0] - ox;  ay = a[1] - oy
    bx = b[0] - ox;  by = b[1] - oy
    cx = c[0] - ox;  cy = c[1] - oy
    d = (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by)) * 2.0
    if d == 0.0:
        return None
    x = ox + ((ax*ax + ay*ay) * (by - cy) + (bx*bx + by*by) * (cy - ay) + (cx*cx + cy*cy) * (ay - by)) / d
    y = oy + ((ax*ax + ay*ay) * (cx - bx) + (bx*bx + by*by) * (ax - cx) + (cx*cx + cy*cy) * (bx - ax)) / d
    ra = hypot(x - a[0], y - a[1])
    rb = hypot(x - b[0], y - b[1])
    rc = hypot(x - c[0], y - c[1])
    return (x, y, max(ra, rb, rc))


def in_circle(p, c):
    """
    Returns True if and only if point p is in circle c.
    """
    return c is not None and hypot(p[0] - c[0], p[1] - c[1]) <= c[2] * (1+1e-14)


def cross_prod(px, py, qx, qy, rx, ry):
    """
    Returns twice the signed area of the triangle defined by the three points.
    """
    return (qx - px) * (ry - py) - (qy - py) * (rx - px)
