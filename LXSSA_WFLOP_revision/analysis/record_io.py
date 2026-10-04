"""Round-trip serialization for future layouts; does not alter archived CSVs."""
import math


def encode_coordinates(layout):
    """Serialize N x 2 finite coordinates with binary64 round-trip precision."""
    points = list(layout)
    if not points:
        raise ValueError("layout must contain at least one turbine")
    out = []
    for point in points:
        if len(point) != 2:
            raise ValueError("each turbine must have two coordinates")
        x, y = map(float, point)
        if not math.isfinite(x) or not math.isfinite(y):
            raise ValueError("coordinates must be finite")
        out.append(f"{x:.17g} {y:.17g}")
    return ";".join(out)


def decode_coordinates(value):
    points = []
    for row in value.split(";"):
        pair = row.split()
        if len(pair) != 2:
            raise ValueError("invalid coordinate pair")
        x, y = map(float, pair)
        if not math.isfinite(x) or not math.isfinite(y):
            raise ValueError("coordinates must be finite")
        points.append((x, y))
    return points
