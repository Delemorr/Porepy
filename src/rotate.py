import math

def rotate(origin, point, angle):
    """
    Поворот точки против часовой стрелки на заданный угол вокруг заданного начала координат.
     Угол следует указывать в радианах.
    """
    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy
