def add(x, y):
    if not isinstance(x, (int, float)) or not isinstance(y, (int, float)):
      raise TypeError("Inputs must be numbers")
    return x + y