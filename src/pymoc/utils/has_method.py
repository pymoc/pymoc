
def has_method(obj, method_name):
  return hasattr(obj, method_name) and callable(getattr(obj, method_name))