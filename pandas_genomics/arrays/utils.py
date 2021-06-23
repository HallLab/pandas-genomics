def required_ploidy(n, return_val):
    """
    Decorator for methods on GenotypeArrays that returns a given value if the ploidy is not n
    """

    def decorator(func):
        import functools

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            self = args[0]
            if self.variant.ploidy != n:
                return return_val
            else:
                return func(*args, **kwargs)

        return wrapper

    return decorator
