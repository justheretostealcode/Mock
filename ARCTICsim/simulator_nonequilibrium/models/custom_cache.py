import functools

CACHE_VALUES = True

CACHED_FUNCS = []

def cache_this(func):
    def identity(*args, **kwargs):
        result = func(*args, **kwargs)
        return result

    if CACHE_VALUES:
        wrapper = functools.cache(func)
        CACHED_FUNCS.append(func)
    else:
        wrapper = identity
    return wrapper



def clear_cache():
    for func in CACHED_FUNCS:
        functools.cache(func).cache_clear()