import functools

CACHE_VALUES = False

CACHED_FUNCS = []

#class CustomCache:

def cache_this(func):
    def identity(*args, **kwargs):
        result = func(*args, **kwargs)
        return result

    if CACHE_VALUES:
        wrapper = functools.lru_cache(func)
        CACHED_FUNCS.append(wrapper)
    else:
        wrapper = identity
    return wrapper


# This function currently cleans the whole cache, regardless of any specific object etc.
def clear_cache():
    for func in CACHED_FUNCS:
        if hasattr(func, "cache_clear"):
            func.cache_clear()
