def combine(*args):
    if args == ():
        return []
    elif len(args) == 1:
        return [[a] for a in args[0]]
    else:
        return [a + [b] for a in combine(*(args[:-1])) for b in args[-1]]
