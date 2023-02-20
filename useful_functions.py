from itertools import tee, islice, chain, takewhile, repeat
import gzip


def previousAndNext(some_iterable):
    # get previous and next item from iterable

    prevs, items, nexts = tee(some_iterable, 3)
    prevs = chain([None], prevs)
    nexts = chain(islice(nexts, 1, None), [None])
    return zip(prevs, items, nexts)


def rawInCount(filename):
    # count total number of lines in a file

    f = open(filename, 'rb')
    bufgen = takewhile(lambda x: x, (f.raw.read(1024*1024) for _ in repeat(None)))
    return sum(buf.count(b'\n') for buf in bufgen)
