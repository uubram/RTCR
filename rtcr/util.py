# Description: Some utility functions for the Recover T Cell Receptor (RTCR)
# pipeline.
import logging
logger = logging.getLogger(__name__)

from itertools import izip, count
from math import sqrt, log, factorial, pow
import re
from time import time
from sys import stdout
from tempfile import mkdtemp
import multiprocessing as mp
import shutil
import warnings

try:
    from scipy import stats
    binom_pmf = stats.binom.pmf
    norm_ppf = stats.norm.ppf
except:
    logger.warning("Unable to load scipy")

    # Binomial coefficient
    def nCk(n,k):
        if k > n:
            return 0
        return factorial(n) / factorial(k) / factorial(n-k)

    # Binomial probability mass function
    def binom_pmf(x,size,prob):
        return nCk(size,x)*pow(prob,x)*pow(1-prob,size-x)

    # Normal variable, percent point function (inverse of cdf)
    # Based on:
    # https://svn.r-project.org/R/trunk/src/nmath/qnorm.c
    def norm_ppf(p):
        assert(0 < p < 1)
        q = p - .5
        if .075 <= p <= .925:
            r = .180625 - q * q
            val = q * (((((((r * 2509.0809287301226727 + \
    33430.575583588128105) * r + 67265.770927008700853) * r + \
    45921.953931549871457) * r + 13731.693765509461125) * r + \
    1971.5909503065514427) * r + 133.14166789178437745) * r + \
    3.387132872796366608) \
    / (((((((r * 5226.495278852854561 + \
    28729.085735721942674) * r + 39307.89580009271061) * r + \
    21213.794301586595867) * r + 5394.1960214247511077) * r + \
    687.1870074920579083) * r + 42.313330701600911252) * r + 1.)
        else: # closer than 0.075 from {0,1} boundary
            # r = min(p, 1-p) < 0.075
            if q > 0:
                r = 1-p
            else:
                r = p
            r = sqrt(-log(r))
            if r <= 5.: # <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11
                r += -1.6
                val = (((((((r * 7.7454501427834140764e-4 + \
    .0227238449892691845833) * r + .24178072517745061177) * \
    r + 1.27045825245236838258) * r + 3.64784832476320460504) * \
    r + 5.7694972214606914055) * r + 4.6303378461565452959) * r + \
    1.42343711074968357734)\
    / (((((((r *\
    1.05075007164441684324e-9 + 5.475938084995344946e-4) * \
    r + .0151986665636164571966) * r + \
    .14810397642748007459) * r + .68976733498510000455) * \
    r + 1.6763848301838038494) * r + 2.05319162663775882187) * r + 1.)
            else: # very close to  0 or 1
                r += -5.
                val = (((((((r * 2.01033439929228813265e-7 + \
    2.71155556874348757815e-5) * r + \
    .0012426609473880784386) * r + .026532189526576123093) * \
    r + .29656057182850489123) * r + \
    1.7848265399172913358) * r + 5.4637849111641143699) * \
    r + 6.6579046435011037772) \
    / (((((((r * \
    2.04426310338993978564e-15 + 1.4215117583164458887e-7)* \
    r + 1.8463183175100546818e-5) * r + \
    7.868691311456132591e-4) * r + .0148753612908506148525) \
    * r + .13692988092273580531) * r + \
    .59983220655588793769) * r + 1.)
            if q < 0.0:
                val = -val
        return val

# TemporaryDirectory class, its code is based on
# https://hg.python.org/cpython/file/3.5/Lib/tempfile.py
class TemporaryDirectory(object):
    """Create and return a temporary directory.  This has the same
    behavior as mkdtemp but can be used as a context manager.  For
    example:

        with TemporaryDirectory() as tmpdir:
            ...

    Upon exiting the context, the directory and everything contained
    in it are removed.
    """

    def __init__(self, suffix = "", prefix = "tmp", dir = None):
        self._closed = False
        self.name = mkdtemp(suffix, prefix, dir)

    def __repr__(self):
        return "<{} {!r}>".format(self.__class__.__name__, self.name)

    def __enter__(self):
        return self.name

    def __exit__(self, exc, value, tb):
        self.cleanup()

    def cleanup(self):
        shutil.rmtree(self.name)

def enum(*names):
    return type("Enum", (object,), dict(izip(names, count())))

# Class that uses two message queues for communication between parent and
# child processes.
class ConnectedProcess(mp.Process):
    def __init__(self, *args, **kwargs):
        super(ConnectedProcess, self).__init__(*args, **kwargs)
        # assuming this is the parent process unless and until the run()
        # method is invoked.
        self._is_parent = True
        self._parent_msgs = mp.Queue() # holds messages produced by parent
        self._child_msgs = mp.Queue() # holds messages produced by child

    def run(self):
        self._is_parent = False
        super(ConnectedProcess, self).run()

    def send(self, msg, block = True):
        if self._is_parent:
            self._parent_msgs.put(msg, block = block)
        else:
            self._child_msgs.put(msg, block = block)

    def recv(self, block = True):
        if self._is_parent:
            return self._child_msgs.get(block = block)
        else:
            return self._parent_msgs.get(block = block)

    def send_nowait(self, msg):
        self.send(msg, block = False)

    def recv_nowait(self):
        return self.recv(block = False)

class Task(object):
    def __init__(self, func, args):
        self.func = func
        self.args = args

    def __call__(self):
        return self.func(*self.args)
 
class ConnectedConsumer(ConnectedProcess):
    def __init__(self, task_queue, result_queue, initializer = None, \
            initargs = ()):
        super(ConnectedConsumer, self).__init__()
        self.task_queue = task_queue
        self.result_queue = result_queue
        self._initializer = initializer
        self._initargs = initargs

    def run(self):
        super(ConnectedConsumer, self).run()

        if not self._initializer is None:
            self._initializer(*self._initargs)

        while True:
            next_task = self.task_queue.get()
            if next_task is None: # sentinel received, shut down
                self.task_queue.task_done()
                break
            result = next_task()
            self.task_queue.task_done()
            self.result_queue.put(result)

class ConnectedConsumerPool(object):
    def __init__(self, n_consumers, initializer = None, initargs = ()):
        self._tasks = mp.JoinableQueue()
        self._ntasks = 0
        self._results = mp.Queue()
        self._state = "OPEN"
        self._done = False

        self._consumers = [ ConnectedConsumer(self._tasks, self._results,
            initializer = initializer, initargs = initargs) \
                    for i in xrange(n_consumers) ]
        for consumer in self._consumers:
            consumer.daemon = True

    def add_task(self, func, args):
        assert(self._state == "OPEN")
        assert(type(args) is tuple)
        self._tasks.put(Task(func, args))
        self._ntasks += 1

    def start(self):
        assert(self._state == "OPEN")
        if self._state == "OPEN":
            self._state = "CLOSE"

        # Add a sentinel to tasks for every consumer
        for i in xrange(len(self._consumers)):
            self._tasks.put(None)

        for consumer in self._consumers:
            consumer.start()

    def join(self):
        assert(self._state in ("CLOSE", "TERMINATE"))
        if self._state == "CLOSE":
            self._tasks.join()
        elif self._state == "TERMINATE":
            for consumer in self._consumers:
                consumer.join()

    def all_tasks_done(self):
        return self._state == "CLOSE" and self._tasks.empty()

    def terminate(self):
        assert(self._state == "CLOSE")
        self._state = "TERMINATE"
        for consumer in self._consumers:
            if consumer.exitcode is None:
                consumer.terminate()
                consumer.join()

    def get_consumers(self):
        assert(self._state == "CLOSE")
        return self._consumers

    def get_results(self):
        assert(self._state == "CLOSE")
        res = []
        while self._ntasks > 0:
            res += [self._results.get()]
            self._ntasks -= 1
        return res
