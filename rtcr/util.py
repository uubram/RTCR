# Description: Some utility functions for the Recover T Cell Receptor (RTCR)
# pipeline.
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

from itertools import izip, count
from math import sqrt, log, factorial, pow
import re
from time import time, sleep
from sys import stdout
import sys
from tempfile import mkdtemp
import threading
import multiprocessing as mp
import shutil
import warnings
import signal

codon2aa = {
    "TTT" : "F",
    "TTC" : "F",
    "TTA" : "L",
    "TTG" : "L",
    "CTT" : "L",
    "CTC" : "L",
    "CTA" : "L",
    "CTG" : "L",
    "ATT" : "I",
    "ATC" : "I",
    "ATA" : "I",
    "ATG" : "M",
    "GTT" : "V",
    "GTC" : "V",
    "GTA" : "V",
    "GTG" : "V",
    "TCT" : "S",
    "TCC" : "S",
    "TCA" : "S",
    "TCG" : "S",
    "CCT" : "P",
    "CCC" : "P",
    "CCA" : "P",
    "CCG" : "P",
    "ACT" : "T",
    "ACC" : "T",
    "ACA" : "T",
    "ACG" : "T",
    "GCT" : "A",
    "GCC" : "A",
    "GCA" : "A",
    "GCG" : "A",
    "TAT" : "Y",
    "TAC" : "Y",
    "TAA" : "*",
    "TAG" : "*",
    "CAT" : "H",
    "CAC" : "H",
    "CAA" : "Q",
    "CAG" : "Q",
    "AAT" : "N",
    "AAC" : "N",
    "AAA" : "K",
    "AAG" : "K",
    "GAT" : "D",
    "GAC" : "D",
    "GAA" : "E",
    "GAG" : "E",
    "TGT" : "C",
    "TGC" : "C",
    "TGA" : "*",
    "TGG" : "W",
    "CGT" : "R",
    "CGC" : "R",
    "CGA" : "R",
    "CGG" : "R",
    "AGT" : "S",
    "AGC" : "S",
    "AGA" : "R",
    "AGG" : "R",
    "GGT" : "G",
    "GGC" : "G",
    "GGA" : "G",
    "GGG" : "G"}

def nt2aa(nt_seq):
    n = int(len(nt_seq) / 3)
    return ''.join([codon2aa.get(nt_seq[3*i:3*i + 3],"?") \
            for i in xrange(n)])

def clone2str(clone, fmt):
    """Turn a clone into a string using provided format.

    :fmt: format string
    """
    nt_seq = clone.seq
    aa_seq = nt2aa(nt_seq)
    count = clone.count
    n_stop_codons = sum([aa == '*' for aa in aa_seq])
    min_phred = min(clone.qual)
    frame = len(clone.seq) % 3
    qualstr = '|'.join(map(str, clone.qual))

    if hasattr(clone, "v"):
        vid = clone.v.allele.name
        jid = clone.j.allele.name
        ve = clone.v.end
        js = clone.j.start
    else:
        vid = None
        jid = None
        ve = 0
        js = 0

    js1 = js + 1 # 1-based start nucleotide of J gene
    return fmt%locals()

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
        # Let parent process handle interrupt signal
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        super(ConnectedProcess, self).run()
        self.__cleanup()

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

    def __cleanup(self):
        logger.debug('Closing message queues')
        self._parent_msgs.close()
        self._child_msgs.close()

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
        self._task_queue = task_queue
        self._result_queue = result_queue
        self._initializer = initializer
        self._initargs = initargs

    def run(self):
        super(ConnectedConsumer, self).run()
        if not self._initializer is None:
            self._initializer(*self._initargs)

        try:
            while True:
                    next_task = self._task_queue.get()
                    if next_task is None: # sentinel received, shut down
                        logger.debug('Received sentinel')
                        break
                    result = next_task()
                    self._result_queue.put(result)
        except:
            logger.error('Exception in consumer',
                    exc_info = sys.exc_info())
            raise
        finally:
            self.__cleanup()

    def __cleanup(self):
        logger.debug('Closing task and result queues')
        self._task_queue.close()
        self._result_queue.close()

class ConnectedConsumerPool(object):
    def __init__(self, n_consumers, initializer = None, initargs = ()):
        # NOTE: making the _tasks a JoinableQueue and then joining it resulted
        # in the main process to hang if it tries to exit while a pool is
        # still running (despite using daemon threads).
        self._tasks = mp.Queue()
        self._results = mp.Queue()
        self._ntasks = 0

        # Without cancel_join_thread the pool can cause the main process to
        # hang if the pool is stopped before it finished sending tasks.
        self._tasks.cancel_join_thread()

        self._accept_tasks = True
        self._started = False
        self._done = False
        self._run = threading.Event()
        self._error = threading.Event()

        self._consumer_handler_interval = 1

        self._consumers = [ ConnectedConsumer(self._tasks, self._results,
            initializer = initializer, initargs = initargs) \
                    for i in xrange(n_consumers) ]
        for consumer in self._consumers:
            consumer.daemon = True

    def _can_continue(self):
        return self._run.is_set() and not self._error.is_set()

    def _handle_consumers(self):
        while self._can_continue():
            logger.debug('Checking consumers')
            for consumer in self._consumers:
                # NOTE: it is important to check if the exitcode is not None
                # BEFORE checking its numerical value. In the reverse order a
                # proper consumer exit can be seen as a bad one, because just
                # before the exit the exitcode is None (and exitcode != 0
                # evaluates to True), and if the consumer then exits properly,
                # the exitcode evaluates to 0 (and exitcode is not None).
                if consumer.exitcode is not None and consumer.exitcode != 0:
                    logger.error('Consumer %s died with exitcode: %s'%(
                        consumer, consumer.exitcode))
                    self._error.set()
            sleep(self._consumer_handler_interval)
        logger.debug('exit')

    def _handle_results(self):
        while self._can_continue() and self._ntasks > 0:
            logger.debug('getting result (%s tasks remaining)'%self._ntasks)
            self._collected_results += [self._results.get()]
            self._ntasks -= 1
        logger.debug('exit')

    def _clean_consumers(self):
        logger.debug('Cleaning up consumers')
        for consumer in self._consumers:
            if consumer.exitcode is not None:
                logger.debug('Cleaning consumer %s'%consumer)
                consumer.join()

    def _terminate_consumers(self):
        logger.debug('Terminating consumers')
        for consumer in self._consumers:
            if consumer.exitcode is None:
                logger.debug('Terminating consumer: %s'%consumer)
                consumer.terminate()

    def _cleanup(self):
        logger.debug('Cleaning up pool')
        self._clean_consumers()
        self._consumer_handler.join()
        # NOTE: Not joining self._result_handler because it should already have
        # been joined before the the clean up. And if not, then joining it now
        # can cause the main process to hang as it will be waiting for new
        # results from the already stopped consumers. 
        self._tasks.close()
        self._results.close()

    def add_task(self, func, args):
        if not self._accept_tasks:
            raise ValueError('Pool does not accept new tasks')
        if not type(args) is tuple:
            raise ValueError('args is not a tuple')
        self._tasks.put(Task(func, args))
        self._ntasks += 1

    def start(self):
        if self._started:
            raise ValueError('Pool already started')
        self._started = True
        self._accept_tasks = False
        self._run.set()

        # Add a sentinel to tasks for every consumer
        for i in xrange(len(self._consumers)):
            self._tasks.put(None)

        for consumer in self._consumers:
            consumer.start()

        self._consumer_handler = threading.Thread(
                target = self._handle_consumers)
        self._consumer_handler.daemon = True 
        self._consumer_handler.name = 'Consumer handler'
        self._consumer_handler.start()

        self._collected_results = []

        self._result_handler = threading.Thread(target = self._handle_results)
        self._result_handler.daemon = True
        self._result_handler.name = 'Result handler'
        self._result_handler.start()

    def join(self):
        logger.debug('Joining result handler')
        self._result_handler.join()

        # All results are in, stop and clean up the pool
        self._run.clear()
        self._cleanup()
        if self._ntasks == 0 and not self._error.is_set():
            self._done = True

    def terminate(self):
        if not self._run.is_set():
            raise ValueError('Pool is not running')

        logger.debug('Terminating pool')
        self._run.clear()
        self._terminate_consumers()
        self._cleanup()

    @property
    def results(self):
        if not self._done:
            raise ValueError('Pool not done yet, results may be incomplete')

        return self._collected_results

    @property
    def task_count(self):
        return self._ntasks
