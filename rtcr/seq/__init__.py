import logging
logger = logging.getLogger(__name__)

try:
    logger.debug("Loading cseq module.")
    from cseq import *
    from seq import Sequence, QSequence, IntervalAnnotation,\
            ConsensusIntervalAnnotation
except:
    logger.warning("Unable to load cseq module, falling back on seq module.")
    from seq import *
