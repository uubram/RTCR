import logging
logger = logging.getLogger(__name__)

try:
    logger.debug("Loading vtrie module.")
    from vtrie import Trie
except:
    logger.warning("Unable to load vtrie module, falling back on trie module.")
    from trie import Trie
