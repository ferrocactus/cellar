import logging

FORMAT = '%(levelname)s:%(name)s: %(message)s'
logging.basicConfig(format=FORMAT)


def setup_logger(name, lvl='INFO'):
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, lvl))
    return logger
