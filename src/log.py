import logging

FORMAT = '%(levelname)s:%(name)s: %(message)s'
#logging.basicConfig(filename='cellar.log', format=FORMAT)
logging.basicConfig(format=FORMAT)


def setup_logger(name, lvl='INFO'):
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, lvl))
    return logger
