import logging
import coloredlogs


def cli_log(name="frisk-cli"):
    # prints log to stdout and also saves to specified log file
    logger = logging.getLogger(name)
    formatter = coloredlogs.ColoredFormatter(fmt="%(asctime)s: %(message)s",
                                             datefmt="%H:%M:%S")
    stderrh = logging.StreamHandler()
    stderrh.setFormatter(formatter)
    logger.addHandler(stderrh)
    logger.setLevel(logging.DEBUG)
    return logger

CLI_LOG = cli_log()
