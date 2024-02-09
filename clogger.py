# Logging, see here: https://stackoverflow.com/questions/3383865/how-to-log-error-to-file-and-not-fail-on-exception
# and here: https://stackoverflow.com/questions/8050775/using-pythons-logging-module-to-log-all-exceptions-and-errors
# and here: https://stackoverflow.com/questions/6234405/logging-uncaught-exceptions-in-python/16993115#16993115

import logging
import sys

logging.basicConfig(
    filename="/tmp/MarzConverter.log",
    level=logging.DEBUG,
    format="%(asctime)s %(levelname)s %(name)s %(message)s",
)

logger = logging.getLogger(__name__)
handler = logging.StreamHandler(stream=sys.stdout)
logger.addHandler(handler)


def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.critical(
        "Uncaught exception, consider reporting the log (/tmp/MarzConverter.log).",
        exc_info=(exc_type, exc_value, exc_traceback),
    )
