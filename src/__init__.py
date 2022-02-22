import logging
from satmut_utils.definitions import LOG_FORMATTER

# Create the root logger for all subpackages, submodules
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Set a console handler
# Logfile handlers will be set in subpackages and submodules to resolve to a user-provided output directory
console_handler = logging.StreamHandler()
console_handler.setFormatter(LOG_FORMATTER)
logger.addHandler(console_handler)
