"""
ExomeFlow — Whole Exome Sequencing analysis pipeline.

Author: Robin Kumar, AIIMS New Delhi
"""

__version__ = "2.2.13"
__author__ = "Robin Kumar"
__email__ = "itsrobintomar@gmail.com"

# Registers the custom logging.Logger.success() method (SUCCESS level,
# between INFO and WARNING) as a side effect of importing this module.
# Every step module calls logger.success(...) unconditionally, so this must
# happen on `import exomeflow.<anything>`, not just when exomeflow.logger
# happens to already be imported by something else first — found via audit:
# `import exomeflow.cnv` alone (which never imports exomeflow.logger) used
# to raise AttributeError the first time it logged a completion message.
from exomeflow import logger as _logger  # noqa: E402,F401
