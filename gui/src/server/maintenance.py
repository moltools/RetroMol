"""Module for running maintenance on items stored in Redis."""

import time
import logging
import os

from routes.session_store import mark_stale_processing_items


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


INTERVAL_SECONDS = int(os.getenv("JOB_WATCHDOG_INTERVAL_SECONDS", "130"))


def main() -> None:
    """
    Main loop for running maintenance tasks.
    """
    logger.info(f"Starting maintenance loop with interval {INTERVAL_SECONDS} seconds")
    while True:
        try:
            updated = mark_stale_processing_items()
            logger.info(f"Marked {updated} stale processing items as timeout error")
        except Exception as e:
            logger.exception(f"Error during maintenance loop: {e}")
        time.sleep(INTERVAL_SECONDS)


if __name__ == "__main__":
    main()