"""This module defines global constants and configuration for the RetroMol package."""

import logging
import os

# Global logger name for RetroMol
LOGGER_NAME = "retromol"
LOGGER_LEVEL = int(os.getenv("LOGGER_LEVEL", logging.INFO))


# Default timeout for running RetroMol on a single molecule
DEFAULT_TIMEOUT_RUN_RETROMOL = 5  # seconds
TIMEOUT_RUN_RETROMOL = int(os.getenv("TIMEOUT", DEFAULT_TIMEOUT_RUN_RETROMOL))


# Default timeout for computing optimal mappings
DEFAULT_TIMEOUT_OPTIMAL_MAPPINGS = 30  # seconds
TIMEOUT_OPTIMAL_MAPPINGS = int(os.getenv("TIMEOUT_OPTIMAL_MAPPINGS", DEFAULT_TIMEOUT_OPTIMAL_MAPPINGS))


# Default timeout for computing linear readouts
DEFAULT_TIMEOUT_LINEAR_READOUT = 30  # seconds
TIMEOUT_LINEAR_READOUT = int(os.getenv("TIMEOUT_LINEAR_READOUT", DEFAULT_TIMEOUT_LINEAR_READOUT))
