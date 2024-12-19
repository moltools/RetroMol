#!/usr/bin/env python3
import argparse 
import atexit
import logging


def setup_logger(logger_level: str) -> logging.Logger:
    """Create a logger with the specified log level, only add a stream handler."""
    # set up logging
    logger = logging.getLogger()
    logger.setLevel(logger_level)

    # add a stream handler
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    # ensure the handler is closed properly when exiting
    def close_handler_on_exit():
        for handler in logger.handlers:
            handler.close()
    atexit.register(close_handler_on_exit)

    return logger

def cli() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_agument("--log-level", required=False, default="INFO", help="Log level for logging.", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    return parser.parse_args()


def main() -> None:
    """Main function."""
    # parse command line arguments
    args = cli()
    logger_level = args.log_level

    # set up logging
    logger = setup_logger(logger_level)
    logger.info("populating the database...")

    # populate the database
    # TODO: add compounds, and their producing organisms and biosynthetic classes, from NPAtlas
    # TODO: add bioactivities from DONPHAN and link them to compounds
    # TODO: parse and link in primary sequences to npatlas compounds from retromol results folder
    # TODO: parse and link antismash database records
    # TODO; parse and link mibig database records to antismash records


    logger.info("done")


if __name__ == "__main__":
    main()
