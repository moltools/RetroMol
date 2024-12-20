import atexit
import logging
import os
import typing as ty

import psycopg2


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


def connect_to_database() -> ty.Tuple[psycopg2._ext.connection, psycopg2._ext.cursor]:
    """Connect to the database and return the connection and cursor."""
    # connect to the database
    conn = psycopg2.connect(
        dbname=os.getenv("DB_NAME", "retromol"),
        user=os.getenv("DB_USER", "davidmeijer"),
        password=os.getenv("DB_PASSWORD", "postgres"),
        host=os.getenv("DB_HOST", "localhost"),
        port=os.getenv("DB_PORT", "5432")
    )

    # create cursor and ensure the connection is closed properly when exiting
    cur = conn.cursor()
    def close_connection_on_exit():
        conn.close()
    atexit.register(close_connection_on_exit)

    return conn, cur
