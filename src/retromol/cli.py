import argparse
import csv
import os
import logging
import threading
import time
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Callable, Any, Dict, Optional, Tuple

from tqdm import tqdm

from retromol.api import run_retromol


class Item:
    """Item for processing."""

    def __init__(self, name: str, smiles: str) -> None:
        """Initialize an Item instance.

        :param name: Item name.
        :type name: str
        :param smiles: Item SMILES.
        :type smiles: str
        """
        # replace forbidden characters in the name
        forbidden_chars = ["<", ">", ":", '"', "/", "\\", "|", "?", "*"]
        for char in forbidden_chars:
            name = name.replace(char, "_")

        self._name = name
        self._smiles = smiles

    @property
    def name(self) -> str:
        """Item name."""
        return self._name
    
    @property
    def smiles(self) -> str:
        """Item SMILES."""
        return self._smiles
    
    def __str__(self) -> str:
        """String representation of the item."""
        return self.name


class TimeoutError(Exception):
    """Custom exception raised when a function exceeds the specified timeout."""
    pass


def setup_logger(item_folder: str, logger_level: str = "DEBUG", log_file_name: str = "run.log") -> logging.Logger:
    """Sets up a logger that writes to a log file in the specified item folder.

    :param item_folder: Path to the item-specific folder where the log file will be stored.
    :type item_folder: str
    :param level: Logging level, defaults to "DEBUG".
    :type level: str, optional
    :param log_file_name: Name of the log file, defaults to "run.log".
    :type log_file_name: str, optional
    :return: A configured logger instance.
    :rtype: logging.Logger
    """
    log_filepath = os.path.join(item_folder, log_file_name)

    # configure the logger
    logger = logging.getLogger(item_folder)  # unique logger per folder
    logger.setLevel(logging.getLevelName(logger_level))

    # remove existing handlers to prevent duplicate logs in the same run
    if logger.hasHandlers(): logger.handlers.clear()

    handler = logging.FileHandler(log_filepath)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    return logger


def run_with_timeout(
    func: Callable[..., Any],
    args: Tuple[Any, ...] = (),
    kwargs: Optional[Dict[str, Any]] = None,
    timeout: int = 5,
    logger: Optional[logging.Logger] = None
) -> Any:
    """Runs a function with a timeout.

    :param func: The function to run.
    :type func: Callable[..., Any]
    :param args: Positional arguments for the function.
    :type args: tuple
    :param kwargs: Keyword arguments for the function.
    :type kwargs: dict, optional
    :param timeout: Maximum time allowed for the function to run, in seconds.
    :type timeout: int
    :param logger: Logger instance to use for logging.
    :type logger: logging.Logger, optional
    :return: The result of the function if completed within timeout.
    :rtype: Any
    :raises TimeoutError: If the function times out.
    :raises Exception: If the function raises any other exception.
    """
    kwargs = kwargs or {}
    result_container: Dict[str, Any] = {"result": None, "error": None}

    def target() -> None:
        try:
            result_container["result"] = func(*args, **kwargs)
        except Exception as e:
            result_container["error"] = e

    thread = threading.Thread(target=target)
    thread.start()
    thread.join(timeout)

    if thread.is_alive():
        if logger: logger.error(f"function '{func.__name__}' timed out after {timeout} second(s)")
        raise TimeoutError(f"function '{func.__name__}' timed out after {timeout} second(s)")

    if result_container["error"]:
        if logger: logger.error(f"error in function '{func.__name__}': {result_container['error']}")
        raise result_container["error"]

    if logger: logger.info(f"function '{func.__name__}' completed successfully")
    return result_container["result"]


def parse_item(item: Item, item_folder: str, logger: logging.Logger) -> None:
    """Simulates processing of an item by sleeping for the specified duration.

    :param item: The item to process.
    :type item: Item
    :param item_folder: Path to the item's specific folder.
    :type item_folder: str
    :param logger: Logger instance to log messages.
    :type logger: logging.Logger
    """
    # log start of task
    start_time = datetime.now()
    if logger: logger.info(f"starting task for item {item.name} at {start_time}")

    # run RetroMol on item
    run_retromol(item.name, item.smiles, logger)

    # log end of task
    end_time = datetime.now()
    if logger: logger.info(f"completed task for item {item} at {end_time}")


def process_item(item: Item, base_output_folder: str, timeout: int = 5) -> Optional[Tuple[str, str]]:
    """Processes an item with a timeout and creates an item-specific folder.

    :param item: The item to process.
    :type item: Item
    :param base_output_folder: Path to the base output folder where item-specific folders will be created.
    :type base_output_folder: str
    :param timeout: Maximum time allowed for the task to run, in seconds.
    :type timeout: int
    :return: A tuple with error type and message if an error occurs, or None if successful.
    :rtype: tuple[str, str] | None
    """
    # Create the item-specific folder
    item_folder = os.path.join(base_output_folder, f"results_{item}")
    os.makedirs(item_folder, exist_ok=True)

    # Set up a logger for this item
    logger = setup_logger(item_folder)

    try:
        # Run the task with a timeout, pass all required arguments
        run_with_timeout(parse_item, args=(item, item_folder, logger), timeout=timeout, logger=logger)
        logger.info(f"item {item} processed successfully")
        return None
    except TimeoutError:
        error_type = "TimeoutError"
        error_message = f"processing item {item} timed out"
        logger.warning(error_message)
        return error_type, error_message
    except Exception as e:
        # get the last traceback
        tb = e.__traceback__
        while tb.tb_next:  # get the last traceback
            tb = tb.tb_next
        logger.error(f"error with item {item}: {e}", exc_info=(type(e), e, tb))

        # return error type and message
        error_type = type(e).__name__
        error_message = str(e)
        logger.error(f"error with item {item}: {error_type} - {error_message}")
        return error_type, error_message


def cli() -> argparse.Namespace:
    """Parse command line arguments.
    
    :return: parsed arguments
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True, help="input file with with IDs and SMILES (csv)")
    parser.add_argument("-o", "--output", type=str, required=True, help="output folder")
    parser.add_argument("-m", "--max-cpus", type=int, default=1, help="maximum number of CPUs to use")
    parser.add_argument("-t", "--timeout", type=int, default=5, help="timeout for each item in seconds")
    parser.add_argument("-l", "--log-level", type=str, default="INFO", help="logging level", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    return parser.parse_args()


def main() -> None:
    """Main entry point for the CLI."""
    start_time = datetime.now()

    # parse CLI arguments
    args = cli()

    # config
    input_file = args.input
    base_output_dir = args.output
    max_cpus = min(os.cpu_count(), args.max_cpus)  # limit to available CPUs
    timeout = args.timeout
    logger_level = args.log_level
    log_file_name = "retromol.log"

    # configure logging
    logger = setup_logger(base_output_dir, logger_level, log_file_name)

    # create output directory if it doesn't exist
    os.makedirs(base_output_dir, exist_ok=True)

    # log input parameters
    logger.info("starting RetroMol with the following parameters:")
    logger.info(f"* input file: {input_file}")
    logger.info(f"* output directory: {base_output_dir}")
    logger.info(f"* maximum CPUs: {max_cpus}")
    logger.info(f"* timeout per item: {timeout} second(s)")

    # read input file
    items = []
    try:
        # parse items from input file
        with open(input_file, "r") as fo:
            fo.readline()  # skip header
            for line in fo:
                name, smiles = line.strip().split(",")
                items.append(Item(name, smiles))

        # check if all names are unique
        names = [item.name for item in items]
        if len(names) != len(set(names)):
            raise ValueError("item names must be unique")
    except Exception as e:
        logger.error(f"error reading input file: {e}")
        print(f"error reading input file: {e}")
        exit(1)
    else:
        logger.info(f"number of items to process: {len(items)}")

    # store results
    results = []

    # use ThreadPoolExecutor to process items in parallel
    with ThreadPoolExecutor(max_workers=max_cpus) as executor:
        future_to_item = {executor.submit(process_item, item, base_output_dir, timeout): item for item in items}

        # progress bar setup
        with tqdm(total=len(future_to_item), desc="processing items", unit="item") as pbar:
            for future in as_completed(future_to_item):
                item = future_to_item[future]
                try:
                    result = future.result()
                    if result:
                        # if there's an error, store error type and message
                        error_type, error_message = result
                        results.append({
                            "item": str(item),
                            "smiles": item.smiles,
                            "status": "failed",
                            "error_type": error_type,
                            "error_message": error_message
                        })
                    else:
                        # if successful, note the success
                        results.append({
                            "item": str(item),
                            "smiles": item.smiles,
                            "status": "success",
                            "error_type": None,
                            "error_message": None
                        })
                except Exception as e:
                    # catch unexpected exceptions from `process_item` itself
                    results.append({
                        "item": str(item),
                        "smiles": item.smiles,
                        "status": "failed",
                        "error_type": "UnexpectedError",
                        "error_message": f"{type(e).__name__}: {str(e)}"
                    })
                finally:
                    # update the progress bar
                    pbar.update(1)

    # sort results by item name
    results = sorted(results, key=lambda x: x["item"])

    # write results to a CSV file
    csv_file_path = os.path.join(base_output_dir, "processing_results.csv")
    with open(csv_file_path, mode="w", newline="") as csvfile:
        fieldnames = ["item", "smiles", "status", "error_type", "error_message"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    end_time = datetime.now()
    msg = f"processing complete\n" \
            f"total runtime: {end_time - start_time}\n" \
            f"results written to: {base_output_dir}\n" \
            f"overview results written to: {csv_file_path}\n" \
            f"log file written to: {os.path.join(base_output_dir, log_file_name)}"
    logger.info(msg)
    print(msg)
    exit(0)


if __name__ == "__main__":
    main()
