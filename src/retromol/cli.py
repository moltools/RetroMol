import argparse
import atexit
import csv
import os
import logging
import threading
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing
multiprocessing.set_start_method('fork')
from typing import Callable, Any, Dict, Optional, Tuple

from tqdm import tqdm

from retromol.api import run_retromol


class Item:
    """Item for processing."""

    def __init__(
        self, 
        name: str, 
        smiles: str,
        preprocessing_rules: Optional[str] = None,
        sequencing_rules: Optional[str] = None,
        motifs: Optional[str] = None,
        only_user_preprocessing_rules: bool = False,
        only_user_sequencing_rules: bool = False,
        only_user_motifs: bool = False
    ) -> None:
        """Initialize an Item instance."""
        # replace forbidden characters in the name
        forbidden_chars = ["<", ">", ":", '"', "/", "\\", "|", "?", "*"]
        for char in forbidden_chars:
            name = name.replace(char, "_")

        self._name = name
        self._smiles = smiles
        self._preprocessing_rules = preprocessing_rules
        self._sequencing_rules = sequencing_rules
        self._motifs = motifs
        self._only_user_preprocessing_rules = only_user_preprocessing_rules
        self._only_user_sequencing_rules = only_user_sequencing_rules
        self._only_user_motifs = only_user_motifs

    @property
    def name(self) -> str:
        """Item name."""
        return self._name
    
    @property
    def smiles(self) -> str:
        """Item SMILES."""
        return self._smiles
    
    @property
    def preprocessing_rules(self) -> Optional[str]:
        """Path to the preprocessing rules file."""
        return self._preprocessing_rules
    
    @property
    def sequencing_rules(self) -> Optional[str]:
        """Path to the sequencing rules file."""
        return self._sequencing_rules
    
    @property
    def motifs(self) -> Optional[str]:
        """Path to the motifs file."""
        return self._motifs
    
    @property
    def only_user_preprocessing_rules(self) -> bool:
        """Whether to use only user-defined preprocessing rules."""
        return self._only_user_preprocessing_rules
    
    @property
    def only_user_sequencing_rules(self) -> bool:
        """Whether to use only user-defined sequencing rules."""
        return self._only_user_sequencing_rules
    
    @property
    def only_user_motifs(self) -> bool: 
        """Whether to use only user-defined motifs."""
        return self._only_user_motifs
    
    def __str__(self) -> str:
        """String representation of the item."""
        return self.name


class TimeoutError(Exception):
    """Custom exception raised when a function exceeds the specified timeout."""
    pass


def setup_logger(item_folder: str, logger_level: str = "DEBUG", log_file_name: str = "run.log", verbose: bool = False) -> logging.Logger:
    """Sets up a logger that writes to a log file in the specified item folder.

    :param item_folder: Path to the item-specific folder where the log file will be stored.
    :type item_folder: str
    :param level: Logging level, defaults to "DEBUG".
    :type level: str, optional
    :param log_file_name: Name of the log file, defaults to "run.log".
    :type log_file_name: str, optional
    :param verbose: Whether to enable verbose logging, defaults to False.
    :type verbose: bool, optional
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

    # if verbose logging is enabled, also log to stdout
    if verbose:
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

    # ensure the handler is closed properly when exiting
    def close_handler_on_exit():
        for handler in logger.handlers:
            handler.close()

    atexit.register(close_handler_on_exit)

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

    thread = threading.Thread(target=target, daemon=True)  # daemon=True makes sure the thread is killed when the main thread exits 
    thread.start()
    thread.join(timeout=float(timeout))

    if thread.is_alive():
        if logger: logger.error(f"function '{func.__name__}' timed out after {timeout} second(s)")
        raise TimeoutError(f"function '{func.__name__}' timed out after {timeout} second(s)")

    if result_container["error"]:
        if logger: logger.error(f"error in function '{func.__name__}': {result_container['error']}")
        raise result_container["error"]

    if logger: logger.info(f"function '{func.__name__}' completed successfully")
    return result_container["result"]


def parse_item(item: Item, item_folder: str, logger: logging.Logger) -> float:
    """Simulates processing of an item by sleeping for the specified duration.

    :param item: The item to process.
    :type item: Item
    :param item_folder: Path to the item's specific folder.
    :type item_folder: str
    :param logger: Logger instance to log messages.
    :type logger: logging.Logger
    :return: The coverage score for the item.
    :rtype: float
    """
    # log start of task
    start_time = datetime.now()
    if logger: logger.info(f"starting task for item {item.name} at {start_time}")

    # run RetroMol on item
    coverage_score = run_retromol(
        item.name, 
        item.smiles, 
        item.preprocessing_rules,
        item.sequencing_rules,
        item.motifs,
        not item.only_user_preprocessing_rules,
        not item.only_user_sequencing_rules,
        not item.only_user_motifs,
        item_folder, 
        logger
    )
    if logger: logger.info(f"coverage score for item {item.name}: {coverage_score}")

    # log end of task
    end_time = datetime.now()
    if logger: logger.info(f"completed task for item {item} at {end_time}")

    return coverage_score


def process_item(item: Item, base_output_folder: str, timeout: int = 5, logger_level: str = "DEBUG", log_file_name: str = "run.log", verbose: bool = False) -> Tuple[Optional[str], Optional[str], Optional[float]]:
    """Processes an item with a timeout and creates an item-specific folder.

    :param item: The item to process.
    :type item: Item
    :param base_output_folder: Path to the base output folder where item-specific folders will be created.
    :type base_output_folder: str
    :param timeout: Maximum time allowed for the task to run, in seconds.
    :type timeout: int
    :param logger_level: Logging level for the logger.
    :type logger_level: str
    :param log_file_name: Name of the log file.
    :type log_file_name: str
    :param verbose: Whether to enable verbose logging.
    :type verbose: bool
    :return: A tuple containing the error type, error message, and coverage score.
    :rtype: Tuple[Optional[str], Optional[str], Optional[float
    """
    # Create the item-specific folder
    item_folder = os.path.join(base_output_folder, f"results_{item}")
    os.makedirs(item_folder, exist_ok=True)

    # Set up a logger for this item
    logger = setup_logger(item_folder, logger_level, log_file_name, verbose)

    try:
        # Run the task with a timeout, pass all required arguments
        coverage_score = run_with_timeout(
            parse_item, 
            args=(
                item, 
                item_folder, 
                logger
            ), 
            timeout=timeout, 
            logger=logger
        )
        logger.info(f"item {item} processed successfully")
        return None, None, coverage_score
    except TimeoutError:
        error_type = "TimeoutError"
        error_message = f"processing item {item} timed out"
        logger.warning(error_message)
        return error_type, error_message, None
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
        return error_type, error_message, None


def cli() -> argparse.Namespace:
    """Parse command line arguments.
    
    :return: parsed arguments
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True, help="input file with with IDs and SMILES (csv)")
    parser.add_argument("-o", "--output", type=str, required=True, help="output folder")
    parser.add_argument("-m", "--max-cpus", type=int, default=1, help="maximum number of CPUs to use")
    parser.add_argument("-t", "--timeout", type=int, default=30, help="timeout for each item in seconds")
    parser.add_argument("-l", "--log-level", type=str, default="INFO", help="logging level", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    parser.add_argument("-v", "--verbose", action="store_true", help="enable verbose logging by also logging to stdout")

    parser.add_argument("-pr", "--preprocessing-rules", type=str, required=False, help="path to json file with preprocessing rules")
    parser.add_argument("-sr", "--sequencing-rules", type=str, required=False, help="path to json file with sequencing rules")
    parser.add_argument("-mo", "--motifs", type=str, required=False, help="path to json file with motifs")
    parser.add_argument("-oupr", "--only-user-preprocessing-rules", action="store_true", help="use only user-defined preprocessing rules")
    parser.add_argument("-ousr", "--only-user-sequencing-rules", action="store_true", help="use only user-defined sequencing rules")
    parser.add_argument("-oumo", "--only-user-motifs", action="store_true", help="use only user-defined motifs")

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
    verbose = args.verbose
    preprocessing_rules = args.preprocessing_rules
    sequencing_rules = args.sequencing_rules
    motifs = args.motifs
    only_user_preprocessing_rules = args.only_user_preprocessing_rules
    only_user_sequencing_rules = args.only_user_sequencing_rules
    only_user_motifs = args.only_user_motifs

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
    logger.info(f"* logging level: {logger_level}")
    logger.info(f"* verbose logging: {'enabled' if verbose else 'disabled'}")
    logger.info(f"* preprocessing rules: {preprocessing_rules}")
    logger.info(f"* sequencing rules: {sequencing_rules}")
    logger.info(f"* motifs: {motifs}")
    logger.info(f"* only user-defined preprocessing rules: {only_user_preprocessing_rules}")
    logger.info(f"* only user-defined sequencing rules: {only_user_sequencing_rules}")
    logger.info(f"* only user-defined motifs: {only_user_motifs}")

    # read input file
    items = []
    try:
        # parse items from input file
        with open(input_file, "r") as fo:
            fo.readline()  # skip header
            for line in fo:
                name, smiles, *_ = line.strip().split(",")
                items.append(
                    Item(
                        name, 
                        smiles,
                        preprocessing_rules,
                        sequencing_rules,
                        motifs,
                        only_user_preprocessing_rules,
                        only_user_sequencing_rules,
                        only_user_motifs
                    )
                )

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
        future_to_item = {executor.submit(process_item, item, base_output_dir, timeout, logger_level, verbose=verbose): item for item in items}

        # progress bar setup
        with tqdm(total=len(future_to_item), desc="processing items", unit="item") as pbar:
            for future in as_completed(future_to_item):
                item = future_to_item[future]
                try:
                    result = future.result()
                    error_type, error_message, coverage_score = result
                    results.append({
                        "item": str(item),
                        "status": "succeeded" if coverage_score is not None else "failed",
                        "coverage_score": coverage_score,
                        "error_type": error_type,
                        "error_message": error_message,
                    })
                except Exception as e:
                    # catch unexpected exceptions from `process_item` itself
                    results.append({
                        "item": str(item),
                        "status": "failed",
                        "coverage_score": None,
                        "error_type": "UnexpectedError",
                        "error_message": f"{type(e).__name__}: {str(e)}",
                    })
                finally:
                    # update the progress bar
                    pbar.update(1)

    # clean thread pool
    executor.shutdown(wait=True)

    # sort results by item name
    results = sorted(results, key=lambda x: x["item"])

    # write results to a CSV file
    csv_file_path = os.path.join(base_output_dir, "processing_results.csv")
    with open(csv_file_path, mode="w", newline="") as csvfile:
        fieldnames = ["item", "status", "coverage_score", "error_type", "error_message"]
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
