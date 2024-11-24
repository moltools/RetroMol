from time import sleep
from random import randint


def run_retromol(smiles: str) -> None:
    """Parse a SMILES string."""
    # sleep random number of seconds between 1 and 7 to simulate different run times
    sleep(randint(1, 7))
