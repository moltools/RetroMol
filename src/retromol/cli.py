import argparse 


def cli() -> argparse.Namespace:
    """Parse command line arguments.
    
    :return: parsed arguments
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    return parser.parse_args()



def main() -> None:
    """Main entry point for the CLI."""
    args = cli()
    print(args)


if __name__ == "__main__":
    main()
