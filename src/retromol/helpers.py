import random 


def get_random_hex_color() -> str:
    return f"#{random.randint(0, 0xFFFFFF):06x}"
