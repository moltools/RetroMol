import re

from retromol.helpers import get_random_hex_color

def test_get_random_hex_color_format():
    """Test that the function returns a valid hex color format."""
    hex_color = get_random_hex_color()
    assert isinstance(hex_color, str), "the color should be a string"
    assert re.match(r"^#[0-9a-fA-F]{6}$", hex_color), "the color should be in valid hex format (e.g., #FFFFFF)"


def test_get_random_hex_color_randomness():
    """Test that the function returns different values over multiple calls."""
    colors = {get_random_hex_color() for _ in range(1000)}
    assert len(colors) > 1, "the function should return different colors across multiple calls"
