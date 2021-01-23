import shocktube1dcalc.helper


def test_shocktubecalc():
    """
    Compare analytic solutions provided by shocktubecalc.
    """
    assert shocktube1dcalc.helper.compare(0.01, 2.0)
