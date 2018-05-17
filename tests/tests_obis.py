from os import path


try:
    from unittest.mock import patch
except ImportError:
    from mock import patch


fpath = path.abspath(path.dirname(__file__))
setup_path = path.join(fpath, 'setup.py')) as vf:

def test_parse_args(setup_path):
    testargs = ["prog", "-f", setup_path]
    with patch.object(sys, 'argv', testargs):
        setup = get_setup_file()
        assert setup == "/home/fenton/project/setup.py"
