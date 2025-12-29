"""Tests for CLI functionality."""

import pytest
from rose2 import cli


def test_main_help():
    """Test that main help runs without error."""
    with pytest.raises(SystemExit) as exc_info:
        cli.main(["--help"])
    assert exc_info.value.code == 0


def test_main_version():
    """Test version display."""
    with pytest.raises(SystemExit) as exc_info:
        cli.main(["--version"])
    # argparse exits with 0 for --version
    assert exc_info.value.code == 0


def test_main_no_args():
    """Test main with no arguments."""
    result = cli.main([])
    assert result == 0  # Should print help and return 0


def test_subcommand_help():
    """Test that subcommand help works."""
    # Test main subcommand help
    with pytest.raises(SystemExit) as exc_info:
        cli.main(["main", "--help"])
    assert exc_info.value.code == 0

    # Test bamToGFF subcommand help
    with pytest.raises(SystemExit) as exc_info:
        cli.main(["bamToGFF", "--help"])
    assert exc_info.value.code == 0

    # Test geneMapper subcommand help
    with pytest.raises(SystemExit) as exc_info:
        cli.main(["geneMapper", "--help"])
    assert exc_info.value.code == 0
