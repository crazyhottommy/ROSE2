"""Tests for rose2.utils module."""

import pytest
from rose2 import utils


class TestLocus:
    """Tests for the Locus class."""

    def test_locus_creation(self):
        """Test basic locus creation."""
        locus = utils.Locus("chr1", 1000, 2000, "+", "test_locus")
        assert locus.chr() == "chr1"
        assert locus.start() == 1000
        assert locus.end() == 2000
        assert locus.sense() == "+"
        assert locus.ID() == "test_locus"

    def test_locus_length(self):
        """Test locus length calculation."""
        locus = utils.Locus("chr1", 1000, 2000, "+")
        assert locus.len() == 1001  # Inclusive of endpoints

    def test_locus_overlaps(self):
        """Test locus overlap detection."""
        locus1 = utils.Locus("chr1", 1000, 2000, "+")
        locus2 = utils.Locus("chr1", 1500, 2500, "+")
        locus3 = utils.Locus("chr1", 3000, 4000, "+")

        assert locus1.overlaps(locus2)
        assert not locus1.overlaps(locus3)

    def test_locus_contains(self):
        """Test locus containment."""
        outer = utils.Locus("chr1", 1000, 3000, "+")
        inner = utils.Locus("chr1", 1500, 2500, "+")
        separate = utils.Locus("chr1", 4000, 5000, "+")

        assert outer.contains(inner)
        assert not outer.contains(separate)
        assert not inner.contains(outer)

    def test_locus_string_representation(self):
        """Test locus string representation."""
        locus = utils.Locus("chr1", 1000, 2000, "+", "test")
        assert str(locus) == "chr1(+):1000-2000"


class TestLocusCollection:
    """Tests for the LocusCollection class."""

    def test_collection_creation(self):
        """Test basic collection creation."""
        loci = [
            utils.Locus("chr1", 1000, 2000, "+", "locus1"),
            utils.Locus("chr1", 3000, 4000, "+", "locus2"),
        ]
        collection = utils.LocusCollection(loci, 500)
        assert len(collection) == 2

    def test_collection_overlap(self):
        """Test getting overlapping loci."""
        loci = [
            utils.Locus("chr1", 1000, 2000, "+", "locus1"),
            utils.Locus("chr1", 1500, 2500, "+", "locus2"),
            utils.Locus("chr1", 5000, 6000, "+", "locus3"),
        ]
        collection = utils.LocusCollection(loci, 500)

        query = utils.Locus("chr1", 1700, 1800, "+")
        overlapping = collection.get_overlap(query, "sense")

        assert len(overlapping) == 2  # Should overlap with locus1 and locus2

    def test_collection_stitch(self):
        """Test stitching nearby loci."""
        loci = [
            utils.Locus("chr1", 1000, 2000, "+", "locus1"),
            utils.Locus("chr1", 2100, 3000, "+", "locus2"),  # 100bp gap
        ]
        collection = utils.LocusCollection(loci, 500)

        # Stitch with 200bp window (should merge)
        stitched = collection.stitch_collection(stitch_window=200)
        assert len(stitched) == 1  # Should be stitched into one

        # Stitch with 50bp window (should not merge)
        collection2 = utils.LocusCollection(loci, 500)
        stitched2 = collection2.stitch_collection(stitch_window=50)
        assert len(stitched2) == 2  # Should remain separate


class TestUtilityFunctions:
    """Tests for utility functions."""

    def test_uniquify(self):
        """Test uniquify function."""
        input_list = [1, 2, 2, 3, 3, 3, 4]
        result = utils.uniquify(input_list)
        assert result == [1, 2, 3, 4]

    def test_uniquify_preserves_order(self):
        """Test that uniquify preserves order."""
        input_list = ["c", "a", "b", "a", "c"]
        result = utils.uniquify(input_list)
        assert result == ["c", "a", "b"]

    def test_order_ascending(self):
        """Test order function in ascending order."""
        input_list = [30, 10, 20]
        result = utils.order(input_list, decreasing=False)
        assert result == [1, 2, 0]  # Indices that would sort the list

    def test_order_descending(self):
        """Test order function in descending order."""
        input_list = [30, 10, 20]
        result = utils.order(input_list, decreasing=True)
        assert result == [0, 2, 1]  # Indices for descending sort

    def test_order_with_none(self):
        """Test order function with None values."""
        input_list = [3, None, 1, None, 2]
        result = utils.order(input_list, none_is_last=True, decreasing=False)
        # Should order: 1, 2, 3, None, None (indices: 2, 4, 0, 1, 3)
        assert result == [2, 4, 0, 1, 3]


class TestFileIO:
    """Tests for file I/O functions."""

    def test_parse_and_unparse_table(self, tmp_path):
        """Test parsing and unparsing tables."""
        # Create test data
        test_table = [
            ["col1", "col2", "col3"],
            ["a", "b", "c"],
            ["d", "e", "f"],
        ]

        # Write table
        output_file = tmp_path / "test_table.txt"
        utils.unparse_table(test_table, str(output_file), "\t")

        # Read it back
        result = utils.parse_table(str(output_file), "\t")

        assert len(result) == 3
        assert result[0] == ["col1", "col2", "col3"]
        assert result[1] == ["a", "b", "c"]
        assert result[2] == ["d", "e", "f"]

    def test_bed_to_gff_conversion(self):
        """Test BED to GFF conversion."""
        bed_data = [
            ["chr1", "1000", "2000", "peak1", "100", "+"],
            ["chr2", "3000", "4000", "peak2", "200", "-"],
        ]

        gff = utils.bed_to_gff(bed_data)

        assert len(gff) == 2
        assert gff[0][0] == "chr1"  # chromosome
        assert gff[0][1] == "peak1"  # name
        assert gff[0][3] == "1000"  # start
        assert gff[0][4] == "2000"  # end
        assert gff[0][6] == "."  # strand placeholder


def test_convert_bitwise_flag():
    """Test SAM flag conversion."""
    assert utils.convert_bitwise_flag(0) == "+"
    assert utils.convert_bitwise_flag(16) == "-"
    assert utils.convert_bitwise_flag("0") == "+"
    assert utils.convert_bitwise_flag("16") == "-"
