#include "unit_test.h"
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/types.h>
#include <tuple>
#include <vector>

/* References:
 * https://bioboot.github.io/bimm143_W20/class-material/nw/
 * https://www.ncbi.nlm.nih.gov/nuccore/NC_045512
 * https://www.ncbi.nlm.nih.gov/datasets/coronavirus/proteins/
 */

using direction = enum : uint8_t { diagonal = 0, horizontal = 1, vertical = 2 };
using NWAlignmentDifference = std::vector<std::vector<int>>;
using NWAlignmentDirection = std::vector<std::vector<direction>>;
using ProteinSequence = std::vector<char>;
using Comparison = std::vector<std::tuple<int, char, char>>;

/*
 * get_id_and_sequence
 *
 * Get an id and sequence (on consecutive lines) of _file_. Store
 * the results in _id_ and _sequence_ and return `true` on success
 * and `false` on failure.
 *
 * input: The file from which to read the id and sequence, as a
 * reference.
 * output: The id read from the file, as a `std::string`.
 * output: The sequence read from the file, as a `std::string`.
 * return: `true` if it was possible to read the id/sequence from
 * the file; false, otherwise.
 */
bool get_id_and_sequence(std::ifstream &file, std::string &id,
                         std::string &sequence) {
  bool valid_sequence = !!std::getline(file, id);
  if (!valid_sequence) {
    return false;
  }
  return !!std::getline(file, sequence);
}

/*
 * open_file
 *
 * Attempt to open a file with a given name.
 *
 * input: The name of the file to open, as a `const` reference to a
 * `std::string`.
 * output: An open (if possible) `std::ifstream` associated with the
 * filename named _filename_.
 * return: `true` if it was possible to open a file named _filename_; false,
 * otherwise.
 */
bool open_file(std::string filename, std::ifstream &file) {
  file.open(filename);
  return file.is_open();
}

/*
 * string_to_protein_sequence
 *
 * Convert a `std::string` whose characters are letters representing
 * different proteins into a `ProteinSequence`.
 *
 * input: The string of characters to convert to a protein sequence, as
 * `const` reference to a `std::string`.
 * return: The converted sequence as a `ProteinSequence`.
 */
ProteinSequence string_to_protein_sequence(const std::string &protein_str) {
  ProteinSequence v(protein_str.length());
  for (int i = 0; i < protein_str.length(); i++) {
    v[i] = protein_str[i];
  }
  return v;
}

/*
 * above_left
 *
 * For a 2D vector of integers return the contents of the cell that is
 * above/left the cell given by row,col, if those two indexes are in bounds.
 * Otherwise, if col is out of bounds, return the literal -1 * row; otherwise,
 * return -1 * col.
 *
 * input: a 2D vector from which to read the values.
 * input: the row relative to which to read.
 * input: the column relative to which to read.
 * return: the value of the cell above/left _row_, _col_ in _difference_vector_,
 * if the row/colum above/left are valid. Otherwise, if col is out of bounds,
 * then return the literal -1 * row; otherwise, return -1 * col.
 */
int above_left(const NWAlignmentDifference &difference_vector, int row,
               int col) {
  int diagonal = 0;
  if (row != 0 && col != 0) {
    diagonal = difference_vector[row - 1][col - 1];
  } else {
    if (col == 0) {
      diagonal = -1 * row;
    } else {
      diagonal = -1 * col;
    }
  }
  return diagonal;
}

/*
 * above
 *
 * For a 2D vector of integers return the contents of the cell that is
 * above the cell given by row,col, if those two indexes are in bounds.
 * Otherwise, if the row above _row_ is out of bounds, then return -1 * (col +
 * 1).
 *
 * input: a 2D vector from which to read the values.
 * input: the row relative to which to read.
 * input: the column relative to which to read.
 * return: the the contents of the cell that is
 * above the cell given by row,col, if those two indexes are in bounds.
 * Otherwise, if the row above _row_ is out of bounds, then return -1 * (col +
 * 1).
 */
int above(const NWAlignmentDifference &difference_vector, int row, int col) {
  if (row == 0) {
    return -1 * (col + 1);
  }
  return difference_vector[row - 1][col];
}

/*
 * left
 *
 * For a 2D vector of integers return the contents of the cell that is
 * left of the cell given by row, col, if that access is in bounds.
 * Otherwise, if the column left of _col_ is out of bounds, then return -1 *
 * (row + 1).
 *
 * input: a 2D vector from which to read the values.
 * input: the row relative to which to read.
 * input: the column relative to which to read.
 * return: the the the contents of the cell that is
 * left of the cell given by row, col, if that access is in bounds.
 * Otherwise, if the column left of _col_ is out of bounds, then return -1 *
 * (row + 1).
 */
int left(const NWAlignmentDifference &difference_vector, int row, int col) {
  if (col == 0) {
    return -1 * (row + 1);
  }
  return difference_vector[row][col - 1];
}

void print_comparison(const Comparison &comparison,
                      const std::string &comparison_id, bool verbose = false) {

  std::cout << comparison_id << ": " << comparison.size();

  if (verbose) {
    for (auto i : comparison) {
      auto index = std::get<0>(i);
      auto first = std::get<1>(i);
      auto second = std::get<2>(i);
      std::cout << "(@" << index << ": " << first << ", " << second << "), ";
    }
  }
  std::cout << "\n";
}

template <typename T>
std::vector<T> reverse_vector(const std::vector<T> &to_reverse) {
  std::vector<T> forward;
  for (auto i : to_reverse) {
    forward.push_back(i);
  }
  return forward;
}

Comparison create_comparison(const ProteinSequence &top,
                             const ProteinSequence &bottom,
                             const NWAlignmentDirection &dirs) {
  int row = bottom.size() - 1;
  int col = top.size() - 1;
  Comparison full_result, diff_result;
  while (row >= 0 && col >= 0) {
    char first{0};
    char second{0};
    if (dirs[row][col] == diagonal) {
      first = bottom[row];
      second = top[col];
      row--;
      col--;
    } else if (dirs[row][col] == horizontal) {
      first = '-';
      second = top[col];
      col--;
    } else {
      first = bottom[row];
      second = '-';
      row--;
    }
    full_result.push_back(std::make_tuple(0, first, second));
  }
  full_result = reverse_vector(full_result);

  auto index = 0;
  for (auto i : full_result) {
    auto first = std::get<1>(i);
    auto second = std::get<2>(i);
    if (first != second) {
      diff_result.push_back(std::make_tuple(index, first, second));
    }
    index++;
  }
  return diff_result;
}

NWAlignmentDirection new_direction_vector(int rows, int cols) {
  std::vector<std::vector<direction>> new_direction_vector(
      rows, std::vector<direction>(cols));
  return new_direction_vector;
}

NWAlignmentDifference new_difference_vector(int rows, int cols) {
  NWAlignmentDifference new_diff_vector(rows, std::vector<int>(cols));
  return new_diff_vector;
}

int s(const char a, const char b) {
  if (a == 'X' || b == 'X') {
    return 0;
  }
  if (a == b) {
    return 1;
  }
  return -1;
}

NWAlignmentDirection needleman_wunsch(const ProteinSequence &top,
                                      const ProteinSequence &bottom) {
  auto bottom_length = bottom.size();
  auto top_length = top.size();
  auto directions = new_direction_vector(bottom_length, top_length);
  auto differences = new_difference_vector(bottom_length, top_length);
  for (int r = 0; r < bottom_length; r++) {
    for (int c = 0; c < top_length; c++) {
      int t, l, d, match, insert, del;
      t = above(differences, r, c);
      l = left(differences, r, c);
      d = above_left(differences, r, c);

      match = d + s(top[c], bottom[r]);
      del = t - 1;
      insert = l - 1;

      if (match > del && match > insert) {
        differences[r][c] = match;
        directions[r][c] = diagonal;
      } else if (del > insert) {
        differences[r][c] = del;
        directions[r][c] = vertical;
      } else {
        differences[r][c] = insert;
        directions[r][c] = horizontal;
      }
    }
  }
  return directions;
}

void run_unit_tests() {
  std::ifstream test_input_file;

  check_result(true, open_file("testing.fasta", test_input_file));
  check_result(true, test_input_file.is_open());

  std::string id, sequence;

  check_result(true, get_id_and_sequence(test_input_file, id, sequence));
  check_result(true, id == std::string{"id1"});
  check_result(true, sequence == std::string{"sequence1"});
  check_result(true, get_id_and_sequence(test_input_file, id, sequence));
  check_result(true, id == std::string{"id2"});
  check_result(true, sequence == std::string{"sequence2"});

  test_input_file.close();

  auto protein_sequence = string_to_protein_sequence(sequence);
  check_result(true,
               std::equal(sequence.begin(), sequence.end(),
                          protein_sequence.begin(), protein_sequence.end()));

  NWAlignmentDifference difference = new_difference_vector(5, 5);
  int counter{0};
  for (auto &row : difference) {
    for (auto &cell : row) {
      cell = counter++;
    }
  }

  check_result(0, above(difference, 1, 0));
  check_result(-2, above(difference, 0, 1));
  check_result(13, above(difference, 3, 3));
  check_result(0, left(difference, 0, 1));
  check_result(17, left(difference, 3, 3));
  check_result(-5, left(difference, 4, 0));
  check_result(12, above_left(difference, 3, 3));
  check_result(-4, above_left(difference, 0, 4));
  check_result(-3, above_left(difference, 3, 0));
}

int main() {

  // run_unit_tests();

  // Declare all the necessary variables up front -- I like doing it this way,
  // but it's a preference.
  std::string reference_file_name{"reference.fasta"};
  std::string comparison_file_name{"comparison.fasta"};
  std::ifstream reference_file, comparison_file;
  std::string reference_sequence{""}, reference_id{""}, comparison_id{""},
      comparison_sequence{""};

  // If it is not possible to open the reference file, print an error and exit.
  if (!open_file(reference_file_name, reference_file)) {
    std::cout << "Oops: Could not open the reference file.\n";
    return 1;
  }
  // If it is not possible to open the comparison file, print an error and exit.
  if (!open_file(comparison_file_name, comparison_file)) {
    std::cout << "Oops: Could not open the comparison file.\n";
    reference_file.close();
    return 1;
  }

  // Attempt to read the reference id and sequence from the reference file. If
  // that fails, print an error and exit.
  if (!get_id_and_sequence(reference_file, reference_id, reference_sequence)) {
    std::cout << "Oops: Could not read the reference sequence from the "
                 "reference file.\n";
    // Don't forget to close the files that we opened!
    reference_file.close();
    comparison_file.close();
    return 1;
  }

  // Convert that reference sequence of characters (which are just characters
  // in biological encoding) in to a proper ProteinSequence.
  auto reference_sequence_vector =
      string_to_protein_sequence(reference_sequence);

  // For as long as it is possible successfully read an id/sequence pair
  // from the comparison file ...
  while (get_id_and_sequence(comparison_file, comparison_id,
                             comparison_sequence)) {
    // Convert that comparison sequence of characters (which are just characters
    // in biological encoding) in to a proper ProteinSequence.
    auto comparison_sequence_vector =
        string_to_protein_sequence(comparison_sequence);

    // Use the needleman_wunsch function to create a direction vector that will
    // help determine the best alignment between the reference and comparison
    // sequence.
    auto dirs =
        needleman_wunsch(reference_sequence_vector, comparison_sequence_vector);

    // With the best alignment in hand, create a comparison between the
    // comparison and the reference sequence.
    auto comparison = create_comparison(reference_sequence_vector,
                                        comparison_sequence_vector, dirs);
    // Now, print out the comparison in the proper format!
    print_comparison(comparison, comparison_id);
  }

  // Don't forget to close the files that we openend!
  reference_file.close();
  comparison_file.close();
}
