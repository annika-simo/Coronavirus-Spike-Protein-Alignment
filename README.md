# Coronavirus-Spike-Protein-Alignment
CS1021 Lab 7

## Introduction
Just a few months ago, most of us were unfamiliar with the technical difference between a pandemic and an endemic and had never seen an N95 mask nor uttered the phrase "social distancing". And yet within weeks of the first reports of a mystery viral [pneumonia](https://www.wsj.com/articles/health-officials-work-to-solve-chinas-mystery-virus-outbreak-11578308757?st=8xu2i33u9axxdta&reflink=desktopwebshare_permalink), we all became experts in the epidemiology of infectious diseases. As of [March 6th, 2022, almost 1 million](https://ourworldindata.org/covid-deaths#what-is-the-cumulative-number-of-confirmed-deaths) people have died in the United States alone from this terrible virus.

In April 2020, researchers published the [first complete genomic sequence](https://www.ncbi.nlm.nih.gov/nuccore/1798174254) of SARS-CoV-2 in the journal [Nature](https://pubmed.ncbi.nlm.nih.gov/32015508/). Because it was the first sample to be sequenced, it is now referred to simply as the reference sequence. Coronaviruses are so-called because of their shape: "[... under electron microscopic examination, each virion is surrounded by a “corona,” or halo. This is due to the presence of viral spike peplomers emanating from each proteinaceous envelope.](https://www.cdc.gov/coronavirus/mers/photos.html)"

[The spike proteins are the means by which the viral cells can infect other cells and cause the disease to multiply within the body](https://www.youtube.com/watch?v=Xuc9D4LVJdg). Proteins are composed of amino acids and the spike protein of SARS-CoV-2 has 1273 amino acids. As of this writing, the United States government hosts a [database](https://www.ncbi.nlm.nih.gov/datasets/coronavirus/proteins/) that contains complete genomic sequences of the spike proteins from over 250,000 samples of the SARS-CoV-2 virus. Long ago, bioinformaticians settled on a standard for representing and storing DNA sequences in digital form. Each of the amino acids that comprise a protein are composed of three nucleic acids (each, in turn, represented as a single letter). A protein is stored as a sequence of those letters. Multiple sequences are stored in the same file according to the FASTA format where each protein comprises two lines: a line containing identification information (it starts with a `>`) and a line containing the sequence of component amino acids.

In this lab we are going to build software that can compare the composition of the spike proteins from each of those sequences with the composition of the spike protein of the reference sequence. To do so, we are going to implement the [Needleman-Wunsch](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm) algorithm. This algorithm was designed in the late 1960s and early 1970s to perform an optimal alignment between two different genetic sequences. The algorithm uses dynamic programming, a topic that we will return to later in the course. The algorithm computes its results by calculating entries in two different specially created 2-dimensional vectors. "Specially created 2-dimensional vectors" is a mouthful. So, we'll shorten the *name of the type* of the two vectors to NWAlignment Direction and NWAlignment Difference. These special 2-dimensional vectors are examples of data structures (*remember the definition of the term data structure?*)

As we discussed in class, you can think of a 2-dimensional vector like a spreadsheet. Columns are vertical and rows are horizontal. In the algorithm, the columns of the "spreadsheet" are the nucleic acids that constitute the Spike protein from the reference sequence and the rows are the nucleic acids that constitute the Spike protein from a comparison sequence. 

The reference protein in the image above is made up of 7 nucleic acids: GCATGCG. The comparison protein is made up of 7 nucleic acids: GATTACA. During the calculation of the alignment, each of the "cells" will hold information about the relationship between an amino acid in the reference protein and an amino acid in the comparison protein. 

For example, the "cell" containing the X will hold some information about the relationship between the 3rd amino acid in the reference protein and the 3rd amino acid in the comparison protein that the algorithm will use to perform its computations. The algorithm addresses each cell in the NWAlignment data structures according to its row and column: the cell containing the X has the index 2,2; the cell above the cell containing the X has the index 1,2; the cell above/left the cell containing the X has the index 1,1; and the cell left of the cell containing the X has the index 2,1.

Accessing the neighbors around a specific cell (above, left and above/left) is normally easy. It is simply a question of modifying the cell's index. However, there are special cases that arise when calculating the surrounding cells on the "edges" of the data structure. Consider the cell labeled Y -- there is no cell above it in the NWAlignment data structure; nor is there a cell above/left. What about the cell labeled Z? There is a cell above, but there is no cell above/left or left. The algorithm specifies particular values in these situations.

The algorithm fills in the cells of the NWAlignment Difference data structure from the top-left cell to bottom-right cell. Once the algorithm calculates values for each of the cells, it proceeds to walk back -- from the bottom-right cell to top-left cell -- filling the cells of an NWAlignment Direction data structure. As it does the backward march, it uses the values in the cells to determine the optimal alignment between the reference and comparison sequence.

In this lab you are going to write *Pangolin*, an application that

1. reads a reference protein sequence from a FASTA-formatted file (named `reference.fasta`),
2. repeatedly reads comparison protein sequences from a second FASTA-formatted file (named `comparison.fasta`),
   1. compares the reference sequence with each comparison sequence,
   2. calculates statistics about the difference between the two sequences, and
   3. prints the statistics to the screen.

You will be implementing the overall execution of the program. While you will *not* be implementing the Needleman-Wunsch algorithm yourself, you will be implementing several helper functions to make the algorithm work properly:

1. `above`
2. `above_left`
3. `left`

To make your job of writing the entire program easier, there are some helpers functions that you will write:

1. `get_id_and_sequence`: Will get the id and list of amino acids from the "next" protein (both as strings) in a FASTA-formatted file; it will return false if there are no more id/sequence pairs to be read and will return true otherwise.
2. `open_file`: Will open a file handle to the file named in the function parameter; it will return false if it could not open the file and true otherwise.
3. `string_to_protein_sequence`: Will convert a string to a sequence of nucleic acids (a ProteinSequence data structure).

To make your life even easier, there are several functions that are already implemented for you:

1. `needleman_wunsch`: Takes a reference protein sequence and a comparison protein sequence as parameters and generates a NWAlignment Direction data structure.
2. `create_comparison`: Takes a reference and comparison ProteinSequence data structure and an NWAlignment Direction data structure as parameters and generates a Comparison data structure.
3. `print_comparison`: Takes a Comparison data structure and the reference sequence *string* as parameters and prints the statistics to the screen (*in the proper format!*).
4. `run_unit_tests`: A function that takes no parameters that will give you an initial impression of whether you implemented your helper functions correctly.
Using the functions that you write and the ones that you are given, pseudo code for your overall solution might look something like this:

```
1. Open the file containing reference protein sequence (reference.fasta) and check whether it opens
2. Open the file containing the comparison protein sequences (comparison.fasta) and check whether it opens
3. Read the reference protein id and sequence from the file containing the reference protein using get_id_and_sequence
4. Convert the reference sequence (which is a string) to a protein sequence (using string_to_protein_sequence)
5. While reading comparison id/sequence pairs from the file containing the comparison proteins (using get_id_and_sequence) succeeds,
- convert the current comparison sequence (which is a string) to a protein sequence (using string_to_protein_sequence)
- call the needleman_wunsch function to generate a NWAlignment Direction data structure
- call create_comparison with the appropriate parameters to generate a Comparison data structure
- call print_comparison using the result of the previous step and the id of the current comparison sequence as parameters 
```

As you are working through this lab, you can confirm that your results are real by comparing with the NIH's implementation online (Links to an external site.).

Good luck!

## Program Design Task
As my dad always said, “If I wait until the last minute, it only takes a minute.” Before you start writing code, please write the pseudocode or draw a flow chart for your implementation of the 6 functions required to complete the implementation of the *Pangolin* biostatistics application.

## Program Design Requirements
Your pseudocode or flow chart must describe the entirety of the process you plan to use to meet the requirements for each helper function. You may choose to write the flow chart or the pseudocode at any level of detail but remember that this is a tool for you! Your pseudocode or flow chart must include a separate, labeled description for each of the 6 helper functions given above. For places where you must handle corner cases, you must describe in your pseudocode how you will handle those.

## Programming Task
Your programming task is to implement the *Pangolin* application, including the 6 helper functions defined above. In this lab you will not prompt the user for input. Your application will output information to the user via the console using the given `print_comparison` function (see above).

Below are the formal specifications for each of the 6 helper functions that you must implement:

| Name	| Parameters |	Return Type |	Semantics |	Tips/Tricks |
| --- | --- | --- | --- | --- |
| `above`	| 1. An NWAlignment Difference named `difference_vector` <br> 2. A row index, named `row` <br> 3. A column index, named `col` | `int`	| If `(row-1)` is greater-than-or-equal-to zero, return the value in the cell at index `(row-1)`, `col`. Otherwise, return `-1 * (col + 1)`. `above` will never be given a `row`/`col` pair that is less than `0` or greater than the bounds of `difference_vector`. | |
| `left` | 1. An NWAlignment Difference named `difference_vector` <br> 2. A row index, named `row` <br> 3. A column index, named `col` | `int` |	If `(col-1)` is greater-than-or-equal-to zero, return the value in the cell at index row, `(col-1)`. Otherwise, return `-1 * (row + 1)`. left will never be given a `row`/`col` pair that is less than `0` or greater than the bounds of `difference_vector`.	| |
| `above_left`	| 1. An NWAlignment Difference named `difference_vector` <br> 2. A row index, named `row` <br> 3. A column index, named `col` | `int`	| If `(row-1)` and `(col-1)` are greater-than-or-equal-to zero, return the value in the cell at index `(row-1), (col-1)`.  Otherwise, if `col` equals `0`, return `-1 * row`. Otherwise, return `-1 * col`. `above_left` will never be given a `row`/`col` pair that is less than `0` or greater than the bounds of `difference_vector`.	| |
| `open_file`	| 1. A `std::string`, named `filename` <br> 2. A `std::ifstream`, named `file`, by reference | `bool`	| If a `std::ifstream` can be opened for the file named `filename`, then `file` is that `std::ifstream` and the function returns `true`. Otherwise, the function returns `false` and `file` is unchanged.	| `std::ifstream` has an [`open()`](https://en.cppreference.com/w/cpp/io/basic_ifstream/open) method. | 
| `get_id_and_sequence`	| 1. A `std::ifstream`, named `file`, by reference <br> 2. A `std::string`, named `id`, by reference <br> 3. A `std::string`, named `sequence`, by reference | `bool`	| If it is possible to read another id/protein sequence from `file`, `id` is that sequence id, `sequence` is that sequence and the function returns `true`. Otherwise, the function returns `false` and `id` and `sequence` are unchanged.	| Consider using the [`getline`](https://en.cppreference.com/w/cpp/string/basic_string/getline) function from the `std` namespace. The function will help you read an entire line of a file into a string. In order to determine whether you successfully read a string from the file using std::getline, try this:<br>`bool successful_read = !!std::getline(`*\<file\>*`, `*\<variable\>*`);` |
| `string_to_protein_sequence` | 1. A string named `protein_str` | `ProteinSequence`	| A `ProteinSequence` that represents that sequence contained in `protein_str`.	| `ProteinSequence` is an alias (thanks to `using`) for a `std::vector` of chars. `std::vector` has a [`push_back`](https://en.cppreference.com/w/cpp/container/vector/push_back) method for appending items. |
  
There is a `reference.fasta` and `comparison.fasta` file included in the skeleton code for this project. If your *Pangolin* application is working correctly, it will generate the following output:
  
```
>QKF30862.1:1-1273: 1
>QKJ89375.1:1-1273: 0
>QKY74774.1:1-1273: 1
>QNN30792.1:1-1273: 2
>QNN30804.1:1-1273: 0
>QNN93353.1:1-1273: 3
>QRR19290.1:1-1273: 2
>QTK30566.1:1-1273: 1
>QYC24546.1:1-1273: 4
>QZL64603.1:1-1273: 5
>QZL65751.1:1-1273: 5
>QZL65763.1:1-1273: 2
>QZL65775.1:1-1273: 2
>QZL65787.1:1-1273: 2
>QZL65799.1:1-1273: 2
```
  
Documenting functions is a very important part of a professional programmer's job. Every programming project has their own preferred format for writing comments that describe functions and their inputs/outputs. For the Pangolin project, the comments for each of the functions you define must conform with the following format:

```cpp
/*
 * <function name>
 *
 * <short description of what the function does>
 *
 * input: <short description of all input parameters>
 * output: <short description of all output parameters
 *          and the return value>
 */
```
  
Programming Requirements
If you are a Windows user, start with the skeleton .zip in this repository. This skeleton provides the starting point for a successful implementation of the *Pangolin* biostatistics application. If you do not use this skeleton code you will not be able to complete this lab.

In addition to the requirements/assumptions stated in the chart in the Programming Task section, there are three other requirements that your Pangolin app must meet:

1. You must read the reference protein sequence from the file `reference.fasta`.
2. You must read the comparison protein sequences from the file `comparison.fasta`.
3. You may use `run_unit_tests` for testing but your final submission must not contain an invocation of that function.
  
## Critical Thinking Task
The Pangolin application is an example of the power of the computer science combined with *big data*. The Internet defines big data as "a field that treats ways to analyze, systematically extract information from, or otherwise deal with data sets that are too large or complex to be dealt with by traditional data-processing application software." Computer scientists write algorithms to make sense of these data sets. The results of their analysis are used throughout society: [predictive policing](https://www.washingtonpost.com/opinions/big-data-may-be-reinforcing-racial-bias-in-the-criminal-justice-system/2017/02/10/d63de518-ee3a-11e6-9973-c5efb7ccfb0d_story.html), [availability of credit](https://www.bostonfed.org/publications/research-department-working-paper/2019/how-magic-a-bullet-is-machine-learning-for-credit-analysis.aspx), [self-driving cars](https://bdtechtalks.com/2020/07/29/self-driving-tesla-car-deep-learning/), and, of course, [advertising](https://www.ibm.com/watson-advertising/thought-leadership/benefits-of-machine-learning-in-advertising).

The utopian vision of universal societal improvement through the application of algorithms to data is clouded by the fact that engineers and developers a) choose the data they analyse and b) write the methods of analysis. These two facts mean that developers have a tremendous responsibility to make sure that their analyses are accurate and equitable. 

Many groups/individuals are studying ways to make sure that the application of algorithms to big data does not disproportionatly harm particular groups of people. Your task is to find and document such a group/individual and summarize their work.

## Critical Thinking Requirement
In 500 words of less, describe a group/individual who is studying ways to make sure that the application of algorithms to big data does not disproportionatly harm particular groups of people. All references to external resources must be properly documented and formatted. The choice of formatting for external references is up to you, but you may find it helpful to consult the Purdue OWL for [information](https://owl.purdue.edu/owl/research_and_citation/apa_style/apa_style_introduction.html). The Purdue OWL also has extensive information on ways to [avoid plagiarism](https://owl.purdue.edu/owl/avoiding_plagiarism/index.html).

## Deliverables
1. The pseudocode describing the algorithms you used to implement the 6 helper functions required to complete the Pangolin biostatistics application in PDF format (named `design.pdf`).
2. The C++ source code for the helper functions required for successful operation of the Pangolin application (named `pangolin.cpp`).
3. The response to the Critical Thinking prompts in PDF format (named `ml.pdf`)

## Related Learning Objectives
1. Writing boolean expressions using relational and logical operators
2. Using if-statements to implement selective program execution
3. Using algorithmic thinking to design programs
4. Writing syntactically correct for/while loops
5. Using/accessing/manipulating 2-dimensional vectors
6. Using methods on objects
