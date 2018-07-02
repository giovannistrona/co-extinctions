"""
READ IN TEXT, LINES or CSV DATA
===============================
@author: Brian Peterson
@since: November 2009
@version: 0.9
@summary: High-level functions to read in text, lines or CSV data with support for iterative stripping & pruning.

Supports ASCII files in the following formats:
    - text (returns string)
    - lines (returns list)
    - CSV (returns list of lists)

Supports the following iterative stripping & pruning options:
    - lines
        - left and/or right 'strip' each line
    - csv_
        - 'strip' each CSV token
        - 'prune' each line where all CSV tokens are empty

"""

import csv
import sys
def text(file_name, fmt='text', sep=',', stop_at_line=None):
    """
    Load Text File
    ==============
    Formats
        - 'text': reads in file as string, contains EOL characters
        - 'lines': reads in file as a list of lines
        - 'csv': reads in a delimited text file as list of lists

    >>> text('test.txt')
    '"  T1","T2  "<a href="file://n1/">\\n1</a>, "1st line\\n2nd line",10/31\\n  3, "3rd line, 10/31"  <a href="file://n4/">\\n4</a>; 4th line; 10/31\\n,,,\\n, ,,\\n\\n'
    >>> text('test.txt', fmt='lines')
    ['"  T1","T2  "<a href="file://n'/">\\n'</a>, '1, "1st line\\n', '2nd line",10/31\\n', '  3, "3rd line, 10/31"  <a href="file://n'/">\\n'</a>, '4; 4th line; 10/31\\n', ',,,\\n', ', ,,\\n', '\\n']
    >>> text('test.txt', fmt='csv')
    [['  T1', 'T2  '], ['1', ' "1st line'], ['2nd line"', '10/31'], ['  3', ' "3rd line', ' 10/31"  '], ['4; 4th line; 10/31'], ['', '', '', ''], ['', ' ', '', ''], []]
    >>> text('test.txt', fmt='csv', sep=';')
    [['  T1,"T2  "'], ['1, "1st line'], ['2nd line",10/31'], ['  3, "3rd line, 10/31"  '], ['4', ' 4th line', ' 10/31'], [',,,'], [', ,,'], []]
    >>> text('test.txt', fmt='lines', stop_at_line=2)
    ['"  T1","T2  "<a href="file://n'/">\\n'</a>, '1, "1st line\\n']

    @param file_name: name of file to load
    @type file_name: String
    @keyword fmt: (string) 'text', 'lines' or 'csv'
    @keyword sep: (string) delimiter in CSV files
    @keyword stop_at_line: (none or integer) line num to stop reading file,
        only applies to 'lines' & 'csv' formats
    @returns: 'text': (string), 'lines': (list), 'csv': (list of lists)
        or none if there is an error reading the file
    @rtype: String: 'text', List: 'lines', List of Lists: 'csv', None if error reading file
    """

    try:
        f = open(file_name)
    except IOError, error_msg:
        sys.stderr.write('load.py: ERROR loading file "%s".\n%s\n' % (file_name, error_msg))
        return
    try:
        if fmt == 'lines':
            if stop_at_line is None:
                result = f.readlines()
            else:
                result = [ f.readline() for line in range(stop_at_line) ]
        elif fmt == 'csv':
            result = list(csv.reader(f, delimiter=sep))
            if stop_at_line is not None:
                result = result[:stop_at_line]
        else:
            result = f.read()
    finally:
        f.close()
    return result

def lines(file_name, strip=False, prune=False, stop_at_line=None):
    """
    Load text file lines into a list
    ================================

    >>> lines('test.txt')
    ['"  T1","T2  "<a href="file://n'/">\\n'</a>, '1, "1st line\\n', '2nd line",10/31\\n', '  3, "3rd line, 10/31"  <a href="file://n'/">\\n'</a>, '4; 4th line; 10/31\\n', ',,,\\n', ', ,,\\n', '\\n']
    >>> lines('test.txt', strip=True)
    ['"  T1","T2  "', '1, "1st line', '2nd line",10/31', '3, "3rd line, 10/31"', '4; 4th line; 10/31', ',,,', ', ,,', '']
    >>> lines('test.txt', strip='left')
    ['"  T1","T2  "<a href="file://n'/">\\n'</a>, '1, "1st line\\n', '2nd line",10/31\\n', '3, "3rd line, 10/31"  <a href="file://n'/">\\n'</a>, '4; 4th line; 10/31\\n', ',,,\\n', ', ,,\\n', '']
    >>> lines('test.txt', strip='right')
    ['"  T1","T2  "', '1, "1st line', '2nd line",10/31', '  3, "3rd line, 10/31"', '4; 4th line; 10/31', ',,,', ', ,,', '']
    >>> lines('test.txt', stop_at_line=2)
    ['"  T1","T2  "<a href="file://n'/">\\n'</a>, '1, "1st line\\n']

    @param file_name: (string) name of file to read
    @keyword strip: (boolean) True: strip whitespace from left & right, 'l' or 'left': strip whitespace from left, 'r' or 'right': strip whitespace from right
    @keyword prune:(boolean) removes blank lines
    @keyword stop_at_line: (int) line number at which to stop loading
    @returns: list of lines without EOL characters
    @rtype: List

    """
    lines = text(file_name, 'lines', stop_at_line=stop_at_line)
    if strip is False:
        return lines
    elif strip is True:
        lines = [ line.strip() for line in lines ]
    elif type(strip) == type(''):
        if len(strip):
            if strip[0] == 'l':
                lines = [ line.lstrip() for line in lines ]
            elif strip[0] == 'r':
                lines = [ line.rstrip() for line in lines ]
    if prune is True:
        lines = [ line for line in lines if line.strip() ]
    return lines 

def csv_(file_name, sep=',', strip= False, prune=False, stop_at_line=None):
    """
    Load Comma Separated Values (CSV) File
    ======================================

    >>> csv_('test.txt')
    [['  T1', 'T2  '], ['1', ' "1st line'], ['2nd line"', '10/31'], ['  3', ' "3rd line', ' 10/31"  '], ['4; 4th line; 10/31'], ['', '', '', ''], ['', ' ', '', ''], []]
    >>> csv_('test.txt', sep=';')
    [['  T1,"T2  "'], ['1, "1st line'], ['2nd line",10/31'], ['  3, "3rd line, 10/31"  '], ['4', ' 4th line', ' 10/31'], [',,,'], [', ,,'], []]
    >>> csv_('test.txt', strip=True)
    [['T1', 'T2'], ['1', '"1st line'], ['2nd line"', '10/31'], ['3', '"3rd line', '10/31"'], ['4; 4th line; 10/31'], ['', '', '', ''], ['', '', '', ''], []]
    >>> csv_('test.txt', prune=True)
    [['  T1', 'T2  '], ['1', ' "1st line'], ['2nd line"', '10/31'], ['  3', ' "3rd line', ' 10/31"  '], ['4; 4th line; 10/31']]
    >>> csv_(['a1,a2','b1, b2'])
    [['a1', 'a2'], ['b1', ' b2']]

    @param file_name: name of file to read or a list of strings
    @type file_name: String
    @keyword sep: (string) separator used by split()
    @type sep: String
    @keyword strip: - (boolean) strip() whitespace from each resulting value in each line
    @keyword prune: - (boolean) remove lines where all resulting values are empty strings
    @returns: list of resulting lines containing a list of separated values
    @rtype: List of Lists

    """
    if type(file_name) == type(''):
        lines = text(file_name, 'csv', sep)
    else:
        lines = list(csv.reader(file_name, delimiter=sep))
    if strip:
        result = []
        for line in lines:
          result.append([ item.strip() for item in line ])
        lines = result
    if prune:
        lines = [ line for line in lines if ''.join(line).strip() ]
    return lines

if __name__ == '__main__':
    import doctest
    doctest.testmod()

