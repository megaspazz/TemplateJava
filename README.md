# TemplateGo
Go Template for Competitive Programming

## Setup

1.  Install the Go programming lanugage.  You can [download it here](https://golang.org/dl/).
2.  Optionally, install a text editor.  I like using [Visual Studio Code](https://code.visualstudio.com/Download).

## How to Use

1.  Put your files in the `src/` folder.  I prefer to put them in subfolders, like `src/codeforces` or `src/topcoder`.
2.  If your program reads from `stdin`, then you can use the template `src/#templates/template.go`, which has wrappers for reading `stdin` and printing to `stdout`.
3.  Put your input in `io/in.txt`.
4.  To build `some_program.go`, run the command `.\run some_program`.  Note that it does not include the `.go` extension.  It will find any file with the matching name in the `src/` folder, so avoid having duplicate names of source files, even if they are in different subdirectories.  This will build and run the source code, putting the executable in the `bin/` folder.
5.  The program output and build errors will be printed in the console, but you can also look at `io/out.txt` and `io/err.txt`.
