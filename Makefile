.PHONY: doc


nothing:
	@- echo
	@- echo "Run Makefile with a command: all/doc/html/latex/clean"
	@- echo

docdir = doc/
srcdir = doc/source/
build_latex = doc/build/latex/
build_html = doc/build/html/
build_singlehtml = doc/build/single_html/

doc:
	@- sphinx-apidoc . --force --full -H QuickReduce -A "Ralf Kotulla" -V "0.9" -R 1234.3 -o $(srcdir)

latex:
	@- sphinx-build -b latex $(srcdir) $(build_latex)
	@- cd $(build_latex); pwd; pdflatex -interaction=nonstopmode QuickReduce.tex; pdflatex -interaction=nonstopmode QuickReduce.tex
	@- cp $(build_latex)/QuickReduce.pdf $(docdir)

html:
	@- sphinx-build -b html $(srcdir) $(build_html)

single_html:
	@- sphinx-build -b singlehtml $(srcdir) $(build_singlehtml)

all: doc html latex single_html

clean:
	rm -rf doc/build/*
