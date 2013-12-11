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


sphinxbuild_options = 
# -D extensions=['sphinxcontrib.napoleon','sphinx.ext.autodoc','sphinx.ext.viewcode',]
sphinx_build = sphinx-build $(sphinxbuild_options)


SVNREV = $(shell svn info | grep Revision | awk '{print $$2;}')

doc:
	@- sphinx-apidoc . --force  -H QuickReduce -R $(SVNREV) -o $(srcdir)

latex: doc
	@- $(sphinx_build) -b latex $(srcdir) $(build_latex)
# 	Run latex twice to get the table of contents right
	@- cd $(build_latex); pwd; pdflatex -interaction=nonstopmode QuickReduce.tex; 
	@- cd $(build_latex); pwd; pdflatex -interaction=nonstopmode QuickReduce.tex; 
	@- cp $(build_latex)/QuickReduce.pdf $(docdir)

html: doc
	@- sphinx-build -b html $(srcdir) $(build_html)

single_html: doc
	@- sphinx-build -b singlehtml $(srcdir) $(build_singlehtml)

all: doc html latex single_html

clean:
	rm -rf doc/build/*
