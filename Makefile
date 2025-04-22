# Define directories
PYTHON_SRC_DIR = python
R_SRC_DIR = R
DOCS_SRC_DIR = docs
DATA_DIR = data

# Define targets
PYTHON_SCRIPTS = $(wildcard $(PYTHON_SRC_DIR)/*.py)
R_SCRIPTS = $(wildcard $(R_SRC_DIR)/*.R)
QUARTO_DOCS = $(wildcard $(DOCS_SRC_DIR)/*.qmd)
HTML_FILES = $(patsubst %.qmd, %.html, $(QUARTO_DOCS))
RUN_R_FLAG = run_R_done.flag

$(info ************************************)
$(info python source directory:   $(PYTHON_SRC_DIR))
$(info R Source director:         $(R_SRC_DIR))
$(info Docs Source director:      $(DOCS_SRC_DIR))
$(info R scripts:                 $(R_SCRIPTS))
$(info Python scripts:            $(PYTHON_SCRIPTS))
$(info Quarto scripts:            $(QUARTO_DOCS))
$(info HTML files:                $(HTML_FILES))
$(info ************************************)

# Default target
all: build_docker run_python run_R render_docs

build_docker:
	@echo "Building Docker image..."
	docker build --tag gcrmn_dev .

# Rule to run Python scripts
run_python: $(PYTHON_SCRIPTS)
	@echo "Running Python scripts..."
	python $^

# Rule to run the R targets pipeline
run_R:
	@echo "Running R targets pipeline..."
	cd R && Rscript -e "targets::tar_make()"

# %.html: %.qmd $(RUN_R_FLAG) $(PYTHON_SRC_DIR)/%.py
%.html: %.qmd
	@echo "Rendering $< to $@..."
	echo "library(quarto); quarto_render(\"$<\")" | R --no-save --no-restore;

# Rule to render all Quarto documents
render_docs: $(HTML_FILES)

# Clean up intermediate files
clean:
	@echo "Cleaning up..."
	rm -rf _targets/ $(RESULTS_DIR)/*.html $(RESULTS_DIR)/*.pdf $(RUN_R_FLAG)
