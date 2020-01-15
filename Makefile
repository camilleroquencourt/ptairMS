PKG_VERSION=$(shell grep '^Version:' DESCRIPTION | sed 's/^Version: //')
ZIPPED_PKG=ptairMS_$(PKG_VERSION).tar.gz

# Set testthat reporter
ifndef TESTTHAT_REPORTER
ifdef VIM
TESTTHAT_REPORTER=summary
else
TESTTHAT_REPORTER=progress
endif
endif

all: R/RcppExports.R
	R CMD SHLIB src/*.cpp src/*.c

check: $(ZIPPED_PKG) R/RcppExports.R
	time R CMD check --no-build-vignettes "$<"

test:
	R -q -e "devtools::test('$(CURDIR)', reporter = c('$(TESTTHAT_REPORTER)', 'fail'))"

$(ZIPPED_PKG) build: doc
	R CMD build .

doc: R/RcppExports.R
	R -q -e "devtools::document('$(CURDIR)')"

R/RcppExports.R: src/*.cpp
	R -q -e "Rcpp::compileAttributes('$(CURDIR)')"

install.deps:
	R -q -e "devtools::install_dev_deps('$(CURDIR)')"

clean:
	$(RM) src/*.o src/*.so src/*.dll

.PHONY: all clean build check doc test
